# Study orchestration around the unchanged estimator kernels.
REVISION_ROOT <- normalizePath('.')
.libPaths(c(file.path(REVISION_ROOT,'library'),.libPaths()))
source('code/common.R')
configuration <- function(mode) {
 stopifnot(mode %in% c('smoke','production'));x<-readRDS('config.rds');x$mode<-mode
 if(mode=='smoke') {
  x$tasks<-x$tasks[x$tasks$scenario%in%x$smoke_scenarios&x$tasks$rep==1L,]
  x$seeds<-x$seeds[x$tasks$id]
  x$execution$verification_chunk_size<-1L
 }
 rownames(x$tasks)<-NULL;x
}
study_grid <- function(X,y,cfg) {
 pr<-softridge:::.prep_design(X,y,TRUE,TRUE);sv<-softridge:::.svd_design(pr$X)
 k<-sv$d[1]^2*cfg$grid$k_ratios;B<-softridge:::ridge_path(pr$X,pr$y,k,sv)
 gamma<-t(vapply(seq_along(k),function(i)cfg$grid$fractions*max(abs(B[,i])),numeric(length(cfg$grid$fractions))))
 list(k=k,gamma=gamma,fraction=cfg$grid$fractions)
}

active_indices <- function(s,cfg) {
 stopifnot(cfg$signal$pattern=='beta_type_1',s$active_count%in%c(8L,40L),s$p>s$active_count)
 as.integer(round(seq(1,s$p,length.out=s$active_count)))
}
true_beta <- function(s,cfg) {
 b<-numeric(s$p);b[active_indices(s,cfg)]<-cfg$signal$amplitude;b
}
covariance <- function(s)s$rho^abs(outer(seq_len(s$p),seq_len(s$p),'-'))
# Dense covariance() is for independent validation only, never data production.
prepare_generation <- function(cfg)invisible(TRUE)
correlated_design <- function(Z,s)softridge:::cpp_ar_design(Z,s$rho)
ar_quad <- function(b,rho)softridge:::cpp_ar_quadratic(as.numeric(b),rho)
# Independent R implementation for validation; no call to the production C++ utility.
ar_quad_reference <- function(b,rho) {
 w<-as.numeric(stats::filter(b,filter=rho,method='recursive',init=0))
 sum(b*b)+2*rho*sum(b[-1L]*head(w,-1L))
}
generate_data <- function(task,cfg) {
 seed<-cfg$seeds[[task$id]];assign('.Random.seed',seed$data,.GlobalEnv)
 s<-cfg$scenarios[cfg$scenarios$scenario==task$scenario,];b<-true_beta(s,cfg)
 X<-correlated_design(matrix(rnorm(s$n*s$p),s$n),s)
 sigma<-sqrt(ar_quad(b,s$rho)/s$snr);y<-drop(X%*%b)+rnorm(s$n,sd=sigma)
 assign('.Random.seed',seed$fold,.GlobalEnv);fid<-sample(rep(seq_len(cfg$nfolds),length.out=s$n))
 list(X=X,y=y,beta=b,sigma=sigma,foldid=fid)
}
classify_fit <- function(f) {
 f$raw_coef<-f$coef
 for(m in c('SCAD','MCP'))if(!isTRUE(f$diagnostics[[m]]$selected_converged))f$status[m]<-'incomplete_solver_path'
 stopifnot(all(unlist(f$diagnostics$glmnet)==0))
 for(m in names(f$status)[f$status!='ok'])f$coef[[m]]<-rep(NA_real_,length(f$coef[[m]]))
 f
}
projection_target <- function(X,beta,q=20L) {
 pr<-softridge:::.prep_design(X,rep(0,nrow(X)),TRUE,TRUE);sv<-softridge:::.svd_design(pr$X);b<-pr$xs*beta
 th<-drop(sv$v%*%crossprod(sv$v,b));s<-sum(beta!=0);S<-order(-abs(th),seq_along(th))[seq_len(s)];truth<-integer(length(beta));truth[S]<-1
 q<-min(q,length(sv$d));thq<-drop(sv$v[,seq_len(q),drop=FALSE]%*%crossprod(sv$v[,seq_len(q),drop=FALSE],b))
 list(theta=th,truth=truth,xs=pr$xs,relative_distance=sqrt(sum((th-b)^2)/sum(b^2)),target_overlap=sum(truth!=0&beta!=0)/s,
 geometry=data.frame(rank=length(sv$d),q=q,lambda_q=sv$d[q]^2,lambda_next=if(q<length(sv$d))sv$d[q+1]^2 else NA_real_,
 theta_norm_sq=sum(th^2),standardized_beta_norm_sq=sum(b^2),target_overlap=sum(truth!=0&beta!=0)/s,
 relative_null_distance=sqrt(sum((th-b)^2)/sum(b^2)),relative_leading_distance=sqrt(sum((thq-b)^2)/sum(b^2))))
}

metrics_for <- function(z,task,cfg,reference=FALSE) {
 f<-z$fit;d<-z$data;s<-cfg$scenarios[cfg$scenarios$scenario==task$scenario,]
 target<-projection_target(d$X,d$beta)
 quad<-if(reference)ar_quad_reference else ar_quad
 do.call(rbind,lapply(names(f$coef),function(m){
  cf<-f$coef[[m]];b<-cf[-1];ok<-f$status[m]=='ok'
  base<-data.frame(id=task$id,stage='simulation',scenario=task$scenario,rep=task$rep,method=m,status=unname(f$status[m]),
   n=nrow(d$X),p=ncol(d$X),active_count=s$active_count,Size=if(ok)sum(b!=0)else NA_real_,time=unname(f$time[m]))
  met<-if(ok)c(PE=cf[1]^2+quad(b-d$beta,s$rho),L2=sqrt(sum((b-d$beta)^2)),
   theta_loss=sum((target$xs*b-target$theta)^2),MCC_beta=unname(support_metrics(b,d$beta)['MCC']),
   MCC_theta=unname(support_metrics(b,target$truth)['MCC']),MCC_theta_top_s=unname(support_metrics(b,target$truth)['MCC']),TPR=mean(b[d$beta!=0]!=0),FPR=mean(b[d$beta==0]!=0)) else
   setNames(rep(NA_real_,8),c('PE','L2','theta_loss','MCC_beta','MCC_theta','MCC_theta_top_s','TPR','FPR'))
  cbind(base,as.data.frame(as.list(met)))
 }))
}
run_task <- function(task,cfg) {
 d<-generate_data(task,cfg);g<-study_grid(d$X,d$y,cfg)
 f<-classify_fit(fit_all(d$X,d$y,foldid=d$foldid,grid=g))
 f<-add_oracles(f,d$X,d$y,d$beta,NULL,g,cfg$oracle_policy)
 z<-list(data=d,fit=f,grid=g);z$res<-metrics_for(z,task,cfg)
 z$geometry<-cbind(id=task$id,scenario=task$scenario,rep=task$rep,projection_target(d$X,d$beta)$geometry)
 z
}
validate_payload <- function(a,task,cfg) {
 z<-a$value;d<-z$data;f<-z$fit;s<-cfg$scenarios[cfg$scenarios$scenario==task$scenario,]
 stopifnot(identical(a$id,task$id),identical(d,generate_data(task,cfg)),identical(f$foldid,d$foldid),
  identical(names(f$coef),cfg$methods),identical(names(f$status),cfg$methods),
  all(f$status%in%c('ok','incomplete_solver_path','not_applicable_singular_support')),
  isTRUE(all.equal(z$grid,study_grid(d$X,d$y,cfg),tolerance=1e-12)),
  abs(d$sigma^2-ar_quad_reference(d$beta,s$rho)/s$snr)<1e-10*(1+d$sigma^2))
 for(m in cfg$methods)stopifnot(length(f$coef[[m]])==s$p+1L,if(f$status[m]=='ok')all(is.finite(f$coef[[m]]))else all(is.na(f$coef[[m]])),
  is.finite(f$time[m]),f$time[m]>=0)
 pr<-softridge:::.prep_design(d$X,d$y,TRUE,TRUE);ss<-svd(pr$X);kp<-ss$d>1e-10*max(ss$d)
 for(m in c('HR','HR_S','HardRidge','HardRidge_S')) {
  dd<-f$diagnostics[[m]];j<-which(dd$cv_error==min(dd$cv_error),arr.ind=TRUE)[1,]
  stopifnot(identical(dim(dd$cv_error),c(20L,41L)),all(is.finite(dd$cv_error)),identical(dd$foldid,d$foldid),dd$k==z$grid$k[j[1]])
  b<-drop(ss$v[,kp,drop=FALSE]%*%(ss$d[kp]/(ss$d[kp]^2+dd$k)*crossprod(ss$u[,kp,drop=FALSE],pr$y)))
  g<-z$grid$fraction[j[2]]*max(abs(b));th<-if(grepl('Hard',m))b*(abs(b)>g)else sign(b)*pmax(abs(b)-g,0)
  cc<-1;if(m%in%c('HR_S','HardRidge_S')){v<-drop(pr$X%*%th);cc<-if(sum(v^2)>0)max(0,sum(v*pr$y)/sum(v^2))else 1}
  th<-cc*th/pr$xs;expected<-c(pr$ym-sum(th*pr$xm),th)
  stopifnot(abs(dd$scale-cc)<1e-7*(1+abs(cc)),max(abs(expected-f$coef[[m]]))<1e-7*(1+max(abs(expected))),
   isTRUE(all.equal(dd$coef,as.numeric(f$coef[[m]]),check.attributes=FALSE)))
 }
 for(m in c('SCAD','MCP')) {
  dd<-f$diagnostics[[m]];available<-which(colSums(is.finite(dd$squared_errors))==nrow(d$X))
  complete<-length(available)==length(dd$lambda)&&all(rowSums(dd$iterations,na.rm=TRUE)<cfg$ncvreg$max_iter)&&!length(dd$warnings)
  stopifnot(identical(dd$omitted_path_indices,setdiff(seq_along(dd$lambda),available)),isTRUE(dd$selected_converged)==complete,(f$status[m]=='ok')==complete)
  if(length(available)) {
   ee<-colMeans(dd$squared_errors[,available,drop=FALSE]);j<-available[which.min(ee)]
   stopifnot(j==dd$index,max(abs(ee-dd$cv_error[available]))<1e-10,isTRUE(all.equal(f$raw_coef[[m]],as.numeric(dd$full_coef[,j]))))
  }
 }
 stopifnot(all(unlist(f$diagnostics$glmnet)==0),identical(f$diagnostics$enet$alpha_grid,cfg$enet_alpha))
 en<-f$diagnostics$enet;stopifnot(en$selected_alpha==en$alpha_grid[which.min(en$alpha_minimum)])
 for(m in c('Ridge','Lasso','ENet')){
  dd<-f$diagnostics$glmnet_surfaces[[m]]
  stopifnot(all(is.finite(dd$cvm)),dd$lambda_min==dd$lambda[which.min(dd$cvm)])
 }
 dd<-f$diagnostics$Oracle_OLS;qr_ref<-qr(cbind(1,d$X[,which(d$beta!=0),drop=FALSE]));available<-qr_ref$rank==sum(d$beta!=0)+1L
 stopifnot(dd$available==available,dd$rank==qr_ref$rank,(f$status['Oracle_OLS']=='ok')==available)
 if(available){cf<-numeric(s$p+1L);cf[c(1L,which(d$beta!=0)+1L)]<-qr.coef(qr_ref,d$y);stopifnot(max(abs(cf-f$coef$Oracle_OLS))<1e-9)}
 stopifnot(identical(z$geometry,cbind(id=task$id,scenario=task$scenario,rep=task$rep,projection_target(d$X,d$beta)$geometry)),
  isTRUE(all.equal(z$res,metrics_for(z,task,cfg,reference=TRUE),tolerance=1e-10)))
 # Independent projection and loss in exactly the same coefficient coordinates.
 xs_ref<-sqrt(colMeans(sweep(d$X,2,colMeans(d$X))^2));b_ref<-xs_ref*d$beta
 theta_ref<-drop(ss$v[,kp,drop=FALSE]%*%crossprod(ss$v[,kp,drop=FALSE],b_ref))
 target_ref<-projection_target(d$X,d$beta)
 stopifnot(max(abs(target_ref$theta-theta_ref))<1e-8*(1+max(abs(theta_ref))),
  identical(z$res$MCC_theta,z$res$MCC_theta_top_s))
 for(m in cfg$methods)if(f$status[m]=='ok') {
  expected<-sum((xs_ref*f$coef[[m]][-1]-theta_ref)^2)
  stopifnot(abs(z$res$theta_loss[z$res$method==m]-expected)<1e-8*(1+expected))
 }
 # Independent confusion-matrix MCC formula, including undefined-denominator policy.
 for(m in cfg$methods)if(f$status[m]=='ok') {
  selected<-f$coef[[m]][-1]!=0;truth<-d$beta!=0
  tp<-sum(selected&truth);tn<-sum(!selected&!truth);fp<-sum(selected&!truth);fn<-sum(!selected&truth)
  den<-sqrt(as.double(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));expected<-if(den==0)NA_real_ else (tp*tn-fp*fn)/den
  stopifnot(isTRUE(all.equal(unname(z$res$MCC_beta[z$res$method==m]),expected,tolerance=1e-12)))
 }
 TRUE
}
