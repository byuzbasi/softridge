# Matched gamma=0 ablation: stored scaled CV scores; no new folds or CV calls.
control_unit <- function(a,cfg,task=NULL) {
 z<-a$value;dat<-z$data;f<-z$fit;dd<-f$diagnostics$HR_S
 E<-dd$cv_error;K<-z$grid$k;fractions<-z$grid$fraction
 stopifnot(is.matrix(E),identical(dim(E),c(20L,41L)),all(is.finite(E)),
  length(K)==20L,length(fractions)==41L,fractions[1]==0,
  identical(dd$foldid,dat$foldid),identical(dd$grid$k,K),
  all(f$status=='ok'),length(unique(dat$foldid))==10L)
 best<-which(E==min(E),arr.ind=TRUE)[1,];j<-which.min(E[,1])
 stopifnot(dd$k==K[best[1]],abs(dd$min_error-min(E))<1e-10*(1+abs(min(E))))
 if(is.null(task))task<-cfg$tasks[cfg$tasks$id==a$id,]
 stopifnot(nrow(task)==1L,identical(task$id,a$id))
 s<-cfg$scenarios[cfg$scenarios$scenario==task$scenario,]
 old<-z$res[z$res$method=='HR_S',];stopifnot(nrow(old)==1L)
 reused<-best[2]==1L
 if(reused) {
  stopifnot(j==best[1],dd$gamma==0)
  cf<-f$coef$HR_S;cc<-dd$scale;loss<-old$theta_loss;pe<-old$PE
 } else {
  pr<-softridge:::.prep_design(dat$X,dat$y,TRUE,TRUE)
  sv<-softridge:::.svd_design(pr$X)
  b<-drop(softridge:::ridge_path(pr$X,pr$y,K[j],sv))
  cc<-softridge:::.ls_scale(pr$X,pr$y,b)[1];u<-cc*b
  cf<-c(pr$ym-sum((u/pr$xs)*pr$xm),u/pr$xs)
  theta<-drop(sv$v%*%crossprod(sv$v,pr$xs*dat$beta))
  loss<-sum((u-theta)^2)
  pe<-cf[1]^2+softridge:::cpp_ar_quadratic(cf[-1]-dat$beta,s$rho)
  pred<-drop(pr$X%*%b);cc_ref<-if(sum(pred^2)>0)max(0,sum(pr$y*pred)/sum(pred^2))else 1
  stopifnot(abs(cc-cc_ref)<1e-10*(1+abs(cc_ref)))
 }
 support<-function(b){
  selected<-b!=0;truth<-dat$beta!=0
  mltools::mcc(TP=sum(selected&truth),FP=sum(selected&!truth),
   TN=sum(!selected&!truth),FN=sum(!selected&truth))
 }
 mcc<-support(cf[-1]);sstr_mcc<-support(f$coef$HR_S[-1])
 stopifnot(all(is.finite(cf)),all(is.finite(c(loss,pe,mcc))),loss>=0,pe>=0,
  min(E)<=E[j,1]+1e-12*(1+abs(E[j,1])))
 if(reused)stopifnot(identical(cf,f$coef$HR_S),loss==old$theta_loss,pe==old$PE,mcc==sstr_mcc)
 row<-data.frame(id=a$id,scenario=task$scenario,rep=task$rep,n=s$n,p=s$p,s=s$active_count,
  rho=s$rho,SNR=s$snr,k_control=K[j],scale_control=cc,control_cv=E[j,1],SSTR_cv=min(E),
  SSTR_zero_threshold=reused,reused_SSTR=reused,control_Size=sum(cf[-1]!=0),
  control_theta_loss=unname(loss),control_PE=unname(pe),control_MCC=unname(mcc),
  SSTR_theta_loss=old$theta_loss,SSTR_PE=old$PE,SSTR_MCC=unname(sstr_mcc),
  theta_gain=unname(loss-old$theta_loss),PE_gain=unname(pe-old$PE),MCC_gain=unname(sstr_mcc-mcc))
 list(row=row,coef=cf,k_index=j)
}

control_summary <- function(d) {
 do.call(rbind,lapply(split(d,d$scenario),function(z){
  base<-z[1,c('scenario','n','p','s','rho','SNR')]
  do.call(rbind,lapply(c('theta_loss','PE','MCC'),function(metric){
   control<-z[[paste0('control_',metric)]];focal<-z[[paste0('SSTR_',metric)]]
   v<-z[[if(metric=='theta_loss')'theta_gain'else paste0(metric,'_gain')]]
   mu<-mean(v);se<-if(nrow(z)>1)sd(v)/sqrt(nrow(z))else NA_real_
   cbind(base,data.frame(metric=metric,n_pairs=nrow(z),control_mean=mean(control),
    control_MCSE=if(nrow(z)>1)sd(control)/sqrt(nrow(z))else NA_real_,
    SSTR_mean=mean(focal),SSTR_MCSE=if(nrow(z)>1)sd(focal)/sqrt(nrow(z))else NA_real_,
    gain_SSTR=mu,paired_MCSE=se,lower=mu-1.96*se,upper=mu+1.96*se,
    zero_threshold_n=sum(z$SSTR_zero_threshold),reused_n=sum(z$reused_SSTR),
    interpretation='positive gain favors SSTR; descriptive MC interval; unadjusted'))
  }))
 }))
}
