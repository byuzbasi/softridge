# One saved replication; two new procedures, with one shared corrected CV surface.
native_mcc <- function(coef,beta) {
 selected<-coef!=0;truth<-beta!=0
 unname(mltools::mcc(TP=sum(selected&truth),FP=sum(selected&!truth),
   TN=sum(!selected&!truth),FN=sum(!selected&truth)))
}

corrected_unit <- function(a,cfg,task=NULL) {
 z<-a$value;dat<-z$data
 if(is.null(task))task<-cfg$tasks[cfg$tasks$id==a$id,]
 stopifnot(nrow(task)==1L,identical(task$id,a$id),all(z$fit$status=='ok'),
   identical(z$fit$diagnostics$HR_S$foldid,dat$foldid))
 s<-cfg$scenarios[cfg$scenarios$scenario==task$scenario,]
 f<-scaleCorrection::corrected_pair(dat$X,dat$y,z$grid,dat$foldid,a$id)
 stopifnot(f$SSTR_c0$cv<=f$Ridge_c0$cv+1e-12*(1+f$Ridge_c0$cv),
   f$Ridge_c0$fraction==0,all(is.finite(f$cv_error)))
 theta<-drop(f$projection$v%*%crossprod(f$projection$v,f$projection$xs*dat$beta))
 row<-data.frame(id=a$id,scenario=task$scenario,rep=task$rep,
   n=s$n,p=s$p,s=s$active_count,rho=s$rho,SNR=s$snr)
 for(m in c('SSTR_c0','Ridge_c0')){
   fit<-f[[m]];cf<-fit$coef
   stopifnot(all(is.finite(cf)),fit$scale>=0,fit$scale<=fit$ls_scale+1e-10*(1+fit$ls_scale))
   values<-c(theta_loss=sum((f$projection$xs*cf[-1]-theta)^2),
     PE=cf[1]^2+softridge:::cpp_ar_quadratic(cf[-1]-dat$beta,s$rho),
     MCC=native_mcc(cf[-1],dat$beta),Size=sum(cf[-1]!=0),k=fit$k,
     k_index=fit$k_index,fraction=fit$fraction,gamma_index=fit$gamma_index,
     gamma=fit$gamma,scale=fit$scale,LS_scale=fit$ls_scale,df=fit$df,cv=fit$cv)
   stopifnot(all(is.finite(values)),values['theta_loss']>=0,values['PE']>=0)
   for(nm in names(values))row[[paste0(m,'_',nm)]]<-unname(values[nm])
 }
 if(f$SSTR_c0$fraction==0)stopifnot(identical(f$SSTR_c0$coef,f$Ridge_c0$coef))
 row$full_sigma2<-tail(f$variance,1)[[1]]$sigma2
 row$variance_floor_count<-sum(vapply(f$variance,`[[`,logical(1),'denominator_floor_used'))
 row$cv_surface_sha256<-digest::digest(f$cv_error,algo='sha256')
 # Store complete scores, selected coefficients and variance diagnostics for audit.
 list(row=row,detail=list(id=a$id,cv_error=f$cv_error,
   coef=list(SSTR_c0=f$SSTR_c0$coef,Ridge_c0=f$Ridge_c0$coef),variance=f$variance))
}

corrected_summary <- function(d) {
 do.call(rbind,lapply(split(d,d$scenario),function(z){
  base<-z[1,c('scenario','n','p','s','rho','SNR')]
  do.call(rbind,lapply(c('theta_loss','PE','MCC'),function(metric){
   control<-z[[paste0('Ridge_c0_',metric)]];focal<-z[[paste0('SSTR_c0_',metric)]]
   v<-if(metric=='MCC')focal-control else control-focal
   mu<-mean(v);se<-if(nrow(z)>1)sd(v)/sqrt(nrow(z))else NA_real_
   cbind(base,data.frame(metric=metric,n_pairs=nrow(z),control_mean=mean(control),
    control_MCSE=if(nrow(z)>1)sd(control)/sqrt(nrow(z))else NA_real_,
    SSTR_mean=mean(focal),SSTR_MCSE=if(nrow(z)>1)sd(focal)/sqrt(nrow(z))else NA_real_,
    gain_SSTR=mu,paired_MCSE=se,lower=mu-1.96*se,upper=mu+1.96*se,
    zero_threshold_n=sum(z$SSTR_c0_fraction==0),zero_scale_n=sum(z$SSTR_c0_scale==0),
    interpretation='positive gain favors SSTR-c0; descriptive MC interval; unadjusted'))
  }))
 }))
}

validate_corrected_details <- function(rows,details) {
 stopifnot(length(details)==nrow(rows))
 for(i in seq_len(nrow(rows))){
  z<-details[[i]];r<-rows[i,]
  stopifnot(z$id==r$id,identical(dim(z$cv_error),c(20L,41L)),
   all(is.finite(z$cv_error)),digest::digest(z$cv_error,algo='sha256')==r$cv_surface_sha256)
  best<-which(z$cv_error==min(z$cv_error),arr.ind=TRUE)[1,]
  stopifnot(r$SSTR_c0_k_index==best[1],r$SSTR_c0_gamma_index==best[2],
    r$Ridge_c0_k_index==which.min(z$cv_error[,1]),
    r$SSTR_c0_cv==min(z$cv_error),r$Ridge_c0_cv==min(z$cv_error[,1]),
    length(z$variance)==11L,all(vapply(z$variance,function(v)is.finite(v$sigma2)&&v$sigma2>=0,logical(1))))
  for(m in c('SSTR_c0','Ridge_c0'))stopifnot(length(z$coef[[m]])==r$p+1,
    all(is.finite(z$coef[[m]])),sum(z$coef[[m]][-1]!=0)==r[[paste0(m,'_Size')]],
    r[[paste0(m,'_scale')]]>=0,
    r[[paste0(m,'_scale')]]<=r[[paste0(m,'_LS_scale')]]+1e-10*(1+r[[paste0(m,'_LS_scale')]]))
 }
 TRUE
}
