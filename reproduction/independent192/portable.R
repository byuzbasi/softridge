# Orchestration and final metric conventions around the frozen study functions.
load_independent <- function(root) {
  setwd(file.path(root,'independent192'))
  source('code/study.R',local=.GlobalEnv)
  source('code/control_core.R',local=.GlobalEnv)
  source('code/corrected_core.R',local=.GlobalEnv)
  source('code/queue.R',local=.GlobalEnv)
  readRDS('config.rds')
}
reproduce_unit <- function(task,cfg) {
  original<-run_task(task,cfg)
  payload<-list(id=task$id,value=original)
  validate_payload(payload,task,cfg)
  control<-control_unit(payload,cfg,task)
  corrected<-corrected_unit(payload,cfg,task)
  stopifnot(validate_corrected_details(corrected$row,list(corrected$detail)))
  data<-original$data
  target<-projection_target(data$X,data$beta)
  s<-cfg$scenarios[cfg$scenarios$scenario==task$scenario,]
  coefficients<-c(original$fit$coef,list(Ridge_LS=control$coef),corrected$detail$coef)
  labels<-c(Ridge='Ridge',Lasso='Lasso',ENet='ENet',SCAD='SCAD',MCP='MCP',
    HardRidge='HTR',HR='STR',HR_S='SSTR-LS',HardRidge_S='SHTR-LS',
    Oracle_OLS='Oracle-OLS',Ridge_LS='LS-scaled Ridge',SSTR_c0='SSTR',Ridge_c0='Ridge-c0')
  rows<-do.call(rbind,lapply(names(coefficients),function(method) {
    cf<-coefficients[[method]]
    data.frame(id=task$id,scenario=task$scenario,rep=task$rep,n=s$n,p=s$p,s=s$active_count,
      rho=s$rho,SNR=s$snr,method=unname(labels[method]),method_key=method,status='ok',
      theta_loss=sum((target$xs*cf[-1]-target$theta)^2),
      PE=cf[1]^2+ar_quad(cf[-1]-data$beta,s$rho),MCC=native_mcc(cf[-1],data$beta),
      Size=sum(cf[-1]!=0),oracle=method=='Oracle_OLS')
  }))
  stopifnot(nrow(rows)==13L,all(is.finite(rows$theta_loss)),all(is.finite(rows$PE)),
    all(is.finite(rows$MCC)),all(rows$MCC>=-1&rows$MCC<=1))
  list(id=task$id,metrics=rows,
    input_sha256=digest::digest(data,algo='sha256'),corrected=corrected$row,
    corrected_cv_sha256=corrected$row$cv_surface_sha256)
}
run_independent <- function(root,args) {
  cfg<-load_independent(root)
  id<-argument(args,'--task')
  full<-'--all'%in%args
  if(is.null(id)&&!full)stop('Choose one --task ID, or explicitly use --all --confirm-full.')
  if(full&&!('--confirm-full'%in%args))stop('A full study requires --confirm-full.')
  stopifnot(is.null(id)||!full)
  tasks<-if(full)cfg$tasks else cfg$tasks[cfg$tasks$id==id,,drop=FALSE]
  if(!nrow(tasks))stop('Unknown task ID; consult independent192/tasks.csv.')
  output<-argument(args,'--output')
  if(is.null(output))stop('Supply a new --output directory.')
  if(!grepl('^(/|[A-Za-z]:)',output))output<-file.path(root,output)
  workers<-as.integer(argument(args,'--workers','1'))
  if(is.na(workers)||workers<1L||workers>55L)stop('workers must be between 1 and 55.')
  if(.Platform$OS.type!='unix')stop('The frozen worker helper requires macOS/Linux; use Linux/WSL on Windows.')
  signature<-digest::digest(list(config=cfg,ids=tasks$id,
    source=read.csv(file.path(root,'source_manifest.csv')),
    versions=current_versions('independent192'),workers=workers,metric_policy='native mltools 0.3.5'),algo='sha256')
  if(dir.exists(output)) {
    if(!('--resume'%in%args))stop('Output exists: refusing overwrite; use --resume only for this same run.')
    stopifnot(readRDS(file.path(output,'definition.rds'))$signature==signature)
  }else {
    dir.create(output,recursive=TRUE)
    atomic_rds(list(signature=signature,ids=tasks$id,workers=workers,
      versions=current_versions('independent192')),file.path(output,'definition.rds'))
  }
  lock<-file.path(output,'.lock')
  if(!dir.create(lock,showWarnings=FALSE))stop('Output is locked by another run. Verify inactivity before removing the lock.')
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  N<-nrow(tasks);B<-min(10L,N);batch<-ceiling(seq_len(N)*B/N)
  manifest_file<-file.path(output,'pack_manifest.csv')
  manifest<-if(file.exists(manifest_file))read.csv(manifest_file)else
    data.frame(batch=integer(),path=character(),bytes=double(),sha256=character(),units=integer())
  stopifnot(identical(as.integer(manifest$batch),seq_len(nrow(manifest))))
  for(j in seq_len(nrow(manifest))) {
    path<-file.path(output,manifest$path[j]);pack<-readRDS(path)
    stopifnot(file.info(path)$size==manifest$bytes[j],digest::digest(file=path,algo='sha256')==manifest$sha256[j],
      pack$signature==signature,identical(pack$ids,tasks$id[batch==j]))
  }
  committed<-sum(manifest$units);done<-committed;running<-0L;failed<-0L;last<-'';error_message<-''
  start<-tick<-Sys.time();initial<-done
  progress<-function(force=FALSE,phase='fitting') {
    now<-Sys.time();elapsed<-as.numeric(difftime(now,start,units='secs'))
    if(!force&&as.numeric(difftime(now,tick,units='secs'))<30)return()
    tick<<-now;rate<-if(done>initial)(done-initial)/max(elapsed,.001)else NA_real_
    eta<-if(done==N)0 else if(is.finite(rate))(N-done)/rate else NA_real_
    state<-list(run_id=basename(output),phase=phase,total=N,validated_completed=done,
      work_unit='one original data/fold seed; all thirteen procedures fitted and metrics validated',
      durably_committed=committed,running=running,failed=failed,pending=max(0,N-done-running-failed),
      percent=100*done/N,elapsed_seconds=elapsed,throughput_per_second=rate,
      computational_eta_seconds=eta,eta_scope='fitting only; final aggregation and verification excluded',
      estimated_finish_utc=if(is.finite(eta))format(now+eta,tz='UTC',usetz=TRUE)else'estimating',
      heartbeat_utc=format(now,tz='UTC',usetz=TRUE),last_completed_unit=last,workflow_error=error_message)
    tmp<-file.path(output,'progress.pending.json');jsonlite::write_json(state,tmp,auto_unbox=TRUE,pretty=TRUE,na='null')
    stopifnot(file.rename(tmp,file.path(output,'progress.json')))
    history<-file.path(output,'progress.tsv')
    write.table(data.frame(timestamp=state$heartbeat_utc,phase=phase,validated=done,committed=committed,
      total=N,running=running,failed=failed,elapsed=elapsed,rate=rate,eta=eta,last=last),history,
      sep='\t',row.names=FALSE,col.names=!file.exists(history),append=file.exists(history),quote=FALSE)
    cat(sprintf('%s %s validated=%d/%d committed=%d running=%d failed=%d ETA fitting=%s seconds\n',
      state$heartbeat_utc,phase,done,N,committed,running,failed,if(is.finite(eta))round(eta)else'estimating'));flush.console()
  }
  progress(TRUE)
  for(j in setdiff(seq_len(B),manifest$batch)) {
    ii<-which(batch==j);buffer<-vector('list',length(ii))
    tryCatch(run_queue(seq_along(ii),workers,
      work=function(k)reproduce_unit(tasks[ii[k],,drop=FALSE],cfg),
      accept=function(k,z,seconds){buffer[[k]]<<-z;done<<-done+1L;last<<-z$id},
      status=function(event,active,pending,i){running<<-length(active);progress()}),
      error=function(e){failed<<-1L;running<<-0L;error_message<<-conditionMessage(e);progress(TRUE,'failed');stop(e)})
    stopifnot(identical(vapply(buffer,`[[`,character(1),'id'),tasks$id[ii]))
    name<-sprintf('part_%02d.rds',j);path<-file.path(output,name)
    if(file.exists(path))stop('Uncommitted pack exists; preserve it and choose a new output directory.')
    atomic_rds(list(signature=signature,ids=tasks$id[ii],units=buffer),path)
    manifest<-rbind(manifest,data.frame(batch=j,path=name,bytes=file.info(path)$size,
      sha256=digest::digest(file=path,algo='sha256'),units=length(ii)))
    csv17(manifest,manifest_file);committed<-done;progress(TRUE,'pack_committed')
  }
  progress(TRUE,'aggregation')
  metrics<-do.call(rbind,lapply(manifest$path,function(name) {
    pack<-readRDS(file.path(output,name));do.call(rbind,lapply(pack$units,`[[`,'metrics'))
  }))
  stopifnot(nrow(metrics)==N*13L,length(unique(metrics$id))==N)
  csv17(metrics,file.path(output,'metrics.csv'))
  atomic_rds(metrics,file.path(output,'metrics.rds'))
  paths<-c('definition.rds','pack_manifest.csv','metrics.csv','metrics.rds',manifest$path)
  csv17(data.frame(path=paths,bytes=file.info(file.path(output,paths))$size,
    sha256=vapply(file.path(output,paths),function(p)digest::digest(file=p,algo='sha256'),character(1))),
    file.path(output,'final_manifest.csv'))
  writeLines(signature,file.path(output,'COMPLETED'))
  source(file.path(root,'verify_run.R'))
  verify_run(output)
  progress(TRUE,'COMPLETED')
  cat('Verified ',N,' saved-seed replications, ',nrow(metrics),' method records.\n',sep='')
}
