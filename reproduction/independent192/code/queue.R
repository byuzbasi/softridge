# Local fork workers compute only; the parent alone validates and writes files.
run_queue <- function(indices, workers, work, accept, status, stop_after_one=FALSE, should_stop=function()FALSE) {
  stopifnot(workers>=1L,workers<=56L,!anyDuplicated(indices))
  active<-list();pending<-indices;accepted<-0L;peak<-0L;stop_requested<-FALSE
  clean<-function() {
    if(length(active)) {
      for(a in active)try(tools::pskill(a$job$pid),silent=TRUE)
      try(parallel::mccollect(lapply(active,`[[`,'job'),wait=FALSE),silent=TRUE)
      active<<-list()
    }
  }
  on.exit(clean(),add=TRUE)
  snapshot<-function(event,i=NA_integer_)status(event,lapply(active,function(a)list(index=a$index,pid=a$job$pid,start=a$start)),pending,i)
  launch_pending<-function() {
    if(isTRUE(should_stop())){stop_requested<<-TRUE;pending<<-integer()}
    while(length(pending)&&length(active)<workers) {
      i<-pending[1];pending<<-pending[-1];start<-Sys.time()
      job<-parallel::mcparallel(tryCatch(list(ok=TRUE,value=work(i)),
        error=function(e)list(ok=FALSE,error=conditionMessage(e))),mc.set.seed=FALSE,silent=TRUE)
      active[[as.character(job$pid)]]<<-list(job=job,index=i,start=start)
      peak<<-max(peak,length(active));snapshot('launch',i)
    }
  }
  while(length(pending)||length(active)) {
    launch_pending()
    if(!length(active))break
    ready<-parallel::mccollect(lapply(active,`[[`,'job'),wait=FALSE)
    if(is.null(ready)){snapshot('heartbeat');Sys.sleep(.1);next}
    # Collected children are already reaped: never send signals to their old PIDs.
    meta<-active[names(ready)];active[names(ready)]<-NULL
    for(pid in names(ready)) {
      ans<-ready[[pid]];a<-meta[[pid]]
      if(!is.list(ans)||!isTRUE(ans$ok)) {
        snapshot('failure',a$index)
        stop(sprintf('Task index %d failed: %s',a$index,
          if(is.list(ans)&&!is.null(ans$error))ans$error else 'worker exited without a result'))
      }
      accept(a$index,ans$value,as.numeric(difftime(Sys.time(),a$start,units='secs')))
      accepted<-accepted+1L;snapshot('accepted',a$index)
      if(stop_after_one){clean();snapshot('stopped',a$index);return(invisible(list(accepted=accepted,peak=peak)))}
      # Keep workers occupied while the parent validates other collected results.
      launch_pending()
    }
  }
  snapshot('drained');invisible(list(accepted=accepted,peak=peak,stop_requested=stop_requested))
}
