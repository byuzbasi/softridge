# Run the reference formulas with their original grids, seeds and two fork streams.
run_reference <- function(root,args) {
  output<-argument(args,'--output')
  if(is.null(output))stop('Supply a new --output directory.')
  if(!grepl('^(/|[A-Za-z]:)',output))output<-file.path(root,output)
  exists<-dir.exists(output)
  if(exists&&!('--resume'%in%args))stop('Output exists: refusing overwrite. Use --resume only for the same recorded run.')
  if(.Platform$OS.type!='unix')
    stop('Reference reproducibility requires macOS/Linux with the original two fork-worker streams.')
  full<-'--all'%in%args
  if(full&&!('--confirm-full'%in%args))stop('Full reference study requires --confirm-full.')
  if('--real-data'%in%args)data_profile(root)
  if(!full&&!('--smoke'%in%args))stop('Choose --smoke or explicitly --all --confirm-full.')
  signature<-digest::digest(list(study='reference22',arguments=args[!args%in%'--resume'],
    source=read.csv(file.path(root,'source_manifest.csv')),versions=current_versions('reference22')),algo='sha256')
  if(exists)stopifnot(readRDS(file.path(output,'definition.rds'))$signature==signature)else {
    dir.create(output,recursive=TRUE)
    atomic_rds(list(signature=signature,versions=current_versions('reference22')),file.path(output,'definition.rds'))
  }
  # Run a copy in a new workspace so immutable supplied results cannot be overwritten.
  work<-file.path(output,'work');dir.create(work,showWarnings=FALSE)
  dir.create(file.path(work,'results'),showWarnings=FALSE);dir.create(file.path(work,'figures'),showWarnings=FALSE)
  if('--real-data'%in%args) {
    dir.create(file.path(output,'data'),showWarnings=FALSE)
    dest<-file.path(output,'data/riboflavin.RData')
    if(!file.exists(dest))stopifnot(file.copy(file.path(root,'reference22/data/riboflavin.RData'),dest))
    stopifnot(digest::digest(file=dest,algo='sha256')==
      digest::digest(file=file.path(root,'reference22/data/riboflavin.RData'),algo='sha256'))
  }
  for(name in c('simulation.R','realdata.R','bands.R'))
    if(!exists)file.copy(file.path(root,'reference22/code',name),file.path(work,name))
  for(name in c('simulation.R','realdata.R','bands.R'))stopifnot(
    digest::digest(file=file.path(work,name),algo='sha256')==
      digest::digest(file=file.path(root,'reference22/code',name),algo='sha256'))
  setwd(work)
  reference<-new.env(parent=.GlobalEnv)
  # The reference run used exactly two fork streams. Do not let a platform's
  # core-detection result change its random-number allocation.
  reference$detectCores<-function(...)2L
  source('simulation.R',local=reference)
  stopifnot(reference$NCORES==2L)
  original_evaluate<-reference$evaluate
  reference$evaluate<-function(cf,beta,Sigma) {
    selected<-cf[-1]!=0;truth<-beta!=0
    c(original_evaluate(cf,beta,Sigma),MCC=unname(mltools::mcc(
      TP=sum(selected&truth),FP=sum(selected&!truth),TN=sum(!selected&!truth),FN=sum(!selected&truth))))
  }
  if('--smoke'%in%args) {
    # An explicitly bounded correctness check, not a scientific performance study.
    X<-matrix(seq_len(12*16),12,16)/100
    X[,seq(2,16,2)]<-sin(X[,seq(2,16,2)])
    y<-seq(-1,1,length.out=12)
    fit<-softridge::hr(X,y,k=1,gamma=0)
    stopifnot(all(is.finite(fit$coef)))
    saveRDS(list(passed=TRUE,scope='fixed-matrix package smoke; no Monte Carlo replications'),
      file.path(output,'smoke_check.rds'))
    cat('Reference package smoke passed; original simulation functions parsed and loaded.\n')
    return(invisible(NULL))
  }
  selected<-strsplit(argument(args,'--studies','2,4,5'),',',fixed=TRUE)[[1]]
  stopifnot(all(selected%in%c('1','2','3','4','5')),!anyDuplicated(selected))
  if('1'%in%selected)reference$study1()
  if('2'%in%selected){z<-reference$study2();print(reference$summarise2(z));reference$plot2(z);reference$tables2(z)}
  if('4'%in%selected){z<-reference$study4();print(reference$summarise4(z));reference$plot4(z);reference$tables4(z)}
  if('5'%in%selected){z<-reference$study5();reference$table5(z)}
  if('3'%in%selected)reference$study3()
  if('--real-data'%in%args)source('realdata.R',local=reference)
  expected<-c('2'=8L,'4'=8L,'5'=6L)
  for(study in intersect(selected,names(expected))) {
    z<-readRDS(file.path(work,'results',paste0('study',study,'.rds')))
    stopifnot(nrow(z)==expected[study]*1500L,all(table(z$scenario)==1500L),
      all(table(z$method)==expected[study]*100L),all(is.finite(z$PE)),all(is.finite(z$MCC)),
      all(z$MCC>=-1&z$MCC<=1))
  }
  saveRDS(current_versions('reference22'),file.path(output,'software_versions.rds'))
  paths<-list.files(work,recursive=TRUE,full.names=TRUE)
  csv17(data.frame(path=substring(paths,nchar(output)+2L),bytes=file.info(paths)$size,
    sha256=vapply(paths,function(p)digest::digest(file=p,algo='sha256'),character(1))),
    file.path(output,'result_manifest.csv'))
  writeLines('completed requested reference analyses',file.path(output,'COMPLETED'))
}
