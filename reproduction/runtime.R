# Portable execution helpers; estimator kernels are loaded from versioned packages.
argument <- function(args, key, default=NULL) {
  at <- match(key,args)
  if(is.na(at))return(default)
  if(at==length(args)||startsWith(args[at+1L],'--'))stop('Missing value for ',key)
  args[at+1L]
}
thread_limits <- function() {
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
             BLIS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
}
private_library <- function(root,study,create=FALSE) {
  path <- file.path(root,'.r-libraries',study)
  if(create)dir.create(path,recursive=TRUE,showWarnings=FALSE)
  if(dir.exists(path)).libPaths(c(path,.libPaths()))
  path
}
dependencies <- function(study) {
  if(study=='reference22')c('stats','parallel','glmnet','ncvreg','digest','mltools') else
    c('stats','parallel','Rcpp','RcppArmadillo','glmnet','ncvreg','mltools','digest','jsonlite')
}
check_sources <- function(root) {
  manifest<-read.csv(file.path(root,'source_manifest.csv'))
  for(i in seq_len(nrow(manifest))) {
    path<-file.path(root,manifest$path[i])
    if(!file.exists(path)||file.info(path)$size!=manifest$bytes[i]||
       digest::digest(file=path,algo='sha256')!=manifest$sha256[i])
      stop('Source checksum mismatch: ',manifest$path[i])
  }
  invisible(TRUE)
}
check_dependencies <- function(study,installed=TRUE) {
  packages <- dependencies(study)
  if(installed)packages<-c(packages,'softridge',if(study=='independent192')'scaleCorrection')
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop('Missing dependencies: ',paste(missing,collapse=', '),
    '. Install these explicitly before continuing; no automatic download is attempted.')
  if(installed) {
    expected<-if(study=='reference22')'0.4.0'else'0.3.0.1'
    stopifnot(as.character(packageVersion('softridge'))==expected)
    if(study=='independent192')stopifnot(as.character(packageVersion('scaleCorrection'))=='0.1.0')
  }
  invisible(TRUE)
}
install_study <- function(root,study) {
  check_dependencies(study,installed=FALSE)
  check_sources(root)
  lib <- private_library(root,study,create=TRUE)
  archives <- if(study=='reference22')'softridge_0.4.0.tar.gz' else
    c('softridge_0.3.0.1.tar.gz','scaleCorrection_0.1.0.tar.gz')
  Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
  for(name in archives) {
    path <- file.path(root,'packages',study,name)
    status<-system2(file.path(R.home('bin'),'R'),c('CMD','INSTALL','--no-multiarch',
      paste0('--library=',shQuote(lib)),shQuote(path)))
    if(status!=0L)stop('Package installation failed: ',name)
  }
  cat('Installed study packages into its separate project library.\n')
}
current_versions <- function(study) {
  packages<-c(dependencies(study),'softridge',if(study=='independent192')'scaleCorrection')
  list(R=as.character(getRversion()),platform=R.version$platform,packages=setNames(vapply(packages,
    function(p)as.character(packageVersion(p)),character(1)),packages))
}
check_recorded_versions <- function(study,allow=FALSE) {
  expected <- if(study=='reference22')c(glmnet='4.1-8',ncvreg='3.16.0.7',mltools='0.3.5') else
    c(Rcpp='1.1.1.1.1',RcppArmadillo='15.2.6.1',glmnet='5.0',ncvreg='3.16.0',mltools='0.3.5',digest='0.6.35',jsonlite='2.0.0')
  actual<-vapply(names(expected),function(p)as.character(packageVersion(p)),character(1))
  expected_R<-if(study=='reference22')'4.3.3'else'4.3.0'
  difference<-c(if(as.character(getRversion())!=expected_R)paste0('R ',expected_R,' -> ',getRversion()),
    paste0(names(expected)[actual!=expected],' ',expected[actual!=expected],' -> ',actual[actual!=expected]))
  if(length(difference)&&!allow)stop('Recorded-environment differences: ',paste(difference,collapse='; '),
    '. Use the recorded environment or explicitly pass --allow-version-differences; new runs remain separate.')
  invisible(difference)
}
atomic_rds <- function(object,path) {
  temporary<-paste0(path,'.pending-',Sys.getpid());saveRDS(object,temporary,version=3)
  if(!file.rename(temporary,path))stop('Atomic write failed')
}
csv17 <- function(object,path) {
  for(name in names(object))if(is.double(object[[name]]))object[[name]]<-vapply(object[[name]],
    function(v)if(is.na(v))NA_character_ else sprintf('%.17g',v),character(1))
  temporary<-paste0(path,'.pending-',Sys.getpid());write.csv(object,temporary,row.names=FALSE)
  if(!file.rename(temporary,path))stop('Atomic CSV write failed')
}
data_profile <- function(root) {
  path<-file.path(root,'reference22/data/riboflavin.RData')
  if(!file.exists(path))stop('Run fetch_riboflavin.R before requesting --real-data.')
  expected<-'ff5c2257a46367992cfb329041ffb5fef77663761811fef8c55043c274d8d2c5'
  if(digest::digest(file=path,algo='sha256')!=expected)stop('Riboflavin checksum mismatch; file retained.')
  envir<-new.env(parent=baseenv());load(path,envir=envir)
  x<-envir$riboflavin
  stopifnot(nrow(x$x)==71L,ncol(x$x)==4088L,length(x$y)==71L,
    length(unique(sub('.*(Fbat[0-9]+).*','\\1',rownames(x$x))))==28L)
  cat('Verified supplied riboflavin dimensions: 71 samples, 4088 genes, 28 batches.\n')
  invisible(TRUE)
}
