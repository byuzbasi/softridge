#!/usr/bin/env Rscript
script<-sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))][1])
if(!file.exists(script))script<-gsub('~+~',' ',script,fixed=TRUE)
root<-dirname(normalizePath(script,mustWork=TRUE))
source(file.path(root,'runtime.R'))
thread_limits()
args<-commandArgs(TRUE)
if(!length(args)||args[1]%in%c('help','--help')) {
  cat(paste(c('Portable reproduction commands (run from STR_SSTR_Software_Run_Code/):',
    'Rscript --vanilla RUNME.R install --study reference22',
    'Rscript --vanilla RUNME.R install --study independent192',
    'Rscript --vanilla RUNME.R check --study reference22',
    'Rscript --vanilla RUNME.R check --study independent192',
    'Rscript --vanilla RUNME.R reference22 --smoke --output new_reference_smoke --allow-version-differences',
    'Rscript --vanilla RUNME.R independent192 --task ars_s001_r001 --output new_unit --allow-version-differences',
    'A full study additionally requires --all --confirm-full; see RUN_GUIDE.txt.',
    'No fitting, installation or network access occurs without an explicit command.'),collapse='\n'),'\n')
  quit(save='no')
}
command<-args[1]
study<-if(command%in%c('reference22','independent192'))command else argument(args,'--study')
if(length(study)!=1L||!study%in%c('reference22','independent192'))stop('Specify --study reference22 or --study independent192.')
invisible(private_library(root,study))
if(command=='install') {install_study(root,study);quit(save='no')}
check_dependencies(study)
check_sources(root)
if(command=='check') {
  print(current_versions(study))
  if(study=='reference22') {
    if(file.exists(file.path(root,'reference22/data/riboflavin.RData')))data_profile(root) else
      cat('Simulation sources checked. For real data, first run fetch_riboflavin.R.\n')
  }
  if(study=='independent192') {
    cfg<-readRDS(file.path(root,'independent192/config.rds'))
    stopifnot(nrow(cfg$scenarios)==192L,nrow(cfg$tasks)==96000L,
      length(cfg$seeds)==96000L,length(cfg$grid$k_ratios)==20L,length(cfg$grid$fractions)==41L)
    cat('Verified 192 scenarios, 96000 task/seed identities and 20 x 41 tuning grid; no fitting.\n')
  }
  quit(save='no')
}
if(command=='verify') {
  source(file.path(root,'verify_run.R'))
  output<-argument(args,'--output');if(is.null(output))stop('Supply --output.')
  if(!grepl('^(/|[A-Za-z]:)',output))output<-file.path(root,output)
  verify_run(output)
  quit(save='no')
}
if(!command%in%c('reference22','independent192'))stop('Unknown command. Use --help.')
check_recorded_versions(study,allow='--allow-version-differences'%in%args)
if(command=='reference22') {
  source(file.path(root,'reference22/portable.R'))
  run_reference(root,args)
}else{
  source(file.path(root,'independent192/portable.R'))
  run_independent(root,args)
}
