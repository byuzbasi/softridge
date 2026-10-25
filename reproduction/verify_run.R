# Validate every committed pack and result against its signature and checksums.
verify_run <- function(output) {
  definition<-readRDS(file.path(output,'definition.rds'))
  manifest<-read.csv(file.path(output,'final_manifest.csv'))
  for(i in seq_len(nrow(manifest))) {
    path<-file.path(output,manifest$path[i])
    stopifnot(file.info(path)$size==manifest$bytes[i],
      digest::digest(file=path,algo='sha256')==manifest$sha256[i])
  }
  packs<-read.csv(file.path(output,'pack_manifest.csv'))
  ids<-character()
  for(i in seq_len(nrow(packs))) {
    pack<-readRDS(file.path(output,packs$path[i]))
    stopifnot(pack$signature==definition$signature,
      identical(pack$ids,vapply(pack$units,`[[`,character(1),'id')),
      length(pack$ids)==packs$units[i])
    for(unit in pack$units) {
      m<-unit$metrics
      stopifnot(nrow(m)==13L,all(m$id==unit$id),!anyDuplicated(m$method_key),
        all(is.finite(m$PE)),all(is.finite(m$theta_loss)),all(is.finite(m$MCC)),
        all(m$PE>=0),all(m$theta_loss>=0),all(m$MCC>=-1&m$MCC<=1))
    }
    ids<-c(ids,pack$ids)
  }
  metrics<-readRDS(file.path(output,'metrics.rds'))
  stopifnot(identical(ids,definition$ids),!anyDuplicated(ids),nrow(metrics)==length(ids)*13L,
    identical(unique(metrics$id),ids),readLines(file.path(output,'COMPLETED'))==definition$signature)
  cat('VERIFIED REPRODUCED RUN: ',length(ids),' units, ',nrow(metrics),
    ' method records; all packs, scientific dimensions and checksums passed.\n',sep='')
  invisible(TRUE)
}
