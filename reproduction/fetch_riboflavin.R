#!/usr/bin/env Rscript
# Explicit optional acquisition; never called automatically by a simulation.
script<-sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))][1])
if(!file.exists(script))script<-gsub('~+~',' ',script,fixed=TRUE)
root<-dirname(normalizePath(script,mustWork=TRUE))
if(!requireNamespace('digest',quietly=TRUE))stop('Install the declared digest dependency first.')
expected<-'ff5c2257a46367992cfb329041ffb5fef77663761811fef8c55043c274d8d2c5'
dest<-file.path(root,'reference22/data/riboflavin.RData')
if(file.exists(dest)) {
  if(digest::digest(file=dest,algo='sha256')!=expected)stop('Existing file differs; refusing overwrite.')
  cat('Verified existing public data; no download needed.\n');quit(save='no')
}
dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
tmp<-tempfile('riboflavin-',tmpdir=dirname(dest))
url<-'https://raw.githubusercontent.com/cran/hdi/b31648aa8110030e868a778a79e6cea70bc5d719/data/riboflavin.RData'
tryCatch({
  download.file(url,tmp,mode='wb',quiet=FALSE)
  if(digest::digest(file=tmp,algo='sha256')!=expected)stop('Downloaded checksum mismatch.')
  if(file.exists(dest))stop('Destination appeared during download; refusing overwrite.')
  if(!file.rename(tmp,dest))stop('Could not publish verified data.')
},finally=if(file.exists(tmp))unlink(tmp))
cat('Verified public riboflavin data from pinned source.\n')
