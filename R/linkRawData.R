checkDirs <- function(b)
{
	return(file.exists(b))
}

verifyDataDirs <- function(map, platform='HumanMethylation450', path='/export/uec-gs1/laird/shared/production/methylation', version='0') 
{ # {{{
  if(length(unique(map$diseaseabr) == 1)) disease = toupper(map$diseaseabr[1])
  else stop("You have multiple (or no) disease abbreviations in your map!")
  raw.path = paste(path, ifelse(grepl('HumanMethylation450',platform,ignore=T),
                                'meth450k', 'meth27k'), 'raw', disease, sep='/')
  dir.create(raw.path, recursive=TRUE)
  tcga.path = paste(path, ifelse(grepl('HumanMethylation450',platform,ignore=T),
                                 'meth450k','meth27k'), 'tcga',disease, sep='/')
  dir.create(tcga.path, recursive=TRUE)
  dirs = list(disease=raw.path, archive=tcga.path)
  return(dirs)
} # }}}

linkRawData <- function(map, unlink.old.files=TRUE, platform='HumanMethylation450', path='/export/uec-gs1/laird/shared/production/methylation', logfile=NULL, idatPath=NULL) 
{ # {{{
  diseasedir = verifyDataDirs(map, platform, path)$disease
  oldwd = getwd()
  setwd(diseasedir)
  stopifnot('barcode' %in% names(map))
  rownames(map) = map$barcode
  platform.path = paste(path, 
                        ifelse(grepl('HumanMethylation450',platform,ignore=T),
                               'meth450k', 'meth27k'), 'raw', sep='/')
  if(is.null(idatPath)){
    check <- sapply(substr(map$barcode,1,10), function(x){checkDirs(paste(platform.path, x, sep='/'))})
    if(!all(check)){
      failed <- paste(unique(names(check[which(check == FALSE)])), collapse=",")
      stop(paste("Barcodes", failed, "were not found"))
    }
  }
  if(unlink.old.files==TRUE) unlink(list.files()) # get your old crap outta here
  if(!is.null(logfile)) unlink(logfile)
  for( beadchip in unique(substr(map$barcode, 1, 10))) {
    which.samples = which(substr(map$barcode, 1, 10) == beadchip)
    beadchip.batch = paste(unique(map$TCGA.BATCH[which.samples]),collapse='/')
    batches = paste('batch', beadchip.batch)
    msg = paste('linking samples from beadchip', beadchip,'(',batches,')')
    if(!is.null(logfile)) cat(msg,"\n",file=logfile,append=T) else message(msg)
    raw.path = ifelse(is.null(idatPath), paste('..', beadchip, sep='/'), idatPath)
    cmd = paste('ln -f ', paste(raw.path, '/', beadchip, '.sdf . 2>&1', sep=''))
    if(!is.null(logfile)) cat('#', cmd, file=logfile, append=T)
    msg = system(command=cmd, intern=TRUE)
    if(!is.null(logfile) & !is.null(msg)) cat(msg,"\n",file=logfile,append=T) 
    for( chip in grep(beadchip, map$barcode, value=T) ) {
      chip.path = ifelse(is.null(idatPath), paste(raw.path, chip, sep='/'), paste(idatPath, beadchip, sep="/"))
      cmd = paste('ln ', chip.path, '_*.idat . 2>&1', sep='')
      # if(!is.null(logfile)) cat('#', cmd, file=logfile, append=T)
      msg = system(command=cmd, intern=TRUE)
      # if(!is.null(logfile)&!is.null(msg)) cat(msg,"\n",file=logfile,append=T)
      # else message(msg)
    }
  }
  disease = unique(map$diseaseabr)
  mapfile = paste(disease, 'mappings', 'csv', sep='.')
  if(unlink.old.files) write.csv(map, file=mapfile)
  msg = paste('Linked files for', disease, 'and wrote map to', mapfile)
  if(!is.null(logfile)) cat(msg, "\n", file=logfile, append=T) else message(msg)
  setwd(oldwd)
} # }}}
