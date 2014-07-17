## build levels 2 & 3 by batch
buildArchive2 <- function(map, base=NULL, platform='HumanMethylation450', parallel=FALSE)
{
  x <- map
  map <- pData(x) 
  stopifnot('TCGA.ID' %in% names(map))
  stopifnot('TCGA.BATCH' %in% names(map))
  sampleNames(x) <- map$TCGA.ID
  if(!('BATCH.ID' %in% names(map))) { # {{{
    map$BATCH.ID <- as.numeric(as.factor(map$TCGA.BATCH))
  } # }}}
  stopifnot('diseaseabr' %in% names(map))
  bs <- unique(map$BATCH.ID)
  bs <- bs[order(bs)]
  
  message('Writing level 2 and level 3 data...')
  
  for(b in bs) {
    packL2(x, b, base, parallel)
    packL3(x, b, base, parallel)
  }
}


packL2 <- function(x, batch.id, base=NULL, parallel=F)
{
  disease <- unique(x$diseaseabr)
  platform <- gsub('k$','',gsub('Illumina','',annotation(x)))
  if(is.null(base)) { 
    message('Base directory not provided.. Writing to current working directory..')
    base <- getwd()
  } 
  
  diseasestub <- paste('jhu-usc.edu', disease, sep='_')
  stopifnot('TCGA.ID' %in% varLabels(x))
  if(!identical(sampleNames(x), x$TCGA.ID)) sampleNames(x) <- x$TCGA.ID
  pkged <- file.path(base, "Level_2")
  message(paste('Creating directory', pkged, '...'))
  if(!file.exists(pkged)) dir.create(pkged)
  setwd(pkged)
  
  stopifnot('TCGA.BATCH' %in% varLabels(x))
  batchname <- unique(x$TCGA.BATCH[ which(x$BATCH.ID == batch.id) ])
  stopifnot(length(batchname)==1)
  subjects <- sampleNames(x)[ which(x$BATCH.ID == batch.id) ]
  
  write.level2 <- function(s) {		
    xs <- which(sampleNames(x) == s)
    message(paste("Writing level 2 data for sample", which(subjects == s),
                  "of", length(subjects), "in TCGA batch", batchname))
    lvl2data <- data.frame( M=methylated(x)[,xs], 
                            U=unmethylated(x)[,xs],
                            P=pvals(x)[,xs] )
    rownames(lvl2data) <- featureNames(x)
    dump.file <- paste(paste(diseasestub,platform,batch.id,'lvl-2',s,'txt',sep='.'))
    headers1 <- paste('Hybridization REF', s, s, s, sep="\t")
    headers2 <- paste('Composite Element REF', 'Methylated_Intensity', 
                      'Unmethylated_Intensity', 'Detection_P_value', sep="\t")
    cat(headers1, "\n", sep='', file=dump.file)
    cat(headers2, "\n", sep='', file=dump.file, append=TRUE)
    write.table(lvl2data, file=dump.file, append=TRUE, quote=FALSE,
                row.names=TRUE, col.names=FALSE, sep="\t")
    return('done')
  }
  
  if(parallel) {
    if(!require(parallel)) require(multicore)
    results <- mclapply(subjects, write.level2)
  } else { 
    results <- lapply(subjects, write.level2)
  }
  gc(,T)
}

packL3 <- function(x, batch.id, base=NULL, parallel=F)
{
  disease <- unique(x$diseaseabr)
  platform <- gsub('k$','',gsub('Illumina','',annotation(x)))
  if(is.null(base)) { 
    message('Base directory not provided.. Writing to current working directory..')
    base <- getwd()
  }
  
  diseasestub <- paste('jhu-usc.edu', disease, sep='_')
  stopifnot('TCGA.ID' %in% varLabels(x))
  if(!identical(sampleNames(x), x$TCGA.ID)) sampleNames(x) <- x$TCGA.ID
  pkged <- file.path(base, "Level_3")
  
  message(paste('Creating directory', pkged, '...'))
  if(!file.exists(pkged)) dir.create(pkged)
  setwd(pkged)
  
  batchname <- unique(x$TCGA.BATCH[ which(x$BATCH.ID == batch.id) ])
  stopifnot(length(batchname)==1)
  subjects <- sampleNames(x)[ which(x$BATCH.ID == batch.id) ]
  
  betas(x) <- methylated(x)/total.intensity(x)
  
  if( !('mask' %in% fvarLabels(x))){
    ifelse(platform == 'HumanMethylation27k', data(probesToMask.27k), data(probesToMask))
    fData(x)$mask <- 0
    fData(x)$mask[ which(featureNames(x) %in% names(toMask))] <- 1
  }
  
  betas(x)[ which(fData(x)$mask==1), ] <- NA
  l3headers <- c('Composite Element REF','Beta_value',
                 'Gene_Symbol','Chromosome','Genomic_Coordinate')
  additional.columns <- l3headers[3:5] # as on the above line
  if(!all(additional.columns %in% fvarLabels(x))) { # {{{
    data(level3.symbols.hg19) # merged from 450k and 27k manifests
    missing.cols <- setdiff(additional.columns, fvarLabels(x))
    fData(x) <- cbind(fData(x), level3.symbols.hg19[featureNames(x), missing.cols])
  }
  
  write.level3 <- function(s) {		
    xs <- which(sampleNames(x) == s)
    message(paste("Writing level 3 data for sample", which(subjects == s),
                  "of", length(subjects)))
    lvl3data <- data.frame(Beta=betas(x)[,xs], 
                           Gene_Symbol=fData(x)[,'Gene_Symbol'],
                           Chromosome=fData(x)[,'Chromosome'],
                           Genomic_Coordinate=as.integer(fData(x)[,'Genomic_Coordinate']),
                           Pval=pvals(x)[,xs])
    rownames(lvl3data) <- featureNames(x)
    lvl3data[ which(lvl3data$Pval > 0.05), 'Beta' ] <- NA
    lvl3data$Pval <- NULL # no longer included in Level 3 output
    lvl3data$Chromosome[which(is.na(lvl3data$Chromosome))] <- NA
    lvl3data$Genomic_Coordinate[which(is.na(lvl3data$Genomic_Coordinate))] <- 0
    dump.file <- paste(paste(diseasestub,platform,batch.id,'lvl-3',s,'txt',sep='.'))
    headers1 <- paste('Hybridization REF', s, s, s, sep="\t")
    headers2 <- paste('Reporter ID', 'Beta_value', 
                      'Gene_Symbol', 'Chromosome', 'Genomic_Coordinate', sep="\t")
    cat(headers1, "\n", sep='', file=dump.file)
    cat(headers2, "\n", sep='', file=dump.file, append=TRUE)
    write.table(lvl3data, file=dump.file, append=TRUE, quote=FALSE,
                row.names=TRUE, col.names=FALSE, sep="\t")
    return('done')
  }
  
  if(parallel) {
    if(!require(parallel)) require(multicore)
    results <- mclapply(subjects, write.level3)
  } else { 
    results <- lapply(subjects, write.level3)
  }
  gc(,T)
}
