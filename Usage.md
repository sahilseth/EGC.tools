Testing USC Methylation pipeline
---
author: "Sahil Seth"
date: "7/9/2014"
output: html_document
---


### 1. Installing the required packages

```r
###### ---------- INSTALL dependencies of EGC.tools
source("http://bioconductor.org/biocLite.R")
biocLite('matrixStats')
biocLite("methylumi")
biocLite("FDb.InfiniumMethylation.hg19")
biocLite("AnnotationForge")
biocLite("IlluminaHumanMethylation450k.db")
biocLite("IlluminaHumanMethylation27k.db")
install.packages("RJSONIO")
## install.packages("RMySQL", type="source") ## removed from dependencies

require(devtools)
## does not work since RMySQL is missing
#install_github(repo = "EGC.tools", username = "uscepigenomecenter", branch = "recovery")

###### ---------- remove reference of MYSQL from NAMESPACE
install_github(repo = "EGC.tools", username = "sahilseth",  ref = "recovery")
```

### 2. Download some sample data from TCGA portal
```
cd ~/tmp/methylation_test

baseurl="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/hnsc/cgcc/jhu-usc.edu"
wget ${baseurl}/humanmethylation450/methylation/jhu-usc.edu_HNSC.HumanMethylation450.Level_1.9.8.0.tar.gz
wget ${baseurl}/humanmethylation450/methylation/jhu-usc.edu_HNSC.HumanMethylation450.aux.1.8.0.tar.gz

tar -zxvf jhu-usc.edu_HNSC.HumanMethylation450.Level_1.9.8.0.tar.gz
tar -zxvf jhu-usc.edu_HNSC.HumanMethylation450.aux.1.8.0.tar.gz
```

### 3. Setup input files and paramters

```r
## setup my WD
basepath = "~/tmp/methylation_test"
levelipath = file.path(basepath, "jhu-usc.edu_HNSC.HumanMethylation450.Level_1.9.8.0")
auxpath = file.path(basepath, "jhu-usc.edu_HNSC.HumanMethylation450.aux.1.8.0")
cohort = "HNSC"
verbose=FALSE
cores = 24
```


### 4. Read mapping file and subset it according to data files

```r
library(EGC.tools, quietly=!verbose, warn.conflicts=verbose)
library(IlluminaHumanMethylation450k.db, quietly=!verbose, warn.conflicts=verbose)
library(parallel)
setwd(basepath)
options('mc.cores' = cores)

###### --------- reading input data
mapping <- read.csv(file=sprintf("%s/%s.mappings.csv", basepath, cohort), stringsAsFactors=FALSE)
idats <-  gsub("(.*)_[Grn|Red]*.idat", "\\1", list.files(levelipath, pattern="idat"))

###### --------- subset rows where idats are available
mapping2 <- mapping[mapping$barcode %in% idats,]
```

### 5. Run the Methylation pipeline
#### 5.1 Reading the data

```r
###### --------- reading idats
TUMOR <- methylumIDAT(mapping2, parallel=TRUE, idatPath=levelipath)
```

```
## 0 HumanMethylation27 samples found
## 19 HumanMethylation450 samples found
```

#### 5.2 Normalization and other processing steps

```r
###### --------- reordering probes
data(probe.ordering)
if(!identical(featureNames(TUMOR), probe.ordering)){
  message("Reordering the probes")
  TUMOR <- TUMOR[match(probe.ordering, featureNames(TUMOR)), ]
  }

###### --------- Background correction
TUMOR <- methylumi.bgcorr(TUMOR)
```

```
## Background mean & SD estimated from 178406 probes
## Background mean & SD estimated from 92596 probes
```

```r
###### --------- Reduce Size of Dataset
TUMOR <- stripMethyLumiSet(TUMOR)

###### --------- Dye-Bias Equalization (Only for 450k data)
TUMOR <- normalizeMethyLumiSet(TUMOR)
```

```
## Normalizing via Illumina controls...
## Using sample number 11 as reference level...
```

#### 5.3 building an archive
```r
###### --------- Generate Level 2 and 3 files
buildArchive2(TUMOR, base = basepath, parallel = TRUE)
```

```
## Writing level 2 and level 3 data...
## Creating directory ~/tmp/methylation_test/Level_2 ...
## Creating directory ~/tmp/methylation_test/Level_3 ...
```

#### 5.4 bathwise folders: Expiremental
- This assumes a unique and strict directory structure in your home folder


```r
mapping <- read.csv(file=sprintf("%s/%s.mappings.csv", auxpath, cohort), stringsAsFactors=FALSE)
###### --------- setting up input paths
outpath = sprintf("%s/tcga/%s", basepath, cohort)
rawpath = sprintf("%s/raw/%s", basepath, cohort) ## raw path
## create the directory skeleton
system(sprintf("mkdir -p %s", dirname(rawpath)))
system(sprintf("rm -r %s; mkdir -p %s", outpath, outpath))
## place auxpath in rawpath
system(sprintf("ln -sf %s %s", auxpath, rawpath))
system(sprintf("ln -sf %s ~/meth450k", basepath))
###### --------- Generate Level 1, 2 and 3 files (along with sdrf)
try(buildArchive(TUMOR, base=basepath))
```

