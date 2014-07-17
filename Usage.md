---
output: html_document
---

```r
source("http://bioconductor.org/biocLite.R")

###### ---------- INSTALL dependencies
biocLite('matrixStats')
biocLite("methylumi")
biocLite("FDb.InfiniumMethylation.hg19")
biocLite("AnnotationForge")
biocLite("IlluminaHumanMethylation450k.db")
install.packages("RJSONIO")
## install.packages("RMySQL", type="source")

###### ---------- INSTALL packages not in 'pacakge description'
install.packages("IlluminaHumanMethylation450k.db")


require(devtools)
## does not work since RMySQL is missing
#install_github(repo = "EGC.tools", username = "uscepigenomecenter", branch = "recovery")

###### ---------- remove reference of MYSQL from NAMESPACE
install_local("~/projects/broad_gdac/EGC.tools")
```

### Download some sample data from TCGA portal
```{download}
cd ~/tmp/methylation_test

wget https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/hnsc/cgcc/jhu-usc.edu/humanmethylation450/methylation/jhu-usc.edu_HNSC.HumanMethylation450.Level_1.9.8.0.tar.gz
wget https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/hnsc/cgcc/jhu-usc.edu/humanmethylation450/methylation/jhu-usc.edu_HNSC.HumanMethylation450.aux.1.8.0.tar.gz

tar -zxvf jhu-usc.edu_HNSC.HumanMethylation450.Level_1.9.8.0.tar.gz
tar -zxvf jhu-usc.edu_HNSC.HumanMethylation450.Level_3.9.8.0.tar.gz

```

#### setup input files

```r
## setup my WD
basepath = "~/tmp/methylation_test"
levelipath = file.path(basepath, "jhu-usc.edu_HNSC.HumanMethylation450.Level_1.9.8.0")
auxpath = file.path(basepath, "jhu-usc.edu_HNSC.HumanMethylation450.aux.1.8.0")
cohort = "HNSC"
```


### 

```r
library(EGC.tools)
library(parallel)
cores=24
setwd(basepath)
options('mc.cores'=cores)

###### --------- setting up input paths
outpath = sprintf("tcga/%s", cohort)
rawpath = sprintf("raw/%s", cohort) ## raw path
## create the directory skeleton
tmp <- sapply(c(outpath, dirname(rawpath)), dir.create, recursive = TRUE, showWarnings = FALSE)
## place auxpath in rawpath
system(sprintf("ln -sf %s %s", auxpath, rawpath))

###### --------- reading input data
mapping <- read.csv(file=sprintf("%s/%s.mappings.csv",auxpath, cohort), stringsAsFactors=FALSE)
idats <-  gsub("(.*)_[Grn|Red]*.idat", "\\1", list.files(levelipath, pattern="idat"))

###### --------- subset rows where idats are available
mapping2 <- mapping[mapping$barcode %in% idats,]

###### --------- reading idats
TUMOR <- methylumIDAT(mapping2, parallel=TRUE, idatPath=levelipath)
```

```
## 0 HumanMethylation27 samples found
## 19 HumanMethylation450 samples found
```

```r
###### --------- reordering probes
data(probe.ordering)
if(!identical(featureNames(TUMOR), probe.ordering)){
  message("Reordering the probes")
  TUMOR <- TUMOR[match(probe.ordering, featureNames(TUMOR)), ]
  }
```

```
## Reordering the probes
```

```r
###### --------- Background correction
TUMOR <- methylumi.bgcorrect(TUMOR)
```

```
## Error: could not find function "methylumi.bgcorrect"
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

```r
###### --------- Generate Level 1, 2 and 3 files
undebug(buildArchive)
```

```
## Warning: argument is not being debugged
```

```r
undebug(makeArchiveDirs)
```

```
## Warning: argument is not being debugged
```

```r
buildArchive(mapping2, base=basepath)
```

```
## Creating archive directories...
```

```
## Error: object 'x' not found
```

