---
title: "R Notebook"
output: github_document
---
#Méthode appliquée à partir du tutoriel : https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
##Amplicon bioinformatics: from raw reads to tables


```r
.cran_packages <- c("ggplot2", "gridExtra", "devtools")
install.packages(.cran_packages) 
```

```
## Installing packages into '/usr/local/lib/R/site-library'
## (as 'lib' is unspecified)
```

```r
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
BiocManager::install(.bioc_packages)
```

```
## 'getOption("repos")' replaces Bioconductor standard repositories, see
## '?repositories' for details
## 
## replacement repositories:
##     CRAN: https://packagemanager.rstudio.com/cran/__linux__/focal/latest
```

```
## Bioconductor version 3.14 (BiocManager 1.30.16), R 4.1.2 (2021-11-01)
```

```
## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
##   re-install: 'dada2' 'phyloseq' 'DECIPHER' 'phangorn'
```

```
## Installation paths not writeable, unable to update packages
##   path: /usr/local/lib/R/library
##   packages:
##     Matrix
```

```r
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
```

```
## Loading required package: ggplot2
```

```
## Loading required package: gridExtra
```

```
## Loading required package: devtools
```

```
## Loading required package: usethis
```

```
## Loading required package: dada2
```

```
## Loading required package: Rcpp
```

```
## Loading required package: phyloseq
```

```
## Loading required package: DECIPHER
```

```
## Loading required package: Biostrings
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following object is masked from 'package:gridExtra':
## 
##     combine
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:phyloseq':
## 
##     distance
```

```
## Loading required package: XVector
```

```
## Loading required package: GenomeInfoDb
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## Loading required package: RSQLite
```

```
## Loading required package: parallel
```

```
## Loading required package: phangorn
```

```
## Loading required package: ape
```

```
## 
## Attaching package: 'ape'
```

```
## The following object is masked from 'package:Biostrings':
## 
##     complement
```

```
##   ggplot2 gridExtra  devtools     dada2  phyloseq  DECIPHER  phangorn 
##      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE
```


```r
set.seed(100)
miseq_path <- "/home/rstudio/MiSeq_SOP"
list.files(miseq_path)
```

```
##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
## [43] "mouse.dpw.metadata"            "mouse.time.design"            
## [45] "stability.batch"               "stability.files"
```

##Filter and Trim


```r
fnFs <- sort(list.files(miseq_path, pattern = "_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern = "_R2_001.fastq"))
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
#file.path permet d'adapter une fonction à tous les languages de programmation (unix/python/R etc.)
```



```r
fnFs[1:3]
```

```
## [1] "/home/rstudio/MiSeq_SOP/F3D0_S188_L001_R1_001.fastq"  
## [2] "/home/rstudio/MiSeq_SOP/F3D1_S189_L001_R1_001.fastq"  
## [3] "/home/rstudio/MiSeq_SOP/F3D141_S207_L001_R1_001.fastq"
```

```r
## [1] "./MiSeq_SOP/F3D0_S188_L001_R1_001.fastq"   "./MiSeq_SOP/F3D1_S189_L001_R1_001.fastq"  
## [3] "./MiSeq_SOP/F3D141_S207_L001_R1_001.fastq"
fnRs[1:3]
```

```
## [1] "/home/rstudio/MiSeq_SOP/F3D0_S188_L001_R2_001.fastq"  
## [2] "/home/rstudio/MiSeq_SOP/F3D1_S189_L001_R2_001.fastq"  
## [3] "/home/rstudio/MiSeq_SOP/F3D141_S207_L001_R2_001.fastq"
```

```r
## [1] "./MiSeq_SOP/F3D0_S188_L001_R2_001.fastq"   "./MiSeq_SOP/F3D1_S189_L001_R2_001.fastq"  
## [3] "./MiSeq_SOP/F3D141_S207_L001_R2_001.fastq"
plotQualityProfile(fnFs[1:2])
```

```
## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
## "none")` instead.
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->



```r
plotQualityProfile(fnRs[1:2])
```

```
## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
## "none")` instead.
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->



```r
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
```

```r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
```

```
##                               reads.in reads.out
## F3D0_S188_L001_R1_001.fastq       7793      7113
## F3D1_S189_L001_R1_001.fastq       5869      5299
## F3D141_S207_L001_R1_001.fastq     5958      5463
## F3D142_S208_L001_R1_001.fastq     3183      2914
## F3D143_S209_L001_R1_001.fastq     3178      2941
## F3D144_S210_L001_R1_001.fastq     4827      4312
```
##Infer sequence variants
###Dereplication

```r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D0_F_filt.fastq.gz
```

```
## Encountered 1979 unique sequences from 7113 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D1_F_filt.fastq.gz
```

```
## Encountered 1639 unique sequences from 5299 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D141_F_filt.fastq.gz
```

```
## Encountered 1477 unique sequences from 5463 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D142_F_filt.fastq.gz
```

```
## Encountered 904 unique sequences from 2914 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D143_F_filt.fastq.gz
```

```
## Encountered 939 unique sequences from 2941 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D144_F_filt.fastq.gz
```

```
## Encountered 1267 unique sequences from 4312 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D145_F_filt.fastq.gz
```

```
## Encountered 1756 unique sequences from 6741 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D146_F_filt.fastq.gz
```

```
## Encountered 1438 unique sequences from 4560 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D147_F_filt.fastq.gz
```

```
## Encountered 3590 unique sequences from 15637 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D148_F_filt.fastq.gz
```

```
## Encountered 2762 unique sequences from 11413 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D149_F_filt.fastq.gz
```

```
## Encountered 3021 unique sequences from 12017 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D150_F_filt.fastq.gz
```

```
## Encountered 1566 unique sequences from 5032 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D2_F_filt.fastq.gz
```

```
## Encountered 3707 unique sequences from 18075 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D3_F_filt.fastq.gz
```

```
## Encountered 1479 unique sequences from 6250 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D5_F_filt.fastq.gz
```

```
## Encountered 1195 unique sequences from 4052 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D6_F_filt.fastq.gz
```

```
## Encountered 1832 unique sequences from 7369 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D7_F_filt.fastq.gz
```

```
## Encountered 1183 unique sequences from 4765 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D8_F_filt.fastq.gz
```

```
## Encountered 1382 unique sequences from 4871 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D9_F_filt.fastq.gz
```

```
## Encountered 1709 unique sequences from 6504 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/Mock_F_filt.fastq.gz
```

```
## Encountered 897 unique sequences from 4314 total sequences read.
```

```r
derepRs <- derepFastq(filtRs, verbose=TRUE)
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D0_R_filt.fastq.gz
```

```
## Encountered 1660 unique sequences from 7113 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D1_R_filt.fastq.gz
```

```
## Encountered 1349 unique sequences from 5299 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D141_R_filt.fastq.gz
```

```
## Encountered 1335 unique sequences from 5463 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D142_R_filt.fastq.gz
```

```
## Encountered 853 unique sequences from 2914 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D143_R_filt.fastq.gz
```

```
## Encountered 880 unique sequences from 2941 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D144_R_filt.fastq.gz
```

```
## Encountered 1286 unique sequences from 4312 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D145_R_filt.fastq.gz
```

```
## Encountered 1803 unique sequences from 6741 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D146_R_filt.fastq.gz
```

```
## Encountered 1265 unique sequences from 4560 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D147_R_filt.fastq.gz
```

```
## Encountered 3414 unique sequences from 15637 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D148_R_filt.fastq.gz
```

```
## Encountered 2522 unique sequences from 11413 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D149_R_filt.fastq.gz
```

```
## Encountered 2771 unique sequences from 12017 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D150_R_filt.fastq.gz
```

```
## Encountered 1415 unique sequences from 5032 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D2_R_filt.fastq.gz
```

```
## Encountered 3290 unique sequences from 18075 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D3_R_filt.fastq.gz
```

```
## Encountered 1390 unique sequences from 6250 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D5_R_filt.fastq.gz
```

```
## Encountered 1134 unique sequences from 4052 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D6_R_filt.fastq.gz
```

```
## Encountered 1635 unique sequences from 7369 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D7_R_filt.fastq.gz
```

```
## Encountered 1084 unique sequences from 4765 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D8_R_filt.fastq.gz
```

```
## Encountered 1161 unique sequences from 4871 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D9_R_filt.fastq.gz
```

```
## Encountered 1502 unique sequences from 6504 total sequences read.
```

```
## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/Mock_R_filt.fastq.gz
```

```
## Encountered 732 unique sequences from 4314 total sequences read.
```

```r
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames
```

```r
errF <- learnErrors(filtFs, multithread=TRUE)
```

```
## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.
```

```r
errR <- learnErrors(filtRs, multithread=TRUE)
```

```
## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.
```

```r
plotErrors(errF)
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

```r
plotErrors(errR)
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

```r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

```
## Sample 1 - 7113 reads in 1979 unique sequences.
## Sample 2 - 5299 reads in 1639 unique sequences.
## Sample 3 - 5463 reads in 1477 unique sequences.
## Sample 4 - 2914 reads in 904 unique sequences.
## Sample 5 - 2941 reads in 939 unique sequences.
## Sample 6 - 4312 reads in 1267 unique sequences.
## Sample 7 - 6741 reads in 1756 unique sequences.
## Sample 8 - 4560 reads in 1438 unique sequences.
## Sample 9 - 15637 reads in 3590 unique sequences.
## Sample 10 - 11413 reads in 2762 unique sequences.
## Sample 11 - 12017 reads in 3021 unique sequences.
## Sample 12 - 5032 reads in 1566 unique sequences.
## Sample 13 - 18075 reads in 3707 unique sequences.
## Sample 14 - 6250 reads in 1479 unique sequences.
## Sample 15 - 4052 reads in 1195 unique sequences.
## Sample 16 - 7369 reads in 1832 unique sequences.
## Sample 17 - 4765 reads in 1183 unique sequences.
## Sample 18 - 4871 reads in 1382 unique sequences.
## Sample 19 - 6504 reads in 1709 unique sequences.
## Sample 20 - 4314 reads in 897 unique sequences.
```

```r
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

```
## Sample 1 - 7113 reads in 1660 unique sequences.
## Sample 2 - 5299 reads in 1349 unique sequences.
## Sample 3 - 5463 reads in 1335 unique sequences.
## Sample 4 - 2914 reads in 853 unique sequences.
## Sample 5 - 2941 reads in 880 unique sequences.
## Sample 6 - 4312 reads in 1286 unique sequences.
## Sample 7 - 6741 reads in 1803 unique sequences.
## Sample 8 - 4560 reads in 1265 unique sequences.
## Sample 9 - 15637 reads in 3414 unique sequences.
## Sample 10 - 11413 reads in 2522 unique sequences.
## Sample 11 - 12017 reads in 2771 unique sequences.
## Sample 12 - 5032 reads in 1415 unique sequences.
## Sample 13 - 18075 reads in 3290 unique sequences.
## Sample 14 - 6250 reads in 1390 unique sequences.
## Sample 15 - 4052 reads in 1134 unique sequences.
## Sample 16 - 7369 reads in 1635 unique sequences.
## Sample 17 - 4765 reads in 1084 unique sequences.
## Sample 18 - 4871 reads in 1161 unique sequences.
## Sample 19 - 6504 reads in 1502 unique sequences.
## Sample 20 - 4314 reads in 732 unique sequences.
```

```r
dadaFs[[1]]
```

```
## dada-class: object describing DADA2 denoising results
## 128 sequence variants were inferred from 1979 input unique sequences.
## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
```
##Construct sequence table and remove chimeras

```r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
```

```r
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))
```

```
## 
## 251 252 253 254 255 
##   1  85 186   5   2
```

```bash
cd ~
wget  https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
```

```
## --2021-12-24 17:40:06--  https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
## Resolving zenodo.org (zenodo.org)... 137.138.76.77
## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
## HTTP request sent, awaiting response... 200 OK
## Length: 137283333 (131M) [application/octet-stream]
## Saving to: ‘silva_nr99_v138.1_train_set.fa.gz.3’
## 
##      0K .......... .......... .......... .......... ..........  0% 12.5M 10s
##     50K .......... .......... .......... .......... ..........  0% 11.7M 11s
##    100K .......... .......... .......... .......... ..........  0% 13.2M 11s
##    150K .......... .......... .......... .......... ..........  0% 11.7M 11s
##    200K .......... .......... .......... .......... ..........  0%  106M 9s
##    250K .......... .......... .......... .......... ..........  0% 97.7M 8s
##    300K .......... .......... .......... .......... ..........  0% 15.9M 8s
##    350K .......... .......... .......... .......... ..........  0% 67.4M 7s
##    400K .......... .......... .......... .......... ..........  0% 63.6M 6s
##    450K .......... .......... .......... .......... ..........  0% 17.6M 6s
##    500K .......... .......... .......... .......... ..........  0% 97.0M 6s
##    550K .......... .......... .......... .......... ..........  0% 35.9M 6s
##    600K .......... .......... .......... .......... ..........  0% 33.6M 6s
##    650K .......... .......... .......... .......... ..........  0% 28.9M 6s
##    700K .......... .......... .......... .......... ..........  0% 93.6M 5s
##    750K .......... .......... .......... .......... ..........  0% 31.2M 5s
##    800K .......... .......... .......... .......... ..........  0% 34.3M 5s
##    850K .......... .......... .......... .......... ..........  0% 24.5M 5s
##    900K .......... .......... .......... .......... ..........  0% 83.2M 5s
##    950K .......... .......... .......... .......... ..........  0% 74.9M 5s
##   1000K .......... .......... .......... .......... ..........  0% 20.1M 5s
##   1050K .......... .......... .......... .......... ..........  0% 44.5M 5s
##   1100K .......... .......... .......... .......... ..........  0% 40.8M 5s
##   1150K .......... .......... .......... .......... ..........  0% 21.9M 5s
##   1200K .......... .......... .......... .......... ..........  0% 95.7M 5s
##   1250K .......... .......... .......... .......... ..........  0% 49.9M 5s
##   1300K .......... .......... .......... .......... ..........  1% 19.3M 5s
##   1350K .......... .......... .......... .......... ..........  1% 77.2M 5s
##   1400K .......... .......... .......... .......... ..........  1% 96.8M 4s
##   1450K .......... .......... .......... .......... ..........  1% 13.6M 5s
##   1500K .......... .......... .......... .......... ..........  1%  101M 4s
##   1550K .......... .......... .......... .......... ..........  1% 15.1M 5s
##   1600K .......... .......... .......... .......... ..........  1% 85.2M 4s
##   1650K .......... .......... .......... .......... ..........  1%  109M 4s
##   1700K .......... .......... .......... .......... ..........  1% 14.5M 5s
##   1750K .......... .......... .......... .......... ..........  1% 76.4M 4s
##   1800K .......... .......... .......... .......... ..........  1% 14.5M 5s
##   1850K .......... .......... .......... .......... ..........  1% 97.9M 4s
##   1900K .......... .......... .......... .......... ..........  1% 98.1M 4s
##   1950K .......... .......... .......... .......... ..........  1% 14.7M 4s
##   2000K .......... .......... .......... .......... ..........  1%  105M 4s
##   2050K .......... .......... .......... .......... ..........  1% 14.7M 5s
##   2100K .......... .......... .......... .......... ..........  1%  102M 4s
##   2150K .......... .......... .......... .......... ..........  1% 91.4M 4s
##   2200K .......... .......... .......... .......... ..........  1% 15.4M 4s
##   2250K .......... .......... .......... .......... ..........  1%  113M 4s
##   2300K .......... .......... .......... .......... ..........  1% 48.9M 4s
##   2350K .......... .......... .......... .......... ..........  1% 27.1M 4s
##   2400K .......... .......... .......... .......... ..........  1%  111M 4s
##   2450K .......... .......... .......... .......... ..........  1% 29.7M 4s
##   2500K .......... .......... .......... .......... ..........  1% 35.8M 4s
##   2550K .......... .......... .......... .......... ..........  1% 25.8M 4s
##   2600K .......... .......... .......... .......... ..........  1% 81.0M 4s
##   2650K .......... .......... .......... .......... ..........  2% 35.2M 4s
##   2700K .......... .......... .......... .......... ..........  2% 35.1M 4s
##   2750K .......... .......... .......... .......... ..........  2% 23.0M 4s
##   2800K .......... .......... .......... .......... ..........  2% 47.6M 4s
##   2850K .......... .......... .......... .......... ..........  2% 76.5M 4s
##   2900K .......... .......... .......... .......... ..........  2% 23.7M 4s
##   2950K .......... .......... .......... .......... ..........  2% 30.8M 4s
##   3000K .......... .......... .......... .......... ..........  2% 27.4M 4s
##   3050K .......... .......... .......... .......... ..........  2% 60.3M 4s
##   3100K .......... .......... .......... .......... ..........  2% 34.7M 4s
##   3150K .......... .......... .......... .......... ..........  2% 26.6M 4s
##   3200K .......... .......... .......... .......... ..........  2% 29.4M 4s
##   3250K .......... .......... .......... .......... ..........  2% 27.8M 4s
##   3300K .......... .......... .......... .......... ..........  2% 65.9M 4s
##   3350K .......... .......... .......... .......... ..........  2% 35.3M 4s
##   3400K .......... .......... .......... .......... ..........  2% 16.3M 4s
##   3450K .......... .......... .......... .......... ..........  2% 67.9M 4s
##   3500K .......... .......... .......... .......... ..........  2% 89.6M 4s
##   3550K .......... .......... .......... .......... ..........  2% 66.2M 4s
##   3600K .......... .......... .......... .......... ..........  2% 28.4M 4s
##   3650K .......... .......... .......... .......... ..........  2% 64.2M 4s
##   3700K .......... .......... .......... .......... ..........  2% 72.1M 4s
##   3750K .......... .......... .......... .......... ..........  2% 80.6M 4s
##   3800K .......... .......... .......... .......... ..........  2% 12.6M 4s
##   3850K .......... .......... .......... .......... ..........  2% 90.6M 4s
##   3900K .......... .......... .......... .......... ..........  2% 91.9M 4s
##   3950K .......... .......... .......... .......... ..........  2% 19.4M 4s
##   4000K .......... .......... .......... .......... ..........  3% 68.0M 4s
##   4050K .......... .......... .......... .......... ..........  3% 61.7M 4s
##   4100K .......... .......... .......... .......... ..........  3% 84.7M 4s
##   4150K .......... .......... .......... .......... ..........  3% 27.5M 4s
##   4200K .......... .......... .......... .......... ..........  3% 79.1M 4s
##   4250K .......... .......... .......... .......... ..........  3% 67.6M 4s
##   4300K .......... .......... .......... .......... ..........  3% 72.2M 4s
##   4350K .......... .......... .......... .......... ..........  3% 27.6M 4s
##   4400K .......... .......... .......... .......... ..........  3% 30.3M 4s
##   4450K .......... .......... .......... .......... ..........  3% 71.6M 4s
##   4500K .......... .......... .......... .......... ..........  3% 90.7M 4s
##   4550K .......... .......... .......... .......... ..........  3% 40.7M 4s
##   4600K .......... .......... .......... .......... ..........  3% 26.5M 4s
##   4650K .......... .......... .......... .......... ..........  3% 82.3M 4s
##   4700K .......... .......... .......... .......... ..........  3% 81.0M 4s
##   4750K .......... .......... .......... .......... ..........  3% 56.0M 4s
##   4800K .......... .......... .......... .......... ..........  3% 21.5M 4s
##   4850K .......... .......... .......... .......... ..........  3% 77.2M 4s
##   4900K .......... .......... .......... .......... ..........  3% 84.6M 4s
##   4950K .......... .......... .......... .......... ..........  3% 66.7M 4s
##   5000K .......... .......... .......... .......... ..........  3% 21.5M 4s
##   5050K .......... .......... .......... .......... ..........  3% 70.0M 4s
##   5100K .......... .......... .......... .......... ..........  3% 91.9M 4s
##   5150K .......... .......... .......... .......... ..........  3% 32.6M 4s
##   5200K .......... .......... .......... .......... ..........  3% 70.7M 4s
##   5250K .......... .......... .......... .......... ..........  3% 96.9M 4s
##   5300K .......... .......... .......... .......... ..........  3% 71.8M 4s
##   5350K .......... .......... .......... .......... ..........  4% 34.0M 4s
##   5400K .......... .......... .......... .......... ..........  4% 80.2M 4s
##   5450K .......... .......... .......... .......... ..........  4% 76.7M 4s
##   5500K .......... .......... .......... .......... ..........  4% 61.7M 4s
##   5550K .......... .......... .......... .......... ..........  4% 33.7M 4s
##   5600K .......... .......... .......... .......... ..........  4% 80.2M 4s
##   5650K .......... .......... .......... .......... ..........  4% 92.2M 4s
##   5700K .......... .......... .......... .......... ..........  4% 44.7M 4s
##   5750K .......... .......... .......... .......... ..........  4% 23.6M 4s
##   5800K .......... .......... .......... .......... ..........  4% 93.6M 3s
##   5850K .......... .......... .......... .......... ..........  4%  104M 3s
##   5900K .......... .......... .......... .......... ..........  4%  104M 3s
##   5950K .......... .......... .......... .......... ..........  4% 19.7M 3s
##   6000K .......... .......... .......... .......... ..........  4% 79.4M 3s
##   6050K .......... .......... .......... .......... ..........  4% 97.1M 3s
##   6100K .......... .......... .......... .......... ..........  4% 30.1M 3s
##   6150K .......... .......... .......... .......... ..........  4% 70.7M 3s
##   6200K .......... .......... .......... .......... ..........  4% 98.9M 3s
##   6250K .......... .......... .......... .......... ..........  4% 72.5M 3s
##   6300K .......... .......... .......... .......... ..........  4% 29.7M 3s
##   6350K .......... .......... .......... .......... ..........  4% 45.3M 3s
##   6400K .......... .......... .......... .......... ..........  4% 66.8M 3s
##   6450K .......... .......... .......... .......... ..........  4%  106M 3s
##   6500K .......... .......... .......... .......... ..........  4% 50.3M 3s
##   6550K .......... .......... .......... .......... ..........  4% 21.2M 3s
##   6600K .......... .......... .......... .......... ..........  4% 88.1M 3s
##   6650K .......... .......... .......... .......... ..........  4% 83.0M 3s
##   6700K .......... .......... .......... .......... ..........  5% 88.7M 3s
##   6750K .......... .......... .......... .......... ..........  5% 18.3M 3s
##   6800K .......... .......... .......... .......... ..........  5% 78.7M 3s
##   6850K .......... .......... .......... .......... ..........  5% 73.9M 3s
##   6900K .......... .......... .......... .......... ..........  5% 91.3M 3s
##   6950K .......... .......... .......... .......... ..........  5% 21.7M 3s
##   7000K .......... .......... .......... .......... ..........  5% 93.0M 3s
##   7050K .......... .......... .......... .......... ..........  5% 89.1M 3s
##   7100K .......... .......... .......... .......... ..........  5% 93.2M 3s
##   7150K .......... .......... .......... .......... ..........  5% 21.7M 3s
##   7200K .......... .......... .......... .......... ..........  5% 82.8M 3s
##   7250K .......... .......... .......... .......... ..........  5%  106M 3s
##   7300K .......... .......... .......... .......... ..........  5% 96.9M 3s
##   7350K .......... .......... .......... .......... ..........  5% 58.0M 3s
##   7400K .......... .......... .......... .......... ..........  5% 29.1M 3s
##   7450K .......... .......... .......... .......... ..........  5% 90.2M 3s
##   7500K .......... .......... .......... .......... ..........  5% 92.1M 3s
##   7550K .......... .......... .......... .......... ..........  5% 50.5M 3s
##   7600K .......... .......... .......... .......... ..........  5% 95.4M 3s
##   7650K .......... .......... .......... .......... ..........  5% 24.1M 3s
##   7700K .......... .......... .......... .......... ..........  5% 96.2M 3s
##   7750K .......... .......... .......... .......... ..........  5% 92.5M 3s
##   7800K .......... .......... .......... .......... ..........  5%  109M 3s
##   7850K .......... .......... .......... .......... ..........  5% 85.6M 3s
##   7900K .......... .......... .......... .......... ..........  5% 27.6M 3s
##   7950K .......... .......... .......... .......... ..........  5% 70.0M 3s
##   8000K .......... .......... .......... .......... ..........  6% 91.6M 3s
##   8050K .......... .......... .......... .......... ..........  6%  102M 3s
##   8100K .......... .......... .......... .......... ..........  6% 81.4M 3s
##   8150K .......... .......... .......... .......... ..........  6% 25.7M 3s
##   8200K .......... .......... .......... .......... ..........  6% 90.7M 3s
##   8250K .......... .......... .......... .......... ..........  6%  114M 3s
##   8300K .......... .......... .......... .......... ..........  6% 93.5M 3s
##   8350K .......... .......... .......... .......... ..........  6% 79.3M 3s
##   8400K .......... .......... .......... .......... ..........  6% 25.3M 3s
##   8450K .......... .......... .......... .......... ..........  6% 64.3M 3s
##   8500K .......... .......... .......... .......... ..........  6% 79.3M 3s
##   8550K .......... .......... .......... .......... ..........  6% 64.5M 3s
##   8600K .......... .......... .......... .......... ..........  6% 89.4M 3s
##   8650K .......... .......... .......... .......... ..........  6% 37.4M 3s
##   8700K .......... .......... .......... .......... ..........  6% 66.4M 3s
##   8750K .......... .......... .......... .......... ..........  6% 72.3M 3s
##   8800K .......... .......... .......... .......... ..........  6% 35.7M 3s
##   8850K .......... .......... .......... .......... ..........  6% 87.8M 3s
##   8900K .......... .......... .......... .......... ..........  6% 69.0M 3s
##   8950K .......... .......... .......... .......... ..........  6% 70.8M 3s
##   9000K .......... .......... .......... .......... ..........  6% 84.5M 3s
##   9050K .......... .......... .......... .......... ..........  6% 94.8M 3s
##   9100K .......... .......... .......... .......... ..........  6% 62.7M 3s
##   9150K .......... .......... .......... .......... ..........  6% 72.7M 3s
##   9200K .......... .......... .......... .......... ..........  6% 66.0M 3s
##   9250K .......... .......... .......... .......... ..........  6% 55.0M 3s
##   9300K .......... .......... .......... .......... ..........  6% 59.1M 3s
##   9350K .......... .......... .......... .......... ..........  7% 63.5M 3s
##   9400K .......... .......... .......... .......... ..........  7% 91.3M 3s
##   9450K .......... .......... .......... .......... ..........  7% 89.5M 3s
##   9500K .......... .......... .......... .......... ..........  7% 47.4M 3s
##   9550K .......... .......... .......... .......... ..........  7% 50.3M 3s
##   9600K .......... .......... .......... .......... ..........  7% 85.3M 3s
##   9650K .......... .......... .......... .......... ..........  7% 91.7M 3s
##   9700K .......... .......... .......... .......... ..........  7% 91.2M 3s
##   9750K .......... .......... .......... .......... ..........  7% 84.2M 3s
##   9800K .......... .......... .......... .......... ..........  7% 57.1M 3s
##   9850K .......... .......... .......... .......... ..........  7% 86.6M 3s
##   9900K .......... .......... .......... .......... ..........  7% 88.0M 3s
##   9950K .......... .......... .......... .......... ..........  7% 78.3M 3s
##  10000K .......... .......... .......... .......... ..........  7% 65.6M 3s
##  10050K .......... .......... .......... .......... ..........  7% 37.4M 3s
##  10100K .......... .......... .......... .......... ..........  7% 86.6M 3s
##  10150K .......... .......... .......... .......... ..........  7% 73.6M 3s
##  10200K .......... .......... .......... .......... ..........  7% 79.4M 3s
##  10250K .......... .......... .......... .......... ..........  7% 93.7M 3s
##  10300K .......... .......... .......... .......... ..........  7% 26.0M 3s
##  10350K .......... .......... .......... .......... ..........  7% 64.1M 3s
##  10400K .......... .......... .......... .......... ..........  7% 89.4M 3s
##  10450K .......... .......... .......... .......... ..........  7% 85.3M 3s
##  10500K .......... .......... .......... .......... ..........  7%  102M 3s
##  10550K .......... .......... .......... .......... ..........  7% 90.8M 3s
##  10600K .......... .......... .......... .......... ..........  7% 65.5M 3s
##  10650K .......... .......... .......... .......... ..........  7% 71.9M 3s
##  10700K .......... .......... .......... .......... ..........  8% 85.8M 3s
##  10750K .......... .......... .......... .......... ..........  8% 78.3M 3s
##  10800K .......... .......... .......... .......... ..........  8% 65.1M 3s
##  10850K .......... .......... .......... .......... ..........  8% 88.1M 3s
##  10900K .......... .......... .......... .......... ..........  8% 50.1M 3s
##  10950K .......... .......... .......... .......... ..........  8% 65.5M 3s
##  11000K .......... .......... .......... .......... ..........  8% 78.5M 3s
##  11050K .......... .......... .......... .......... ..........  8% 37.6M 3s
##  11100K .......... .......... .......... .......... ..........  8% 83.0M 3s
##  11150K .......... .......... .......... .......... ..........  8% 54.0M 3s
##  11200K .......... .......... .......... .......... ..........  8%  104M 3s
##  11250K .......... .......... .......... .......... ..........  8% 84.5M 3s
##  11300K .......... .......... .......... .......... ..........  8%  105M 3s
##  11350K .......... .......... .......... .......... ..........  8% 22.1M 3s
##  11400K .......... .......... .......... .......... ..........  8% 72.1M 3s
##  11450K .......... .......... .......... .......... ..........  8% 81.3M 3s
##  11500K .......... .......... .......... .......... ..........  8% 61.9M 3s
##  11550K .......... .......... .......... .......... ..........  8% 85.7M 3s
##  11600K .......... .......... .......... .......... ..........  8% 23.1M 3s
##  11650K .......... .......... .......... .......... ..........  8% 78.3M 3s
##  11700K .......... .......... .......... .......... ..........  8% 72.8M 3s
##  11750K .......... .......... .......... .......... ..........  8% 73.6M 3s
##  11800K .......... .......... .......... .......... ..........  8% 90.1M 3s
##  11850K .......... .......... .......... .......... ..........  8% 89.2M 3s
##  11900K .......... .......... .......... .......... ..........  8% 70.8M 3s
##  11950K .......... .......... .......... .......... ..........  8% 68.5M 3s
##  12000K .......... .......... .......... .......... ..........  8% 84.1M 3s
##  12050K .......... .......... .......... .......... ..........  9% 64.6M 3s
##  12100K .......... .......... .......... .......... ..........  9% 53.8M 3s
##  12150K .......... .......... .......... .......... ..........  9% 54.5M 3s
##  12200K .......... .......... .......... .......... ..........  9% 86.0M 3s
##  12250K .......... .......... .......... .......... ..........  9% 17.8M 3s
##  12300K .......... .......... .......... .......... ..........  9%  102M 3s
##  12350K .......... .......... .......... .......... ..........  9% 77.6M 3s
##  12400K .......... .......... .......... .......... ..........  9%  106M 3s
##  12450K .......... .......... .......... .......... ..........  9% 12.9M 3s
##  12500K .......... .......... .......... .......... ..........  9% 54.7M 3s
##  12550K .......... .......... .......... .......... ..........  9% 85.9M 3s
##  12600K .......... .......... .......... .......... ..........  9%  107M 3s
##  12650K .......... .......... .......... .......... ..........  9%  104M 3s
##  12700K .......... .......... .......... .......... ..........  9% 87.9M 3s
##  12750K .......... .......... .......... .......... ..........  9% 12.0M 3s
##  12800K .......... .......... .......... .......... ..........  9%  103M 3s
##  12850K .......... .......... .......... .......... ..........  9%  105M 3s
##  12900K .......... .......... .......... .......... ..........  9%  119M 3s
##  12950K .......... .......... .......... .......... ..........  9% 89.0M 3s
##  13000K .......... .......... .......... .......... ..........  9% 23.5M 3s
##  13050K .......... .......... .......... .......... ..........  9% 99.7M 3s
##  13100K .......... .......... .......... .......... ..........  9% 86.2M 3s
##  13150K .......... .......... .......... .......... ..........  9% 63.3M 3s
##  13200K .......... .......... .......... .......... ..........  9%  110M 3s
##  13250K .......... .......... .......... .......... ..........  9%  102M 3s
##  13300K .......... .......... .......... .......... ..........  9% 68.7M 3s
##  13350K .......... .......... .......... .......... ..........  9% 34.4M 3s
##  13400K .......... .......... .......... .......... .......... 10% 57.9M 3s
##  13450K .......... .......... .......... .......... .......... 10% 79.1M 3s
##  13500K .......... .......... .......... .......... .......... 10% 85.8M 3s
##  13550K .......... .......... .......... .......... .......... 10% 86.1M 3s
##  13600K .......... .......... .......... .......... .......... 10% 57.7M 3s
##  13650K .......... .......... .......... .......... .......... 10% 52.6M 3s
##  13700K .......... .......... .......... .......... .......... 10% 63.8M 3s
##  13750K .......... .......... .......... .......... .......... 10% 48.2M 3s
##  13800K .......... .......... .......... .......... .......... 10% 41.4M 3s
##  13850K .......... .......... .......... .......... .......... 10% 62.6M 3s
##  13900K .......... .......... .......... .......... .......... 10% 22.0M 3s
##  13950K .......... .......... .......... .......... .......... 10% 53.3M 3s
##  14000K .......... .......... .......... .......... .......... 10% 84.6M 3s
##  14050K .......... .......... .......... .......... .......... 10% 81.5M 3s
##  14100K .......... .......... .......... .......... .......... 10% 55.6M 3s
##  14150K .......... .......... .......... .......... .......... 10% 64.9M 3s
##  14200K .......... .......... .......... .......... .......... 10% 55.7M 3s
##  14250K .......... .......... .......... .......... .......... 10% 69.5M 3s
##  14300K .......... .......... .......... .......... .......... 10% 79.1M 3s
##  14350K .......... .......... .......... .......... .......... 10% 57.7M 3s
##  14400K .......... .......... .......... .......... .......... 10% 90.2M 3s
##  14450K .......... .......... .......... .......... .......... 10% 63.3M 3s
##  14500K .......... .......... .......... .......... .......... 10% 20.0M 3s
##  14550K .......... .......... .......... .......... .......... 10% 61.2M 3s
##  14600K .......... .......... .......... .......... .......... 10% 87.5M 3s
##  14650K .......... .......... .......... .......... .......... 10% 94.1M 3s
##  14700K .......... .......... .......... .......... .......... 11% 75.3M 3s
##  14750K .......... .......... .......... .......... .......... 11% 76.3M 3s
##  14800K .......... .......... .......... .......... .......... 11% 52.8M 3s
##  14850K .......... .......... .......... .......... .......... 11% 63.2M 3s
##  14900K .......... .......... .......... .......... .......... 11% 90.7M 3s
##  14950K .......... .......... .......... .......... .......... 11% 78.9M 3s
##  15000K .......... .......... .......... .......... .......... 11% 37.3M 3s
##  15050K .......... .......... .......... .......... .......... 11% 77.5M 3s
##  15100K .......... .......... .......... .......... .......... 11% 63.3M 2s
##  15150K .......... .......... .......... .......... .......... 11% 66.1M 2s
##  15200K .......... .......... .......... .......... .......... 11% 93.4M 2s
##  15250K .......... .......... .......... .......... .......... 11% 86.6M 2s
##  15300K .......... .......... .......... .......... .......... 11% 72.0M 2s
##  15350K .......... .......... .......... .......... .......... 11% 55.0M 2s
##  15400K .......... .......... .......... .......... .......... 11% 76.1M 2s
##  15450K .......... .......... .......... .......... .......... 11% 73.6M 2s
##  15500K .......... .......... .......... .......... .......... 11%  104M 2s
##  15550K .......... .......... .......... .......... .......... 11% 58.5M 2s
##  15600K .......... .......... .......... .......... .......... 11% 22.8M 2s
##  15650K .......... .......... .......... .......... .......... 11% 94.2M 2s
##  15700K .......... .......... .......... .......... .......... 11% 99.5M 2s
##  15750K .......... .......... .......... .......... .......... 11% 86.9M 2s
##  15800K .......... .......... .......... .......... .......... 11%  104M 2s
##  15850K .......... .......... .......... .......... .......... 11% 19.1M 2s
##  15900K .......... .......... .......... .......... .......... 11% 61.6M 2s
##  15950K .......... .......... .......... .......... .......... 11% 57.4M 2s
##  16000K .......... .......... .......... .......... .......... 11%  107M 2s
##  16050K .......... .......... .......... .......... .......... 12%  101M 2s
##  16100K .......... .......... .......... .......... .......... 12% 96.4M 2s
##  16150K .......... .......... .......... .......... .......... 12% 56.9M 2s
##  16200K .......... .......... .......... .......... .......... 12% 85.7M 2s
##  16250K .......... .......... .......... .......... .......... 12% 88.5M 2s
##  16300K .......... .......... .......... .......... .......... 12%  100M 2s
##  16350K .......... .......... .......... .......... .......... 12% 69.6M 2s
##  16400K .......... .......... .......... .......... .......... 12% 95.2M 2s
##  16450K .......... .......... .......... .......... .......... 12%  104M 2s
##  16500K .......... .......... .......... .......... .......... 12% 47.0M 2s
##  16550K .......... .......... .......... .......... .......... 12% 71.0M 2s
##  16600K .......... .......... .......... .......... .......... 12%  103M 2s
##  16650K .......... .......... .......... .......... .......... 12% 83.2M 2s
##  16700K .......... .......... .......... .......... .......... 12% 85.0M 2s
##  16750K .......... .......... .......... .......... .......... 12% 78.0M 2s
##  16800K .......... .......... .......... .......... .......... 12% 79.8M 2s
##  16850K .......... .......... .......... .......... .......... 12% 81.1M 2s
##  16900K .......... .......... .......... .......... .......... 12% 92.5M 2s
##  16950K .......... .......... .......... .......... .......... 12% 97.1M 2s
##  17000K .......... .......... .......... .......... .......... 12% 44.8M 2s
##  17050K .......... .......... .......... .......... .......... 12% 89.1M 2s
##  17100K .......... .......... .......... .......... .......... 12% 99.2M 2s
##  17150K .......... .......... .......... .......... .......... 12% 77.0M 2s
##  17200K .......... .......... .......... .......... .......... 12% 85.7M 2s
##  17250K .......... .......... .......... .......... .......... 12% 29.3M 2s
##  17300K .......... .......... .......... .......... .......... 12% 79.3M 2s
##  17350K .......... .......... .......... .......... .......... 12% 80.2M 2s
##  17400K .......... .......... .......... .......... .......... 13%  105M 2s
##  17450K .......... .......... .......... .......... .......... 13% 72.7M 2s
##  17500K .......... .......... .......... .......... .......... 13% 96.0M 2s
##  17550K .......... .......... .......... .......... .......... 13% 51.5M 2s
##  17600K .......... .......... .......... .......... .......... 13% 89.3M 2s
##  17650K .......... .......... .......... .......... .......... 13% 92.7M 2s
##  17700K .......... .......... .......... .......... .......... 13% 88.7M 2s
##  17750K .......... .......... .......... .......... .......... 13% 40.3M 2s
##  17800K .......... .......... .......... .......... .......... 13% 85.8M 2s
##  17850K .......... .......... .......... .......... .......... 13% 95.3M 2s
##  17900K .......... .......... .......... .......... .......... 13% 19.7M 2s
##  17950K .......... .......... .......... .......... .......... 13% 68.2M 2s
##  18000K .......... .......... .......... .......... .......... 13% 92.9M 2s
##  18050K .......... .......... .......... .......... .......... 13%  107M 2s
##  18100K .......... .......... .......... .......... .......... 13%  107M 2s
##  18150K .......... .......... .......... .......... .......... 13% 84.8M 2s
##  18200K .......... .......... .......... .......... .......... 13% 22.8M 2s
##  18250K .......... .......... .......... .......... .......... 13%  100M 2s
##  18300K .......... .......... .......... .......... .......... 13% 39.5M 2s
##  18350K .......... .......... .......... .......... .......... 13% 75.4M 2s
##  18400K .......... .......... .......... .......... .......... 13%  118M 2s
##  18450K .......... .......... .......... .......... .......... 13% 82.8M 2s
##  18500K .......... .......... .......... .......... .......... 13% 93.1M 2s
##  18550K .......... .......... .......... .......... .......... 13% 93.8M 2s
##  18600K .......... .......... .......... .......... .......... 13% 38.2M 2s
##  18650K .......... .......... .......... .......... .......... 13% 68.6M 2s
##  18700K .......... .......... .......... .......... .......... 13% 95.4M 2s
##  18750K .......... .......... .......... .......... .......... 14% 80.9M 2s
##  18800K .......... .......... .......... .......... .......... 14% 95.1M 2s
##  18850K .......... .......... .......... .......... .......... 14%  118M 2s
##  18900K .......... .......... .......... .......... .......... 14% 72.9M 2s
##  18950K .......... .......... .......... .......... .......... 14% 69.6M 2s
##  19000K .......... .......... .......... .......... .......... 14% 85.3M 2s
##  19050K .......... .......... .......... .......... .......... 14% 80.5M 2s
##  19100K .......... .......... .......... .......... .......... 14% 52.4M 2s
##  19150K .......... .......... .......... .......... .......... 14% 55.6M 2s
##  19200K .......... .......... .......... .......... .......... 14% 68.4M 2s
##  19250K .......... .......... .......... .......... .......... 14%  103M 2s
##  19300K .......... .......... .......... .......... .......... 14% 97.8M 2s
##  19350K .......... .......... .......... .......... .......... 14% 88.7M 2s
##  19400K .......... .......... .......... .......... .......... 14% 39.2M 2s
##  19450K .......... .......... .......... .......... .......... 14% 88.1M 2s
##  19500K .......... .......... .......... .......... .......... 14% 69.2M 2s
##  19550K .......... .......... .......... .......... .......... 14% 79.3M 2s
##  19600K .......... .......... .......... .......... .......... 14% 77.3M 2s
##  19650K .......... .......... .......... .......... .......... 14% 92.4M 2s
##  19700K .......... .......... .......... .......... .......... 14%  108M 2s
##  19750K .......... .......... .......... .......... .......... 14%  103M 2s
##  19800K .......... .......... .......... .......... .......... 14% 52.0M 2s
##  19850K .......... .......... .......... .......... .......... 14% 91.5M 2s
##  19900K .......... .......... .......... .......... .......... 14% 68.2M 2s
##  19950K .......... .......... .......... .......... .......... 14% 80.7M 2s
##  20000K .......... .......... .......... .......... .......... 14% 87.1M 2s
##  20050K .......... .......... .......... .......... .......... 14% 90.6M 2s
##  20100K .......... .......... .......... .......... .......... 15% 39.9M 2s
##  20150K .......... .......... .......... .......... .......... 15% 76.3M 2s
##  20200K .......... .......... .......... .......... .......... 15%  101M 2s
##  20250K .......... .......... .......... .......... .......... 15% 79.5M 2s
##  20300K .......... .......... .......... .......... .......... 15% 71.3M 2s
##  20350K .......... .......... .......... .......... .......... 15% 62.4M 2s
##  20400K .......... .......... .......... .......... .......... 15% 57.0M 2s
##  20450K .......... .......... .......... .......... .......... 15% 98.9M 2s
##  20500K .......... .......... .......... .......... .......... 15% 74.4M 2s
##  20550K .......... .......... .......... .......... .......... 15% 56.6M 2s
##  20600K .......... .......... .......... .......... .......... 15% 73.7M 2s
##  20650K .......... .......... .......... .......... .......... 15% 94.1M 2s
##  20700K .......... .......... .......... .......... .......... 15% 77.3M 2s
##  20750K .......... .......... .......... .......... .......... 15% 71.3M 2s
##  20800K .......... .......... .......... .......... .......... 15% 67.2M 2s
##  20850K .......... .......... .......... .......... .......... 15% 77.2M 2s
##  20900K .......... .......... .......... .......... .......... 15% 76.6M 2s
##  20950K .......... .......... .......... .......... .......... 15% 81.8M 2s
##  21000K .......... .......... .......... .......... .......... 15% 63.8M 2s
##  21050K .......... .......... .......... .......... .......... 15% 81.4M 2s
##  21100K .......... .......... .......... .......... .......... 15% 82.2M 2s
##  21150K .......... .......... .......... .......... .......... 15% 64.7M 2s
##  21200K .......... .......... .......... .......... .......... 15% 96.5M 2s
##  21250K .......... .......... .......... .......... .......... 15% 71.8M 2s
##  21300K .......... .......... .......... .......... .......... 15% 84.2M 2s
##  21350K .......... .......... .......... .......... .......... 15% 64.4M 2s
##  21400K .......... .......... .......... .......... .......... 15% 86.2M 2s
##  21450K .......... .......... .......... .......... .......... 16% 64.5M 2s
##  21500K .......... .......... .......... .......... .......... 16% 76.9M 2s
##  21550K .......... .......... .......... .......... .......... 16% 72.3M 2s
##  21600K .......... .......... .......... .......... .......... 16% 81.9M 2s
##  21650K .......... .......... .......... .......... .......... 16% 82.1M 2s
##  21700K .......... .......... .......... .......... .......... 16% 77.8M 2s
##  21750K .......... .......... .......... .......... .......... 16% 80.8M 2s
##  21800K .......... .......... .......... .......... .......... 16% 80.1M 2s
##  21850K .......... .......... .......... .......... .......... 16%  110M 2s
##  21900K .......... .......... .......... .......... .......... 16% 90.7M 2s
##  21950K .......... .......... .......... .......... .......... 16% 66.9M 2s
##  22000K .......... .......... .......... .......... .......... 16%  106M 2s
##  22050K .......... .......... .......... .......... .......... 16% 46.0M 2s
##  22100K .......... .......... .......... .......... .......... 16% 68.1M 2s
##  22150K .......... .......... .......... .......... .......... 16% 78.6M 2s
##  22200K .......... .......... .......... .......... .......... 16%  100M 2s
##  22250K .......... .......... .......... .......... .......... 16%  115M 2s
##  22300K .......... .......... .......... .......... .......... 16%  102M 2s
##  22350K .......... .......... .......... .......... .......... 16% 86.7M 2s
##  22400K .......... .......... .......... .......... .......... 16% 73.1M 2s
##  22450K .......... .......... .......... .......... .......... 16% 83.8M 2s
##  22500K .......... .......... .......... .......... .......... 16% 85.5M 2s
##  22550K .......... .......... .......... .......... .......... 16% 72.4M 2s
##  22600K .......... .......... .......... .......... .......... 16% 94.8M 2s
##  22650K .......... .......... .......... .......... .......... 16% 83.4M 2s
##  22700K .......... .......... .......... .......... .......... 16% 79.3M 2s
##  22750K .......... .......... .......... .......... .......... 17% 80.4M 2s
##  22800K .......... .......... .......... .......... .......... 17% 80.9M 2s
##  22850K .......... .......... .......... .......... .......... 17% 95.9M 2s
##  22900K .......... .......... .......... .......... .......... 17%  113M 2s
##  22950K .......... .......... .......... .......... .......... 17% 63.8M 2s
##  23000K .......... .......... .......... .......... .......... 17% 78.6M 2s
##  23050K .......... .......... .......... .......... .......... 17%  100M 2s
##  23100K .......... .......... .......... .......... .......... 17% 84.3M 2s
##  23150K .......... .......... .......... .......... .......... 17% 79.9M 2s
##  23200K .......... .......... .......... .......... .......... 17% 81.4M 2s
##  23250K .......... .......... .......... .......... .......... 17% 93.1M 2s
##  23300K .......... .......... .......... .......... .......... 17%  105M 2s
##  23350K .......... .......... .......... .......... .......... 17% 75.8M 2s
##  23400K .......... .......... .......... .......... .......... 17% 87.9M 2s
##  23450K .......... .......... .......... .......... .......... 17% 98.0M 2s
##  23500K .......... .......... .......... .......... .......... 17%  108M 2s
##  23550K .......... .......... .......... .......... .......... 17% 12.9M 2s
##  23600K .......... .......... .......... .......... .......... 17%  109M 2s
##  23650K .......... .......... .......... .......... .......... 17% 86.3M 2s
##  23700K .......... .......... .......... .......... .......... 17%  107M 2s
##  23750K .......... .......... .......... .......... .......... 17%  110M 2s
##  23800K .......... .......... .......... .......... .......... 17%  134M 2s
##  23850K .......... .......... .......... .......... .......... 17% 41.2M 2s
##  23900K .......... .......... .......... .......... .......... 17%  107M 2s
##  23950K .......... .......... .......... .......... .......... 17% 68.6M 2s
##  24000K .......... .......... .......... .......... .......... 17% 79.6M 2s
##  24050K .......... .......... .......... .......... .......... 17% 86.4M 2s
##  24100K .......... .......... .......... .......... .......... 18% 98.2M 2s
##  24150K .......... .......... .......... .......... .......... 18% 93.8M 2s
##  24200K .......... .......... .......... .......... .......... 18%  105M 2s
##  24250K .......... .......... .......... .......... .......... 18% 95.5M 2s
##  24300K .......... .......... .......... .......... .......... 18% 5.64M 2s
##  24350K .......... .......... .......... .......... .......... 18% 87.4M 2s
##  24400K .......... .......... .......... .......... .......... 18% 99.1M 2s
##  24450K .......... .......... .......... .......... .......... 18% 98.1M 2s
##  24500K .......... .......... .......... .......... .......... 18% 87.8M 2s
##  24550K .......... .......... .......... .......... .......... 18%  123M 2s
##  24600K .......... .......... .......... .......... .......... 18%  115M 2s
##  24650K .......... .......... .......... .......... .......... 18% 74.4M 2s
##  24700K .......... .......... .......... .......... .......... 18% 14.8M 2s
##  24750K .......... .......... .......... .......... .......... 18%  102M 2s
##  24800K .......... .......... .......... .......... .......... 18%  105M 2s
##  24850K .......... .......... .......... .......... .......... 18%  135M 2s
##  24900K .......... .......... .......... .......... .......... 18%  131M 2s
##  24950K .......... .......... .......... .......... .......... 18%  119M 2s
##  25000K .......... .......... .......... .......... .......... 18% 51.8M 2s
##  25050K .......... .......... .......... .......... .......... 18% 35.7M 2s
##  25100K .......... .......... .......... .......... .......... 18% 56.1M 2s
##  25150K .......... .......... .......... .......... .......... 18% 22.6M 2s
##  25200K .......... .......... .......... .......... .......... 18%  108M 2s
##  25250K .......... .......... .......... .......... .......... 18%  102M 2s
##  25300K .......... .......... .......... .......... .......... 18%  131M 2s
##  25350K .......... .......... .......... .......... .......... 18%  122M 2s
##  25400K .......... .......... .......... .......... .......... 18%  136M 2s
##  25450K .......... .......... .......... .......... .......... 19% 36.5M 2s
##  25500K .......... .......... .......... .......... .......... 19% 37.5M 2s
##  25550K .......... .......... .......... .......... .......... 19% 97.4M 2s
##  25600K .......... .......... .......... .......... .......... 19%  100M 2s
##  25650K .......... .......... .......... .......... .......... 19% 54.7M 2s
##  25700K .......... .......... .......... .......... .......... 19% 38.3M 2s
##  25750K .......... .......... .......... .......... .......... 19% 32.8M 2s
##  25800K .......... .......... .......... .......... .......... 19% 7.53M 2s
##  25850K .......... .......... .......... .......... .......... 19% 57.9M 2s
##  25900K .......... .......... .......... .......... .......... 19% 95.1M 2s
##  25950K .......... .......... .......... .......... .......... 19% 67.8M 2s
##  26000K .......... .......... .......... .......... .......... 19%  106M 2s
##  26050K .......... .......... .......... .......... .......... 19% 96.1M 2s
##  26100K .......... .......... .......... .......... .......... 19%  101M 2s
##  26150K .......... .......... .......... .......... .......... 19% 92.9M 2s
##  26200K .......... .......... .......... .......... .......... 19% 56.1M 2s
##  26250K .......... .......... .......... .......... .......... 19% 87.7M 2s
##  26300K .......... .......... .......... .......... .......... 19% 84.2M 2s
##  26350K .......... .......... .......... .......... .......... 19% 51.1M 2s
##  26400K .......... .......... .......... .......... .......... 19% 97.7M 2s
##  26450K .......... .......... .......... .......... .......... 19% 53.0M 2s
##  26500K .......... .......... .......... .......... .......... 19% 24.4M 2s
##  26550K .......... .......... .......... .......... .......... 19% 92.3M 2s
##  26600K .......... .......... .......... .......... .......... 19% 97.7M 2s
##  26650K .......... .......... .......... .......... .......... 19% 17.2M 2s
##  26700K .......... .......... .......... .......... .......... 19% 79.1M 2s
##  26750K .......... .......... .......... .......... .......... 19% 86.8M 2s
##  26800K .......... .......... .......... .......... .......... 20%  101M 2s
##  26850K .......... .......... .......... .......... .......... 20% 82.7M 2s
##  26900K .......... .......... .......... .......... .......... 20%  114M 2s
##  26950K .......... .......... .......... .......... .......... 20% 48.7M 2s
##  27000K .......... .......... .......... .......... .......... 20% 79.2M 2s
##  27050K .......... .......... .......... .......... .......... 20% 49.1M 2s
##  27100K .......... .......... .......... .......... .......... 20% 88.9M 2s
##  27150K .......... .......... .......... .......... .......... 20% 40.6M 2s
##  27200K .......... .......... .......... .......... .......... 20%  102M 2s
##  27250K .......... .......... .......... .......... .......... 20% 92.2M 2s
##  27300K .......... .......... .......... .......... .......... 20% 95.0M 2s
##  27350K .......... .......... .......... .......... .......... 20% 49.3M 2s
##  27400K .......... .......... .......... .......... .......... 20%  112M 2s
##  27450K .......... .......... .......... .......... .......... 20% 27.4M 2s
##  27500K .......... .......... .......... .......... .......... 20% 77.5M 2s
##  27550K .......... .......... .......... .......... .......... 20% 52.6M 2s
##  27600K .......... .......... .......... .......... .......... 20% 88.6M 2s
##  27650K .......... .......... .......... .......... .......... 20% 85.1M 2s
##  27700K .......... .......... .......... .......... .......... 20%  103M 2s
##  27750K .......... .......... .......... .......... .......... 20% 84.6M 2s
##  27800K .......... .......... .......... .......... .......... 20% 19.5M 2s
##  27850K .......... .......... .......... .......... .......... 20% 65.7M 2s
##  27900K .......... .......... .......... .......... .......... 20%  118M 2s
##  27950K .......... .......... .......... .......... .......... 20% 84.7M 2s
##  28000K .......... .......... .......... .......... .......... 20%  119M 2s
##  28050K .......... .......... .......... .......... .......... 20%  122M 2s
##  28100K .......... .......... .......... .......... .......... 20% 96.2M 2s
##  28150K .......... .......... .......... .......... .......... 21% 15.7M 2s
##  28200K .......... .......... .......... .......... .......... 21%  105M 2s
##  28250K .......... .......... .......... .......... .......... 21% 93.0M 2s
##  28300K .......... .......... .......... .......... .......... 21% 99.6M 2s
##  28350K .......... .......... .......... .......... .......... 21% 72.2M 2s
##  28400K .......... .......... .......... .......... .......... 21%  119M 2s
##  28450K .......... .......... .......... .......... .......... 21%  109M 2s
##  28500K .......... .......... .......... .......... .......... 21% 11.6M 2s
##  28550K .......... .......... .......... .......... .......... 21% 93.3M 2s
##  28600K .......... .......... .......... .......... .......... 21%  113M 2s
##  28650K .......... .......... .......... .......... .......... 21%  125M 2s
##  28700K .......... .......... .......... .......... .......... 21%  119M 2s
##  28750K .......... .......... .......... .......... .......... 21% 96.1M 2s
##  28800K .......... .......... .......... .......... .......... 21% 21.8M 2s
##  28850K .......... .......... .......... .......... .......... 21% 75.3M 2s
##  28900K .......... .......... .......... .......... .......... 21% 55.3M 2s
##  28950K .......... .......... .......... .......... .......... 21% 92.1M 2s
##  29000K .......... .......... .......... .......... .......... 21% 72.7M 2s
##  29050K .......... .......... .......... .......... .......... 21%  103M 2s
##  29100K .......... .......... .......... .......... .......... 21%  116M 2s
##  29150K .......... .......... .......... .......... .......... 21%  107M 2s
##  29200K .......... .......... .......... .......... .......... 21% 5.04M 2s
##  29250K .......... .......... .......... .......... .......... 21%  117M 2s
##  29300K .......... .......... .......... .......... .......... 21%  125M 2s
##  29350K .......... .......... .......... .......... .......... 21%  104M 2s
##  29400K .......... .......... .......... .......... .......... 21%  119M 2s
##  29450K .......... .......... .......... .......... .......... 22%  130M 2s
##  29500K .......... .......... .......... .......... .......... 22%  117M 2s
##  29550K .......... .......... .......... .......... .......... 22% 27.2M 2s
##  29600K .......... .......... .......... .......... .......... 22%  105M 2s
##  29650K .......... .......... .......... .......... .......... 22% 48.6M 2s
##  29700K .......... .......... .......... .......... .......... 22% 76.1M 2s
##  29750K .......... .......... .......... .......... .......... 22% 73.2M 2s
##  29800K .......... .......... .......... .......... .......... 22% 91.1M 2s
##  29850K .......... .......... .......... .......... .......... 22%  115M 2s
##  29900K .......... .......... .......... .......... .......... 22%  125M 2s
##  29950K .......... .......... .......... .......... .......... 22% 40.5M 2s
##  30000K .......... .......... .......... .......... .......... 22%  101M 2s
##  30050K .......... .......... .......... .......... .......... 22%  132M 2s
##  30100K .......... .......... .......... .......... .......... 22% 71.7M 2s
##  30150K .......... .......... .......... .......... .......... 22% 69.9M 2s
##  30200K .......... .......... .......... .......... .......... 22% 87.6M 2s
##  30250K .......... .......... .......... .......... .......... 22%  100M 2s
##  30300K .......... .......... .......... .......... .......... 22%  108M 2s
##  30350K .......... .......... .......... .......... .......... 22% 74.0M 2s
##  30400K .......... .......... .......... .......... .......... 22% 33.6M 2s
##  30450K .......... .......... .......... .......... .......... 22% 95.2M 2s
##  30500K .......... .......... .......... .......... .......... 22%  129M 2s
##  30550K .......... .......... .......... .......... .......... 22% 30.9M 2s
##  30600K .......... .......... .......... .......... .......... 22% 95.3M 2s
##  30650K .......... .......... .......... .......... .......... 22%  109M 2s
##  30700K .......... .......... .......... .......... .......... 22% 62.9M 2s
##  30750K .......... .......... .......... .......... .......... 22% 79.6M 2s
##  30800K .......... .......... .......... .......... .......... 23% 97.1M 2s
##  30850K .......... .......... .......... .......... .......... 23% 19.9M 2s
##  30900K .......... .......... .......... .......... .......... 23% 68.6M 2s
##  30950K .......... .......... .......... .......... .......... 23% 61.9M 2s
##  31000K .......... .......... .......... .......... .......... 23% 84.5M 2s
##  31050K .......... .......... .......... .......... .......... 23% 68.4M 2s
##  31100K .......... .......... .......... .......... .......... 23%  117M 2s
##  31150K .......... .......... .......... .......... .......... 23% 73.6M 2s
##  31200K .......... .......... .......... .......... .......... 23% 75.9M 2s
##  31250K .......... .......... .......... .......... .......... 23%  128M 2s
##  31300K .......... .......... .......... .......... .......... 23% 45.9M 2s
##  31350K .......... .......... .......... .......... .......... 23% 22.3M 2s
##  31400K .......... .......... .......... .......... .......... 23% 79.7M 2s
##  31450K .......... .......... .......... .......... .......... 23%  121M 2s
##  31500K .......... .......... .......... .......... .......... 23%  128M 2s
##  31550K .......... .......... .......... .......... .......... 23% 93.0M 2s
##  31600K .......... .......... .......... .......... .......... 23%  109M 2s
##  31650K .......... .......... .......... .......... .......... 23% 91.4M 2s
##  31700K .......... .......... .......... .......... .......... 23% 38.8M 2s
##  31750K .......... .......... .......... .......... .......... 23% 89.1M 2s
##  31800K .......... .......... .......... .......... .......... 23%  132M 2s
##  31850K .......... .......... .......... .......... .......... 23% 14.9M 2s
##  31900K .......... .......... .......... .......... .......... 23%  118M 2s
##  31950K .......... .......... .......... .......... .......... 23%  101M 2s
##  32000K .......... .......... .......... .......... .......... 23%  126M 2s
##  32050K .......... .......... .......... .......... .......... 23%  129M 2s
##  32100K .......... .......... .......... .......... .......... 23% 47.1M 2s
##  32150K .......... .......... .......... .......... .......... 24% 46.5M 2s
##  32200K .......... .......... .......... .......... .......... 24% 92.6M 2s
##  32250K .......... .......... .......... .......... .......... 24% 55.3M 2s
##  32300K .......... .......... .......... .......... .......... 24% 87.6M 2s
##  32350K .......... .......... .......... .......... .......... 24% 76.0M 2s
##  32400K .......... .......... .......... .......... .......... 24%  109M 2s
##  32450K .......... .......... .......... .......... .......... 24% 96.0M 2s
##  32500K .......... .......... .......... .......... .......... 24% 78.8M 2s
##  32550K .......... .......... .......... .......... .......... 24% 30.4M 2s
##  32600K .......... .......... .......... .......... .......... 24% 73.4M 2s
##  32650K .......... .......... .......... .......... .......... 24%  108M 2s
##  32700K .......... .......... .......... .......... .......... 24% 46.1M 2s
##  32750K .......... .......... .......... .......... .......... 24% 73.2M 2s
##  32800K .......... .......... .......... .......... .......... 24%  102M 2s
##  32850K .......... .......... .......... .......... .......... 24%  100M 2s
##  32900K .......... .......... .......... .......... .......... 24%  103M 2s
##  32950K .......... .......... .......... .......... .......... 24% 90.0M 2s
##  33000K .......... .......... .......... .......... .......... 24% 23.3M 2s
##  33050K .......... .......... .......... .......... .......... 24% 19.4M 2s
##  33100K .......... .......... .......... .......... .......... 24%  102M 2s
##  33150K .......... .......... .......... .......... .......... 24% 94.1M 2s
##  33200K .......... .......... .......... .......... .......... 24%  115M 2s
##  33250K .......... .......... .......... .......... .......... 24%  117M 2s
##  33300K .......... .......... .......... .......... .......... 24%  108M 2s
##  33350K .......... .......... .......... .......... .......... 24% 60.4M 2s
##  33400K .......... .......... .......... .......... .......... 24% 44.0M 2s
##  33450K .......... .......... .......... .......... .......... 24% 45.5M 2s
##  33500K .......... .......... .......... .......... .......... 25% 84.3M 2s
##  33550K .......... .......... .......... .......... .......... 25% 80.6M 2s
##  33600K .......... .......... .......... .......... .......... 25% 95.0M 2s
##  33650K .......... .......... .......... .......... .......... 25%  118M 2s
##  33700K .......... .......... .......... .......... .......... 25% 94.8M 2s
##  33750K .......... .......... .......... .......... .......... 25% 82.1M 2s
##  33800K .......... .......... .......... .......... .......... 25% 88.8M 2s
##  33850K .......... .......... .......... .......... .......... 25% 85.5M 2s
##  33900K .......... .......... .......... .......... .......... 25% 89.8M 2s
##  33950K .......... .......... .......... .......... .......... 25% 78.9M 2s
##  34000K .......... .......... .......... .......... .......... 25% 69.0M 2s
##  34050K .......... .......... .......... .......... .......... 25% 85.1M 2s
##  34100K .......... .......... .......... .......... .......... 25% 85.3M 2s
##  34150K .......... .......... .......... .......... .......... 25% 72.7M 2s
##  34200K .......... .......... .......... .......... .......... 25% 94.6M 2s
##  34250K .......... .......... .......... .......... .......... 25% 73.7M 2s
##  34300K .......... .......... .......... .......... .......... 25%  104M 2s
##  34350K .......... .......... .......... .......... .......... 25% 88.3M 2s
##  34400K .......... .......... .......... .......... .......... 25% 81.4M 2s
##  34450K .......... .......... .......... .......... .......... 25% 64.2M 2s
##  34500K .......... .......... .......... .......... .......... 25%  119M 2s
##  34550K .......... .......... .......... .......... .......... 25% 85.2M 2s
##  34600K .......... .......... .......... .......... .......... 25% 66.9M 2s
##  34650K .......... .......... .......... .......... .......... 25% 75.7M 2s
##  34700K .......... .......... .......... .......... .......... 25%  108M 2s
##  34750K .......... .......... .......... .......... .......... 25% 70.7M 2s
##  34800K .......... .......... .......... .......... .......... 25% 88.4M 2s
##  34850K .......... .......... .......... .......... .......... 26%  107M 2s
##  34900K .......... .......... .......... .......... .......... 26%  105M 2s
##  34950K .......... .......... .......... .......... .......... 26% 51.3M 2s
##  35000K .......... .......... .......... .......... .......... 26%  107M 2s
##  35050K .......... .......... .......... .......... .......... 26% 91.7M 2s
##  35100K .......... .......... .......... .......... .......... 26% 56.9M 2s
##  35150K .......... .......... .......... .......... .......... 26% 92.9M 2s
##  35200K .......... .......... .......... .......... .......... 26%  103M 2s
##  35250K .......... .......... .......... .......... .......... 26%  111M 2s
##  35300K .......... .......... .......... .......... .......... 26%  100M 2s
##  35350K .......... .......... .......... .......... .......... 26% 67.9M 2s
##  35400K .......... .......... .......... .......... .......... 26% 90.3M 2s
##  35450K .......... .......... .......... .......... .......... 26% 90.8M 2s
##  35500K .......... .......... .......... .......... .......... 26% 74.8M 2s
##  35550K .......... .......... .......... .......... .......... 26% 70.6M 2s
##  35600K .......... .......... .......... .......... .......... 26% 97.3M 2s
##  35650K .......... .......... .......... .......... .......... 26% 74.8M 2s
##  35700K .......... .......... .......... .......... .......... 26%  109M 2s
##  35750K .......... .......... .......... .......... .......... 26%  103M 2s
##  35800K .......... .......... .......... .......... .......... 26% 99.2M 2s
##  35850K .......... .......... .......... .......... .......... 26% 62.7M 2s
##  35900K .......... .......... .......... .......... .......... 26% 14.2M 2s
##  35950K .......... .......... .......... .......... .......... 26% 82.8M 2s
##  36000K .......... .......... .......... .......... .......... 26%  131M 2s
##  36050K .......... .......... .......... .......... .......... 26%  117M 2s
##  36100K .......... .......... .......... .......... .......... 26%  137M 2s
##  36150K .......... .......... .......... .......... .......... 27%  114M 2s
##  36200K .......... .......... .......... .......... .......... 27%  124M 2s
##  36250K .......... .......... .......... .......... .......... 27% 18.2M 2s
##  36300K .......... .......... .......... .......... .......... 27%  102M 2s
##  36350K .......... .......... .......... .......... .......... 27% 85.6M 2s
##  36400K .......... .......... .......... .......... .......... 27%  116M 2s
##  36450K .......... .......... .......... .......... .......... 27%  139M 2s
##  36500K .......... .......... .......... .......... .......... 27%  135M 2s
##  36550K .......... .......... .......... .......... .......... 27%  112M 2s
##  36600K .......... .......... .......... .......... .......... 27% 37.7M 2s
##  36650K .......... .......... .......... .......... .......... 27% 89.5M 2s
##  36700K .......... .......... .......... .......... .......... 27%  130M 2s
##  36750K .......... .......... .......... .......... .......... 27% 86.4M 2s
##  36800K .......... .......... .......... .......... .......... 27%  105M 2s
##  36850K .......... .......... .......... .......... .......... 27%  128M 2s
##  36900K .......... .......... .......... .......... .......... 27%  129M 2s
##  36950K .......... .......... .......... .......... .......... 27% 48.3M 2s
##  37000K .......... .......... .......... .......... .......... 27% 58.4M 2s
##  37050K .......... .......... .......... .......... .......... 27% 47.4M 2s
##  37100K .......... .......... .......... .......... .......... 27% 95.3M 2s
##  37150K .......... .......... .......... .......... .......... 27% 27.3M 2s
##  37200K .......... .......... .......... .......... .......... 27%  105M 2s
##  37250K .......... .......... .......... .......... .......... 27%  134M 2s
##  37300K .......... .......... .......... .......... .......... 27%  123M 2s
##  37350K .......... .......... .......... .......... .......... 27%  114M 2s
##  37400K .......... .......... .......... .......... .......... 27%  131M 2s
##  37450K .......... .......... .......... .......... .......... 27% 28.9M 2s
##  37500K .......... .......... .......... .......... .......... 28%  115M 2s
##  37550K .......... .......... .......... .......... .......... 28% 98.7M 2s
##  37600K .......... .......... .......... .......... .......... 28%  132M 2s
##  37650K .......... .......... .......... .......... .......... 28% 62.2M 2s
##  37700K .......... .......... .......... .......... .......... 28%  129M 2s
##  37750K .......... .......... .......... .......... .......... 28% 32.8M 2s
##  37800K .......... .......... .......... .......... .......... 28%  103M 2s
##  37850K .......... .......... .......... .......... .......... 28% 48.2M 2s
##  37900K .......... .......... .......... .......... .......... 28% 97.9M 2s
##  37950K .......... .......... .......... .......... .......... 28% 42.7M 2s
##  38000K .......... .......... .......... .......... .......... 28%  100M 2s
##  38050K .......... .......... .......... .......... .......... 28% 45.2M 2s
##  38100K .......... .......... .......... .......... .......... 28% 87.2M 2s
##  38150K .......... .......... .......... .......... .......... 28% 93.3M 2s
##  38200K .......... .......... .......... .......... .......... 28%  112M 2s
##  38250K .......... .......... .......... .......... .......... 28% 40.6M 2s
##  38300K .......... .......... .......... .......... .......... 28% 84.3M 2s
##  38350K .......... .......... .......... .......... .......... 28% 96.0M 2s
##  38400K .......... .......... .......... .......... .......... 28% 54.3M 2s
##  38450K .......... .......... .......... .......... .......... 28%  105M 2s
##  38500K .......... .......... .......... .......... .......... 28% 96.0M 2s
##  38550K .......... .......... .......... .......... .......... 28% 86.4M 2s
##  38600K .......... .......... .......... .......... .......... 28% 47.2M 2s
##  38650K .......... .......... .......... .......... .......... 28% 91.9M 2s
##  38700K .......... .......... .......... .......... .......... 28% 61.7M 2s
##  38750K .......... .......... .......... .......... .......... 28% 58.0M 2s
##  38800K .......... .......... .......... .......... .......... 28%  119M 2s
##  38850K .......... .......... .......... .......... .......... 29% 76.4M 2s
##  38900K .......... .......... .......... .......... .......... 29% 57.2M 2s
##  38950K .......... .......... .......... .......... .......... 29% 69.8M 2s
##  39000K .......... .......... .......... .......... .......... 29% 92.3M 2s
##  39050K .......... .......... .......... .......... .......... 29% 43.4M 2s
##  39100K .......... .......... .......... .......... .......... 29% 98.0M 2s
##  39150K .......... .......... .......... .......... .......... 29% 81.2M 2s
##  39200K .......... .......... .......... .......... .......... 29% 94.9M 2s
##  39250K .......... .......... .......... .......... .......... 29%  107M 2s
##  39300K .......... .......... .......... .......... .......... 29% 81.8M 2s
##  39350K .......... .......... .......... .......... .......... 29% 38.8M 2s
##  39400K .......... .......... .......... .......... .......... 29%  105M 2s
##  39450K .......... .......... .......... .......... .......... 29%  127M 2s
##  39500K .......... .......... .......... .......... .......... 29% 71.4M 2s
##  39550K .......... .......... .......... .......... .......... 29% 84.3M 2s
##  39600K .......... .......... .......... .......... .......... 29% 71.3M 2s
##  39650K .......... .......... .......... .......... .......... 29% 55.7M 2s
##  39700K .......... .......... .......... .......... .......... 29% 76.7M 2s
##  39750K .......... .......... .......... .......... .......... 29% 94.2M 2s
##  39800K .......... .......... .......... .......... .......... 29% 46.3M 2s
##  39850K .......... .......... .......... .......... .......... 29%  106M 2s
##  39900K .......... .......... .......... .......... .......... 29%  116M 2s
##  39950K .......... .......... .......... .......... .......... 29% 72.4M 2s
##  40000K .......... .......... .......... .......... .......... 29%  102M 2s
##  40050K .......... .......... .......... .......... .......... 29% 90.7M 2s
##  40100K .......... .......... .......... .......... .......... 29% 58.9M 2s
##  40150K .......... .......... .......... .......... .......... 29% 96.9M 2s
##  40200K .......... .......... .......... .......... .......... 30% 69.2M 2s
##  40250K .......... .......... .......... .......... .......... 30% 95.5M 2s
##  40300K .......... .......... .......... .......... .......... 30% 77.2M 2s
##  40350K .......... .......... .......... .......... .......... 30% 77.7M 2s
##  40400K .......... .......... .......... .......... .......... 30%  112M 2s
##  40450K .......... .......... .......... .......... .......... 30%  103M 2s
##  40500K .......... .......... .......... .......... .......... 30% 66.3M 2s
##  40550K .......... .......... .......... .......... .......... 30% 87.3M 2s
##  40600K .......... .......... .......... .......... .......... 30% 43.3M 2s
##  40650K .......... .......... .......... .......... .......... 30% 89.6M 2s
##  40700K .......... .......... .......... .......... .......... 30%  123M 2s
##  40750K .......... .......... .......... .......... .......... 30% 90.4M 2s
##  40800K .......... .......... .......... .......... .......... 30%  114M 2s
##  40850K .......... .......... .......... .......... .......... 30% 82.9M 2s
##  40900K .......... .......... .......... .......... .......... 30%  104M 2s
##  40950K .......... .......... .......... .......... .......... 30% 76.5M 2s
##  41000K .......... .......... .......... .......... .......... 30% 74.3M 2s
##  41050K .......... .......... .......... .......... .......... 30% 97.4M 2s
##  41100K .......... .......... .......... .......... .......... 30%  108M 2s
##  41150K .......... .......... .......... .......... .......... 30% 64.4M 2s
##  41200K .......... .......... .......... .......... .......... 30% 85.0M 2s
##  41250K .......... .......... .......... .......... .......... 30% 93.4M 2s
##  41300K .......... .......... .......... .......... .......... 30% 68.6M 2s
##  41350K .......... .......... .......... .......... .......... 30% 94.3M 2s
##  41400K .......... .......... .......... .......... .......... 30%  101M 2s
##  41450K .......... .......... .......... .......... .......... 30% 75.3M 2s
##  41500K .......... .......... .......... .......... .......... 30% 99.0M 2s
##  41550K .......... .......... .......... .......... .......... 31% 64.6M 2s
##  41600K .......... .......... .......... .......... .......... 31%  111M 2s
##  41650K .......... .......... .......... .......... .......... 31% 82.3M 2s
##  41700K .......... .......... .......... .......... .......... 31%  111M 2s
##  41750K .......... .......... .......... .......... .......... 31% 71.1M 2s
##  41800K .......... .......... .......... .......... .......... 31%  122M 2s
##  41850K .......... .......... .......... .......... .......... 31% 99.3M 2s
##  41900K .......... .......... .......... .......... .......... 31% 80.1M 2s
##  41950K .......... .......... .......... .......... .......... 31% 81.0M 2s
##  42000K .......... .......... .......... .......... .......... 31% 85.0M 2s
##  42050K .......... .......... .......... .......... .......... 31% 81.0M 2s
##  42100K .......... .......... .......... .......... .......... 31%  110M 2s
##  42150K .......... .......... .......... .......... .......... 31% 86.7M 2s
##  42200K .......... .......... .......... .......... .......... 31% 95.4M 2s
##  42250K .......... .......... .......... .......... .......... 31% 84.4M 2s
##  42300K .......... .......... .......... .......... .......... 31% 95.0M 2s
##  42350K .......... .......... .......... .......... .......... 31% 82.2M 2s
##  42400K .......... .......... .......... .......... .......... 31% 85.8M 2s
##  42450K .......... .......... .......... .......... .......... 31% 97.9M 2s
##  42500K .......... .......... .......... .......... .......... 31% 85.6M 2s
##  42550K .......... .......... .......... .......... .......... 31% 91.9M 2s
##  42600K .......... .......... .......... .......... .......... 31% 21.6M 2s
##  42650K .......... .......... .......... .......... .......... 31%  127M 2s
##  42700K .......... .......... .......... .......... .......... 31%  139M 2s
##  42750K .......... .......... .......... .......... .......... 31% 52.4M 2s
##  42800K .......... .......... .......... .......... .......... 31% 8.11M 2s
##  42850K .......... .......... .......... .......... .......... 31% 79.5M 2s
##  42900K .......... .......... .......... .......... .......... 32%  122M 2s
##  42950K .......... .......... .......... .......... .......... 32% 94.2M 2s
##  43000K .......... .......... .......... .......... .......... 32% 89.9M 2s
##  43050K .......... .......... .......... .......... .......... 32%  116M 2s
##  43100K .......... .......... .......... .......... .......... 32%  107M 2s
##  43150K .......... .......... .......... .......... .......... 32%  115M 2s
##  43200K .......... .......... .......... .......... .......... 32%  140M 2s
##  43250K .......... .......... .......... .......... .......... 32%  109M 2s
##  43300K .......... .......... .......... .......... .......... 32%  111M 2s
##  43350K .......... .......... .......... .......... .......... 32%  109M 2s
##  43400K .......... .......... .......... .......... .......... 32%  116M 2s
##  43450K .......... .......... .......... .......... .......... 32%  139M 2s
##  43500K .......... .......... .......... .......... .......... 32%  123M 2s
##  43550K .......... .......... .......... .......... .......... 32%  101M 2s
##  43600K .......... .......... .......... .......... .......... 32%  134M 2s
##  43650K .......... .......... .......... .......... .......... 32%  127M 2s
##  43700K .......... .......... .......... .......... .......... 32%  136M 2s
##  43750K .......... .......... .......... .......... .......... 32%  103M 2s
##  43800K .......... .......... .......... .......... .......... 32%  140M 2s
##  43850K .......... .......... .......... .......... .......... 32%  129M 2s
##  43900K .......... .......... .......... .......... .......... 32%  138M 2s
##  43950K .......... .......... .......... .......... .......... 32% 54.6M 2s
##  44000K .......... .......... .......... .......... .......... 32%  123M 2s
##  44050K .......... .......... .......... .......... .......... 32%  142M 2s
##  44100K .......... .......... .......... .......... .......... 32% 94.5M 2s
##  44150K .......... .......... .......... .......... .......... 32% 14.2M 2s
##  44200K .......... .......... .......... .......... .......... 33%  107M 2s
##  44250K .......... .......... .......... .......... .......... 33%  107M 2s
##  44300K .......... .......... .......... .......... .......... 33%  117M 2s
##  44350K .......... .......... .......... .......... .......... 33%  101M 2s
##  44400K .......... .......... .......... .......... .......... 33% 24.1M 2s
##  44450K .......... .......... .......... .......... .......... 33% 75.2M 2s
##  44500K .......... .......... .......... .......... .......... 33%  114M 2s
##  44550K .......... .......... .......... .......... .......... 33% 99.4M 2s
##  44600K .......... .......... .......... .......... .......... 33%  100M 2s
##  44650K .......... .......... .......... .......... .......... 33%  119M 2s
##  44700K .......... .......... .......... .......... .......... 33% 42.5M 2s
##  44750K .......... .......... .......... .......... .......... 33% 30.2M 2s
##  44800K .......... .......... .......... .......... .......... 33%  126M 2s
##  44850K .......... .......... .......... .......... .......... 33%  134M 2s
##  44900K .......... .......... .......... .......... .......... 33%  142M 2s
##  44950K .......... .......... .......... .......... .......... 33% 24.0M 2s
##  45000K .......... .......... .......... .......... .......... 33%  113M 2s
##  45050K .......... .......... .......... .......... .......... 33%  115M 2s
##  45100K .......... .......... .......... .......... .......... 33%  115M 2s
##  45150K .......... .......... .......... .......... .......... 33%  120M 2s
##  45200K .......... .......... .......... .......... .......... 33% 32.7M 2s
##  45250K .......... .......... .......... .......... .......... 33% 30.0M 2s
##  45300K .......... .......... .......... .......... .......... 33%  120M 2s
##  45350K .......... .......... .......... .......... .......... 33%  118M 2s
##  45400K .......... .......... .......... .......... .......... 33%  114M 1s
##  45450K .......... .......... .......... .......... .......... 33%  126M 1s
##  45500K .......... .......... .......... .......... .......... 33% 26.0M 1s
##  45550K .......... .......... .......... .......... .......... 34% 73.4M 1s
##  45600K .......... .......... .......... .......... .......... 34%  107M 1s
##  45650K .......... .......... .......... .......... .......... 34%  131M 1s
##  45700K .......... .......... .......... .......... .......... 34%  125M 1s
##  45750K .......... .......... .......... .......... .......... 34% 22.2M 1s
##  45800K .......... .......... .......... .......... .......... 34%  110M 1s
##  45850K .......... .......... .......... .......... .......... 34%  112M 1s
##  45900K .......... .......... .......... .......... .......... 34%  120M 1s
##  45950K .......... .......... .......... .......... .......... 34%  107M 1s
##  46000K .......... .......... .......... .......... .......... 34% 36.7M 1s
##  46050K .......... .......... .......... .......... .......... 34% 86.2M 1s
##  46100K .......... .......... .......... .......... .......... 34% 94.1M 1s
##  46150K .......... .......... .......... .......... .......... 34%  103M 1s
##  46200K .......... .......... .......... .......... .......... 34%  114M 1s
##  46250K .......... .......... .......... .......... .......... 34%  137M 1s
##  46300K .......... .......... .......... .......... .......... 34% 11.6M 1s
##  46350K .......... .......... .......... .......... .......... 34% 83.7M 1s
##  46400K .......... .......... .......... .......... .......... 34% 93.5M 1s
##  46450K .......... .......... .......... .......... .......... 34%  127M 1s
##  46500K .......... .......... .......... .......... .......... 34% 16.6M 1s
##  46550K .......... .......... .......... .......... .......... 34% 91.0M 1s
##  46600K .......... .......... .......... .......... .......... 34% 80.6M 1s
##  46650K .......... .......... .......... .......... .......... 34% 82.7M 1s
##  46700K .......... .......... .......... .......... .......... 34%  111M 1s
##  46750K .......... .......... .......... .......... .......... 34% 25.0M 1s
##  46800K .......... .......... .......... .......... .......... 34% 61.0M 1s
##  46850K .......... .......... .......... .......... .......... 34% 67.1M 1s
##  46900K .......... .......... .......... .......... .......... 35% 98.9M 1s
##  46950K .......... .......... .......... .......... .......... 35% 64.8M 1s
##  47000K .......... .......... .......... .......... .......... 35% 98.3M 1s
##  47050K .......... .......... .......... .......... .......... 35% 38.2M 1s
##  47100K .......... .......... .......... .......... .......... 35% 94.3M 1s
##  47150K .......... .......... .......... .......... .......... 35% 38.6M 1s
##  47200K .......... .......... .......... .......... .......... 35%  102M 1s
##  47250K .......... .......... .......... .......... .......... 35%  106M 1s
##  47300K .......... .......... .......... .......... .......... 35% 31.2M 1s
##  47350K .......... .......... .......... .......... .......... 35%  107M 1s
##  47400K .......... .......... .......... .......... .......... 35%  106M 1s
##  47450K .......... .......... .......... .......... .......... 35%  119M 1s
##  47500K .......... .......... .......... .......... .......... 35%  114M 1s
##  47550K .......... .......... .......... .......... .......... 35% 27.0M 1s
##  47600K .......... .......... .......... .......... .......... 35% 74.5M 1s
##  47650K .......... .......... .......... .......... .......... 35% 24.2M 1s
##  47700K .......... .......... .......... .......... .......... 35% 98.1M 1s
##  47750K .......... .......... .......... .......... .......... 35% 90.4M 1s
##  47800K .......... .......... .......... .......... .......... 35% 42.7M 1s
##  47850K .......... .......... .......... .......... .......... 35% 94.8M 1s
##  47900K .......... .......... .......... .......... .......... 35%  125M 1s
##  47950K .......... .......... .......... .......... .......... 35% 57.3M 1s
##  48000K .......... .......... .......... .......... .......... 35% 58.1M 1s
##  48050K .......... .......... .......... .......... .......... 35% 88.2M 1s
##  48100K .......... .......... .......... .......... .......... 35% 78.1M 1s
##  48150K .......... .......... .......... .......... .......... 35%  104M 1s
##  48200K .......... .......... .......... .......... .......... 35%  108M 1s
##  48250K .......... .......... .......... .......... .......... 36% 50.2M 1s
##  48300K .......... .......... .......... .......... .......... 36% 32.1M 1s
##  48350K .......... .......... .......... .......... .......... 36% 41.8M 1s
##  48400K .......... .......... .......... .......... .......... 36% 76.5M 1s
##  48450K .......... .......... .......... .......... .......... 36%  122M 1s
##  48500K .......... .......... .......... .......... .......... 36% 27.8M 1s
##  48550K .......... .......... .......... .......... .......... 36%  105M 1s
##  48600K .......... .......... .......... .......... .......... 36% 76.9M 1s
##  48650K .......... .......... .......... .......... .......... 36%  116M 1s
##  48700K .......... .......... .......... .......... .......... 36%  107M 1s
##  48750K .......... .......... .......... .......... .......... 36% 40.2M 1s
##  48800K .......... .......... .......... .......... .......... 36%  102M 1s
##  48850K .......... .......... .......... .......... .......... 36% 65.8M 1s
##  48900K .......... .......... .......... .......... .......... 36% 97.7M 1s
##  48950K .......... .......... .......... .......... .......... 36% 57.6M 1s
##  49000K .......... .......... .......... .......... .......... 36% 22.1M 1s
##  49050K .......... .......... .......... .......... .......... 36% 95.4M 1s
##  49100K .......... .......... .......... .......... .......... 36%  128M 1s
##  49150K .......... .......... .......... .......... .......... 36% 84.5M 1s
##  49200K .......... .......... .......... .......... .......... 36%  122M 1s
##  49250K .......... .......... .......... .......... .......... 36% 42.0M 1s
##  49300K .......... .......... .......... .......... .......... 36% 83.7M 1s
##  49350K .......... .......... .......... .......... .......... 36% 89.0M 1s
##  49400K .......... .......... .......... .......... .......... 36% 36.7M 1s
##  49450K .......... .......... .......... .......... .......... 36%  113M 1s
##  49500K .......... .......... .......... .......... .......... 36% 30.8M 1s
##  49550K .......... .......... .......... .......... .......... 36% 69.9M 1s
##  49600K .......... .......... .......... .......... .......... 37%  112M 1s
##  49650K .......... .......... .......... .......... .......... 37% 44.1M 1s
##  49700K .......... .......... .......... .......... .......... 37% 99.2M 1s
##  49750K .......... .......... .......... .......... .......... 37% 94.7M 1s
##  49800K .......... .......... .......... .......... .......... 37% 57.3M 1s
##  49850K .......... .......... .......... .......... .......... 37% 80.2M 1s
##  49900K .......... .......... .......... .......... .......... 37% 22.9M 1s
##  49950K .......... .......... .......... .......... .......... 37% 97.6M 1s
##  50000K .......... .......... .......... .......... .......... 37% 82.3M 1s
##  50050K .......... .......... .......... .......... .......... 37% 78.9M 1s
##  50100K .......... .......... .......... .......... .......... 37% 95.3M 1s
##  50150K .......... .......... .......... .......... .......... 37% 26.3M 1s
##  50200K .......... .......... .......... .......... .......... 37%  114M 1s
##  50250K .......... .......... .......... .......... .......... 37% 95.2M 1s
##  50300K .......... .......... .......... .......... .......... 37%  105M 1s
##  50350K .......... .......... .......... .......... .......... 37% 89.3M 1s
##  50400K .......... .......... .......... .......... .......... 37% 18.9M 1s
##  50450K .......... .......... .......... .......... .......... 37% 81.5M 1s
##  50500K .......... .......... .......... .......... .......... 37%  104M 1s
##  50550K .......... .......... .......... .......... .......... 37% 91.0M 1s
##  50600K .......... .......... .......... .......... .......... 37%  125M 1s
##  50650K .......... .......... .......... .......... .......... 37% 21.5M 1s
##  50700K .......... .......... .......... .......... .......... 37% 90.5M 1s
##  50750K .......... .......... .......... .......... .......... 37% 91.9M 1s
##  50800K .......... .......... .......... .......... .......... 37% 80.7M 1s
##  50850K .......... .......... .......... .......... .......... 37%  109M 1s
##  50900K .......... .......... .......... .......... .......... 38% 23.7M 1s
##  50950K .......... .......... .......... .......... .......... 38% 96.3M 1s
##  51000K .......... .......... .......... .......... .......... 38%  122M 1s
##  51050K .......... .......... .......... .......... .......... 38% 59.8M 1s
##  51100K .......... .......... .......... .......... .......... 38%  112M 1s
##  51150K .......... .......... .......... .......... .......... 38% 24.8M 1s
##  51200K .......... .......... .......... .......... .......... 38% 96.1M 1s
##  51250K .......... .......... .......... .......... .......... 38%  120M 1s
##  51300K .......... .......... .......... .......... .......... 38% 53.2M 1s
##  51350K .......... .......... .......... .......... .......... 38%  106M 1s
##  51400K .......... .......... .......... .......... .......... 38% 27.7M 1s
##  51450K .......... .......... .......... .......... .......... 38% 87.7M 1s
##  51500K .......... .......... .......... .......... .......... 38%  118M 1s
##  51550K .......... .......... .......... .......... .......... 38% 87.1M 1s
##  51600K .......... .......... .......... .......... .......... 38%  101M 1s
##  51650K .......... .......... .......... .......... .......... 38% 29.2M 1s
##  51700K .......... .......... .......... .......... .......... 38% 80.5M 1s
##  51750K .......... .......... .......... .......... .......... 38%  102M 1s
##  51800K .......... .......... .......... .......... .......... 38%  105M 1s
##  51850K .......... .......... .......... .......... .......... 38%  120M 1s
##  51900K .......... .......... .......... .......... .......... 38%  112M 1s
##  51950K .......... .......... .......... .......... .......... 38% 24.9M 1s
##  52000K .......... .......... .......... .......... .......... 38% 95.8M 1s
##  52050K .......... .......... .......... .......... .......... 38%  108M 1s
##  52100K .......... .......... .......... .......... .......... 38%  119M 1s
##  52150K .......... .......... .......... .......... .......... 38% 94.2M 1s
##  52200K .......... .......... .......... .......... .......... 38% 37.5M 1s
##  52250K .......... .......... .......... .......... .......... 39% 83.4M 1s
##  52300K .......... .......... .......... .......... .......... 39% 70.1M 1s
##  52350K .......... .......... .......... .......... .......... 39% 89.1M 1s
##  52400K .......... .......... .......... .......... .......... 39%  107M 1s
##  52450K .......... .......... .......... .......... .......... 39% 42.3M 1s
##  52500K .......... .......... .......... .......... .......... 39%  109M 1s
##  52550K .......... .......... .......... .......... .......... 39% 29.8M 1s
##  52600K .......... .......... .......... .......... .......... 39%  106M 1s
##  52650K .......... .......... .......... .......... .......... 39%  115M 1s
##  52700K .......... .......... .......... .......... .......... 39% 88.7M 1s
##  52750K .......... .......... .......... .......... .......... 39% 83.2M 1s
##  52800K .......... .......... .......... .......... .......... 39% 23.0M 1s
##  52850K .......... .......... .......... .......... .......... 39%  109M 1s
##  52900K .......... .......... .......... .......... .......... 39% 96.7M 1s
##  52950K .......... .......... .......... .......... .......... 39%  108M 1s
##  53000K .......... .......... .......... .......... .......... 39%  102M 1s
##  53050K .......... .......... .......... .......... .......... 39% 26.5M 1s
##  53100K .......... .......... .......... .......... .......... 39%  120M 1s
##  53150K .......... .......... .......... .......... .......... 39%  100M 1s
##  53200K .......... .......... .......... .......... .......... 39%  121M 1s
##  53250K .......... .......... .......... .......... .......... 39%  126M 1s
##  53300K .......... .......... .......... .......... .......... 39% 37.6M 1s
##  53350K .......... .......... .......... .......... .......... 39% 56.7M 1s
##  53400K .......... .......... .......... .......... .......... 39%  104M 1s
##  53450K .......... .......... .......... .......... .......... 39%  107M 1s
##  53500K .......... .......... .......... .......... .......... 39%  113M 1s
##  53550K .......... .......... .......... .......... .......... 39% 24.5M 1s
##  53600K .......... .......... .......... .......... .......... 40% 74.3M 1s
##  53650K .......... .......... .......... .......... .......... 40% 97.1M 1s
##  53700K .......... .......... .......... .......... .......... 40% 97.6M 1s
##  53750K .......... .......... .......... .......... .......... 40% 95.6M 1s
##  53800K .......... .......... .......... .......... .......... 40%  118M 1s
##  53850K .......... .......... .......... .......... .......... 40% 35.1M 1s
##  53900K .......... .......... .......... .......... .......... 40% 31.7M 1s
##  53950K .......... .......... .......... .......... .......... 40% 81.5M 1s
##  54000K .......... .......... .......... .......... .......... 40%  114M 1s
##  54050K .......... .......... .......... .......... .......... 40%  113M 1s
##  54100K .......... .......... .......... .......... .......... 40%  102M 1s
##  54150K .......... .......... .......... .......... .......... 40% 24.7M 1s
##  54200K .......... .......... .......... .......... .......... 40% 87.7M 1s
##  54250K .......... .......... .......... .......... .......... 40% 99.9M 1s
##  54300K .......... .......... .......... .......... .......... 40% 90.7M 1s
##  54350K .......... .......... .......... .......... .......... 40%  104M 1s
##  54400K .......... .......... .......... .......... .......... 40% 19.6M 1s
##  54450K .......... .......... .......... .......... .......... 40% 85.6M 1s
##  54500K .......... .......... .......... .......... .......... 40% 96.0M 1s
##  54550K .......... .......... .......... .......... .......... 40% 93.8M 1s
##  54600K .......... .......... .......... .......... .......... 40%  114M 1s
##  54650K .......... .......... .......... .......... .......... 40% 24.4M 1s
##  54700K .......... .......... .......... .......... .......... 40% 83.3M 1s
##  54750K .......... .......... .......... .......... .......... 40% 81.2M 1s
##  54800K .......... .......... .......... .......... .......... 40% 92.3M 1s
##  54850K .......... .......... .......... .......... .......... 40%  121M 1s
##  54900K .......... .......... .......... .......... .......... 40% 25.3M 1s
##  54950K .......... .......... .......... .......... .......... 41% 90.6M 1s
##  55000K .......... .......... .......... .......... .......... 41%  113M 1s
##  55050K .......... .......... .......... .......... .......... 41% 57.2M 1s
##  55100K .......... .......... .......... .......... .......... 41%  113M 1s
##  55150K .......... .......... .......... .......... .......... 41% 19.6M 1s
##  55200K .......... .......... .......... .......... .......... 41%  113M 1s
##  55250K .......... .......... .......... .......... .......... 41%  110M 1s
##  55300K .......... .......... .......... .......... .......... 41% 62.0M 1s
##  55350K .......... .......... .......... .......... .......... 41%  102M 1s
##  55400K .......... .......... .......... .......... .......... 41%  114M 1s
##  55450K .......... .......... .......... .......... .......... 41% 46.1M 1s
##  55500K .......... .......... .......... .......... .......... 41%  113M 1s
##  55550K .......... .......... .......... .......... .......... 41% 51.4M 1s
##  55600K .......... .......... .......... .......... .......... 41%  105M 1s
##  55650K .......... .......... .......... .......... .......... 41% 31.7M 1s
##  55700K .......... .......... .......... .......... .......... 41% 93.0M 1s
##  55750K .......... .......... .......... .......... .......... 41% 97.6M 1s
##  55800K .......... .......... .......... .......... .......... 41%  108M 1s
##  55850K .......... .......... .......... .......... .......... 41%  104M 1s
##  55900K .......... .......... .......... .......... .......... 41%  117M 1s
##  55950K .......... .......... .......... .......... .......... 41% 33.7M 1s
##  56000K .......... .......... .......... .......... .......... 41%  106M 1s
##  56050K .......... .......... .......... .......... .......... 41% 38.5M 1s
##  56100K .......... .......... .......... .......... .......... 41%  108M 1s
##  56150K .......... .......... .......... .......... .......... 41%  101M 1s
##  56200K .......... .......... .......... .......... .......... 41% 79.7M 1s
##  56250K .......... .......... .......... .......... .......... 41%  117M 1s
##  56300K .......... .......... .......... .......... .......... 42% 35.8M 1s
##  56350K .......... .......... .......... .......... .......... 42% 79.4M 1s
##  56400K .......... .......... .......... .......... .......... 42%  113M 1s
##  56450K .......... .......... .......... .......... .......... 42% 64.6M 1s
##  56500K .......... .......... .......... .......... .......... 42%  104M 1s
##  56550K .......... .......... .......... .......... .......... 42% 36.8M 1s
##  56600K .......... .......... .......... .......... .......... 42% 93.3M 1s
##  56650K .......... .......... .......... .......... .......... 42%  112M 1s
##  56700K .......... .......... .......... .......... .......... 42% 28.7M 1s
##  56750K .......... .......... .......... .......... .......... 42% 80.2M 1s
##  56800K .......... .......... .......... .......... .......... 42%  111M 1s
##  56850K .......... .......... .......... .......... .......... 42% 97.0M 1s
##  56900K .......... .......... .......... .......... .......... 42% 95.8M 1s
##  56950K .......... .......... .......... .......... .......... 42% 15.1M 1s
##  57000K .......... .......... .......... .......... .......... 42% 86.5M 1s
##  57050K .......... .......... .......... .......... .......... 42%  110M 1s
##  57100K .......... .......... .......... .......... .......... 42%  118M 1s
##  57150K .......... .......... .......... .......... .......... 42% 89.6M 1s
##  57200K .......... .......... .......... .......... .......... 42% 32.7M 1s
##  57250K .......... .......... .......... .......... .......... 42% 92.3M 1s
##  57300K .......... .......... .......... .......... .......... 42%  113M 1s
##  57350K .......... .......... .......... .......... .......... 42% 61.8M 1s
##  57400K .......... .......... .......... .......... .......... 42%  110M 1s
##  57450K .......... .......... .......... .......... .......... 42% 41.1M 1s
##  57500K .......... .......... .......... .......... .......... 42%  102M 1s
##  57550K .......... .......... .......... .......... .......... 42% 44.1M 1s
##  57600K .......... .......... .......... .......... .......... 43% 85.6M 1s
##  57650K .......... .......... .......... .......... .......... 43%  119M 1s
##  57700K .......... .......... .......... .......... .......... 43%  117M 1s
##  57750K .......... .......... .......... .......... .......... 43% 34.4M 1s
##  57800K .......... .......... .......... .......... .......... 43%  101M 1s
##  57850K .......... .......... .......... .......... .......... 43% 44.2M 1s
##  57900K .......... .......... .......... .......... .......... 43%  104M 1s
##  57950K .......... .......... .......... .......... .......... 43% 99.1M 1s
##  58000K .......... .......... .......... .......... .......... 43% 80.2M 1s
##  58050K .......... .......... .......... .......... .......... 43%  101M 1s
##  58100K .......... .......... .......... .......... .......... 43% 20.3M 1s
##  58150K .......... .......... .......... .......... .......... 43%  103M 1s
##  58200K .......... .......... .......... .......... .......... 43%  113M 1s
##  58250K .......... .......... .......... .......... .......... 43%  115M 1s
##  58300K .......... .......... .......... .......... .......... 43%  126M 1s
##  58350K .......... .......... .......... .......... .......... 43% 15.2M 1s
##  58400K .......... .......... .......... .......... .......... 43% 47.0M 1s
##  58450K .......... .......... .......... .......... .......... 43% 92.0M 1s
##  58500K .......... .......... .......... .......... .......... 43% 97.0M 1s
##  58550K .......... .......... .......... .......... .......... 43% 99.9M 1s
##  58600K .......... .......... .......... .......... .......... 43% 34.5M 1s
##  58650K .......... .......... .......... .......... .......... 43%  106M 1s
##  58700K .......... .......... .......... .......... .......... 43%  124M 1s
##  58750K .......... .......... .......... .......... .......... 43% 39.1M 1s
##  58800K .......... .......... .......... .......... .......... 43%  124M 1s
##  58850K .......... .......... .......... .......... .......... 43% 35.4M 1s
##  58900K .......... .......... .......... .......... .......... 43%  112M 1s
##  58950K .......... .......... .......... .......... .......... 44%  114M 1s
##  59000K .......... .......... .......... .......... .......... 44% 34.9M 1s
##  59050K .......... .......... .......... .......... .......... 44%  124M 1s
##  59100K .......... .......... .......... .......... .......... 44% 32.0M 1s
##  59150K .......... .......... .......... .......... .......... 44%  104M 1s
##  59200K .......... .......... .......... .......... .......... 44% 40.5M 1s
##  59250K .......... .......... .......... .......... .......... 44% 98.1M 1s
##  59300K .......... .......... .......... .......... .......... 44%  128M 1s
##  59350K .......... .......... .......... .......... .......... 44% 40.5M 1s
##  59400K .......... .......... .......... .......... .......... 44%  102M 1s
##  59450K .......... .......... .......... .......... .......... 44%  131M 1s
##  59500K .......... .......... .......... .......... .......... 44% 42.3M 1s
##  59550K .......... .......... .......... .......... .......... 44% 96.6M 1s
##  59600K .......... .......... .......... .......... .......... 44% 38.6M 1s
##  59650K .......... .......... .......... .......... .......... 44%  104M 1s
##  59700K .......... .......... .......... .......... .......... 44%  131M 1s
##  59750K .......... .......... .......... .......... .......... 44% 29.5M 1s
##  59800K .......... .......... .......... .......... .......... 44%  127M 1s
##  59850K .......... .......... .......... .......... .......... 44% 51.6M 1s
##  59900K .......... .......... .......... .......... .......... 44%  122M 1s
##  59950K .......... .......... .......... .......... .......... 44% 70.7M 1s
##  60000K .......... .......... .......... .......... .......... 44% 34.9M 1s
##  60050K .......... .......... .......... .......... .......... 44%  125M 1s
##  60100K .......... .......... .......... .......... .......... 44% 29.3M 1s
##  60150K .......... .......... .......... .......... .......... 44% 87.6M 1s
##  60200K .......... .......... .......... .......... .......... 44%  121M 1s
##  60250K .......... .......... .......... .......... .......... 44% 67.5M 1s
##  60300K .......... .......... .......... .......... .......... 45%  120M 1s
##  60350K .......... .......... .......... .......... .......... 45% 30.9M 1s
##  60400K .......... .......... .......... .......... .......... 45% 73.8M 1s
##  60450K .......... .......... .......... .......... .......... 45%  131M 1s
##  60500K .......... .......... .......... .......... .......... 45%  108M 1s
##  60550K .......... .......... .......... .......... .......... 45%  101M 1s
##  60600K .......... .......... .......... .......... .......... 45%  135M 1s
##  60650K .......... .......... .......... .......... .......... 45% 20.5M 1s
##  60700K .......... .......... .......... .......... .......... 45%  133M 1s
##  60750K .......... .......... .......... .......... .......... 45% 59.8M 1s
##  60800K .......... .......... .......... .......... .......... 45%  115M 1s
##  60850K .......... .......... .......... .......... .......... 45%  113M 1s
##  60900K .......... .......... .......... .......... .......... 45% 29.1M 1s
##  60950K .......... .......... .......... .......... .......... 45% 85.1M 1s
##  61000K .......... .......... .......... .......... .......... 45% 42.9M 1s
##  61050K .......... .......... .......... .......... .......... 45%  104M 1s
##  61100K .......... .......... .......... .......... .......... 45%  101M 1s
##  61150K .......... .......... .......... .......... .......... 45% 65.9M 1s
##  61200K .......... .......... .......... .......... .......... 45% 57.2M 1s
##  61250K .......... .......... .......... .......... .......... 45% 27.0M 1s
##  61300K .......... .......... .......... .......... .......... 45%  111M 1s
##  61350K .......... .......... .......... .......... .......... 45% 7.02M 1s
##  61400K .......... .......... .......... .......... .......... 45%  115M 1s
##  61450K .......... .......... .......... .......... .......... 45%  131M 1s
##  61500K .......... .......... .......... .......... .......... 45%  126M 1s
##  61550K .......... .......... .......... .......... .......... 45%  114M 1s
##  61600K .......... .......... .......... .......... .......... 45%  132M 1s
##  61650K .......... .......... .......... .......... .......... 46% 6.02M 1s
##  61700K .......... .......... .......... .......... .......... 46%  116M 1s
##  61750K .......... .......... .......... .......... .......... 46%  120M 1s
##  61800K .......... .......... .......... .......... .......... 46%  138M 1s
##  61850K .......... .......... .......... .......... .......... 46%  125M 1s
##  61900K .......... .......... .......... .......... .......... 46%  141M 1s
##  61950K .......... .......... .......... .......... .......... 46%  113M 1s
##  62000K .......... .......... .......... .......... .......... 46%  135M 1s
##  62050K .......... .......... .......... .......... .......... 46%  143M 1s
##  62100K .......... .......... .......... .......... .......... 46% 22.2M 1s
##  62150K .......... .......... .......... .......... .......... 46%  104M 1s
##  62200K .......... .......... .......... .......... .......... 46%  136M 1s
##  62250K .......... .......... .......... .......... .......... 46%  138M 1s
##  62300K .......... .......... .......... .......... .......... 46% 17.2M 1s
##  62350K .......... .......... .......... .......... .......... 46% 68.8M 1s
##  62400K .......... .......... .......... .......... .......... 46%  114M 1s
##  62450K .......... .......... .......... .......... .......... 46% 18.9M 1s
##  62500K .......... .......... .......... .......... .......... 46% 94.8M 1s
##  62550K .......... .......... .......... .......... .......... 46%  110M 1s
##  62600K .......... .......... .......... .......... .......... 46%  123M 1s
##  62650K .......... .......... .......... .......... .......... 46% 17.1M 1s
##  62700K .......... .......... .......... .......... .......... 46% 94.9M 1s
##  62750K .......... .......... .......... .......... .......... 46%  107M 1s
##  62800K .......... .......... .......... .......... .......... 46%  123M 1s
##  62850K .......... .......... .......... .......... .......... 46% 27.4M 1s
##  62900K .......... .......... .......... .......... .......... 46% 52.6M 1s
##  62950K .......... .......... .......... .......... .......... 46% 88.1M 1s
##  63000K .......... .......... .......... .......... .......... 47% 68.7M 1s
##  63050K .......... .......... .......... .......... .......... 47% 27.5M 1s
##  63100K .......... .......... .......... .......... .......... 47% 66.1M 1s
##  63150K .......... .......... .......... .......... .......... 47% 51.8M 1s
##  63200K .......... .......... .......... .......... .......... 47% 49.4M 1s
##  63250K .......... .......... .......... .......... .......... 47% 97.4M 1s
##  63300K .......... .......... .......... .......... .......... 47% 27.9M 1s
##  63350K .......... .......... .......... .......... .......... 47% 53.7M 1s
##  63400K .......... .......... .......... .......... .......... 47% 71.7M 1s
##  63450K .......... .......... .......... .......... .......... 47%  103M 1s
##  63500K .......... .......... .......... .......... .......... 47% 16.1M 1s
##  63550K .......... .......... .......... .......... .......... 47% 90.5M 1s
##  63600K .......... .......... .......... .......... .......... 47%  141M 1s
##  63650K .......... .......... .......... .......... .......... 47% 16.9M 1s
##  63700K .......... .......... .......... .......... .......... 47% 92.5M 1s
##  63750K .......... .......... .......... .......... .......... 47%  108M 1s
##  63800K .......... .......... .......... .......... .......... 47%  129M 1s
##  63850K .......... .......... .......... .......... .......... 47% 17.0M 1s
##  63900K .......... .......... .......... .......... .......... 47% 98.9M 1s
##  63950K .......... .......... .......... .......... .......... 47%  113M 1s
##  64000K .......... .......... .......... .......... .......... 47%  132M 1s
##  64050K .......... .......... .......... .......... .......... 47% 14.6M 1s
##  64100K .......... .......... .......... .......... .......... 47%  108M 1s
##  64150K .......... .......... .......... .......... .......... 47% 3.43M 1s
##  64200K .......... .......... .......... .......... .......... 47%  126M 1s
##  64250K .......... .......... .......... .......... .......... 47% 16.5M 1s
##  64300K .......... .......... .......... .......... .......... 47% 70.9M 1s
##  64350K .......... .......... .......... .......... .......... 48%  119M 1s
##  64400K .......... .......... .......... .......... .......... 48%  140M 1s
##  64450K .......... .......... .......... .......... .......... 48%  127M 1s
##  64500K .......... .......... .......... .......... .......... 48% 36.2M 1s
##  64550K .......... .......... .......... .......... .......... 48% 40.6M 1s
##  64600K .......... .......... .......... .......... .......... 48%  105M 1s
##  64650K .......... .......... .......... .......... .......... 48%  136M 1s
##  64700K .......... .......... .......... .......... .......... 48% 15.9M 1s
##  64750K .......... .......... .......... .......... .......... 48%  103M 1s
##  64800K .......... .......... .......... .......... .......... 48% 90.1M 1s
##  64850K .......... .......... .......... .......... .......... 48% 18.6M 1s
##  64900K .......... .......... .......... .......... .......... 48%  130M 1s
##  64950K .......... .......... .......... .......... .......... 48% 89.3M 1s
##  65000K .......... .......... .......... .......... .......... 48% 19.6M 1s
##  65050K .......... .......... .......... .......... .......... 48%  111M 1s
##  65100K .......... .......... .......... .......... .......... 48% 78.9M 1s
##  65150K .......... .......... .......... .......... .......... 48% 20.7M 1s
##  65200K .......... .......... .......... .......... .......... 48% 46.2M 1s
##  65250K .......... .......... .......... .......... .......... 48% 17.3M 1s
##  65300K .......... .......... .......... .......... .......... 48% 93.8M 1s
##  65350K .......... .......... .......... .......... .......... 48%  109M 1s
##  65400K .......... .......... .......... .......... .......... 48% 16.1M 1s
##  65450K .......... .......... .......... .......... .......... 48%  132M 1s
##  65500K .......... .......... .......... .......... .......... 48% 59.7M 1s
##  65550K .......... .......... .......... .......... .......... 48% 26.1M 1s
##  65600K .......... .......... .......... .......... .......... 48%  133M 1s
##  65650K .......... .......... .......... .......... .......... 49% 30.4M 1s
##  65700K .......... .......... .......... .......... .......... 49% 26.3M 1s
##  65750K .......... .......... .......... .......... .......... 49% 35.3M 1s
##  65800K .......... .......... .......... .......... .......... 49%  103M 1s
##  65850K .......... .......... .......... .......... .......... 49% 27.5M 1s
##  65900K .......... .......... .......... .......... .......... 49% 26.9M 1s
##  65950K .......... .......... .......... .......... .......... 49% 35.4M 1s
##  66000K .......... .......... .......... .......... .......... 49% 21.3M 1s
##  66050K .......... .......... .......... .......... .......... 49%  114M 1s
##  66100K .......... .......... .......... .......... .......... 49% 38.8M 1s
##  66150K .......... .......... .......... .......... .......... 49% 20.2M 1s
##  66200K .......... .......... .......... .......... .......... 49% 41.6M 1s
##  66250K .......... .......... .......... .......... .......... 49%  120M 1s
##  66300K .......... .......... .......... .......... .......... 49% 28.3M 1s
##  66350K .......... .......... .......... .......... .......... 49% 25.3M 1s
##  66400K .......... .......... .......... .......... .......... 49% 58.4M 1s
##  66450K .......... .......... .......... .......... .......... 49% 16.5M 1s
##  66500K .......... .......... .......... .......... .......... 49% 94.5M 1s
##  66550K .......... .......... .......... .......... .......... 49%  101M 1s
##  66600K .......... .......... .......... .......... .......... 49% 18.1M 1s
##  66650K .......... .......... .......... .......... .......... 49% 78.1M 1s
##  66700K .......... .......... .......... .......... .......... 49%  116M 1s
##  66750K .......... .......... .......... .......... .......... 49% 19.5M 1s
##  66800K .......... .......... .......... .......... .......... 49% 59.3M 1s
##  66850K .......... .......... .......... .......... .......... 49% 21.3M 1s
##  66900K .......... .......... .......... .......... .......... 49% 86.0M 1s
##  66950K .......... .......... .......... .......... .......... 49% 64.1M 1s
##  67000K .......... .......... .......... .......... .......... 50% 29.6M 1s
##  67050K .......... .......... .......... .......... .......... 50% 35.4M 1s
##  67100K .......... .......... .......... .......... .......... 50% 54.7M 1s
##  67150K .......... .......... .......... .......... .......... 50% 17.3M 1s
##  67200K .......... .......... .......... .......... .......... 50%  115M 1s
##  67250K .......... .......... .......... .......... .......... 50% 42.9M 1s
##  67300K .......... .......... .......... .......... .......... 50% 13.9M 1s
##  67350K .......... .......... .......... .......... .......... 50% 89.1M 1s
##  67400K .......... .......... .......... .......... .......... 50% 14.6M 1s
##  67450K .......... .......... .......... .......... .......... 50% 80.8M 1s
##  67500K .......... .......... .......... .......... .......... 50%  137M 1s
##  67550K .......... .......... .......... .......... .......... 50% 13.8M 1s
##  67600K .......... .......... .......... .......... .......... 50%  130M 1s
##  67650K .......... .......... .......... .......... .......... 50% 14.1M 1s
##  67700K .......... .......... .......... .......... .......... 50% 87.4M 1s
##  67750K .......... .......... .......... .......... .......... 50% 94.2M 1s
##  67800K .......... .......... .......... .......... .......... 50% 13.1M 1s
##  67850K .......... .......... .......... .......... .......... 50%  125M 1s
##  67900K .......... .......... .......... .......... .......... 50%  131M 1s
##  67950K .......... .......... .......... .......... .......... 50% 15.0M 1s
##  68000K .......... .......... .......... .......... .......... 50% 93.5M 1s
##  68050K .......... .......... .......... .......... .......... 50%  130M 1s
##  68100K .......... .......... .......... .......... .......... 50% 17.8M 1s
##  68150K .......... .......... .......... .......... .......... 50% 92.3M 1s
##  68200K .......... .......... .......... .......... .......... 50% 94.0M 1s
##  68250K .......... .......... .......... .......... .......... 50% 17.8M 1s
##  68300K .......... .......... .......... .......... .......... 50% 63.1M 1s
##  68350K .......... .......... .......... .......... .......... 51% 15.0M 1s
##  68400K .......... .......... .......... .......... .......... 51% 99.8M 1s
##  68450K .......... .......... .......... .......... .......... 51%  134M 1s
##  68500K .......... .......... .......... .......... .......... 51% 15.2M 1s
##  68550K .......... .......... .......... .......... .......... 51%  102M 1s
##  68600K .......... .......... .......... .......... .......... 51%  125M 1s
##  68650K .......... .......... .......... .......... .......... 51% 12.9M 1s
##  68700K .......... .......... .......... .......... .......... 51%  124M 1s
##  68750K .......... .......... .......... .......... .......... 51%  117M 1s
##  68800K .......... .......... .......... .......... .......... 51% 13.1M 1s
##  68850K .......... .......... .......... .......... .......... 51%  128M 1s
##  68900K .......... .......... .......... .......... .......... 51%  127M 1s
##  68950K .......... .......... .......... .......... .......... 51% 17.6M 1s
##  69000K .......... .......... .......... .......... .......... 51%  113M 1s
##  69050K .......... .......... .......... .......... .......... 51%  139M 1s
##  69100K .......... .......... .......... .......... .......... 51% 14.6M 1s
##  69150K .......... .......... .......... .......... .......... 51%  103M 1s
##  69200K .......... .......... .......... .......... .......... 51% 14.1M 1s
##  69250K .......... .......... .......... .......... .......... 51% 98.4M 1s
##  69300K .......... .......... .......... .......... .......... 51%  127M 1s
##  69350K .......... .......... .......... .......... .......... 51% 16.8M 1s
##  69400K .......... .......... .......... .......... .......... 51% 74.1M 1s
##  69450K .......... .......... .......... .......... .......... 51% 25.0M 1s
##  69500K .......... .......... .......... .......... .......... 51%  103M 1s
##  69550K .......... .......... .......... .......... .......... 51% 18.7M 1s
##  69600K .......... .......... .......... .......... .......... 51% 84.3M 1s
##  69650K .......... .......... .......... .......... .......... 51% 27.5M 1s
##  69700K .......... .......... .......... .......... .......... 52% 30.8M 1s
##  69750K .......... .......... .......... .......... .......... 52% 23.5M 1s
##  69800K .......... .......... .......... .......... .......... 52%  131M 1s
##  69850K .......... .......... .......... .......... .......... 52% 29.4M 1s
##  69900K .......... .......... .......... .......... .......... 52% 25.9M 1s
##  69950K .......... .......... .......... .......... .......... 52% 37.3M 1s
##  70000K .......... .......... .......... .......... .......... 52% 27.0M 1s
##  70050K .......... .......... .......... .......... .......... 52% 94.6M 1s
##  70100K .......... .......... .......... .......... .......... 52% 33.4M 1s
##  70150K .......... .......... .......... .......... .......... 52% 24.7M 1s
##  70200K .......... .......... .......... .......... .......... 52% 29.4M 1s
##  70250K .......... .......... .......... .......... .......... 52%  127M 1s
##  70300K .......... .......... .......... .......... .......... 52% 25.2M 1s
##  70350K .......... .......... .......... .......... .......... 52% 29.1M 1s
##  70400K .......... .......... .......... .......... .......... 52% 25.3M 1s
##  70450K .......... .......... .......... .......... .......... 52% 28.5M 1s
##  70500K .......... .......... .......... .......... .......... 52%  101M 1s
##  70550K .......... .......... .......... .......... .......... 52% 25.8M 1s
##  70600K .......... .......... .......... .......... .......... 52% 29.5M 1s
##  70650K .......... .......... .......... .......... .......... 52% 26.6M 1s
##  70700K .......... .......... .......... .......... .......... 52%  124M 1s
##  70750K .......... .......... .......... .......... .......... 52% 30.4M 1s
##  70800K .......... .......... .......... .......... .......... 52% 24.7M 1s
##  70850K .......... .......... .......... .......... .......... 52% 34.2M 1s
##  70900K .......... .......... .......... .......... .......... 52% 22.4M 1s
##  70950K .......... .......... .......... .......... .......... 52%  109M 1s
##  71000K .......... .......... .......... .......... .......... 52% 35.6M 1s
##  71050K .......... .......... .......... .......... .......... 53% 20.6M 1s
##  71100K .......... .......... .......... .......... .......... 53% 56.2M 1s
##  71150K .......... .......... .......... .......... .......... 53% 38.1M 1s
##  71200K .......... .......... .......... .......... .......... 53% 73.5M 1s
##  71250K .......... .......... .......... .......... .......... 53% 22.3M 1s
##  71300K .......... .......... .......... .......... .......... 53% 46.0M 1s
##  71350K .......... .......... .......... .......... .......... 53% 20.9M 1s
##  71400K .......... .......... .......... .......... .......... 53%  105M 1s
##  71450K .......... .......... .......... .......... .......... 53% 37.9M 1s
##  71500K .......... .......... .......... .......... .......... 53% 23.2M 1s
##  71550K .......... .......... .......... .......... .......... 53% 28.9M 1s
##  71600K .......... .......... .......... .......... .......... 53% 27.8M 1s
##  71650K .......... .......... .......... .......... .......... 53% 36.9M 1s
##  71700K .......... .......... .......... .......... .......... 53% 93.0M 1s
##  71750K .......... .......... .......... .......... .......... 53% 21.0M 1s
##  71800K .......... .......... .......... .......... .......... 53% 45.2M 1s
##  71850K .......... .......... .......... .......... .......... 53% 19.5M 1s
##  71900K .......... .......... .......... .......... .......... 53%  127M 1s
##  71950K .......... .......... .......... .......... .......... 53% 24.4M 1s
##  72000K .......... .......... .......... .......... .......... 53% 34.0M 1s
##  72050K .......... .......... .......... .......... .......... 53% 22.5M 1s
##  72100K .......... .......... .......... .......... .......... 53% 34.9M 1s
##  72150K .......... .......... .......... .......... .......... 53% 99.8M 1s
##  72200K .......... .......... .......... .......... .......... 53% 24.0M 1s
##  72250K .......... .......... .......... .......... .......... 53% 32.6M 1s
##  72300K .......... .......... .......... .......... .......... 53% 25.3M 1s
##  72350K .......... .......... .......... .......... .......... 54% 31.7M 1s
##  72400K .......... .......... .......... .......... .......... 54%  109M 1s
##  72450K .......... .......... .......... .......... .......... 54% 28.8M 1s
##  72500K .......... .......... .......... .......... .......... 54% 23.0M 1s
##  72550K .......... .......... .......... .......... .......... 54% 38.7M 1s
##  72600K .......... .......... .......... .......... .......... 54%  131M 1s
##  72650K .......... .......... .......... .......... .......... 54% 26.5M 1s
##  72700K .......... .......... .......... .......... .......... 54% 33.5M 1s
##  72750K .......... .......... .......... .......... .......... 54% 26.8M 1s
##  72800K .......... .......... .......... .......... .......... 54% 28.9M 1s
##  72850K .......... .......... .......... .......... .......... 54%  115M 1s
##  72900K .......... .......... .......... .......... .......... 54% 27.1M 1s
##  72950K .......... .......... .......... .......... .......... 54% 21.1M 1s
##  73000K .......... .......... .......... .......... .......... 54% 40.3M 1s
##  73050K .......... .......... .......... .......... .......... 54%  116M 1s
##  73100K .......... .......... .......... .......... .......... 54% 21.5M 1s
##  73150K .......... .......... .......... .......... .......... 54%  114M 1s
##  73200K .......... .......... .......... .......... .......... 54% 54.0M 1s
##  73250K .......... .......... .......... .......... .......... 54% 21.3M 1s
##  73300K .......... .......... .......... .......... .......... 54% 96.5M 1s
##  73350K .......... .......... .......... .......... .......... 54%  114M 1s
##  73400K .......... .......... .......... .......... .......... 54%  124M 1s
##  73450K .......... .......... .......... .......... .......... 54% 17.3M 1s
##  73500K .......... .......... .......... .......... .......... 54%  106M 1s
##  73550K .......... .......... .......... .......... .......... 54% 94.2M 1s
##  73600K .......... .......... .......... .......... .......... 54% 79.1M 1s
##  73650K .......... .......... .......... .......... .......... 54% 17.3M 1s
##  73700K .......... .......... .......... .......... .......... 55%  120M 1s
##  73750K .......... .......... .......... .......... .......... 55%  114M 1s
##  73800K .......... .......... .......... .......... .......... 55%  122M 1s
##  73850K .......... .......... .......... .......... .......... 55% 16.0M 1s
##  73900K .......... .......... .......... .......... .......... 55% 95.3M 1s
##  73950K .......... .......... .......... .......... .......... 55% 82.0M 1s
##  74000K .......... .......... .......... .......... .......... 55% 21.2M 1s
##  74050K .......... .......... .......... .......... .......... 55%  120M 1s
##  74100K .......... .......... .......... .......... .......... 55% 51.7M 1s
##  74150K .......... .......... .......... .......... .......... 55% 91.8M 1s
##  74200K .......... .......... .......... .......... .......... 55% 20.5M 1s
##  74250K .......... .......... .......... .......... .......... 55% 85.8M 1s
##  74300K .......... .......... .......... .......... .......... 55% 41.8M 1s
##  74350K .......... .......... .......... .......... .......... 55% 94.8M 1s
##  74400K .......... .......... .......... .......... .......... 55% 19.5M 1s
##  74450K .......... .......... .......... .......... .......... 55% 40.9M 1s
##  74500K .......... .......... .......... .......... .......... 55% 99.6M 1s
##  74550K .......... .......... .......... .......... .......... 55%  104M 1s
##  74600K .......... .......... .......... .......... .......... 55% 39.7M 1s
##  74650K .......... .......... .......... .......... .......... 55% 21.9M 1s
##  74700K .......... .......... .......... .......... .......... 55%  109M 1s
##  74750K .......... .......... .......... .......... .......... 55% 93.8M 1s
##  74800K .......... .......... .......... .......... .......... 55% 48.2M 1s
##  74850K .......... .......... .......... .......... .......... 55% 29.3M 1s
##  74900K .......... .......... .......... .......... .......... 55% 48.1M 1s
##  74950K .......... .......... .......... .......... .......... 55% 95.1M 1s
##  75000K .......... .......... .......... .......... .......... 55% 38.5M 1s
##  75050K .......... .......... .......... .......... .......... 56% 37.2M 1s
##  75100K .......... .......... .......... .......... .......... 56% 54.2M 1s
##  75150K .......... .......... .......... .......... .......... 56% 41.1M 1s
##  75200K .......... .......... .......... .......... .......... 56% 99.3M 1s
##  75250K .......... .......... .......... .......... .......... 56% 18.1M 1s
##  75300K .......... .......... .......... .......... .......... 56%  109M 1s
##  75350K .......... .......... .......... .......... .......... 56% 97.5M 1s
##  75400K .......... .......... .......... .......... .......... 56%  140M 1s
##  75450K .......... .......... .......... .......... .......... 56% 14.5M 1s
##  75500K .......... .......... .......... .......... .......... 56%  110M 1s
##  75550K .......... .......... .......... .......... .......... 56%  108M 1s
##  75600K .......... .......... .......... .......... .......... 56% 14.7M 1s
##  75650K .......... .......... .......... .......... .......... 56%  121M 1s
##  75700K .......... .......... .......... .......... .......... 56%  112M 1s
##  75750K .......... .......... .......... .......... .......... 56%  100M 1s
##  75800K .......... .......... .......... .......... .......... 56% 18.3M 1s
##  75850K .......... .......... .......... .......... .......... 56%  103M 1s
##  75900K .......... .......... .......... .......... .......... 56% 86.7M 1s
##  75950K .......... .......... .......... .......... .......... 56% 67.9M 1s
##  76000K .......... .......... .......... .......... .......... 56% 21.2M 1s
##  76050K .......... .......... .......... .......... .......... 56% 82.6M 1s
##  76100K .......... .......... .......... .......... .......... 56%  106M 1s
##  76150K .......... .......... .......... .......... .......... 56% 93.2M 1s
##  76200K .......... .......... .......... .......... .......... 56% 16.5M 1s
##  76250K .......... .......... .......... .......... .......... 56%  126M 1s
##  76300K .......... .......... .......... .......... .......... 56%  114M 1s
##  76350K .......... .......... .......... .......... .......... 56% 14.0M 1s
##  76400K .......... .......... .......... .......... .......... 57% 92.4M 1s
##  76450K .......... .......... .......... .......... .......... 57%  124M 1s
##  76500K .......... .......... .......... .......... .......... 57%  141M 1s
##  76550K .......... .......... .......... .......... .......... 57% 16.3M 1s
##  76600K .......... .......... .......... .......... .......... 57% 60.4M 1s
##  76650K .......... .......... .......... .......... .......... 57% 59.8M 1s
##  76700K .......... .......... .......... .......... .......... 57% 78.0M 1s
##  76750K .......... .......... .......... .......... .......... 57% 24.6M 1s
##  76800K .......... .......... .......... .......... .......... 57% 29.6M 1s
##  76850K .......... .......... .......... .......... .......... 57% 62.6M 1s
##  76900K .......... .......... .......... .......... .......... 57% 77.4M 1s
##  76950K .......... .......... .......... .......... .......... 57% 62.8M 1s
##  77000K .......... .......... .......... .......... .......... 57% 77.0M 1s
##  77050K .......... .......... .......... .......... .......... 57% 23.8M 1s
##  77100K .......... .......... .......... .......... .......... 57% 74.1M 1s
##  77150K .......... .......... .......... .......... .......... 57% 64.1M 1s
##  77200K .......... .......... .......... .......... .......... 57% 52.0M 1s
##  77250K .......... .......... .......... .......... .......... 57% 58.1M 1s
##  77300K .......... .......... .......... .......... .......... 57% 69.0M 1s
##  77350K .......... .......... .......... .......... .......... 57% 63.5M 1s
##  77400K .......... .......... .......... .......... .......... 57% 20.6M 1s
##  77450K .......... .......... .......... .......... .......... 57% 64.2M 1s
##  77500K .......... .......... .......... .......... .......... 57% 78.5M 1s
##  77550K .......... .......... .......... .......... .......... 57% 69.3M 1s
##  77600K .......... .......... .......... .......... .......... 57% 29.9M 1s
##  77650K .......... .......... .......... .......... .......... 57% 47.5M 1s
##  77700K .......... .......... .......... .......... .......... 57% 72.6M 1s
##  77750K .......... .......... .......... .......... .......... 58% 74.4M 1s
##  77800K .......... .......... .......... .......... .......... 58% 36.0M 1s
##  77850K .......... .......... .......... .......... .......... 58% 36.4M 1s
##  77900K .......... .......... .......... .......... .......... 58% 65.6M 1s
##  77950K .......... .......... .......... .......... .......... 58% 59.5M 1s
##  78000K .......... .......... .......... .......... .......... 58% 58.9M 1s
##  78050K .......... .......... .......... .......... .......... 58% 55.9M 1s
##  78100K .......... .......... .......... .......... .......... 58% 33.4M 1s
##  78150K .......... .......... .......... .......... .......... 58% 57.3M 1s
##  78200K .......... .......... .......... .......... .......... 58% 45.0M 1s
##  78250K .......... .......... .......... .......... .......... 58% 78.7M 1s
##  78300K .......... .......... .......... .......... .......... 58% 30.9M 1s
##  78350K .......... .......... .......... .......... .......... 58% 58.0M 1s
##  78400K .......... .......... .......... .......... .......... 58% 41.3M 1s
##  78450K .......... .......... .......... .......... .......... 58% 19.5M 1s
##  78500K .......... .......... .......... .......... .......... 58% 60.1M 1s
##  78550K .......... .......... .......... .......... .......... 58% 26.8M 1s
##  78600K .......... .......... .......... .......... .......... 58% 70.9M 1s
##  78650K .......... .......... .......... .......... .......... 58% 44.3M 1s
##  78700K .......... .......... .......... .......... .......... 58% 61.0M 1s
##  78750K .......... .......... .......... .......... .......... 58% 67.7M 1s
##  78800K .......... .......... .......... .......... .......... 58% 72.5M 1s
##  78850K .......... .......... .......... .......... .......... 58% 21.5M 1s
##  78900K .......... .......... .......... .......... .......... 58% 67.5M 1s
##  78950K .......... .......... .......... .......... .......... 58% 55.4M 1s
##  79000K .......... .......... .......... .......... .......... 58% 83.5M 1s
##  79050K .......... .......... .......... .......... .......... 59% 20.8M 1s
##  79100K .......... .......... .......... .......... .......... 59% 66.4M 1s
##  79150K .......... .......... .......... .......... .......... 59% 57.2M 1s
##  79200K .......... .......... .......... .......... .......... 59% 76.9M 1s
##  79250K .......... .......... .......... .......... .......... 59% 42.8M 1s
##  79300K .......... .......... .......... .......... .......... 59% 69.9M 1s
##  79350K .......... .......... .......... .......... .......... 59% 69.5M 1s
##  79400K .......... .......... .......... .......... .......... 59% 72.2M 1s
##  79450K .......... .......... .......... .......... .......... 59% 22.5M 1s
##  79500K .......... .......... .......... .......... .......... 59% 24.5M 1s
##  79550K .......... .......... .......... .......... .......... 59% 62.1M 1s
##  79600K .......... .......... .......... .......... .......... 59% 51.9M 1s
##  79650K .......... .......... .......... .......... .......... 59% 66.0M 1s
##  79700K .......... .......... .......... .......... .......... 59% 50.9M 1s
##  79750K .......... .......... .......... .......... .......... 59% 59.0M 1s
##  79800K .......... .......... .......... .......... .......... 59% 26.6M 1s
##  79850K .......... .......... .......... .......... .......... 59% 73.6M 1s
##  79900K .......... .......... .......... .......... .......... 59% 85.2M 1s
##  79950K .......... .......... .......... .......... .......... 59% 20.9M 1s
##  80000K .......... .......... .......... .......... .......... 59% 40.6M 1s
##  80050K .......... .......... .......... .......... .......... 59% 89.9M 1s
##  80100K .......... .......... .......... .......... .......... 59% 18.8M 1s
##  80150K .......... .......... .......... .......... .......... 59% 29.1M 1s
##  80200K .......... .......... .......... .......... .......... 59% 73.1M 1s
##  80250K .......... .......... .......... .......... .......... 59% 90.2M 1s
##  80300K .......... .......... .......... .......... .......... 59% 81.5M 1s
##  80350K .......... .......... .......... .......... .......... 59% 16.4M 1s
##  80400K .......... .......... .......... .......... .......... 60% 84.7M 1s
##  80450K .......... .......... .......... .......... .......... 60% 62.2M 1s
##  80500K .......... .......... .......... .......... .......... 60% 21.1M 1s
##  80550K .......... .......... .......... .......... .......... 60% 27.7M 1s
##  80600K .......... .......... .......... .......... .......... 60% 61.5M 1s
##  80650K .......... .......... .......... .......... .......... 60% 82.9M 1s
##  80700K .......... .......... .......... .......... .......... 60% 63.7M 1s
##  80750K .......... .......... .......... .......... .......... 60% 26.7M 1s
##  80800K .......... .......... .......... .......... .......... 60% 78.3M 1s
##  80850K .......... .......... .......... .......... .......... 60% 85.0M 1s
##  80900K .......... .......... .......... .......... .......... 60% 63.0M 1s
##  80950K .......... .......... .......... .......... .......... 60% 21.0M 1s
##  81000K .......... .......... .......... .......... .......... 60% 68.9M 1s
##  81050K .......... .......... .......... .......... .......... 60% 80.6M 1s
##  81100K .......... .......... .......... .......... .......... 60% 94.0M 1s
##  81150K .......... .......... .......... .......... .......... 60% 28.5M 1s
##  81200K .......... .......... .......... .......... .......... 60% 86.6M 1s
##  81250K .......... .......... .......... .......... .......... 60% 11.7M 1s
##  81300K .......... .......... .......... .......... .......... 60% 32.0M 1s
##  81350K .......... .......... .......... .......... .......... 60% 78.6M 1s
##  81400K .......... .......... .......... .......... .......... 60% 88.5M 1s
##  81450K .......... .......... .......... .......... .......... 60% 63.8M 1s
##  81500K .......... .......... .......... .......... .......... 60% 91.2M 1s
##  81550K .......... .......... .......... .......... .......... 60% 25.4M 1s
##  81600K .......... .......... .......... .......... .......... 60% 89.2M 1s
##  81650K .......... .......... .......... .......... .......... 60% 68.6M 1s
##  81700K .......... .......... .......... .......... .......... 60% 19.9M 1s
##  81750K .......... .......... .......... .......... .......... 61% 32.4M 1s
##  81800K .......... .......... .......... .......... .......... 61% 87.3M 1s
##  81850K .......... .......... .......... .......... .......... 61% 90.9M 1s
##  81900K .......... .......... .......... .......... .......... 61% 24.0M 1s
##  81950K .......... .......... .......... .......... .......... 61% 43.3M 1s
##  82000K .......... .......... .......... .......... .......... 61% 95.3M 1s
##  82050K .......... .......... .......... .......... .......... 61% 94.5M 1s
##  82100K .......... .......... .......... .......... .......... 61% 20.2M 1s
##  82150K .......... .......... .......... .......... .......... 61% 46.9M 1s
##  82200K .......... .......... .......... .......... .......... 61% 94.3M 1s
##  82250K .......... .......... .......... .......... .......... 61% 93.2M 1s
##  82300K .......... .......... .......... .......... .......... 61% 22.7M 1s
##  82350K .......... .......... .......... .......... .......... 61% 49.6M 1s
##  82400K .......... .......... .......... .......... .......... 61% 92.3M 1s
##  82450K .......... .......... .......... .......... .......... 61% 33.3M 1s
##  82500K .......... .......... .......... .......... .......... 61% 46.8M 1s
##  82550K .......... .......... .......... .......... .......... 61% 50.3M 1s
##  82600K .......... .......... .......... .......... .......... 61%  101M 1s
##  82650K .......... .......... .......... .......... .......... 61% 36.5M 1s
##  82700K .......... .......... .......... .......... .......... 61% 31.8M 1s
##  82750K .......... .......... .......... .......... .......... 61% 73.2M 1s
##  82800K .......... .......... .......... .......... .......... 61% 96.7M 1s
##  82850K .......... .......... .......... .......... .......... 61% 40.5M 1s
##  82900K .......... .......... .......... .......... .......... 61% 54.3M 1s
##  82950K .......... .......... .......... .......... .......... 61% 60.0M 1s
##  83000K .......... .......... .......... .......... .......... 61% 94.7M 1s
##  83050K .......... .......... .......... .......... .......... 61% 23.2M 1s
##  83100K .......... .......... .......... .......... .......... 62% 45.9M 1s
##  83150K .......... .......... .......... .......... .......... 62% 78.3M 1s
##  83200K .......... .......... .......... .......... .......... 62% 92.1M 1s
##  83250K .......... .......... .......... .......... .......... 62% 26.5M 1s
##  83300K .......... .......... .......... .......... .......... 62% 47.7M 1s
##  83350K .......... .......... .......... .......... .......... 62% 83.3M 1s
##  83400K .......... .......... .......... .......... .......... 62% 93.8M 1s
##  83450K .......... .......... .......... .......... .......... 62% 60.2M 1s
##  83500K .......... .......... .......... .......... .......... 62% 62.4M 1s
##  83550K .......... .......... .......... .......... .......... 62% 41.6M 1s
##  83600K .......... .......... .......... .......... .......... 62% 39.9M 1s
##  83650K .......... .......... .......... .......... .......... 62% 77.8M 1s
##  83700K .......... .......... .......... .......... .......... 62% 51.4M 1s
##  83750K .......... .......... .......... .......... .......... 62% 61.8M 1s
##  83800K .......... .......... .......... .......... .......... 62% 30.6M 1s
##  83850K .......... .......... .......... .......... .......... 62% 93.6M 1s
##  83900K .......... .......... .......... .......... .......... 62% 31.2M 1s
##  83950K .......... .......... .......... .......... .......... 62% 48.2M 1s
##  84000K .......... .......... .......... .......... .......... 62% 71.3M 1s
##  84050K .......... .......... .......... .......... .......... 62% 40.1M 1s
##  84100K .......... .......... .......... .......... .......... 62% 85.0M 1s
##  84150K .......... .......... .......... .......... .......... 62% 39.6M 1s
##  84200K .......... .......... .......... .......... .......... 62% 58.2M 1s
##  84250K .......... .......... .......... .......... .......... 62% 39.2M 1s
##  84300K .......... .......... .......... .......... .......... 62% 89.6M 1s
##  84350K .......... .......... .......... .......... .......... 62% 30.4M 1s
##  84400K .......... .......... .......... .......... .......... 62% 85.5M 1s
##  84450K .......... .......... .......... .......... .......... 63% 32.6M 1s
##  84500K .......... .......... .......... .......... .......... 63% 32.0M 1s
##  84550K .......... .......... .......... .......... .......... 63% 84.8M 1s
##  84600K .......... .......... .......... .......... .......... 63% 89.7M 1s
##  84650K .......... .......... .......... .......... .......... 63% 52.2M 1s
##  84700K .......... .......... .......... .......... .......... 63% 26.3M 1s
##  84750K .......... .......... .......... .......... .......... 63% 71.8M 1s
##  84800K .......... .......... .......... .......... .......... 63% 80.2M 1s
##  84850K .......... .......... .......... .......... .......... 63%  105M 1s
##  84900K .......... .......... .......... .......... .......... 63% 32.2M 1s
##  84950K .......... .......... .......... .......... .......... 63% 48.1M 1s
##  85000K .......... .......... .......... .......... .......... 63% 76.6M 1s
##  85050K .......... .......... .......... .......... .......... 63% 98.5M 1s
##  85100K .......... .......... .......... .......... .......... 63% 12.4M 1s
##  85150K .......... .......... .......... .......... .......... 63% 84.4M 1s
##  85200K .......... .......... .......... .......... .......... 63% 99.4M 1s
##  85250K .......... .......... .......... .......... .......... 63% 42.8M 1s
##  85300K .......... .......... .......... .......... .......... 63% 25.1M 1s
##  85350K .......... .......... .......... .......... .......... 63% 77.7M 1s
##  85400K .......... .......... .......... .......... .......... 63% 65.2M 1s
##  85450K .......... .......... .......... .......... .......... 63% 84.9M 1s
##  85500K .......... .......... .......... .......... .......... 63% 94.9M 1s
##  85550K .......... .......... .......... .......... .......... 63% 18.4M 1s
##  85600K .......... .......... .......... .......... .......... 63% 99.1M 1s
##  85650K .......... .......... .......... .......... .......... 63% 99.4M 1s
##  85700K .......... .......... .......... .......... .......... 63% 18.4M 1s
##  85750K .......... .......... .......... .......... .......... 63% 50.4M 1s
##  85800K .......... .......... .......... .......... .......... 64% 96.4M 1s
##  85850K .......... .......... .......... .......... .......... 64% 96.7M 1s
##  85900K .......... .......... .......... .......... .......... 64% 28.6M 1s
##  85950K .......... .......... .......... .......... .......... 64% 39.5M 1s
##  86000K .......... .......... .......... .......... .......... 64% 97.1M 1s
##  86050K .......... .......... .......... .......... .......... 64% 90.4M 1s
##  86100K .......... .......... .......... .......... .......... 64% 30.5M 1s
##  86150K .......... .......... .......... .......... .......... 64% 43.5M 1s
##  86200K .......... .......... .......... .......... .......... 64% 60.4M 1s
##  86250K .......... .......... .......... .......... .......... 64%  101M 1s
##  86300K .......... .......... .......... .......... .......... 64% 38.4M 1s
##  86350K .......... .......... .......... .......... .......... 64% 50.5M 1s
##  86400K .......... .......... .......... .......... .......... 64% 57.2M 1s
##  86450K .......... .......... .......... .......... .......... 64% 14.8M 1s
##  86500K .......... .......... .......... .......... .......... 64% 88.3M 1s
##  86550K .......... .......... .......... .......... .......... 64% 84.0M 1s
##  86600K .......... .......... .......... .......... .......... 64% 97.4M 1s
##  86650K .......... .......... .......... .......... .......... 64% 19.6M 1s
##  86700K .......... .......... .......... .......... .......... 64% 36.4M 1s
##  86750K .......... .......... .......... .......... .......... 64% 53.9M 1s
##  86800K .......... .......... .......... .......... .......... 64% 79.4M 1s
##  86850K .......... .......... .......... .......... .......... 64% 99.1M 1s
##  86900K .......... .......... .......... .......... .......... 64% 31.5M 1s
##  86950K .......... .......... .......... .......... .......... 64% 41.8M 1s
##  87000K .......... .......... .......... .......... .......... 64% 94.2M 1s
##  87050K .......... .......... .......... .......... .......... 64%  107M 1s
##  87100K .......... .......... .......... .......... .......... 65% 37.2M 1s
##  87150K .......... .......... .......... .......... .......... 65% 30.9M 1s
##  87200K .......... .......... .......... .......... .......... 65% 80.2M 1s
##  87250K .......... .......... .......... .......... .......... 65% 81.7M 1s
##  87300K .......... .......... .......... .......... .......... 65% 60.5M 1s
##  87350K .......... .......... .......... .......... .......... 65% 29.2M 1s
##  87400K .......... .......... .......... .......... .......... 65% 40.8M 1s
##  87450K .......... .......... .......... .......... .......... 65% 51.2M 1s
##  87500K .......... .......... .......... .......... .......... 65% 80.5M 1s
##  87550K .......... .......... .......... .......... .......... 65% 43.9M 1s
##  87600K .......... .......... .......... .......... .......... 65% 31.4M 1s
##  87650K .......... .......... .......... .......... .......... 65% 38.5M 1s
##  87700K .......... .......... .......... .......... .......... 65% 81.9M 1s
##  87750K .......... .......... .......... .......... .......... 65% 57.2M 1s
##  87800K .......... .......... .......... .......... .......... 65% 19.6M 1s
##  87850K .......... .......... .......... .......... .......... 65% 47.4M 1s
##  87900K .......... .......... .......... .......... .......... 65% 86.4M 1s
##  87950K .......... .......... .......... .......... .......... 65% 75.4M 1s
##  88000K .......... .......... .......... .......... .......... 65% 86.6M 1s
##  88050K .......... .......... .......... .......... .......... 65% 22.6M 1s
##  88100K .......... .......... .......... .......... .......... 65% 45.1M 1s
##  88150K .......... .......... .......... .......... .......... 65% 73.4M 1s
##  88200K .......... .......... .......... .......... .......... 65% 90.8M 1s
##  88250K .......... .......... .......... .......... .......... 65% 67.5M 1s
##  88300K .......... .......... .......... .......... .......... 65% 69.0M 1s
##  88350K .......... .......... .......... .......... .......... 65% 47.2M 1s
##  88400K .......... .......... .......... .......... .......... 65% 92.7M 1s
##  88450K .......... .......... .......... .......... .......... 66% 44.9M 1s
##  88500K .......... .......... .......... .......... .......... 66% 45.5M 1s
##  88550K .......... .......... .......... .......... .......... 66% 43.6M 1s
##  88600K .......... .......... .......... .......... .......... 66% 94.2M 1s
##  88650K .......... .......... .......... .......... .......... 66% 38.2M 1s
##  88700K .......... .......... .......... .......... .......... 66% 72.4M 1s
##  88750K .......... .......... .......... .......... .......... 66% 43.8M 1s
##  88800K .......... .......... .......... .......... .......... 66% 37.2M 1s
##  88850K .......... .......... .......... .......... .......... 66% 87.1M 1s
##  88900K .......... .......... .......... .......... .......... 66% 24.6M 1s
##  88950K .......... .......... .......... .......... .......... 66% 82.9M 1s
##  89000K .......... .......... .......... .......... .......... 66% 66.7M 1s
##  89050K .......... .......... .......... .......... .......... 66%  102M 1s
##  89100K .......... .......... .......... .......... .......... 66% 17.1M 1s
##  89150K .......... .......... .......... .......... .......... 66% 68.4M 1s
##  89200K .......... .......... .......... .......... .......... 66%  106M 1s
##  89250K .......... .......... .......... .......... .......... 66% 16.3M 1s
##  89300K .......... .......... .......... .......... .......... 66% 60.5M 1s
##  89350K .......... .......... .......... .......... .......... 66% 40.5M 1s
##  89400K .......... .......... .......... .......... .......... 66% 85.8M 1s
##  89450K .......... .......... .......... .......... .......... 66% 78.9M 1s
##  89500K .......... .......... .......... .......... .......... 66%  100M 1s
##  89550K .......... .......... .......... .......... .......... 66% 47.0M 1s
##  89600K .......... .......... .......... .......... .......... 66% 30.0M 1s
##  89650K .......... .......... .......... .......... .......... 66% 72.3M 1s
##  89700K .......... .......... .......... .......... .......... 66% 21.9M 1s
##  89750K .......... .......... .......... .......... .......... 66% 42.6M 1s
##  89800K .......... .......... .......... .......... .......... 67% 87.5M 1s
##  89850K .......... .......... .......... .......... .......... 67% 99.9M 1s
##  89900K .......... .......... .......... .......... .......... 67% 39.6M 1s
##  89950K .......... .......... .......... .......... .......... 67% 54.4M 1s
##  90000K .......... .......... .......... .......... .......... 67% 62.8M 1s
##  90050K .......... .......... .......... .......... .......... 67% 42.8M 1s
##  90100K .......... .......... .......... .......... .......... 67% 99.5M 1s
##  90150K .......... .......... .......... .......... .......... 67% 55.6M 1s
##  90200K .......... .......... .......... .......... .......... 67% 99.0M 1s
##  90250K .......... .......... .......... .......... .......... 67% 56.2M 1s
##  90300K .......... .......... .......... .......... .......... 67% 38.9M 1s
##  90350K .......... .......... .......... .......... .......... 67% 47.5M 1s
##  90400K .......... .......... .......... .......... .......... 67% 54.6M 1s
##  90450K .......... .......... .......... .......... .......... 67% 96.2M 1s
##  90500K .......... .......... .......... .......... .......... 67% 24.6M 1s
##  90550K .......... .......... .......... .......... .......... 67% 86.4M 1s
##  90600K .......... .......... .......... .......... .......... 67% 82.8M 1s
##  90650K .......... .......... .......... .......... .......... 67%  100M 1s
##  90700K .......... .......... .......... .......... .......... 67% 27.2M 1s
##  90750K .......... .......... .......... .......... .......... 67% 59.7M 1s
##  90800K .......... .......... .......... .......... .......... 67% 80.8M 1s
##  90850K .......... .......... .......... .......... .......... 67% 34.1M 1s
##  90900K .......... .......... .......... .......... .......... 67% 44.1M 1s
##  90950K .......... .......... .......... .......... .......... 67% 52.8M 1s
##  91000K .......... .......... .......... .......... .......... 67% 62.1M 1s
##  91050K .......... .......... .......... .......... .......... 67% 76.8M 1s
##  91100K .......... .......... .......... .......... .......... 67% 22.4M 1s
##  91150K .......... .......... .......... .......... .......... 68% 78.1M 1s
##  91200K .......... .......... .......... .......... .......... 68% 72.5M 1s
##  91250K .......... .......... .......... .......... .......... 68%  102M 1s
##  91300K .......... .......... .......... .......... .......... 68% 33.3M 1s
##  91350K .......... .......... .......... .......... .......... 68% 38.8M 1s
##  91400K .......... .......... .......... .......... .......... 68% 96.9M 1s
##  91450K .......... .......... .......... .......... .......... 68%  102M 1s
##  91500K .......... .......... .......... .......... .......... 68% 34.9M 1s
##  91550K .......... .......... .......... .......... .......... 68% 30.9M 1s
##  91600K .......... .......... .......... .......... .......... 68% 93.7M 1s
##  91650K .......... .......... .......... .......... .......... 68% 96.9M 1s
##  91700K .......... .......... .......... .......... .......... 68% 30.9M 1s
##  91750K .......... .......... .......... .......... .......... 68% 80.0M 1s
##  91800K .......... .......... .......... .......... .......... 68% 78.8M 1s
##  91850K .......... .......... .......... .......... .......... 68% 95.5M 1s
##  91900K .......... .......... .......... .......... .......... 68% 18.6M 1s
##  91950K .......... .......... .......... .......... .......... 68% 49.4M 1s
##  92000K .......... .......... .......... .......... .......... 68% 62.3M 1s
##  92050K .......... .......... .......... .......... .......... 68%  110M 1s
##  92100K .......... .......... .......... .......... .......... 68%  112M 1s
##  92150K .......... .......... .......... .......... .......... 68% 56.4M 1s
##  92200K .......... .......... .......... .......... .......... 68% 32.6M 1s
##  92250K .......... .......... .......... .......... .......... 68% 77.0M 1s
##  92300K .......... .......... .......... .......... .......... 68%  111M 1s
##  92350K .......... .......... .......... .......... .......... 68% 93.7M 1s
##  92400K .......... .......... .......... .......... .......... 68% 26.8M 1s
##  92450K .......... .......... .......... .......... .......... 68% 39.4M 1s
##  92500K .......... .......... .......... .......... .......... 69% 98.3M 1s
##  92550K .......... .......... .......... .......... .......... 69% 98.7M 1s
##  92600K .......... .......... .......... .......... .......... 69% 89.0M 1s
##  92650K .......... .......... .......... .......... .......... 69% 58.1M 1s
##  92700K .......... .......... .......... .......... .......... 69% 34.3M 1s
##  92750K .......... .......... .......... .......... .......... 69% 34.9M 1s
##  92800K .......... .......... .......... .......... .......... 69% 54.3M 1s
##  92850K .......... .......... .......... .......... .......... 69% 99.8M 1s
##  92900K .......... .......... .......... .......... .......... 69%  113M 1s
##  92950K .......... .......... .......... .......... .......... 69% 75.8M 1s
##  93000K .......... .......... .......... .......... .......... 69% 28.6M 1s
##  93050K .......... .......... .......... .......... .......... 69% 58.9M 1s
##  93100K .......... .......... .......... .......... .......... 69% 63.8M 1s
##  93150K .......... .......... .......... .......... .......... 69% 87.9M 1s
##  93200K .......... .......... .......... .......... .......... 69%  114M 1s
##  93250K .......... .......... .......... .......... .......... 69% 42.4M 1s
##  93300K .......... .......... .......... .......... .......... 69% 49.6M 1s
##  93350K .......... .......... .......... .......... .......... 69% 42.0M 1s
##  93400K .......... .......... .......... .......... .......... 69%  101M 1s
##  93450K .......... .......... .......... .......... .......... 69%  112M 1s
##  93500K .......... .......... .......... .......... .......... 69%  116M 1s
##  93550K .......... .......... .......... .......... .......... 69% 32.7M 1s
##  93600K .......... .......... .......... .......... .......... 69% 64.7M 1s
##  93650K .......... .......... .......... .......... .......... 69% 72.5M 1s
##  93700K .......... .......... .......... .......... .......... 69%  111M 1s
##  93750K .......... .......... .......... .......... .......... 69%  100M 1s
##  93800K .......... .......... .......... .......... .......... 70% 39.2M 1s
##  93850K .......... .......... .......... .......... .......... 70% 45.2M 1s
##  93900K .......... .......... .......... .......... .......... 70% 63.8M 1s
##  93950K .......... .......... .......... .......... .......... 70% 94.1M 1s
##  94000K .......... .......... .......... .......... .......... 70%  115M 1s
##  94050K .......... .......... .......... .......... .......... 70% 57.9M 1s
##  94100K .......... .......... .......... .......... .......... 70% 76.9M 1s
##  94150K .......... .......... .......... .......... .......... 70% 62.7M 1s
##  94200K .......... .......... .......... .......... .......... 70% 85.9M 1s
##  94250K .......... .......... .......... .......... .......... 70%  116M 1s
##  94300K .......... .......... .......... .......... .......... 70% 51.0M 1s
##  94350K .......... .......... .......... .......... .......... 70% 57.7M 1s
##  94400K .......... .......... .......... .......... .......... 70% 66.3M 1s
##  94450K .......... .......... .......... .......... .......... 70% 67.9M 1s
##  94500K .......... .......... .......... .......... .......... 70%  102M 1s
##  94550K .......... .......... .......... .......... .......... 70% 20.3M 1s
##  94600K .......... .......... .......... .......... .......... 70% 68.2M 1s
##  94650K .......... .......... .......... .......... .......... 70%  117M 1s
##  94700K .......... .......... .......... .......... .......... 70% 93.3M 1s
##  94750K .......... .......... .......... .......... .......... 70%  103M 1s
##  94800K .......... .......... .......... .......... .......... 70% 20.9M 1s
##  94850K .......... .......... .......... .......... .......... 70% 74.5M 1s
##  94900K .......... .......... .......... .......... .......... 70% 41.7M 1s
##  94950K .......... .......... .......... .......... .......... 70% 97.0M 1s
##  95000K .......... .......... .......... .......... .......... 70%  123M 1s
##  95050K .......... .......... .......... .......... .......... 70% 38.9M 1s
##  95100K .......... .......... .......... .......... .......... 70% 28.9M 1s
##  95150K .......... .......... .......... .......... .......... 71% 64.9M 1s
##  95200K .......... .......... .......... .......... .......... 71% 52.7M 1s
##  95250K .......... .......... .......... .......... .......... 71% 55.1M 1s
##  95300K .......... .......... .......... .......... .......... 71% 92.9M 1s
##  95350K .......... .......... .......... .......... .......... 71% 61.0M 1s
##  95400K .......... .......... .......... .......... .......... 71%  122M 1s
##  95450K .......... .......... .......... .......... .......... 71% 25.7M 1s
##  95500K .......... .......... .......... .......... .......... 71% 63.7M 1s
##  95550K .......... .......... .......... .......... .......... 71% 52.3M 1s
##  95600K .......... .......... .......... .......... .......... 71%  114M 1s
##  95650K .......... .......... .......... .......... .......... 71%  118M 1s
##  95700K .......... .......... .......... .......... .......... 71% 53.4M 1s
##  95750K .......... .......... .......... .......... .......... 71% 56.7M 1s
##  95800K .......... .......... .......... .......... .......... 71% 65.9M 1s
##  95850K .......... .......... .......... .......... .......... 71% 64.6M 1s
##  95900K .......... .......... .......... .......... .......... 71%  117M 1s
##  95950K .......... .......... .......... .......... .......... 71% 55.5M 1s
##  96000K .......... .......... .......... .......... .......... 71% 13.2M 1s
##  96050K .......... .......... .......... .......... .......... 71%  115M 1s
##  96100K .......... .......... .......... .......... .......... 71%  120M 1s
##  96150K .......... .......... .......... .......... .......... 71% 97.7M 1s
##  96200K .......... .......... .......... .......... .......... 71%  118M 1s
##  96250K .......... .......... .......... .......... .......... 71% 22.2M 1s
##  96300K .......... .......... .......... .......... .......... 71% 95.4M 1s
##  96350K .......... .......... .......... .......... .......... 71% 36.1M 1s
##  96400K .......... .......... .......... .......... .......... 71%  120M 1s
##  96450K .......... .......... .......... .......... .......... 71%  113M 1s
##  96500K .......... .......... .......... .......... .......... 72% 86.2M 1s
##  96550K .......... .......... .......... .......... .......... 72% 84.1M 1s
##  96600K .......... .......... .......... .......... .......... 72% 16.7M 1s
##  96650K .......... .......... .......... .......... .......... 72% 41.9M 1s
##  96700K .......... .......... .......... .......... .......... 72% 64.9M 1s
##  96750K .......... .......... .......... .......... .......... 72% 94.7M 1s
##  96800K .......... .......... .......... .......... .......... 72%  115M 1s
##  96850K .......... .......... .......... .......... .......... 72% 95.6M 1s
##  96900K .......... .......... .......... .......... .......... 72% 21.7M 1s
##  96950K .......... .......... .......... .......... .......... 72% 48.2M 1s
##  97000K .......... .......... .......... .......... .......... 72% 62.7M 1s
##  97050K .......... .......... .......... .......... .......... 72%  104M 1s
##  97100K .......... .......... .......... .......... .......... 72%  125M 1s
##  97150K .......... .......... .......... .......... .......... 72% 38.9M 1s
##  97200K .......... .......... .......... .......... .......... 72% 36.1M 1s
##  97250K .......... .......... .......... .......... .......... 72% 87.3M 1s
##  97300K .......... .......... .......... .......... .......... 72% 50.4M 1s
##  97350K .......... .......... .......... .......... .......... 72% 65.2M 1s
##  97400K .......... .......... .......... .......... .......... 72%  127M 1s
##  97450K .......... .......... .......... .......... .......... 72% 42.6M 1s
##  97500K .......... .......... .......... .......... .......... 72%  120M 1s
##  97550K .......... .......... .......... .......... .......... 72% 20.9M 1s
##  97600K .......... .......... .......... .......... .......... 72%  121M 1s
##  97650K .......... .......... .......... .......... .......... 72% 61.0M 1s
##  97700K .......... .......... .......... .......... .......... 72%  114M 1s
##  97750K .......... .......... .......... .......... .......... 72%  112M 1s
##  97800K .......... .......... .......... .......... .......... 72% 31.3M 1s
##  97850K .......... .......... .......... .......... .......... 73% 98.7M 1s
##  97900K .......... .......... .......... .......... .......... 73% 69.7M 1s
##  97950K .......... .......... .......... .......... .......... 73% 35.0M 1s
##  98000K .......... .......... .......... .......... .......... 73% 81.8M 1s
##  98050K .......... .......... .......... .......... .......... 73% 49.0M 1s
##  98100K .......... .......... .......... .......... .......... 73% 87.7M 1s
##  98150K .......... .......... .......... .......... .......... 73% 38.8M 1s
##  98200K .......... .......... .......... .......... .......... 73% 54.0M 1s
##  98250K .......... .......... .......... .......... .......... 73% 64.3M 1s
##  98300K .......... .......... .......... .......... .......... 73%  109M 1s
##  98350K .......... .......... .......... .......... .......... 73% 59.1M 1s
##  98400K .......... .......... .......... .......... .......... 73% 70.9M 1s
##  98450K .......... .......... .......... .......... .......... 73%  114M 1s
##  98500K .......... .......... .......... .......... .......... 73% 68.9M 1s
##  98550K .......... .......... .......... .......... .......... 73% 54.6M 1s
##  98600K .......... .......... .......... .......... .......... 73% 18.8M 1s
##  98650K .......... .......... .......... .......... .......... 73% 72.8M 1s
##  98700K .......... .......... .......... .......... .......... 73%  108M 1s
##  98750K .......... .......... .......... .......... .......... 73% 55.0M 1s
##  98800K .......... .......... .......... .......... .......... 73%  129M 1s
##  98850K .......... .......... .......... .......... .......... 73% 42.6M 1s
##  98900K .......... .......... .......... .......... .......... 73%  120M 1s
##  98950K .......... .......... .......... .......... .......... 73% 16.3M 1s
##  99000K .......... .......... .......... .......... .......... 73%  119M 1s
##  99050K .......... .......... .......... .......... .......... 73%  135M 1s
##  99100K .......... .......... .......... .......... .......... 73%  105M 1s
##  99150K .......... .......... .......... .......... .......... 73%  112M 1s
##  99200K .......... .......... .......... .......... .......... 74% 23.3M 1s
##  99250K .......... .......... .......... .......... .......... 74% 46.6M 1s
##  99300K .......... .......... .......... .......... .......... 74% 78.7M 1s
##  99350K .......... .......... .......... .......... .......... 74% 62.5M 1s
##  99400K .......... .......... .......... .......... .......... 74%  116M 1s
##  99450K .......... .......... .......... .......... .......... 74% 57.6M 1s
##  99500K .......... .......... .......... .......... .......... 74% 23.9M 1s
##  99550K .......... .......... .......... .......... .......... 74% 78.2M 1s
##  99600K .......... .......... .......... .......... .......... 74% 67.5M 1s
##  99650K .......... .......... .......... .......... .......... 74% 69.6M 1s
##  99700K .......... .......... .......... .......... .......... 74% 98.2M 1s
##  99750K .......... .......... .......... .......... .......... 74% 75.2M 1s
##  99800K .......... .......... .......... .......... .......... 74% 74.6M 1s
##  99850K .......... .......... .......... .......... .......... 74%  121M 1s
##  99900K .......... .......... .......... .......... .......... 74% 25.8M 1s
##  99950K .......... .......... .......... .......... .......... 74%  106M 1s
## 100000K .......... .......... .......... .......... .......... 74% 80.5M 1s
## 100050K .......... .......... .......... .......... .......... 74% 97.4M 1s
## 100100K .......... .......... .......... .......... .......... 74%  135M 1s
## 100150K .......... .......... .......... .......... .......... 74% 13.5M 1s
## 100200K .......... .......... .......... .......... .......... 74% 58.8M 1s
## 100250K .......... .......... .......... .......... .......... 74% 49.4M 1s
## 100300K .......... .......... .......... .......... .......... 74%  128M 1s
## 100350K .......... .......... .......... .......... .......... 74%  113M 1s
## 100400K .......... .......... .......... .......... .......... 74% 50.2M 1s
## 100450K .......... .......... .......... .......... .......... 74% 72.7M 1s
## 100500K .......... .......... .......... .......... .......... 75% 75.7M 1s
## 100550K .......... .......... .......... .......... .......... 75% 72.7M 1s
## 100600K .......... .......... .......... .......... .......... 75%  124M 1s
## 100650K .......... .......... .......... .......... .......... 75% 40.9M 1s
## 100700K .......... .......... .......... .......... .......... 75% 34.0M 1s
## 100750K .......... .......... .......... .......... .......... 75%  112M 1s
## 100800K .......... .......... .......... .......... .......... 75%  101M 1s
## 100850K .......... .......... .......... .......... .......... 75%  136M 1s
## 100900K .......... .......... .......... .......... .......... 75% 45.3M 1s
## 100950K .......... .......... .......... .......... .......... 75% 81.4M 1s
## 101000K .......... .......... .......... .......... .......... 75% 21.5M 1s
## 101050K .......... .......... .......... .......... .......... 75% 63.9M 1s
## 101100K .......... .......... .......... .......... .......... 75% 75.3M 1s
## 101150K .......... .......... .......... .......... .......... 75%  113M 1s
## 101200K .......... .......... .......... .......... .......... 75% 25.9M 1s
## 101250K .......... .......... .......... .......... .......... 75% 65.3M 1s
## 101300K .......... .......... .......... .......... .......... 75% 86.8M 1s
## 101350K .......... .......... .......... .......... .......... 75% 64.4M 1s
## 101400K .......... .......... .......... .......... .......... 75%  132M 1s
## 101450K .......... .......... .......... .......... .......... 75% 44.3M 1s
## 101500K .......... .......... .......... .......... .......... 75% 83.0M 1s
## 101550K .......... .......... .......... .......... .......... 75% 68.2M 1s
## 101600K .......... .......... .......... .......... .......... 75%  131M 1s
## 101650K .......... .......... .......... .......... .......... 75% 70.6M 1s
## 101700K .......... .......... .......... .......... .......... 75% 15.5M 1s
## 101750K .......... .......... .......... .......... .......... 75% 74.7M 1s
## 101800K .......... .......... .......... .......... .......... 75%  131M 1s
## 101850K .......... .......... .......... .......... .......... 76%  141M 1s
## 101900K .......... .......... .......... .......... .......... 76%  125M 1s
## 101950K .......... .......... .......... .......... .......... 76% 29.6M 1s
## 102000K .......... .......... .......... .......... .......... 76% 57.0M 1s
## 102050K .......... .......... .......... .......... .......... 76% 59.2M 1s
## 102100K .......... .......... .......... .......... .......... 76% 66.1M 1s
## 102150K .......... .......... .......... .......... .......... 76%  117M 1s
## 102200K .......... .......... .......... .......... .......... 76%  107M 1s
## 102250K .......... .......... .......... .......... .......... 76%  138M 1s
## 102300K .......... .......... .......... .......... .......... 76% 34.5M 1s
## 102350K .......... .......... .......... .......... .......... 76% 37.1M 1s
## 102400K .......... .......... .......... .......... .......... 76% 70.5M 1s
## 102450K .......... .......... .......... .......... .......... 76% 67.4M 1s
## 102500K .......... .......... .......... .......... .......... 76% 70.6M 1s
## 102550K .......... .......... .......... .......... .......... 76%  107M 1s
## 102600K .......... .......... .......... .......... .......... 76% 27.7M 1s
## 102650K .......... .......... .......... .......... .......... 76%  110M 1s
## 102700K .......... .......... .......... .......... .......... 76%  132M 1s
## 102750K .......... .......... .......... .......... .......... 76% 31.5M 1s
## 102800K .......... .......... .......... .......... .......... 76%  131M 1s
## 102850K .......... .......... .......... .......... .......... 76% 51.1M 1s
## 102900K .......... .......... .......... .......... .......... 76% 77.6M 1s
## 102950K .......... .......... .......... .......... .......... 76% 95.0M 1s
## 103000K .......... .......... .......... .......... .......... 76% 30.7M 1s
## 103050K .......... .......... .......... .......... .......... 76%  134M 1s
## 103100K .......... .......... .......... .......... .......... 76%  112M 1s
## 103150K .......... .......... .......... .......... .......... 76% 33.9M 1s
## 103200K .......... .......... .......... .......... .......... 77% 51.7M 1s
## 103250K .......... .......... .......... .......... .......... 77% 73.1M 1s
## 103300K .......... .......... .......... .......... .......... 77% 66.6M 1s
## 103350K .......... .......... .......... .......... .......... 77% 37.7M 1s
## 103400K .......... .......... .......... .......... .......... 77% 49.6M 1s
## 103450K .......... .......... .......... .......... .......... 77% 86.8M 1s
## 103500K .......... .......... .......... .......... .......... 77% 79.3M 1s
## 103550K .......... .......... .......... .......... .......... 77% 30.2M 1s
## 103600K .......... .......... .......... .......... .......... 77% 81.6M 1s
## 103650K .......... .......... .......... .......... .......... 77% 59.1M 1s
## 103700K .......... .......... .......... .......... .......... 77%  100M 1s
## 103750K .......... .......... .......... .......... .......... 77% 48.8M 1s
## 103800K .......... .......... .......... .......... .......... 77% 25.5M 1s
## 103850K .......... .......... .......... .......... .......... 77% 53.0M 1s
## 103900K .......... .......... .......... .......... .......... 77% 31.9M 1s
## 103950K .......... .......... .......... .......... .......... 77% 54.6M 1s
## 104000K .......... .......... .......... .......... .......... 77% 90.4M 1s
## 104050K .......... .......... .......... .......... .......... 77% 63.5M 1s
## 104100K .......... .......... .......... .......... .......... 77% 75.1M 1s
## 104150K .......... .......... .......... .......... .......... 77% 81.5M 1s
## 104200K .......... .......... .......... .......... .......... 77% 30.3M 1s
## 104250K .......... .......... .......... .......... .......... 77% 45.5M 1s
## 104300K .......... .......... .......... .......... .......... 77% 71.0M 1s
## 104350K .......... .......... .......... .......... .......... 77% 81.7M 1s
## 104400K .......... .......... .......... .......... .......... 77% 27.6M 1s
## 104450K .......... .......... .......... .......... .......... 77% 67.4M 1s
## 104500K .......... .......... .......... .......... .......... 77% 32.9M 1s
## 104550K .......... .......... .......... .......... .......... 78% 61.4M 1s
## 104600K .......... .......... .......... .......... .......... 78% 82.2M 1s
## 104650K .......... .......... .......... .......... .......... 78% 79.2M 1s
## 104700K .......... .......... .......... .......... .......... 78% 14.2M 1s
## 104750K .......... .......... .......... .......... .......... 78% 24.1M 1s
## 104800K .......... .......... .......... .......... .......... 78% 75.5M 1s
## 104850K .......... .......... .......... .......... .......... 78% 62.4M 1s
## 104900K .......... .......... .......... .......... .......... 78% 80.3M 1s
## 104950K .......... .......... .......... .......... .......... 78% 69.5M 1s
## 105000K .......... .......... .......... .......... .......... 78% 74.5M 1s
## 105050K .......... .......... .......... .......... .......... 78% 17.7M 1s
## 105100K .......... .......... .......... .......... .......... 78% 67.8M 1s
## 105150K .......... .......... .......... .......... .......... 78% 68.0M 1s
## 105200K .......... .......... .......... .......... .......... 78% 48.5M 1s
## 105250K .......... .......... .......... .......... .......... 78% 72.2M 1s
## 105300K .......... .......... .......... .......... .......... 78% 90.5M 1s
## 105350K .......... .......... .......... .......... .......... 78% 73.5M 1s
## 105400K .......... .......... .......... .......... .......... 78% 43.0M 1s
## 105450K .......... .......... .......... .......... .......... 78% 26.4M 1s
## 105500K .......... .......... .......... .......... .......... 78% 62.3M 1s
## 105550K .......... .......... .......... .......... .......... 78% 24.9M 1s
## 105600K .......... .......... .......... .......... .......... 78% 39.3M 1s
## 105650K .......... .......... .......... .......... .......... 78% 37.3M 1s
## 105700K .......... .......... .......... .......... .......... 78% 73.7M 1s
## 105750K .......... .......... .......... .......... .......... 78% 52.9M 1s
## 105800K .......... .......... .......... .......... .......... 78% 46.0M 1s
## 105850K .......... .......... .......... .......... .......... 78% 48.9M 1s
## 105900K .......... .......... .......... .......... .......... 79% 46.3M 1s
## 105950K .......... .......... .......... .......... .......... 79% 33.8M 1s
## 106000K .......... .......... .......... .......... .......... 79% 79.8M 1s
## 106050K .......... .......... .......... .......... .......... 79% 69.9M 1s
## 106100K .......... .......... .......... .......... .......... 79% 47.4M 1s
## 106150K .......... .......... .......... .......... .......... 79% 62.1M 1s
## 106200K .......... .......... .......... .......... .......... 79% 69.1M 1s
## 106250K .......... .......... .......... .......... .......... 79% 52.7M 1s
## 106300K .......... .......... .......... .......... .......... 79% 48.4M 1s
## 106350K .......... .......... .......... .......... .......... 79% 43.0M 1s
## 106400K .......... .......... .......... .......... .......... 79% 79.6M 1s
## 106450K .......... .......... .......... .......... .......... 79% 91.2M 1s
## 106500K .......... .......... .......... .......... .......... 79% 59.6M 1s
## 106550K .......... .......... .......... .......... .......... 79% 60.9M 1s
## 106600K .......... .......... .......... .......... .......... 79% 67.6M 1s
## 106650K .......... .......... .......... .......... .......... 79% 86.3M 1s
## 106700K .......... .......... .......... .......... .......... 79% 87.3M 1s
## 106750K .......... .......... .......... .......... .......... 79% 77.6M 1s
## 106800K .......... .......... .......... .......... .......... 79% 94.4M 1s
## 106850K .......... .......... .......... .......... .......... 79% 75.4M 1s
## 106900K .......... .......... .......... .......... .......... 79% 64.9M 1s
## 106950K .......... .......... .......... .......... .......... 79% 35.1M 1s
## 107000K .......... .......... .......... .......... .......... 79% 91.4M 1s
## 107050K .......... .......... .......... .......... .......... 79% 91.1M 1s
## 107100K .......... .......... .......... .......... .......... 79% 66.7M 1s
## 107150K .......... .......... .......... .......... .......... 79% 52.2M 1s
## 107200K .......... .......... .......... .......... .......... 79% 56.4M 1s
## 107250K .......... .......... .......... .......... .......... 80% 46.2M 1s
## 107300K .......... .......... .......... .......... .......... 80% 46.2M 1s
## 107350K .......... .......... .......... .......... .......... 80% 68.1M 1s
## 107400K .......... .......... .......... .......... .......... 80% 75.0M 1s
## 107450K .......... .......... .......... .......... .......... 80% 56.8M 1s
## 107500K .......... .......... .......... .......... .......... 80% 69.0M 1s
## 107550K .......... .......... .......... .......... .......... 80% 45.2M 1s
## 107600K .......... .......... .......... .......... .......... 80% 43.0M 0s
## 107650K .......... .......... .......... .......... .......... 80% 61.7M 0s
## 107700K .......... .......... .......... .......... .......... 80% 39.6M 0s
## 107750K .......... .......... .......... .......... .......... 80% 85.6M 0s
## 107800K .......... .......... .......... .......... .......... 80% 72.3M 0s
## 107850K .......... .......... .......... .......... .......... 80% 88.2M 0s
## 107900K .......... .......... .......... .......... .......... 80% 61.9M 0s
## 107950K .......... .......... .......... .......... .......... 80% 26.1M 0s
## 108000K .......... .......... .......... .......... .......... 80% 86.7M 0s
## 108050K .......... .......... .......... .......... .......... 80% 52.7M 0s
## 108100K .......... .......... .......... .......... .......... 80% 98.2M 0s
## 108150K .......... .......... .......... .......... .......... 80% 84.0M 0s
## 108200K .......... .......... .......... .......... .......... 80%  106M 0s
## 108250K .......... .......... .......... .......... .......... 80% 78.0M 0s
## 108300K .......... .......... .......... .......... .......... 80% 74.6M 0s
## 108350K .......... .......... .......... .......... .......... 80% 84.0M 0s
## 108400K .......... .......... .......... .......... .......... 80% 54.6M 0s
## 108450K .......... .......... .......... .......... .......... 80% 96.7M 0s
## 108500K .......... .......... .......... .......... .......... 80% 73.8M 0s
## 108550K .......... .......... .......... .......... .......... 81% 68.2M 0s
## 108600K .......... .......... .......... .......... .......... 81% 48.8M 0s
## 108650K .......... .......... .......... .......... .......... 81%  106M 0s
## 108700K .......... .......... .......... .......... .......... 81% 72.6M 0s
## 108750K .......... .......... .......... .......... .......... 81% 80.8M 0s
## 108800K .......... .......... .......... .......... .......... 81%  104M 0s
## 108850K .......... .......... .......... .......... .......... 81%  103M 0s
## 108900K .......... .......... .......... .......... .......... 81%  102M 0s
## 108950K .......... .......... .......... .......... .......... 81% 71.2M 0s
## 109000K .......... .......... .......... .......... .......... 81%  111M 0s
## 109050K .......... .......... .......... .......... .......... 81% 67.7M 0s
## 109100K .......... .......... .......... .......... .......... 81%  105M 0s
## 109150K .......... .......... .......... .......... .......... 81% 90.9M 0s
## 109200K .......... .......... .......... .......... .......... 81% 89.9M 0s
## 109250K .......... .......... .......... .......... .......... 81%  104M 0s
## 109300K .......... .......... .......... .......... .......... 81%  114M 0s
## 109350K .......... .......... .......... .......... .......... 81% 30.1M 0s
## 109400K .......... .......... .......... .......... .......... 81% 60.1M 0s
## 109450K .......... .......... .......... .......... .......... 81% 29.8M 0s
## 109500K .......... .......... .......... .......... .......... 81% 62.0M 0s
## 109550K .......... .......... .......... .......... .......... 81% 89.4M 0s
## 109600K .......... .......... .......... .......... .......... 81%  107M 0s
## 109650K .......... .......... .......... .......... .......... 81% 83.2M 0s
## 109700K .......... .......... .......... .......... .......... 81%  106M 0s
## 109750K .......... .......... .......... .......... .......... 81% 29.2M 0s
## 109800K .......... .......... .......... .......... .......... 81% 39.6M 0s
## 109850K .......... .......... .......... .......... .......... 81% 60.9M 0s
## 109900K .......... .......... .......... .......... .......... 82% 45.3M 0s
## 109950K .......... .......... .......... .......... .......... 82% 61.8M 0s
## 110000K .......... .......... .......... .......... .......... 82% 58.9M 0s
## 110050K .......... .......... .......... .......... .......... 82% 41.4M 0s
## 110100K .......... .......... .......... .......... .......... 82% 46.0M 0s
## 110150K .......... .......... .......... .......... .......... 82% 72.3M 0s
## 110200K .......... .......... .......... .......... .......... 82% 81.1M 0s
## 110250K .......... .......... .......... .......... .......... 82% 30.8M 0s
## 110300K .......... .......... .......... .......... .......... 82% 49.9M 0s
## 110350K .......... .......... .......... .......... .......... 82% 74.3M 0s
## 110400K .......... .......... .......... .......... .......... 82% 65.4M 0s
## 110450K .......... .......... .......... .......... .......... 82% 85.6M 0s
## 110500K .......... .......... .......... .......... .......... 82% 53.7M 0s
## 110550K .......... .......... .......... .......... .......... 82% 55.4M 0s
## 110600K .......... .......... .......... .......... .......... 82% 86.0M 0s
## 110650K .......... .......... .......... .......... .......... 82% 77.4M 0s
## 110700K .......... .......... .......... .......... .......... 82% 53.0M 0s
## 110750K .......... .......... .......... .......... .......... 82% 46.4M 0s
## 110800K .......... .......... .......... .......... .......... 82% 35.2M 0s
## 110850K .......... .......... .......... .......... .......... 82% 86.5M 0s
## 110900K .......... .......... .......... .......... .......... 82% 53.3M 0s
## 110950K .......... .......... .......... .......... .......... 82% 73.8M 0s
## 111000K .......... .......... .......... .......... .......... 82% 89.0M 0s
## 111050K .......... .......... .......... .......... .......... 82% 98.3M 0s
## 111100K .......... .......... .......... .......... .......... 82% 60.0M 0s
## 111150K .......... .......... .......... .......... .......... 82% 63.7M 0s
## 111200K .......... .......... .......... .......... .......... 82% 90.2M 0s
## 111250K .......... .......... .......... .......... .......... 83% 63.8M 0s
## 111300K .......... .......... .......... .......... .......... 83% 87.0M 0s
## 111350K .......... .......... .......... .......... .......... 83% 80.3M 0s
## 111400K .......... .......... .......... .......... .......... 83% 15.4M 0s
## 111450K .......... .......... .......... .......... .......... 83% 66.3M 0s
## 111500K .......... .......... .......... .......... .......... 83% 94.7M 0s
## 111550K .......... .......... .......... .......... .......... 83% 78.2M 0s
## 111600K .......... .......... .......... .......... .......... 83% 49.2M 0s
## 111650K .......... .......... .......... .......... .......... 83% 71.4M 0s
## 111700K .......... .......... .......... .......... .......... 83% 61.1M 0s
## 111750K .......... .......... .......... .......... .......... 83% 58.6M 0s
## 111800K .......... .......... .......... .......... .......... 83% 77.9M 0s
## 111850K .......... .......... .......... .......... .......... 83% 97.3M 0s
## 111900K .......... .......... .......... .......... .......... 83% 83.8M 0s
## 111950K .......... .......... .......... .......... .......... 83% 81.1M 0s
## 112000K .......... .......... .......... .......... .......... 83% 35.6M 0s
## 112050K .......... .......... .......... .......... .......... 83% 51.3M 0s
## 112100K .......... .......... .......... .......... .......... 83% 95.6M 0s
## 112150K .......... .......... .......... .......... .......... 83% 66.1M 0s
## 112200K .......... .......... .......... .......... .......... 83% 91.3M 0s
## 112250K .......... .......... .......... .......... .......... 83% 25.1M 0s
## 112300K .......... .......... .......... .......... .......... 83% 46.8M 0s
## 112350K .......... .......... .......... .......... .......... 83% 73.9M 0s
## 112400K .......... .......... .......... .......... .......... 83% 99.5M 0s
## 112450K .......... .......... .......... .......... .......... 83% 88.3M 0s
## 112500K .......... .......... .......... .......... .......... 83% 77.1M 0s
## 112550K .......... .......... .......... .......... .......... 83% 61.5M 0s
## 112600K .......... .......... .......... .......... .......... 84% 49.3M 0s
## 112650K .......... .......... .......... .......... .......... 84% 58.6M 0s
## 112700K .......... .......... .......... .......... .......... 84% 91.8M 0s
## 112750K .......... .......... .......... .......... .......... 84% 69.5M 0s
## 112800K .......... .......... .......... .......... .......... 84% 94.3M 0s
## 112850K .......... .......... .......... .......... .......... 84% 26.4M 0s
## 112900K .......... .......... .......... .......... .......... 84% 88.3M 0s
## 112950K .......... .......... .......... .......... .......... 84% 77.0M 0s
## 113000K .......... .......... .......... .......... .......... 84% 79.8M 0s
## 113050K .......... .......... .......... .......... .......... 84% 75.5M 0s
## 113100K .......... .......... .......... .......... .......... 84% 27.7M 0s
## 113150K .......... .......... .......... .......... .......... 84% 70.8M 0s
## 113200K .......... .......... .......... .......... .......... 84% 76.0M 0s
## 113250K .......... .......... .......... .......... .......... 84% 97.3M 0s
## 113300K .......... .......... .......... .......... .......... 84% 95.3M 0s
## 113350K .......... .......... .......... .......... .......... 84% 45.6M 0s
## 113400K .......... .......... .......... .......... .......... 84% 60.4M 0s
## 113450K .......... .......... .......... .......... .......... 84% 55.4M 0s
## 113500K .......... .......... .......... .......... .......... 84% 99.9M 0s
## 113550K .......... .......... .......... .......... .......... 84% 81.9M 0s
## 113600K .......... .......... .......... .......... .......... 84% 48.0M 0s
## 113650K .......... .......... .......... .......... .......... 84% 61.5M 0s
## 113700K .......... .......... .......... .......... .......... 84% 95.4M 0s
## 113750K .......... .......... .......... .......... .......... 84% 87.9M 0s
## 113800K .......... .......... .......... .......... .......... 84% 19.1M 0s
## 113850K .......... .......... .......... .......... .......... 84% 38.2M 0s
## 113900K .......... .......... .......... .......... .......... 84% 91.0M 0s
## 113950K .......... .......... .......... .......... .......... 85% 83.3M 0s
## 114000K .......... .......... .......... .......... .......... 85%  104M 0s
## 114050K .......... .......... .......... .......... .......... 85% 11.3M 0s
## 114100K .......... .......... .......... .......... .......... 85% 67.9M 0s
## 114150K .......... .......... .......... .......... .......... 85% 55.4M 0s
## 114200K .......... .......... .......... .......... .......... 85% 90.6M 0s
## 114250K .......... .......... .......... .......... .......... 85%  107M 0s
## 114300K .......... .......... .......... .......... .......... 85% 56.5M 0s
## 114350K .......... .......... .......... .......... .......... 85% 5.57M 0s
## 114400K .......... .......... .......... .......... .......... 85% 44.8M 0s
## 114450K .......... .......... .......... .......... .......... 85%  101M 0s
## 114500K .......... .......... .......... .......... .......... 85%  111M 0s
## 114550K .......... .......... .......... .......... .......... 85% 95.3M 0s
## 114600K .......... .......... .......... .......... .......... 85% 34.4M 0s
## 114650K .......... .......... .......... .......... .......... 85% 69.8M 0s
## 114700K .......... .......... .......... .......... .......... 85% 44.5M 0s
## 114750K .......... .......... .......... .......... .......... 85% 89.5M 0s
## 114800K .......... .......... .......... .......... .......... 85% 60.6M 0s
## 114850K .......... .......... .......... .......... .......... 85% 66.4M 0s
## 114900K .......... .......... .......... .......... .......... 85% 78.4M 0s
## 114950K .......... .......... .......... .......... .......... 85% 61.2M 0s
## 115000K .......... .......... .......... .......... .......... 85% 96.0M 0s
## 115050K .......... .......... .......... .......... .......... 85%  104M 0s
## 115100K .......... .......... .......... .......... .......... 85% 75.2M 0s
## 115150K .......... .......... .......... .......... .......... 85% 89.9M 0s
## 115200K .......... .......... .......... .......... .......... 85%  109M 0s
## 115250K .......... .......... .......... .......... .......... 86% 64.9M 0s
## 115300K .......... .......... .......... .......... .......... 86%  102M 0s
## 115350K .......... .......... .......... .......... .......... 86% 94.0M 0s
## 115400K .......... .......... .......... .......... .......... 86% 23.5M 0s
## 115450K .......... .......... .......... .......... .......... 86% 75.7M 0s
## 115500K .......... .......... .......... .......... .......... 86% 86.4M 0s
## 115550K .......... .......... .......... .......... .......... 86% 37.2M 0s
## 115600K .......... .......... .......... .......... .......... 86% 67.3M 0s
## 115650K .......... .......... .......... .......... .......... 86% 63.7M 0s
## 115700K .......... .......... .......... .......... .......... 86%  111M 0s
## 115750K .......... .......... .......... .......... .......... 86%  101M 0s
## 115800K .......... .......... .......... .......... .......... 86% 87.6M 0s
## 115850K .......... .......... .......... .......... .......... 86%  106M 0s
## 115900K .......... .......... .......... .......... .......... 86% 60.8M 0s
## 115950K .......... .......... .......... .......... .......... 86% 54.3M 0s
## 116000K .......... .......... .......... .......... .......... 86% 59.9M 0s
## 116050K .......... .......... .......... .......... .......... 86% 94.9M 0s
## 116100K .......... .......... .......... .......... .......... 86%  105M 0s
## 116150K .......... .......... .......... .......... .......... 86% 75.7M 0s
## 116200K .......... .......... .......... .......... .......... 86%  101M 0s
## 116250K .......... .......... .......... .......... .......... 86% 60.1M 0s
## 116300K .......... .......... .......... .......... .......... 86% 72.6M 0s
## 116350K .......... .......... .......... .......... .......... 86% 47.9M 0s
## 116400K .......... .......... .......... .......... .......... 86% 58.1M 0s
## 116450K .......... .......... .......... .......... .......... 86% 66.3M 0s
## 116500K .......... .......... .......... .......... .......... 86% 67.0M 0s
## 116550K .......... .......... .......... .......... .......... 86% 93.3M 0s
## 116600K .......... .......... .......... .......... .......... 87% 97.2M 0s
## 116650K .......... .......... .......... .......... .......... 87% 60.7M 0s
## 116700K .......... .......... .......... .......... .......... 87% 64.8M 0s
## 116750K .......... .......... .......... .......... .......... 87% 49.0M 0s
## 116800K .......... .......... .......... .......... .......... 87%  111M 0s
## 116850K .......... .......... .......... .......... .......... 87%  117M 0s
## 116900K .......... .......... .......... .......... .......... 87%  112M 0s
## 116950K .......... .......... .......... .......... .......... 87% 73.2M 0s
## 117000K .......... .......... .......... .......... .......... 87% 69.5M 0s
## 117050K .......... .......... .......... .......... .......... 87% 96.4M 0s
## 117100K .......... .......... .......... .......... .......... 87% 39.3M 0s
## 117150K .......... .......... .......... .......... .......... 87% 95.9M 0s
## 117200K .......... .......... .......... .......... .......... 87%  114M 0s
## 117250K .......... .......... .......... .......... .......... 87% 95.0M 0s
## 117300K .......... .......... .......... .......... .......... 87%  112M 0s
## 117350K .......... .......... .......... .......... .......... 87% 26.2M 0s
## 117400K .......... .......... .......... .......... .......... 87% 80.8M 0s
## 117450K .......... .......... .......... .......... .......... 87% 73.4M 0s
## 117500K .......... .......... .......... .......... .......... 87% 48.2M 0s
## 117550K .......... .......... .......... .......... .......... 87% 56.9M 0s
## 117600K .......... .......... .......... .......... .......... 87%  119M 0s
## 117650K .......... .......... .......... .......... .......... 87% 89.7M 0s
## 117700K .......... .......... .......... .......... .......... 87%  118M 0s
## 117750K .......... .......... .......... .......... .......... 87% 61.3M 0s
## 117800K .......... .......... .......... .......... .......... 87% 85.4M 0s
## 117850K .......... .......... .......... .......... .......... 87% 58.3M 0s
## 117900K .......... .......... .......... .......... .......... 87% 59.4M 0s
## 117950K .......... .......... .......... .......... .......... 88% 86.1M 0s
## 118000K .......... .......... .......... .......... .......... 88% 78.9M 0s
## 118050K .......... .......... .......... .......... .......... 88% 66.4M 0s
## 118100K .......... .......... .......... .......... .......... 88% 60.8M 0s
## 118150K .......... .......... .......... .......... .......... 88% 55.8M 0s
## 118200K .......... .......... .......... .......... .......... 88%  108M 0s
## 118250K .......... .......... .......... .......... .......... 88%  115M 0s
## 118300K .......... .......... .......... .......... .......... 88%  119M 0s
## 118350K .......... .......... .......... .......... .......... 88% 49.1M 0s
## 118400K .......... .......... .......... .......... .......... 88%  109M 0s
## 118450K .......... .......... .......... .......... .......... 88% 36.1M 0s
## 118500K .......... .......... .......... .......... .......... 88% 86.6M 0s
## 118550K .......... .......... .......... .......... .......... 88% 75.2M 0s
## 118600K .......... .......... .......... .......... .......... 88% 68.8M 0s
## 118650K .......... .......... .......... .......... .......... 88% 82.0M 0s
## 118700K .......... .......... .......... .......... .......... 88% 42.0M 0s
## 118750K .......... .......... .......... .......... .......... 88% 45.5M 0s
## 118800K .......... .......... .......... .......... .......... 88% 60.4M 0s
## 118850K .......... .......... .......... .......... .......... 88% 77.5M 0s
## 118900K .......... .......... .......... .......... .......... 88% 67.3M 0s
## 118950K .......... .......... .......... .......... .......... 88% 63.6M 0s
## 119000K .......... .......... .......... .......... .......... 88% 72.5M 0s
## 119050K .......... .......... .......... .......... .......... 88% 64.4M 0s
## 119100K .......... .......... .......... .......... .......... 88% 68.6M 0s
## 119150K .......... .......... .......... .......... .......... 88% 66.3M 0s
## 119200K .......... .......... .......... .......... .......... 88% 79.6M 0s
## 119250K .......... .......... .......... .......... .......... 88% 80.0M 0s
## 119300K .......... .......... .......... .......... .......... 89% 74.5M 0s
## 119350K .......... .......... .......... .......... .......... 89% 64.1M 0s
## 119400K .......... .......... .......... .......... .......... 89% 63.8M 0s
## 119450K .......... .......... .......... .......... .......... 89% 86.8M 0s
## 119500K .......... .......... .......... .......... .......... 89% 84.5M 0s
## 119550K .......... .......... .......... .......... .......... 89% 66.4M 0s
## 119600K .......... .......... .......... .......... .......... 89% 77.9M 0s
## 119650K .......... .......... .......... .......... .......... 89% 73.1M 0s
## 119700K .......... .......... .......... .......... .......... 89% 76.4M 0s
## 119750K .......... .......... .......... .......... .......... 89% 68.0M 0s
## 119800K .......... .......... .......... .......... .......... 89% 66.3M 0s
## 119850K .......... .......... .......... .......... .......... 89% 79.3M 0s
## 119900K .......... .......... .......... .......... .......... 89% 84.3M 0s
## 119950K .......... .......... .......... .......... .......... 89% 68.5M 0s
## 120000K .......... .......... .......... .......... .......... 89% 78.1M 0s
## 120050K .......... .......... .......... .......... .......... 89% 85.5M 0s
## 120100K .......... .......... .......... .......... .......... 89% 83.1M 0s
## 120150K .......... .......... .......... .......... .......... 89% 77.6M 0s
## 120200K .......... .......... .......... .......... .......... 89% 90.3M 0s
## 120250K .......... .......... .......... .......... .......... 89% 86.7M 0s
## 120300K .......... .......... .......... .......... .......... 89% 80.2M 0s
## 120350K .......... .......... .......... .......... .......... 89% 25.9M 0s
## 120400K .......... .......... .......... .......... .......... 89% 64.5M 0s
## 120450K .......... .......... .......... .......... .......... 89% 75.0M 0s
## 120500K .......... .......... .......... .......... .......... 89% 88.0M 0s
## 120550K .......... .......... .......... .......... .......... 89% 77.4M 0s
## 120600K .......... .......... .......... .......... .......... 89% 94.6M 0s
## 120650K .......... .......... .......... .......... .......... 90% 38.8M 0s
## 120700K .......... .......... .......... .......... .......... 90% 84.2M 0s
## 120750K .......... .......... .......... .......... .......... 90% 66.6M 0s
## 120800K .......... .......... .......... .......... .......... 90% 81.9M 0s
## 120850K .......... .......... .......... .......... .......... 90% 88.6M 0s
## 120900K .......... .......... .......... .......... .......... 90% 89.5M 0s
## 120950K .......... .......... .......... .......... .......... 90% 67.1M 0s
## 121000K .......... .......... .......... .......... .......... 90% 40.4M 0s
## 121050K .......... .......... .......... .......... .......... 90% 71.3M 0s
## 121100K .......... .......... .......... .......... .......... 90% 78.5M 0s
## 121150K .......... .......... .......... .......... .......... 90% 82.1M 0s
## 121200K .......... .......... .......... .......... .......... 90% 97.1M 0s
## 121250K .......... .......... .......... .......... .......... 90% 96.6M 0s
## 121300K .......... .......... .......... .......... .......... 90% 51.5M 0s
## 121350K .......... .......... .......... .......... .......... 90% 68.5M 0s
## 121400K .......... .......... .......... .......... .......... 90% 94.4M 0s
## 121450K .......... .......... .......... .......... .......... 90% 95.0M 0s
## 121500K .......... .......... .......... .......... .......... 90% 91.1M 0s
## 121550K .......... .......... .......... .......... .......... 90% 80.8M 0s
## 121600K .......... .......... .......... .......... .......... 90% 80.5M 0s
## 121650K .......... .......... .......... .......... .......... 90% 91.3M 0s
## 121700K .......... .......... .......... .......... .......... 90% 90.2M 0s
## 121750K .......... .......... .......... .......... .......... 90% 25.3M 0s
## 121800K .......... .......... .......... .......... .......... 90% 87.5M 0s
## 121850K .......... .......... .......... .......... .......... 90% 86.2M 0s
## 121900K .......... .......... .......... .......... .......... 90% 79.2M 0s
## 121950K .......... .......... .......... .......... .......... 91% 83.9M 0s
## 122000K .......... .......... .......... .......... .......... 91% 75.2M 0s
## 122050K .......... .......... .......... .......... .......... 91% 76.9M 0s
## 122100K .......... .......... .......... .......... .......... 91% 83.6M 0s
## 122150K .......... .......... .......... .......... .......... 91% 82.5M 0s
## 122200K .......... .......... .......... .......... .......... 91% 97.7M 0s
## 122250K .......... .......... .......... .......... .......... 91% 86.7M 0s
## 122300K .......... .......... .......... .......... .......... 91%  103M 0s
## 122350K .......... .......... .......... .......... .......... 91% 52.9M 0s
## 122400K .......... .......... .......... .......... .......... 91% 39.8M 0s
## 122450K .......... .......... .......... .......... .......... 91% 89.3M 0s
## 122500K .......... .......... .......... .......... .......... 91% 79.8M 0s
## 122550K .......... .......... .......... .......... .......... 91% 90.0M 0s
## 122600K .......... .......... .......... .......... .......... 91%  106M 0s
## 122650K .......... .......... .......... .......... .......... 91% 97.6M 0s
## 122700K .......... .......... .......... .......... .......... 91%  105M 0s
## 122750K .......... .......... .......... .......... .......... 91% 19.3M 0s
## 122800K .......... .......... .......... .......... .......... 91%  103M 0s
## 122850K .......... .......... .......... .......... .......... 91% 96.7M 0s
## 122900K .......... .......... .......... .......... .......... 91%  103M 0s
## 122950K .......... .......... .......... .......... .......... 91% 88.9M 0s
## 123000K .......... .......... .......... .......... .......... 91% 35.0M 0s
## 123050K .......... .......... .......... .......... .......... 91% 78.3M 0s
## 123100K .......... .......... .......... .......... .......... 91% 93.3M 0s
## 123150K .......... .......... .......... .......... .......... 91% 60.1M 0s
## 123200K .......... .......... .......... .......... .......... 91% 92.3M 0s
## 123250K .......... .......... .......... .......... .......... 91%  101M 0s
## 123300K .......... .......... .......... .......... .......... 92% 89.9M 0s
## 123350K .......... .......... .......... .......... .......... 92% 91.4M 0s
## 123400K .......... .......... .......... .......... .......... 92% 62.4M 0s
## 123450K .......... .......... .......... .......... .......... 92% 75.8M 0s
## 123500K .......... .......... .......... .......... .......... 92% 44.4M 0s
## 123550K .......... .......... .......... .......... .......... 92% 47.9M 0s
## 123600K .......... .......... .......... .......... .......... 92% 67.2M 0s
## 123650K .......... .......... .......... .......... .......... 92% 37.0M 0s
## 123700K .......... .......... .......... .......... .......... 92% 66.0M 0s
## 123750K .......... .......... .......... .......... .......... 92% 45.5M 0s
## 123800K .......... .......... .......... .......... .......... 92% 62.6M 0s
## 123850K .......... .......... .......... .......... .......... 92%  105M 0s
## 123900K .......... .......... .......... .......... .......... 92%  106M 0s
## 123950K .......... .......... .......... .......... .......... 92% 74.3M 0s
## 124000K .......... .......... .......... .......... .......... 92% 98.9M 0s
## 124050K .......... .......... .......... .......... .......... 92%  107M 0s
## 124100K .......... .......... .......... .......... .......... 92%  103M 0s
## 124150K .......... .......... .......... .......... .......... 92% 93.1M 0s
## 124200K .......... .......... .......... .......... .......... 92%  102M 0s
## 124250K .......... .......... .......... .......... .......... 92% 22.6M 0s
## 124300K .......... .......... .......... .......... .......... 92% 94.8M 0s
## 124350K .......... .......... .......... .......... .......... 92% 84.0M 0s
## 124400K .......... .......... .......... .......... .......... 92%  105M 0s
## 124450K .......... .......... .......... .......... .......... 92%  110M 0s
## 124500K .......... .......... .......... .......... .......... 92%  105M 0s
## 124550K .......... .......... .......... .......... .......... 92% 93.0M 0s
## 124600K .......... .......... .......... .......... .......... 92% 98.8M 0s
## 124650K .......... .......... .......... .......... .......... 93% 42.5M 0s
## 124700K .......... .......... .......... .......... .......... 93% 89.3M 0s
## 124750K .......... .......... .......... .......... .......... 93% 80.2M 0s
## 124800K .......... .......... .......... .......... .......... 93%  114M 0s
## 124850K .......... .......... .......... .......... .......... 93%  115M 0s
## 124900K .......... .......... .......... .......... .......... 93%  113M 0s
## 124950K .......... .......... .......... .......... .......... 93%  104M 0s
## 125000K .......... .......... .......... .......... .......... 93% 71.1M 0s
## 125050K .......... .......... .......... .......... .......... 93% 27.1M 0s
## 125100K .......... .......... .......... .......... .......... 93%  112M 0s
## 125150K .......... .......... .......... .......... .......... 93% 91.0M 0s
## 125200K .......... .......... .......... .......... .......... 93%  116M 0s
## 125250K .......... .......... .......... .......... .......... 93%  104M 0s
## 125300K .......... .......... .......... .......... .......... 93%  109M 0s
## 125350K .......... .......... .......... .......... .......... 93% 36.0M 0s
## 125400K .......... .......... .......... .......... .......... 93% 24.6M 0s
## 125450K .......... .......... .......... .......... .......... 93% 73.5M 0s
## 125500K .......... .......... .......... .......... .......... 93% 56.2M 0s
## 125550K .......... .......... .......... .......... .......... 93% 64.6M 0s
## 125600K .......... .......... .......... .......... .......... 93% 76.1M 0s
## 125650K .......... .......... .......... .......... .......... 93% 71.4M 0s
## 125700K .......... .......... .......... .......... .......... 93% 63.1M 0s
## 125750K .......... .......... .......... .......... .......... 93% 63.9M 0s
## 125800K .......... .......... .......... .......... .......... 93% 68.2M 0s
## 125850K .......... .......... .......... .......... .......... 93% 72.9M 0s
## 125900K .......... .......... .......... .......... .......... 93% 69.4M 0s
## 125950K .......... .......... .......... .......... .......... 93% 53.3M 0s
## 126000K .......... .......... .......... .......... .......... 94% 76.6M 0s
## 126050K .......... .......... .......... .......... .......... 94% 78.2M 0s
## 126100K .......... .......... .......... .......... .......... 94% 62.4M 0s
## 126150K .......... .......... .......... .......... .......... 94% 55.4M 0s
## 126200K .......... .......... .......... .......... .......... 94% 74.4M 0s
## 126250K .......... .......... .......... .......... .......... 94% 74.8M 0s
## 126300K .......... .......... .......... .......... .......... 94% 80.5M 0s
## 126350K .......... .......... .......... .......... .......... 94% 61.6M 0s
## 126400K .......... .......... .......... .......... .......... 94% 79.3M 0s
## 126450K .......... .......... .......... .......... .......... 94% 39.7M 0s
## 126500K .......... .......... .......... .......... .......... 94% 80.2M 0s
## 126550K .......... .......... .......... .......... .......... 94% 69.9M 0s
## 126600K .......... .......... .......... .......... .......... 94% 70.9M 0s
## 126650K .......... .......... .......... .......... .......... 94% 63.5M 0s
## 126700K .......... .......... .......... .......... .......... 94% 67.9M 0s
## 126750K .......... .......... .......... .......... .......... 94% 64.7M 0s
## 126800K .......... .......... .......... .......... .......... 94% 69.5M 0s
## 126850K .......... .......... .......... .......... .......... 94% 78.2M 0s
## 126900K .......... .......... .......... .......... .......... 94% 81.9M 0s
## 126950K .......... .......... .......... .......... .......... 94% 69.0M 0s
## 127000K .......... .......... .......... .......... .......... 94% 87.7M 0s
## 127050K .......... .......... .......... .......... .......... 94% 87.2M 0s
## 127100K .......... .......... .......... .......... .......... 94% 86.2M 0s
## 127150K .......... .......... .......... .......... .......... 94% 64.4M 0s
## 127200K .......... .......... .......... .......... .......... 94% 74.0M 0s
## 127250K .......... .......... .......... .......... .......... 94% 76.2M 0s
## 127300K .......... .......... .......... .......... .......... 94% 87.4M 0s
## 127350K .......... .......... .......... .......... .......... 95% 72.8M 0s
## 127400K .......... .......... .......... .......... .......... 95% 91.4M 0s
## 127450K .......... .......... .......... .......... .......... 95% 69.1M 0s
## 127500K .......... .......... .......... .......... .......... 95% 61.3M 0s
## 127550K .......... .......... .......... .......... .......... 95% 73.7M 0s
## 127600K .......... .......... .......... .......... .......... 95% 77.6M 0s
## 127650K .......... .......... .......... .......... .......... 95% 76.2M 0s
## 127700K .......... .......... .......... .......... .......... 95% 94.1M 0s
## 127750K .......... .......... .......... .......... .......... 95% 53.0M 0s
## 127800K .......... .......... .......... .......... .......... 95% 71.5M 0s
## 127850K .......... .......... .......... .......... .......... 95% 92.2M 0s
## 127900K .......... .......... .......... .......... .......... 95% 91.3M 0s
## 127950K .......... .......... .......... .......... .......... 95% 78.4M 0s
## 128000K .......... .......... .......... .......... .......... 95% 18.7M 0s
## 128050K .......... .......... .......... .......... .......... 95% 64.2M 0s
## 128100K .......... .......... .......... .......... .......... 95% 84.1M 0s
## 128150K .......... .......... .......... .......... .......... 95% 75.6M 0s
## 128200K .......... .......... .......... .......... .......... 95% 86.4M 0s
## 128250K .......... .......... .......... .......... .......... 95% 30.5M 0s
## 128300K .......... .......... .......... .......... .......... 95% 67.1M 0s
## 128350K .......... .......... .......... .......... .......... 95% 65.3M 0s
## 128400K .......... .......... .......... .......... .......... 95% 43.2M 0s
## 128450K .......... .......... .......... .......... .......... 95% 66.8M 0s
## 128500K .......... .......... .......... .......... .......... 95% 90.8M 0s
## 128550K .......... .......... .......... .......... .......... 95% 64.6M 0s
## 128600K .......... .......... .......... .......... .......... 95% 87.3M 0s
## 128650K .......... .......... .......... .......... .......... 95% 74.6M 0s
## 128700K .......... .......... .......... .......... .......... 96% 86.7M 0s
## 128750K .......... .......... .......... .......... .......... 96% 74.4M 0s
## 128800K .......... .......... .......... .......... .......... 96% 91.5M 0s
## 128850K .......... .......... .......... .......... .......... 96% 67.8M 0s
## 128900K .......... .......... .......... .......... .......... 96% 64.0M 0s
## 128950K .......... .......... .......... .......... .......... 96% 52.9M 0s
## 129000K .......... .......... .......... .......... .......... 96% 61.4M 0s
## 129050K .......... .......... .......... .......... .......... 96% 82.7M 0s
## 129100K .......... .......... .......... .......... .......... 96% 93.1M 0s
## 129150K .......... .......... .......... .......... .......... 96% 66.7M 0s
## 129200K .......... .......... .......... .......... .......... 96% 75.4M 0s
## 129250K .......... .......... .......... .......... .......... 96% 91.0M 0s
## 129300K .......... .......... .......... .......... .......... 96% 54.8M 0s
## 129350K .......... .......... .......... .......... .......... 96% 65.0M 0s
## 129400K .......... .......... .......... .......... .......... 96% 85.5M 0s
## 129450K .......... .......... .......... .......... .......... 96% 70.8M 0s
## 129500K .......... .......... .......... .......... .......... 96% 77.5M 0s
## 129550K .......... .......... .......... .......... .......... 96% 72.7M 0s
## 129600K .......... .......... .......... .......... .......... 96% 46.7M 0s
## 129650K .......... .......... .......... .......... .......... 96% 70.5M 0s
## 129700K .......... .......... .......... .......... .......... 96% 73.3M 0s
## 129750K .......... .......... .......... .......... .......... 96% 75.4M 0s
## 129800K .......... .......... .......... .......... .......... 96% 85.1M 0s
## 129850K .......... .......... .......... .......... .......... 96% 93.3M 0s
## 129900K .......... .......... .......... .......... .......... 96% 53.5M 0s
## 129950K .......... .......... .......... .......... .......... 96% 78.1M 0s
## 130000K .......... .......... .......... .......... .......... 97% 20.0M 0s
## 130050K .......... .......... .......... .......... .......... 97% 89.9M 0s
## 130100K .......... .......... .......... .......... .......... 97% 85.9M 0s
## 130150K .......... .......... .......... .......... .......... 97% 71.0M 0s
## 130200K .......... .......... .......... .......... .......... 97% 93.8M 0s
## 130250K .......... .......... .......... .......... .......... 97% 99.6M 0s
## 130300K .......... .......... .......... .......... .......... 97% 16.2M 0s
## 130350K .......... .......... .......... .......... .......... 97% 79.0M 0s
## 130400K .......... .......... .......... .......... .......... 97%  101M 0s
## 130450K .......... .......... .......... .......... .......... 97% 95.2M 0s
## 130500K .......... .......... .......... .......... .......... 97%  105M 0s
## 130550K .......... .......... .......... .......... .......... 97% 75.8M 0s
## 130600K .......... .......... .......... .......... .......... 97% 78.6M 0s
## 130650K .......... .......... .......... .......... .......... 97%  100M 0s
## 130700K .......... .......... .......... .......... .......... 97% 30.5M 0s
## 130750K .......... .......... .......... .......... .......... 97% 60.0M 0s
## 130800K .......... .......... .......... .......... .......... 97% 82.6M 0s
## 130850K .......... .......... .......... .......... .......... 97% 77.7M 0s
## 130900K .......... .......... .......... .......... .......... 97%  102M 0s
## 130950K .......... .......... .......... .......... .......... 97% 83.3M 0s
## 131000K .......... .......... .......... .......... .......... 97%  103M 0s
## 131050K .......... .......... .......... .......... .......... 97% 71.4M 0s
## 131100K .......... .......... .......... .......... .......... 97% 64.8M 0s
## 131150K .......... .......... .......... .......... .......... 97% 64.7M 0s
## 131200K .......... .......... .......... .......... .......... 97% 77.3M 0s
## 131250K .......... .......... .......... .......... .......... 97% 83.1M 0s
## 131300K .......... .......... .......... .......... .......... 97% 48.3M 0s
## 131350K .......... .......... .......... .......... .......... 98% 66.6M 0s
## 131400K .......... .......... .......... .......... .......... 98% 73.5M 0s
## 131450K .......... .......... .......... .......... .......... 98% 46.3M 0s
## 131500K .......... .......... .......... .......... .......... 98% 99.7M 0s
## 131550K .......... .......... .......... .......... .......... 98% 42.1M 0s
## 131600K .......... .......... .......... .......... .......... 98% 67.9M 0s
## 131650K .......... .......... .......... .......... .......... 98% 88.5M 0s
## 131700K .......... .......... .......... .......... .......... 98% 74.6M 0s
## 131750K .......... .......... .......... .......... .......... 98% 67.8M 0s
## 131800K .......... .......... .......... .......... .......... 98% 74.4M 0s
## 131850K .......... .......... .......... .......... .......... 98% 80.8M 0s
## 131900K .......... .......... .......... .......... .......... 98% 91.7M 0s
## 131950K .......... .......... .......... .......... .......... 98% 64.3M 0s
## 132000K .......... .......... .......... .......... .......... 98% 85.6M 0s
## 132050K .......... .......... .......... .......... .......... 98% 75.5M 0s
## 132100K .......... .......... .......... .......... .......... 98% 81.2M 0s
## 132150K .......... .......... .......... .......... .......... 98% 80.1M 0s
## 132200K .......... .......... .......... .......... .......... 98% 94.1M 0s
## 132250K .......... .......... .......... .......... .......... 98% 85.5M 0s
## 132300K .......... .......... .......... .......... .......... 98% 56.2M 0s
## 132350K .......... .......... .......... .......... .......... 98% 26.2M 0s
## 132400K .......... .......... .......... .......... .......... 98% 91.6M 0s
## 132450K .......... .......... .......... .......... .......... 98%  102M 0s
## 132500K .......... .......... .......... .......... .......... 98%  101M 0s
## 132550K .......... .......... .......... .......... .......... 98% 87.6M 0s
## 132600K .......... .......... .......... .......... .......... 98% 56.0M 0s
## 132650K .......... .......... .......... .......... .......... 98% 92.2M 0s
## 132700K .......... .......... .......... .......... .......... 99% 27.4M 0s
## 132750K .......... .......... .......... .......... .......... 99% 63.1M 0s
## 132800K .......... .......... .......... .......... .......... 99% 88.6M 0s
## 132850K .......... .......... .......... .......... .......... 99% 97.4M 0s
## 132900K .......... .......... .......... .......... .......... 99%  106M 0s
## 132950K .......... .......... .......... .......... .......... 99% 92.2M 0s
## 133000K .......... .......... .......... .......... .......... 99% 58.7M 0s
## 133050K .......... .......... .......... .......... .......... 99% 61.4M 0s
## 133100K .......... .......... .......... .......... .......... 99% 65.3M 0s
## 133150K .......... .......... .......... .......... .......... 99% 66.3M 0s
## 133200K .......... .......... .......... .......... .......... 99% 86.9M 0s
## 133250K .......... .......... .......... .......... .......... 99% 99.3M 0s
## 133300K .......... .......... .......... .......... .......... 99% 98.8M 0s
## 133350K .......... .......... .......... .......... .......... 99% 20.6M 0s
## 133400K .......... .......... .......... .......... .......... 99% 25.5M 0s
## 133450K .......... .......... .......... .......... .......... 99% 68.2M 0s
## 133500K .......... .......... .......... .......... .......... 99% 70.6M 0s
## 133550K .......... .......... .......... .......... .......... 99% 68.1M 0s
## 133600K .......... .......... .......... .......... .......... 99% 83.5M 0s
## 133650K .......... .......... .......... .......... .......... 99%  109M 0s
## 133700K .......... .......... .......... .......... .......... 99% 99.5M 0s
## 133750K .......... .......... .......... .......... .......... 99% 81.2M 0s
## 133800K .......... .......... .......... .......... .......... 99% 80.1M 0s
## 133850K .......... .......... .......... .......... .......... 99% 71.7M 0s
## 133900K .......... .......... .......... .......... .......... 99% 66.5M 0s
## 133950K .......... .......... .......... .......... .......... 99% 63.4M 0s
## 134000K .......... .......... .......... .......... .......... 99%  103M 0s
## 134050K .......... .....                                      100%  103M=2.4s
## 
## 2021-12-24 17:40:09 (53.7 MB/s) - ‘silva_nr99_v138.1_train_set.fa.gz.3’ saved [137283333/137283333]
```

```r
seqtabNoC <- removeBimeraDenovo(seqtabAll)
```
##Assign taxonomy

```r
fastaRef <-"/home/rstudio/silva_nr99_v138.1_train_set.fa.gz"
taxTab<-assignTaxonomy(seqtabNoC, refFasta=fastaRef, multithread=TRUE)
unname(head(taxTab))
```

```
##      [,1]       [,2]           [,3]          [,4]            [,5]            
## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
##      [,6]         
## [1,] NA           
## [2,] NA           
## [3,] NA           
## [4,] NA           
## [5,] "Bacteroides"
## [6,] NA
```
##Construct phylogenetic tree

```r
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
```

```r
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
```
##Combine data into a phyloseq object

```r
samdf <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/MIMARKS_Data_combined.csv",header=TRUE)
samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21)
samdf <- samdf[!duplicated(samdf$SampleID),] # Remove dupicate entries for reverse reads
rownames(seqtabAll) <- gsub("124", "125", rownames(seqtabAll)) # Fix discrepancy
all(rownames(seqtabAll) %in% samdf$SampleID) # TRUE
```

```
## [1] TRUE
```

```r
rownames(samdf) <- samdf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment",
"host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
"diet", "family_relationship", "genotype", "SampleID") 
samdf <- samdf[rownames(seqtabAll), keep.cols]
```

```r
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxTab),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 218 taxa and 19 samples ]
## sample_data() Sample Data:       [ 19 samples by 14 sample variables ]
## tax_table()   Taxonomy Table:    [ 218 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 218 tips and 216 internal nodes ]
```
#Using phyloseq
##Loading the data

```r
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds")
ps = readRDS(ps_connect)
ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]
```
##Filtering
###Taxonomic Filtering

```r
rank_names(ps)
```

```
## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"
```

```r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

```r
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

```r
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

```
##                         Phylum         1     2
## 1               Actinobacteria 120.15385  1562
## 2                Bacteroidetes 265.52174  6107
## 3  Candidatus_Saccharibacteria 280.00000   280
## 4    Cyanobacteria/Chloroplast  64.25000   257
## 5          Deinococcus-Thermus  52.00000    52
## 6                   Firmicutes 179.24771 58614
## 7                 Fusobacteria   2.00000     2
## 8               Proteobacteria  59.09091   650
## 9                  Tenericutes 234.00000   234
## 10             Verrucomicrobia 104.00000   104
```

```r
# Define phyla to filter
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]
```
###Prevalence Filtering

```r
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

```r
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

```
## [1] 18
```

```r
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```
##Agglomerate taxa

```r
# How many genera would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
```

```
## [1] 49
```

```r
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
```

```r
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

```r
multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
```
##Abundance value transformation

```r
# group plots together
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

```r
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
```

```r
# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
```

```r
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->
##Subset by taxonomy

```r
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->
###Preprocessing

```r
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

```r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

```r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

```
## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
## GCAAGCGTTATCCGGAATGACTGGGCGTAAAGGGTGCGTAGGTGGTTTGGCAAGTTGGTAGCGTAATTCCGGGGCTCAACCTCGGCGCTACTACCAAAACTGCTGGACTTGAGTGCAGGAGGGGTGAATGGAATTCCTAGTGTAGCGGTGGAATGCGTAGATATTAGGAAGAACACCAGCGGCGAAGGCGATTCACTGGACTGTAACTGACACTGAGGCACGAAAGCGTGGGGAG
## -- in the phylogenetic tree in the data you provided.
```

```r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

```r
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->
#Different Ordination Projections

```r
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)
```

```r
which(!rowSums(otu_table(ps)) > 1000)
```

```
## F5D145 M1D149   M1D9 M2D125  M2D19 M3D148 M3D149   M3D3   M3D5   M3D8 
##     69    185    200    204    218    243    244    252    256    260
```

```r
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
```

```r
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

```r
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

```r
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

```r
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
```

```
## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
## GCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCTTGATAAGTCTGAAGTGAAAGGCCAAGGCTTAACCATGGAACTGCTTTGGAAACTATGAGGCTAGAGTGCTGGAGAGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTACTGGACAGAAACTGACGTTGAGGCTCGAAAGCGTGGGGAG
## -- in the phylogenetic tree in the data you provided.
```

```r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->
##Why are the ordination plots so far from square?
###PCA on ranks

```r
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
```

```r
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1
```

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:Biostrings':
## 
##     collapse, intersect, setdiff, setequal, union
```

```
## The following object is masked from 'package:GenomeInfoDb':
## 
##     intersect
```

```
## The following object is masked from 'package:XVector':
## 
##     slice
```

```
## The following objects are masked from 'package:IRanges':
## 
##     collapse, desc, intersect, setdiff, slice, union
```

```
## The following objects are masked from 'package:S4Vectors':
## 
##     first, intersect, rename, setdiff, setequal, union
```

```
## The following objects are masked from 'package:BiocGenerics':
## 
##     combine, intersect, setdiff, union
```

```
## The following object is masked from 'package:gridExtra':
## 
##     combine
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

```
## Joining, by = c("Var1", "Var2")
```

```r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

```
## Joining, by = c("Var1", "Var2")
```

```r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")
sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

```r
library(ade4)
```

```
## 
## Attaching package: 'ade4'
```

```
## The following object is masked from 'package:Biostrings':
## 
##     score
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     score
```

```r
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
```

```
## Joining, by = "SampleID"
```

```r
col_scores <- col_scores %>%
  left_join(tax)
```

```
## Joining, by = "seq"
```

```r
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->
###Canonical correspondence

```r
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)
```

```r
library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))
```

```
## Joining, by = "SampleID"
```

```r
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)
```

```
## Joining, by = "otu_id"
```

```r
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
            size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
        y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

```
## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps

## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

```r
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

```r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

```r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

```
## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
## GCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGTGGCATGGCAAGCCAGAAGTGAAAACCCGGGGCTTAACCCCGCGGATTGCTTTTGGAACTGTCAGGCTGGAGTGCAGGAGGGGCAGGCGGAATTCCTGGTGTAGCGGTGAAATGCGTAGATATCAGGAGGAACACCGGTGGCGAAGGCGGCCTGCTGGACTGTAACTGACACTGAGGCTCGAAAGCGTGGGGA
## -- in the phylogenetic tree in the data you provided.
```

```r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](cc_Richard_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->
