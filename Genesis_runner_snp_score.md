Genesis
================
Rutuja Gupte
2024-07-03

``` r
library(GWASTools)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, table,
    ##     tapply, union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

``` r
library(SNPRelate)
```

    ## Loading required package: gdsfmt

    ## SNPRelate -- supported by Streaming SIMD Extensions 2 (SSE2)

``` r
library(GENESIS)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::combine()    masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ✖ ggplot2::Position() masks BiocGenerics::Position(), base::Position()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
traits <- read.table("traits.txt", header=FALSE)
colnames(traits) <- c("sample.id", "pheno")
pheno <- traits$pheno
mydat <- data.frame(scanID = traits$sample.id,
                    pheno = pheno)
head(mydat)
```

    ##   scanID pheno
    ## 1  33-16 64.75
    ## 2  38-11 92.25
    ## 3   4226 65.50
    ## 4   4722 81.13
    ## 5   A188 27.50
    ## 6  A214N 65.00

``` r
scanAnnot <- ScanAnnotationDataFrame(mydat)
scanAnnot
```

    ## An object of class 'ScanAnnotationDataFrame'
    ##   scans: 1 2 ... 271 (271 total)
    ##   varLabels: scanID pheno
    ##   varMetadata: labelDescription

``` r
snpgdsVCF2GDS("processed.vcf.gz", "snpgds/file_score.gds")
```

    ## Start file conversion from VCF to SNP GDS ...
    ## Method: extracting biallelic SNPs
    ## Number of samples: 271
    ## Parsing "processed.vcf.gz" ...
    ##  import 2828 variants.
    ## + genotype   { Bit2 271x2828, 187.1K } *
    ## Optimize the access efficiency ...
    ## Clean up the fragments of GDS file:
    ##     open the file 'snpgds/file_score.gds' (208.1K)
    ##     # of fragments: 47
    ##     save to 'snpgds/file_score.gds.tmp'
    ##     rename 'snpgds/file_score.gds.tmp' (207.8K, reduced: 324B)
    ##     # of fragments: 20

``` r
geno <- GdsGenotypeReader(filename = "snpgds/file_score.gds")
genoData <- GenotypeData(geno)
```

``` r
nullmod <- fitNullModel(scanAnnot, outcome = "pheno", 
                        family = "gaussian")
```

``` r
genoIterator <- GenotypeBlockIterator(genoData, snpBlock=5000)
assoc <- assocTestSingle(genoIterator, null.model = nullmod,
                         BPPARAM = BiocParallel::SerialParam())
```

    ## Using 1 CPU cores

``` r
write.table(assoc, "genesis_snp_score.txt", sep="\t", row.names=FALSE, col.names=TRUE)
```

``` r
close(geno)
```
