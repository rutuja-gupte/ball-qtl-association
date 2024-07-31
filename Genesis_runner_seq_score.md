Genesis
================
Rutuja Gupte
2024-07-03

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(SeqArray)
```

    ## Loading required package: gdsfmt
    ## 
    ## Attaching package: 'SeqArray'
    ## 
    ## The following object is masked from 'package:stringr':
    ## 
    ##     fixed

``` r
library(Biobase)
```

    ## Loading required package: BiocGenerics
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, table,
    ##     tapply, union, unique, unsplit, which.max, which.min
    ## 
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

``` r
library(SeqVarTools)
library(GENESIS)
```

``` r
seqVCF2GDS("processed.vcf.gz", "seqgds/file_score.gds")
```

    ## Wed Jul 31 15:22:38 2024
    ## Variant Call Format (VCF) Import:
    ##     file:
    ##         processed.vcf.gz (207.1K)
    ##     file format: VCFv4.0
    ##     genome reference: <unknown>
    ##     # of sets of chromosomes (ploidy): 2
    ##     # of samples: 271
    ##     genotype field: GT
    ##     genotype storage: bit2
    ##     compression method: LZMA_RA
    ##     # of samples: 271
    ##     INFO: NS,DP,AF
    ##     FORMAT: AD,DP,GQ,PL
    ## Output:
    ##     seqgds/file_score.gds
    ##     [Progress Info: file_score.gds.progress]
    ## Parsing 'processed.vcf.gz':
    ## + genotype/data   { Bit2 2x271x2828 LZMA_ra(3.14%), 60.0K }
    ## Digests:
    ##     sample.id  [md5: 8fe7d22d1aac449aeb539c3a9c9de00d]
    ##     variant.id  [md5: cbad95ff2b3b24a58b8220b9b1b91e7d]
    ##     position  [md5: 5ad5486d85e9a6bf06ee23a1000eb8c7]
    ##     chromosome  [md5: 35cbdc679bf0faa4738495040cf55d88]
    ##     allele  [md5: 25ce460fcb43681e677f6368240bac88]
    ##     genotype  [md5: 96880de3abd46604929a60f4d9f1dc6d]
    ##     phase  [md5: 15b5fc8412f4967383b711df110b4ccc]
    ##     annotation/id  [md5: 0c4c4343b6c4251cf06d3b3ce0fd3b10]
    ##     annotation/qual  [md5: 7b12f5d3e16220588a82d5c91ff3249f]
    ##     annotation/filter  [md5: 0113515ab32a71180d580ecbb3577ee8]
    ##     annotation/info/NS  [md5: 226aacb5d2c62a8158059f288606ffe4]
    ##     annotation/info/DP  [md5: 226aacb5d2c62a8158059f288606ffe4]
    ##     annotation/info/AF  [md5: 6e6dcb39f4b68c5a7dc6ad670c6ff366]
    ##     annotation/format/AD  [md5: 6e6dcb39f4b68c5a7dc6ad670c6ff366]
    ##     annotation/format/DP  [md5: 6e6dcb39f4b68c5a7dc6ad670c6ff366]
    ##     annotation/format/GQ  [md5: 6e6dcb39f4b68c5a7dc6ad670c6ff366]
    ##     annotation/format/PL  [md5: 6e6dcb39f4b68c5a7dc6ad670c6ff366]
    ## Done.
    ## Wed Jul 31 15:22:38 2024
    ## Optimize the access efficiency ...
    ## Clean up the fragments of GDS file:
    ##     open the file 'seqgds/file_score.gds' (129.2K)
    ##     # of fragments: 159
    ##     save to 'seqgds/file_score.gds.tmp'
    ##     rename 'seqgds/file_score.gds.tmp' (128.2K, reduced: 984B)
    ##     # of fragments: 77
    ## Wed Jul 31 15:22:38 2024

``` r
gds <- seqOpen("seqgds/file_score.gds")
gds
```

    ## Object of class "SeqVarGDSClass"
    ## File: C:\Users\RGupte\OneDrive - Ball Horticultural Company\Association_Documentation\seqgds\file_score.gds (128.2K)
    ## +    [  ] *
    ## |--+ description   [  ] *
    ## |--+ sample.id   { Str8 271 LZMA_ra(54.7%), 829B } *
    ## |--+ variant.id   { Int32 2828 LZMA_ra(13.6%), 1.5K } *
    ## |--+ position   { Int32 2828 LZMA_ra(65.8%), 7.3K } *
    ## |--+ chromosome   { Str8 2828 LZMA_ra(2.36%), 145B } *
    ## |--+ allele   { Str8 2828 LZMA_ra(15.5%), 1.7K } *
    ## |--+ genotype   [  ] *
    ## |  |--+ data   { Bit2 2x271x2828 LZMA_ra(26.2%), 98.2K } *
    ## |  |--+ extra.index   { Int32 3x0 LZMA_ra, 18B } *
    ## |  \--+ extra   { Int16 0 LZMA_ra, 18B }
    ## |--+ phase   [  ]
    ## |  |--+ data   { Bit1 271x2828 LZMA_ra(0.17%), 169B } *
    ## |  |--+ extra.index   { Int32 3x0 LZMA_ra, 18B } *
    ## |  \--+ extra   { Bit1 0 LZMA_ra, 18B }
    ## |--+ annotation   [  ]
    ## |  |--+ id   { Str8 2828 LZMA_ra(23.2%), 7.0K } *
    ## |  |--+ qual   { Float32 2828 LZMA_ra(1.15%), 137B } *
    ## |  |--+ filter   { Int32,factor 2828 LZMA_ra(1.15%), 137B } *
    ## |  |--+ info   [  ]
    ## |  |  |--+ NS   { Int32 2828 LZMA_ra(1.15%), 137B } *
    ## |  |  |--+ DP   { Int32 2828 LZMA_ra(1.15%), 137B } *
    ## |  |  \--+ AF   { Float32 0 LZMA_ra, 18B } *
    ## |  \--+ format   [  ]
    ## |     |--+ AD   [  ] *
    ## |     |  \--+ data   { VL_Int 271x0 LZMA_ra, 18B } *
    ## |     |--+ DP   [  ] *
    ## |     |  \--+ data   { VL_Int 271x0 LZMA_ra, 18B } *
    ## |     |--+ GQ   [  ] *
    ## |     |  \--+ data   { Float32 271x0 LZMA_ra, 18B } *
    ## |     \--+ PL   [  ] *
    ## |        \--+ data   { Float32 271x0 LZMA_ra, 18B } *
    ## \--+ sample.annotation   [  ]

``` r
traits <- read.table("traits.txt", header=FALSE)
colnames(traits) <- c("sample.id", "pheno")

pheno <- traits$pheno
annot <- data.frame(sample.id = traits$sample.id,
                    outcome = pheno)
annot <- annot[match(seqGetData(gds, "sample.id"), annot$sample.id),]
row.names(annot) <- NULL
head(annot)
```

    ##   sample.id outcome
    ## 1     33-16   64.75
    ## 2     38-11   92.25
    ## 3      4226   65.50
    ## 4      4722   81.13
    ## 5      A188   27.50
    ## 6     A214N   65.00

``` r
metadata <- data.frame(labelDescription=c("sample id", 
                                          "phenotype"),
                       row.names=names(annot))
annot <- AnnotatedDataFrame(annot, metadata)
annot
```

    ## An object of class 'AnnotatedDataFrame'
    ##   rowNames: 1 2 ... 271 (271 total)
    ##   varLabels: sample.id outcome
    ##   varMetadata: labelDescription

``` r
all.equal(annot$sample.id, seqGetData(gds, "sample.id"))
```

    ## [1] TRUE

``` r
seqData <- SeqVarData(gds, sampleData=annot)
```

``` r
# fit the null model
nullmod <- fitNullModel(seqData, 
                        family="gaussian",
                        outcome="outcome")
```

``` r
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, nullmod, 
                         verbose=FALSE,
                         BPPARAM=BiocParallel::SerialParam())
```

``` r
write.table(assoc, "genesis_seq_score.txt", sep="\t", row.names=FALSE, col.names=TRUE)
```

``` r
seqClose(gds)
```
