QTL Association
================
Rutuja Gupte

This is a guide to using 3 different softwares for linkage mapping of
Quantitative Trait Loci (QTL). The original guides of all the softwares
are linked at the bottom.

The file formatting guide here is based in R using the vcfR package. I
am assuming that R is installed. If not, please visit [The Comprehensive
R Archive Network](https://cran.rstudio.com/) and find the appropriate
version of R. Here, I am using the R package vcfR to preprocess the
files to generate the appropriate files for all softwares. Some of these
softwares require less preprocessing than others but I am following the
same routine to clean the data to be consistent across softwares. So,
here are the instructions for vcfR installation.

## Installations

Let us begin with the installations. There will be a lot of file
processing using vcfR. I also like to use tidyverse for managing my
data. But you do not have to use tidyverse. If these libraries are not
installed. Run the following lines of code on an R console. vcfR may
take a longer time to install.

    install.packages("tidyverse")
    install.packages("vcfR")

### Troubleshooting

This will work with Windows but may not work with Linux. As a side note,
some R packages need a C and a Fortran compiler to compile some
packages. I already had the GCC compiler and had to download a Fortran
compiler.

    sudo apt install --reinstall gcc-12
    sudo apt-get install gfortran

vcfR also requires ape and vegan packages other than the standard CRAN
repositories. For Linux, if install.packages is not working. Start by
downloading the source code for all 3 packages (.tar.gz files) in the
appropriate folder with the other packages.

The source links for the packages can be found on their official CRAN
website. 1.
[ape](%22https://cran.r-project.org/web/packages/ape/index.html%22) 2.
[vegan](%22https://cran.r-project.org/web/packages/vegan/index.html%22)
3. [vcfR](%22https://cran.r-project.org/web/packages/vcfR/index.html%22)

Then install the 3 packages in this order using R code that looks like
this.

    install.packages("/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4/ape_5.8.tar.gz",lib="/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4", repos = NULL, type="source")
    install.packages("/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4/vegan_2.6-6.1.tar.gz",lib="/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4", repos = NULL, type="source")
    install.packages("/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4/vcfR_1.15.0.tar.gz",lib="/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4", repos = NULL, type="source")

## Preprocessing

These steps are common for all the softwares. Before beginnging the
specific file processing for each software please go through these
steps. I have included notes about which of the steps are strictly
necessary for which tools, so you may be able to skip some things as
needed.

An extra step here is that I am using a sample VCF dataset from pinfsc50
to demonstrate the steps. You do not need to install this package.

``` r
library(pinfsc50)
library(tidyverse)
library(vcfR)
```

Reading the data and displaying some basic information

``` r
vcf_raw <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
vcf_raw <- read.vcfR(vcf_raw)
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 29
    ##   header_line: 30
    ##   variant count: 22031
    ##   column count: 27
    ## Meta line 29 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 22031
    ##   Character matrix gt cols: 27
    ##   skip: 0
    ##   nrows: 22031
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant 11000Processed variant 12000Processed variant 13000Processed variant 14000Processed variant 15000Processed variant 16000Processed variant 17000Processed variant 18000Processed variant 19000Processed variant 20000Processed variant 21000Processed variant 22000Processed variant: 22031
    ## All variants processed

``` r
# # This part is for actual data files instead of pulling data from the package.
# vcf_raw <- read.vcfR("GC1_round1/GC1_round1_gt.vcf")
# vcf_raw

gt_raw <- extract.gt(vcf_raw, element="GT")
dim(gt_raw)
```

    ## [1] 22031    18

``` r
gt_raw[1:5,1:4]
```

    ##                      BL2009P4_us23 DDR7602 IN2009T1_us22 LBUS5
    ## Supercontig_1.50_41  "1|1"         "1|1"   "1|1"         "1|1"
    ## Supercontig_1.50_136 "0|0"         "0|0"   "0|0"         "0|0"
    ## Supercontig_1.50_254 "0|0"         "0|0"   "0|0"         "0|0"
    ## Supercontig_1.50_275 "0|0"         "0|0"   "0|0"         "0|0"
    ## Supercontig_1.50_386 "0|0"         "0|0"   "0|0"         "0|0"

### Removing scaffolds

This is not necessary for this dataset but may be required for real
datasets which may have some scaffolds that are not aligned to any
chromosome. I am also changing the chromosome names to numbers. This is
only required for GModel but it is easier if we do it for all 3.
Chromosome name format really does not matter for GAPIT and for PLINK
there is an extra flag that needs to be passed to allow for non-integer
chromosomes.

**Additional side note about vcfR:** A vcfR object contains meta, fix
and gt, representing the metadata, fixed columns and the variable
genotype columns respectively, which can be accessed using @. Also note
that `vcf@gt` and `extract.gt(vcf, element="GT")` are different. The
former contains multiple kinds of information about the genotypes, has
integer rownames and includes a format column. The later includes only
the genotypes and is named by the chromosome and position.

``` r
rm(gt_raw)
vcf <- vcf_raw

gt <- extract.gt(vcf, element="GT")

# Checking the chromosome names to spot any contigs
scaffolds <- unique(vcf@fix[,"CHROM"])
table(vcf@fix[,"CHROM"])
```

    ## 
    ## Supercontig_1.50 
    ##            22031

``` r
# Needs to be modified based on how the dataset labels their chromosomes
chromosomes <- scaffolds[str_detect(scaffolds, "^Supercontig")]
chromosomes
```

    ## [1] "Supercontig_1.50"

``` r
bool_ser <- vcf@fix[, "CHROM"] %in% chromosomes

# filtering the rows of both the gt and fix part of the data
vcf@gt <- vcf@gt[bool_ser,]
vcf@fix <- vcf@fix[bool_ser,]

rm(bool_ser)

# Sorting the chromosomes for numbering (sometimes causes more harm than good)
# chromosomes <- sort(chromosomes)

# Here I am assuming that we have all the chromosomes starting from 1. Needs additional modification if that is not the case.
chromosomes
```

    ## [1] "Supercontig_1.50"

``` r
for (val in 0:(length(chromosomes))){
  vcf@fix[vcf@fix[,"CHROM"] == chromosomes[val], "CHROM"] <- val
}
head(vcf@fix)
```

    ##      CHROM POS   ID REF  ALT QUAL      FILTER
    ## [1,] "1"   "41"  NA "AT" "A" "4784.43" NA    
    ## [2,] "1"   "136" NA "A"  "C" "550.27"  NA    
    ## [3,] "1"   "254" NA "T"  "G" "774.44"  NA    
    ## [4,] "1"   "275" NA "A"  "G" "714.53"  NA    
    ## [5,] "1"   "386" NA "T"  "G" "876.55"  NA    
    ## [6,] "1"   "462" NA "T"  "G" "1301.07" NA    
    ##      INFO                                                                                                                                                                                                
    ## [1,] "AC=32;AF=1.00;AN=32;DP=174;FS=0.000;InbreedingCoeff=-0.0224;MLEAC=32;MLEAF=1.00;MQ=51.30;MQ0=0;QD=27.50;SOR=4.103"                                                                                 
    ## [2,] "AC=2;AF=0.059;AN=34;BaseQRankSum=-0.116;ClippingRankSum=-0.831;DP=390;FS=0.000;InbreedingCoeff=-0.0292;MLEAC=2;MLEAF=0.059;MQ=52.83;MQ0=0;MQRankSum=3.872;QD=11.01;ReadPosRankSum=2.829;SOR=0.632" 
    ## [3,] "AC=3;AF=0.088;AN=34;BaseQRankSum=-2.565;ClippingRankSum=0.268;DP=514;FS=1.169;InbreedingCoeff=0.5463;MLEAC=2;MLEAF=0.059;MQ=56.79;MQ0=0;MQRankSum=-7.878;QD=16.48;ReadPosRankSum=1.300;SOR=0.804"  
    ## [4,] "AC=3;AF=0.088;AN=34;BaseQRankSum=-3.812;ClippingRankSum=-0.084;DP=514;FS=0.000;InbreedingCoeff=0.5586;MLEAC=3;MLEAF=0.088;MQ=57.07;MQ0=0;MQRankSum=-6.942;QD=15.88;ReadPosRankSum=-0.670;SOR=0.765"
    ## [5,] "AC=3;AF=0.094;AN=32;BaseQRankSum=-4.806;ClippingRankSum=0.793;DP=509;FS=2.356;InbreedingCoeff=0.5896;MLEAC=3;MLEAF=0.094;MQ=57.40;MQ0=0;MQRankSum=-0.200;QD=15.38;ReadPosRankSum=-0.290;SOR=0.876" 
    ## [6,] "AC=3;AF=0.088;AN=34;BaseQRankSum=-4.788;ClippingRankSum=0.096;DP=508;FS=0.000;InbreedingCoeff=0.5423;MLEAC=3;MLEAF=0.088;MQ=58.89;MQ0=0;MQRankSum=-1.160;QD=17.58;ReadPosRankSum=-0.467;SOR=0.581"

``` r
# copying vcf into vcf_raw as a backup before the next step
vcf_raw <- vcf
```

### Removing missing values

There is no one way to do this. The sample dataset here is good enough
that we can remove all rows with missing genotypes and still have a
sizable dataset left over. But that may not always be the case. If it is
not possible to remove all rows with missing values, a viable
alternative is to remove some bad samples and then remove all rows with
missing values. Ultimately, we need to end up with a dataset that has no
missing values. GModel does not support missing values. GAPIT and PLINK
do support missing values but may have different methods of processing
the missing values during analysis leading to different results.

Starting by quantifying and plotting missingness across samples in order
to decide on the best strategy to filter the missing values. Followed by
actually removing the bad samples with an arbitrary threshold.

``` r
# This is here just to make it easier to rerun
vcf <- vcf_raw
dim(vcf@gt)
```

    ## [1] 22031    19

``` r
gt <- extract.gt(vcf, element="GT")

# First looking at missingness across samples to identify any particularly bad samples.
miss_sample <- apply(gt, 2, function(r)mean(is.na(r)))
hist(miss_sample)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Setting an arbitrary threshold
threshold = 0.2
miss_sample[miss_sample>threshold]
```

    ##     P7722     t30-4 
    ## 0.2191004 0.2801961

``` r
# Updating the vcf
vcf@gt <- vcf@gt[, c("FORMAT", colnames(gt)[apply(gt, 2, function(r)mean(is.na(r)) < threshold)])]
gt <- extract.gt(vcf, element="GT")
dim(gt)
```

    ## [1] 22031    16

Now quantifying and plotting missingness across variants. Then removing
all rows with missing values. If there are too many variants with
missing values, this code can be updated to set a tolerance for some
amount of missingness or the previous chunk can be modified to remove
more bad samples.

``` r
# Looking at missingness across variants.   
miss_var <- apply(gt, 1, function(r)mean(is.na(r)))
hist(miss_var)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# Number of rows that will be leftover after removing the missing rows.
sum(miss_var == 0)
```

    ## [1] 16648

``` r
# Updating the vcf
vcf@gt <- vcf@gt[apply(gt, 1, function(r)(mean(is.na(r)) == 0)),]
vcf@fix <- vcf@fix[apply(gt, 1, function(r)(mean(is.na(r)) == 0)),]

gt <- extract.gt(vcf, element="GT")

# Sanity check to see if it actually worked
head(apply(gt, 2, function(r)mean(is.na(r))))
```

    ## BL2009P4_us23       DDR7602 IN2009T1_us22         LBUS5       NL07434 
    ##             0             0             0             0             0 
    ##        P10127 
    ##             0

``` r
# Another sanity check to see if the dimensions are as expected.
dim(gt)
```

    ## [1] 16648    16

``` r
# This is just a quick check to see if there are any genotypes where missing values are given by '.' which is a common practice in the VCF format.
sum(str_detect(gt, "\\."), na.rm = TRUE)
```

    ## [1] 0

Now looking at the updated VCF. At this step, it should say 0 percent
missing data or you may have missed something.

``` r
vcf
```

    ## ***** Object of Class vcfR *****
    ## 16 samples
    ## 1 CHROMs
    ## 16,648 variants
    ## Object size: 17.4 Mb
    ## 0 percent missing data
    ## *****        *****         *****

### Keeping only biallelic loci

Biallelic loci are separated by ‘,’. So, I am trying to detect ‘,’ in
the ALT column and removing rows with multiple alleles listed.

``` r
# Checking for NAs in the ALT column
alt <- vcf@fix[,"ALT"]
vcf@gt <- vcf@gt[!is.na(alt),]
vcf@fix <- vcf@fix[!is.na(alt),]
gt <- extract.gt(vcf, element="GT")
dim(gt)
```

    ## [1] 16648    16

``` r
alt <- vcf@fix[,"ALT"]
vcf@gt <- vcf@gt[!str_detect(alt,","),]
vcf@fix <- vcf@fix[!str_detect(alt,","),]
gt <- extract.gt(vcf, element="GT")

vcf
```

    ## ***** Object of Class vcfR *****
    ## 16 samples
    ## 1 CHROMs
    ## 16,375 variants
    ## Object size: 16.9 Mb
    ## 0 percent missing data
    ## *****        *****         *****

### Keeping only SNPs

These softwares are mainly designed to work for SNPs

``` r
vcf <- extract.indels(vcf, return.indels = FALSE)
vcf
```

    ## ***** Object of Class vcfR *****
    ## 16 samples
    ## 1 CHROMs
    ## 14,643 variants
    ## Object size: 15.3 Mb
    ## 0 percent missing data
    ## *****        *****         *****

This object can be saved as a vcf file that can be used to make the
other models.

*Important realization:* write.vcf writes a gzipped file by default. I
had been using the wrong extension all along.

``` r
write.vcf(vcf, "sample/processed_vcf.vcf.gz")
```

## GModel

## GAPIT

## PLINK
