PLINK
================
Rutuja Gupte

## Installation

    wget "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20240625.zip"
    unzip plink2_linux_x86_64_20240625.zip
    cd plink2_linux_x86_64_20240625
    chmod 755 plink2

## File

I am assuming that the vcf file has already been preprocessed according
to the requirements of this function. The main assumptions being made
here are:  
1. Chromosomes are numbers and can be converted to integers  
2. There are no missing values  
3. The phenotype file has 2 columns of names and no headers

## Actual Software

This is my favorite one because it is entirely command-line based and
does not need any extra processing other than the common processing we
have done so far. Even within the common pre-processing, if you do not
want to do some of that, it is totally fine because that can be
accounted for using flags.

The main command that I would use to run PLINK is:

    ./plink2 --vcf gt5382_processed.vcf.gz --double-id --pheno gt5382_traits.txt --glm allow-no-covars --pca 5 --out gt5382 --autosome-num 8

For multiple phenotypes in 2.0: The same stuff works. No modifications
needed yay!

    ./plink2 --vcf TASSEL_samples_processed.vcf.gz --double-id --pheno TASSEL_samples_traits.txt --glm allow-no-covars --pca 5 --out TASSEL_samples --autosome-num 10

The beauty of it all is that it works in both PowerShell and Linux. Make
sure all the appropriate files are in the directory containing the
executable PLINK file.

Move the .PHENO1.glm file into the working directory to visualize the
data and make Manhattan plots. We will require the R package qqman.

``` r
plink <- read.table("gt5382/gt5382.PHENO1.glm.linear", header = TRUE)
colnames(plink) <- c("CHROM",   "POS",  "ID",   "REF",  "ALT",  "PROVISIONAL_REF",  "A1",   "OMITTED",  "A1_FREQ",  "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P",    "ERRCODE")
plink <- drop_na(plink)
plink$CHR <- plink$CHROM
plink$BP <- plink$POS
plink$SNP <- paste(plink$CHROM, "_", plink$POS, sep="")
plink$method <- "PLINK"
manhattan(plink)
```

![](PLINK_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
