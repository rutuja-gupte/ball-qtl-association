GModel
================
Rutuja Gupte

``` r
library(tidyverse)
library(vcfR)
```

This guide is for GModel2. The main difference between GModel and
GModel2 is the GModel needs prior estimates of the proportion of
variation that is due to genetic effects. Thus, GModel2 takes
considerably longer to run. In fact, among the softwares being compared
here, this is the slowest. However, it is very easy to download and easy
to use but requires files following a highly specific format.

## Installation

Now time to download GModel. Open a terminal.

    wget "https://bernardo-group.org/wp-content/uploads/2021/01/GModel.zip"
    unzip GModel.zip
    cd GModel_release/
    chmod 755 GModel2.exe

Make sure that GModel2.exe is executable and test it.

    ./GModel2.exe

There will be a prompt asking for the parameter file. The GModel
download folder contains sample files for both GModel and GModel2. To
check if your installation has working, pass `ParmsGModel2.csv` as the
parameter file. This file should be a part of the folder that we
downloaded from the website. Again, warning, it can take more than an
hour to run.

## Preparing the files

I am assuming that the vcf file has already been preprocessed according
to the requirements of this function. The main assumptions being made
here are:  
1. Chromosomes are numbers and can be converted to integers  
2. There are no missing values  
3. The phenotype file has 2 columns of names and no headers

The sample phenotype file is “traits.csv”

This function should generate 3 files. The last file “f4.csv” can be
passed as an argument to GModel. “f4.csv” should be updated based on the
number of phenotypes being studied.

``` r
run_gmodel <- function(vcf.name, pheno.name, f1 = "gmodel1.csv", f2 = "gmodel2.csv", f3 = "gnodel3.csv"){
  vcf <- read.vcfR(vcf.name)
  traits <- read.table(pheno.name)
  gt <- extract.gt(vcf, element="GT")
  
  # Marker names and chromosomes file
  t1 <- vcf@fix[,c("POS", "CHROM")]
  t1 <- data.frame(t1)
  t1$Markers <- rownames(gt)
  t1$CHROM <- as.numeric(t1$CHROM)
  t1$POS <- as.numeric(t1$POS)
  t1 <- arrange(t1, CHROM, POS)
  t1 <- t1 %>% select(Markers, CHROM)
  write.table(t1, f1, sep=",", col.names = FALSE, row.names = FALSE)
  
  # SNP Marker data
  # No missing data
  t2 <- gt
  # Allowing for both phased and unphased files
  t2[t2 == "0|0"] <- "1"
  t2[t2 == "0|1"] <- "0"
  t2[t2 == "1|0"] <- "0"
  t2[t2 == "1|1"] <- "-1"
  t2[t2 == "0/0"] <- "1"
  t2[t2 == "0/1"] <- "0"
  t2[t2 == "1/0"] <- "0"
  t2[t2 == "1/1"] <- "-1"
  t2 <- apply(t2, 1:2, as.numeric)
  # Now getting the right order
  t2 <- data.frame(t2)
  t2$chr <- as.numeric(vcf@fix[, "CHROM"])
  t2$pos <- as.numeric(vcf@fix[, "POS"])
  t2 <- arrange(t2, chr, pos)
  t2 <- t2 %>% select(-chr, -pos)
  write.table(t2, f2, sep=",", col.names = FALSE, row.names = FALSE)
  
  # Phenotypic data
  t3 <- data.frame(var1 = str_split_i(colnames(gt), "\\.", 1))
  t3 <- left_join(t3, traits, by = c("var1"="V1"))
  t3 <- t3[, c(3:ncol(t3))]
  t3 <- data.frame(t3) %>% mutate(
    across(everything(), ~replace_na(.x, 0))
  )
  write.table(t3, f3, sep=",", row.names = FALSE)
}
```

Testing the function here

``` r
run_gmodel("sample/processed_vcf.vcf.gz", "sample/traits.txt", "sample/gmodel1.csv", "sample/gmodel2.csv", "sample/gmodel3.csv" )
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 29
    ##   header_line: 30
    ##   variant count: 14643
    ##   column count: 25
    ## Meta line 29 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 14643
    ##   Character matrix gt cols: 25
    ##   skip: 0
    ##   nrows: 14643
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant 11000Processed variant 12000Processed variant 13000Processed variant 14000Processed variant: 14643
    ## All variants processed
