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
run_gmodel <- function(vcf.name, pheno.name, f1 = "gmodel1.csv", f2 = "gmodel2.csv", f3 = "gmodel3.csv"){
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
run_gmodel("gt5382/gt5382_processed.vcf.gz", "gt5382/gt5382_traits.txt", "gt5382/gmodel1.csv", "gt5382/gmodel2.csv", "gt5382/gmodel3.csv" )
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 5633
    ##   header_line: 5634
    ##   variant count: 1526
    ##   column count: 198
    ## Meta line 1000 read in.Meta line 2000 read in.Meta line 3000 read in.Meta line 4000 read in.Meta line 5000 read in.Meta line 5633 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 1526
    ##   Character matrix gt cols: 198
    ##   skip: 0
    ##   nrows: 1526
    ##   row_num: 0
    ## Processed variant 1000Processed variant: 1526
    ## All variants processed

Double check f4.csv here to make sure it has the right number of traits
(especially important for multiple traits)
