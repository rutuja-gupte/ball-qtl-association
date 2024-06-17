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
library(tidyverse)
library(vcfR)
```

Reading the data and displaying some basic information

``` r
# This part is for actual data files instead of pulling data from the package.
vcf_raw <- read.vcfR("gt5382/5382--4.vcf.gz")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 5633
    ##   header_line: 5634
    ##   variant count: 11845
    ##   column count: 278
    ## Meta line 1000 read in.Meta line 2000 read in.Meta line 3000 read in.Meta line 4000 read in.Meta line 5000 read in.Meta line 5633 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 11845
    ##   Character matrix gt cols: 278
    ##   skip: 0
    ##   nrows: 11845
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant 11000Processed variant: 11845
    ## All variants processed

``` r
vcf_raw
```

    ## ***** Object of Class vcfR *****
    ## 269 samples
    ## 51 CHROMs
    ## 11,845 variants
    ## Object size: 76.6 Mb
    ## 0 percent missing data
    ## *****        *****         *****

``` r
gt_raw <- extract.gt(vcf_raw, element="GT")
dim(gt_raw)
```

    ## [1] 11845   269

``` r
gt_raw[1:5,1:4]
```

    ##                      BHC_145802_P015_WH09 BHC_145802_P015_WH08
    ## Chr1_RagTag_46992_1  NA                   NA                  
    ## Chr1_RagTag_101865_2 NA                   NA                  
    ## Chr1_RagTag_106556_3 NA                   "0/1"               
    ## Chr1_RagTag_111155_4 "1/1"                "0/1"               
    ## Chr1_RagTag_114979_5 NA                   "0/0"               
    ##                      BHC_145802_P015_WH07 BHC_145802_P015_WH05
    ## Chr1_RagTag_46992_1  NA                   NA                  
    ## Chr1_RagTag_101865_2 NA                   NA                  
    ## Chr1_RagTag_106556_3 NA                   NA                  
    ## Chr1_RagTag_111155_4 NA                   "0/1"               
    ## Chr1_RagTag_114979_5 NA                   "0/0"

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
    ##                    Chr1_RagTag                    Chr2_RagTag 
    ##                           1738                           1274 
    ##                    Chr3_RagTag                    Chr4_RagTag 
    ##                           1072                           1778 
    ##                    Chr5_RagTag                    Chr6_RagTag 
    ##                           1440                           2116 
    ##                    Chr7_RagTag                    Chr8_RagTag 
    ##                            975                           1263 
    ##  scaffold_13_2405393_2456755_+  scaffold_13_2501113_2557509_+ 
    ##                              1                              1 
    ## scaffold_2_43610542_44010732_+ scaffold_2_44010733_44200631_+ 
    ##                              2                              7 
    ## scaffold_2_44378654_44665405_+   scaffold_4_8821050_8835591_+ 
    ##                              1                              1 
    ## scaffold_5_50378355_50389848_+ scaffold_5_50870056_51027793_+ 
    ##                              2                              1 
    ## scaffold_5_66654508_66671430_+ scaffold_5_66676892_66738772_+ 
    ##                              1                              1 
    ## scaffold_5_66738773_66748559_+ scaffold_5_66748560_66760450_+ 
    ##                              4                              1 
    ## scaffold_5_66760451_66773884_+ scaffold_5_66773885_66778886_+ 
    ##                             10                              1 
    ## scaffold_5_66788642_66834494_+ scaffold_5_66834495_66841691_+ 
    ##                             15                              1 
    ## scaffold_5_66841692_66864802_+ scaffold_5_66896087_66946511_+ 
    ##                              3                              4 
    ## scaffold_5_66981881_67018507_+ scaffold_5_67018508_67065375_+ 
    ##                             10                             18 
    ## scaffold_5_67065376_67084796_+ scaffold_5_67152005_67179424_+ 
    ##                              3                              6 
    ## scaffold_5_67179425_67201789_+ scaffold_5_67218797_67238136_+ 
    ##                              5                              3 
    ## scaffold_5_67238137_67243095_+ scaffold_5_67243096_67249609_+ 
    ##                              2                              5 
    ## scaffold_5_67362151_67397328_+ scaffold_5_67397329_67422125_+ 
    ##                              3                              6 
    ## scaffold_5_67422126_67437538_+ scaffold_5_67447924_67536504_+ 
    ##                              1                              8 
    ## scaffold_5_67536505_67557225_+ scaffold_5_67638663_67679663_+ 
    ##                              8                              5 
    ## scaffold_5_67691885_67710111_+ scaffold_5_67710112_67718854_+ 
    ##                              2                              3 
    ## scaffold_5_67736142_67809078_+ scaffold_5_68004996_68017440_+ 
    ##                             23                              1 
    ## scaffold_5_68017441_68056512_+ scaffold_5_68056513_68064199_+ 
    ##                              3                              1 
    ## scaffold_5_68064200_68084460_+ scaffold_5_68123145_68227850_+ 
    ##                              4                              4 
    ## scaffold_5_68425210_68511892_+ scaffold_5_68609885_68632876_+ 
    ##                              3                              1 
    ## scaffold_6_67387877_67436584_+ 
    ##                              4

``` r
# Needs to be modified based on how the dataset labels their chromosomes
chromosomes <- scaffolds[str_detect(scaffolds, "^Chr")]
chromosomes
```

    ## [1] "Chr1_RagTag" "Chr2_RagTag" "Chr3_RagTag" "Chr4_RagTag" "Chr5_RagTag"
    ## [6] "Chr6_RagTag" "Chr7_RagTag" "Chr8_RagTag"

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

    ## [1] "Chr1_RagTag" "Chr2_RagTag" "Chr3_RagTag" "Chr4_RagTag" "Chr5_RagTag"
    ## [6] "Chr6_RagTag" "Chr7_RagTag" "Chr8_RagTag"

``` r
for (val in 1:(length(chromosomes))){
  vcf@fix[vcf@fix[,"CHROM"] == chromosomes[val], "CHROM"] <- val
}
head(vcf@fix)
```

    ##      CHROM POS      ID REF ALT QUAL        FILTER
    ## [1,] "1"   "46992"  NA "G" "A" "11593.6"   NA    
    ## [2,] "1"   "101865" NA "A" "G" "0.0118675" NA    
    ## [3,] "1"   "106556" NA "C" "T" "28887.8"   NA    
    ## [4,] "1"   "111155" NA "C" "T" "105687"    NA    
    ## [5,] "1"   "114979" NA "T" NA  "0"         NA    
    ## [6,] "1"   "198530" NA "T" "C" "0"         NA    
    ##      INFO                                                                                                                                                                                                                                                                                                                                                                                                                                       
    ## [1,] "AB=0.502674;ABP=3.05675;AC=464;AF=0.246023;AN=1886;AO=817;CIGAR=1X;DP=3393;DPB=3393;DPRA=1.103;EPP=1751.14;EPPR=5536.04;GTI=99;LEN=1;MEANALT=1.00323;MQM=60;MQMR=60;NS=1136;NUMALT=1;ODDS=0.0340914;PAIRED=0.995104;PAIREDR=0.998056;PAO=0;PQA=0;PQR=0;PRO=0;QA=29818;QR=93356;RO=2572;RPL=3;RPP=1751.14;RPPR=5536.04;RPR=814;RUN=1;SAF=817;SAP=1777.1;SAR=0;SRF=2572;SRP=5588.04;SRR=0;TYPE=snp;technology.ILLUMINA=1"                   
    ## [2,] "AB=0.25;ABP=18.2106;AC=6;AF=0.00309598;AN=1938;AO=7;CIGAR=1X;DP=4267;DPB=4267;DPRA=1.06016;EPP=5.80219;EPPR=8809.29;GTI=5;LEN=1;MEANALT=1;MQM=47;MQMR=59.7996;NS=1136;NUMALT=1;ODDS=6.12006;PAIRED=0.571429;PAIREDR=0.996946;PAO=0;PQA=0;PQR=0;PRO=0;QA=247;QR=155249;RO=4257;RPL=6;RPP=10.7656;RPPR=8607;RPR=1;RUN=1;SAF=6;SAP=10.7656;SAR=1;SRF=4219;SRP=8919.85;SRR=38;TYPE=snp;technology.ILLUMINA=1"                                 
    ## [3,] "AB=0.489855;ABP=4.86077;AC=521;AF=0.244601;AN=2130;AO=1785;CIGAR=1X;DP=7772;DPB=7772;DPRA=1.01202;EPP=3758.44;EPPR=11201.9;GTI=53;LEN=1;MEANALT=1.01279;MQM=59.9457;MQMR=59.5413;NS=1136;NUMALT=1;ODDS=0.0382684;PAIRED=0.993838;PAIREDR=0.984199;PAO=0;PQA=0;PQR=0;PRO=0;QA=62819;QR=214231;RO=5949;RPL=15;RPP=3749.9;RPPR=9763.23;RPR=1770;RUN=1;SAF=1776;SAP=3801.31;SAR=9;SRF=5761;SRP=11339.8;SRR=188;TYPE=snp;technology.ILLUMINA=1"
    ## [4,] "AB=0.501043;ABP=3.06918;AC=572;AF=0.258123;AN=2216;AO=5504;CIGAR=1X;DP=23243;DPB=23243;DPRA=0.902641;EPP=11920.1;EPPR=38191;GTI=12;LEN=1;MEANALT=1.03596;MQM=60;MQMR=59.9994;NS=1136;NUMALT=1;ODDS=0.322958;PAIRED=0.997638;PAIREDR=0.997233;PAO=0;PQA=0;PQR=0;PRO=0;QA=200305;QR=648116;RO=17706;RPL=5494;RPP=11868.1;RPPR=38225.6;RPR=10;RUN=1;SAF=5498;SAP=11902.7;SAR=6;SRF=17698;SRP=38381.6;SRR=8;TYPE=snp;technology.ILLUMINA=1"   
    ## [5,] "DP=23570;DPB=23570;EPPR=49541.8;GTI=0;MQMR=60;NS=1136;NUMALT=0;ODDS=0;PAIREDR=0.994219;PQR=0;PRO=0;QR=841006;RO=23524;RPPR=12054.2"                                                                                                                                                                                                                                                                                                       
    ## [6,] "AB=0.5;ABP=3.0103;AC=16;AF=0.0536913;AN=298;AO=11;CIGAR=1X;DP=295;DPB=295;DPRA=1.94196;EPP=4.78696;EPPR=117.971;GTI=48;LEN=1;MEANALT=1;MQM=55.6364;MQMR=57.7022;NS=1136;NUMALT=1;ODDS=1848.11;PAIRED=1;PAIREDR=0.963235;PAO=0;PQA=0;PQR=0;PRO=0;QA=403;QR=9358;RO=272;RPL=0;RPP=26.8965;RPPR=35.7101;RPR=11;RUN=1;SAF=4;SAP=4.78696;SAR=7;SRF=74;SRP=125.762;SRR=198;TYPE=snp;technology.ILLUMINA=1"

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

    ## [1] 11656   270

``` r
gt <- extract.gt(vcf, element="GT")

# First looking at missingness across samples to identify any particularly bad samples.
miss_sample <- apply(gt, 2, function(r)mean(is.na(r)))
hist(miss_sample)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Setting an arbitrary threshold
threshold = 0.6
miss_sample[miss_sample>threshold]
```

    ## BHC_145802_P015_WH07 BHC_145802_P015_WG11 BHC_145802_P015_WG09 
    ##            0.7686170            0.6935484            0.6680679 
    ## BHC_145802_P015_WG08 BHC_145802_P015_WG07 BHC_145802_P015_WG05 
    ##            0.6546843            0.8402539            0.7161119 
    ## BHC_145802_P015_WE12 BHC_145802_P015_WE07 BHC_145802_P015_WD12 
    ##            0.6303191            0.6868566            0.6285175 
    ## BHC_145802_P015_WD11 BHC_145802_P015_WD10 BHC_145802_P015_WD08 
    ##            0.6027797            0.6545127            0.6624056 
    ## BHC_145802_P015_WC12 BHC_145802_P015_WC09 BHC_145802_P015_WE10 
    ##            0.6347804            0.7050446            0.6491935 
    ## BHC_145802_P015_WB09 BHC_145802_P015_WC11 BHC_145802_P015_WB06 
    ##            0.6496225            0.6380405            0.7108785 
    ## BHC_145802_P015_WA12 BHC_145802_P015_WA07 BHC_145802_P015_WA10 
    ##            0.6040666            0.6268874            0.7238332 
    ## BHC_145802_P013_WH02 BHC_145802_P013_WH01 BHC_145802_P013_WG02 
    ##            0.6619767            0.6460192            0.7101064 
    ## BHC_145802_P013_WF01 BHC_145802_P013_WD03 BHC_145802_P013_WD02 
    ##            0.6092141            0.6568291            0.6727865 
    ## BHC_145802_P013_WC03 BHC_145802_P015_WF10 BHC_145802_P013_WA03 
    ##            0.6839396            0.6287749            0.6723576 
    ## BHC_145802_P013_WA02 BHC_145802_P012_WH08 BHC_145802_P012_WG10 
    ##            0.6796500            0.6129032            0.6300618 
    ## BHC_145802_P012_WF11 BHC_145802_P015_WB12 BHC_145802_P012_WD11 
    ##            0.6502231            0.6962080            0.6763040 
    ## BHC_145802_P008_WA11 BHC_145802_P010_WD02 BHC_145802_P013_WE01 
    ##            0.9914207            0.9990563            0.6538264 
    ## BHC_145802_P015_WF05 BHC_145802_P015_WE09 BHC_145802_P005_WD03 
    ##            0.6902883            0.7036719            0.6926905 
    ## BHC_145802_P013_WC02 BHC_145802_P013_WB02 BHC_145802_P013_WB03 
    ##            0.6430165            0.7086479            0.6727865 
    ## BHC_145802_P015_WF12 BHC_145802_P015_WH11 BHC_145802_P015_WF11 
    ##            0.6607756            0.6022649            0.6884866 
    ## BHC_145802_P008_WC12 BHC_145802_P013_WF02 BHC_145802_P012_WH12 
    ##            0.9969973            0.6924331            0.6369252 
    ## BHC_145802_P015_WB08 BHC_145802_P015_WD06 BHC_145802_P015_WG10 
    ##            0.6264585            0.6422443            0.6465340 
    ## BHC_145802_P015_WA09 BHC_145802_P009_WG09 BHC_145802_P015_WF08 
    ##            0.6964653            0.7267502            0.6541695 
    ## BHC_145802_P015_WC10 BHC_145802_P015_WE11 BHC_145802_P015_WB10 
    ##            0.7265786            0.6660089            0.6783631 
    ## BHC_145802_P009_WA01 BHC_145802_P015_WA11 BHC_145802_P009_WD11 
    ##            0.6576012            0.6388126            0.9977694 
    ## BHC_145802_P009_WE09 BHC_145802_P009_WF12 BHC_145802_P009_WF01 
    ##            0.6164207            0.9987989            0.6177934 
    ## BHC_145802_P015_WF09 BHC_145802_P015_WA06 BHC_145802_P015_WB11 
    ##            0.6237131            0.6411290            0.6565717 
    ## BHC_145802_P010_WD04 BHC_145802_P010_WF01 BHC_145802_P005_WA03 
    ##            0.6100721            0.9981984            0.7143102 
    ## BHC_145802_P015_WB07 BHC_145802_P013_WE02 BHC_145802_P012_WC11 
    ##            0.6884008            0.6419012            0.7425360 
    ## BHC_145802_P012_WA11 BHC_145802_P012_WB10 BHC_145802_P012_WB11 
    ##            0.6402711            0.6106726            0.6609472 
    ## BHC_145802_P012_WB12 BHC_145802_P015_WC08 
    ##            0.7213452            0.6322066

``` r
# Updating the vcf
vcf@gt <- vcf@gt[, c("FORMAT", colnames(gt)[apply(gt, 2, function(r)mean(is.na(r)) < threshold)])]
gt <- extract.gt(vcf, element="GT")
dim(gt)
```

    ## [1] 11656   189

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

    ## [1] 2072

``` r
# Updating the vcf
vcf@gt <- vcf@gt[apply(gt, 1, function(r)(mean(is.na(r)) == 0)),]
vcf@fix <- vcf@fix[apply(gt, 1, function(r)(mean(is.na(r)) == 0)),]

gt <- extract.gt(vcf, element="GT")

# Sanity check to see if it actually worked
head(apply(gt, 2, function(r)mean(is.na(r))))
```

    ## BHC_145802_P015_WH09 BHC_145802_P015_WH08 BHC_145802_P015_WH05 
    ##                    0                    0                    0 
    ## BHC_145802_P015_WF06 BHC_145802_P015_WE08 BHC_145802_P015_WC07 
    ##                    0                    0                    0

``` r
# Another sanity check to see if the dimensions are as expected.
dim(gt)
```

    ## [1] 2072  189

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
    ## 189 samples
    ## 8 CHROMs
    ## 2,072 variants
    ## Object size: 24.8 Mb
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

    ## [1] 1526  189

``` r
alt <- vcf@fix[,"ALT"]
vcf@gt <- vcf@gt[!str_detect(alt,","),]
vcf@fix <- vcf@fix[!str_detect(alt,","),]
gt <- extract.gt(vcf, element="GT")

vcf
```

    ## ***** Object of Class vcfR *****
    ## 189 samples
    ## 8 CHROMs
    ## 1,526 variants
    ## Object size: 21.3 Mb
    ## 0 percent missing data
    ## *****        *****         *****

### Keeping only SNPs

These softwares are mainly designed to work for SNPs

``` r
vcf <- extract.indels(vcf, return.indels = FALSE)
vcf
```

    ## ***** Object of Class vcfR *****
    ## 189 samples
    ## 8 CHROMs
    ## 1,526 variants
    ## Object size: 21.3 Mb
    ## 0 percent missing data
    ## *****        *****         *****

### Modifying the phenotype file

Also converting all phenotypes to integers because PLINK needs all
integers.

``` r
# Also removing samples for which we do not have phenotype information
traits <- read.table("gt5382/5382--4.pheno_1100.txt")
traits <- drop_na(traits)
vcf@gt <- vcf@gt[, c("FORMAT", colnames(gt)[colnames(gt) %in% traits[,1]])]
gt <- extract.gt(vcf, element="GT")
dim(gt)
```

    ## [1] 1526  189

``` r
traits <- traits[traits[,1] %in% colnames(gt),]
traits <- traits %>% mutate(across(3:ncol(traits), round))
```

This object can be saved as a vcf file that can be used to make the
other models.

*Important realization:* write.vcf writes a gzipped file by default. I
had been using the wrong extension all along.

``` r
write.vcf(vcf, "gt5382/gt5382_processed.vcf.gz")
write.table(traits, "gt5382/gt5382_traits.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
```

## GModel

## GAPIT

## PLINK
