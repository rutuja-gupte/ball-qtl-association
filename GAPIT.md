GAPIT
================
Rutuja Gupte

GAPIT is available as an R package that is quite easy to download.
Please make sure you have a C/C++ compiler available.

    source("http://zzlab.net/GAPIT/GAPIT.library.R") 
    source("http://zzlab.net/GAPIT/gapit_functions.txt") 

This was the installation method that worked best for me. But if this
does not work, try the other method

    install.packages("devtools") 
    devtools::install_github("jiabowang/GAPIT3",force=TRUE) 
    library(GAPIT3) 

This may or may not work but the first which is why I prefer the first
way.

The file processing here does not require creating separate files but I
found that it was easier to do that way. This function is written to
prepare all the required files and immediately run GAPIT since the
entire process can be done in R. It will generate multiple files and
will give warnings. (But that is mainly a version and deprecation issue)

I am assuming that the vcf file has already been preprocessed according
to the requirements of this function. The main assumptions being made
here are:  
1. Chromosomes are numbers and can be converted to integers  
2. There are no missing values  
3. The phenotype file has 2 columns of names and no headers

*Important note:* GAPIT does support missing values

``` r
run_gapit <- function(vcf.name, pheno.name, f1 = "gapit1.txt", f2 = "gapit2.txt", f3 = "gapit3.txt"){
  vcf <- read.vcfR(vcf.name)
  traits <- read.table(pheno.name)
  gt <- extract.gt(vcf, element="GT")
  
  # Phenotype data
  t1 <- data.frame(var1 = colnames(gt))
  t1 <- left_join(t1, traits, by = c("var1"="V1"))
  t1 <- drop_na(t1)
  t1 <- t1[, c(2:ncol(t1))]
  t1 <- data.frame(t1) %>% mutate(
    across(everything(), ~replace_na(.x, 0))
  )
  write.table(t1, f1, sep="\t", row.names=FALSE, quote=FALSE)
  myY <- read.table(f1, head = TRUE)
  
  # Genotype data
  t2 <- gt
  # Allowing for both phased and unphased files
  t2[t2 == "0|0"] <- "0"
  t2[t2 == "0|1"] <- "1"
  t2[t2 == "1|0"] <- "1"
  t2[t2 == "1|1"] <- "2"
  t2[t2 == "0/0"] <- "0"
  t2[t2 == "0/1"] <- "1"
  t2[t2 == "1/0"] <- "1"
  t2[t2 == "1/1"] <- "2"
  t2 <- apply(t2, 1:2, as.numeric)
  t2 <- t(t2)
  t2 <- data.frame(t2)
  t2$taxa <- rownames(t2)
  t2 <- t2 %>% relocate(taxa)
  # column_names <- colnames(t2)
  # column_names[1] <- "taxa"
  write.table(t2, f2, sep="\t", row.names=FALSE, quote=FALSE)
  myGD <- read.table(f2, head = TRUE)
  
  # Genetic map
  t3 <- vcf@fix[,c("POS", "CHROM")]
  t3 <- data.frame(t3)
  t3$Markers <- rownames(gt)
  t3$CHROM <- as.numeric(t3$CHROM)
  t3$POS <- as.numeric(t3$POS)
  t3 <- t3 %>% relocate(Markers, CHROM, POS)
  write.table(t3, f3, sep="\t", row.names=FALSE, quote=FALSE)
  myGM <- read.table(f3, head = TRUE)
  
  myGAPIT=GAPIT(
    Y=myY[,], #first column is ID
    GD=myGD,
    GM=myGM,
    PCA.total=3,
    model=c(
      "GLM",
      "SUPER",
      "FarmCPU", 
      "BLINK"),
    kinship.algorithm = "VanRaden",
    kinship.cluster = "average",
    Multiple_analysis=TRUE,)
}
```

``` r
run_gapit("gt5382/gt5382_processed.vcf.gz", "gt5382/gt5382_traits.txt", "gt5382/gapit1.txt", "gt5382/gapit2.txt", "gt5382/gapit3.txt" )
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
    ## [1] "--------------------- Welcome to GAPIT ----------------------------"
    ## [1] "GLM"     "SUPER"   "FarmCPU" "BLINK"  
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 1 model in all."
    ## [1] "GLM"
    ## [1] "GAPIT.DP in process..."
    ## [1] "GAPIT will filter marker with MAF setting !!"
    ## [1] "The markers will be filtered by SNP.MAF: 0"
    ## maf_index
    ## TRUE 
    ## 1526 
    ## [1] "Calculating kinship..."
    ## [1] "Number of individuals and SNPs are  189  and  1526"
    ## [1] "Calculating kinship with VanRaden method..."
    ## [1] "substracting P..."
    ## [1] "Getting X'X..."
    ## [1] "Adjusting..."
    ## [1] "Calculating kinship with VanRaden method: done"
    ## [1] "kinship calculated"
    ## [1] "Creating heat map for kinship..."

    ## [1] "Kinship heat map created"
    ## [1] "Adding IDs to kinship..."
    ## [1] "Writing kinship to file..."
    ## [1] "Kinship save as file"
    ## [1] "Kinship created!"
    ## [1] "Calling prcomp..."
    ## [1] "Creating PCA graphs..."

    ## [1] "Joining taxa..."
    ## [1] "Exporting PCs..."
    ## [1] "PC created"
    ## [1] "Filting marker for GAPIT.Genotype.View function ..."

    ## [1] "GAPIT.Genotype.View . pdfs generate.successfully!"
    ## [1] 189   4
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V3"
    ## [1] "GAPIT.Phenotype.View in press..."

    ## [1] "GAPIT.Phenotype.View output pdf has been generated successfully!"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "VanRaden"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  189  individuals in phenotype data."
    ## [1] "There are  1526  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 189 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 189   4
    ## [1] "GAPIT.IC accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT.SS in process..."
    ## [1] "GAPIT will be into GWAS approach..."
    ## [1] "VanRaden"
    ## [1] "The GAPIT would go into Main..."
    ## [1] "------------Examining data (QC)------------------------------------------"
    ## [1] "Try to group from and to were set to 1"
    ## [1] "------------Examining data (QC) done-------------------------------------"
    ## [1] "-------------------Sandwich burger and dressing------------------------"

    ## The upper bound of groups (group.to) is not sufficient. both boundries were set to a and GLM is performed!

    ## The upper bound of groups (group.to) is not sufficient. both boundries were set to their minimum and GLM is performed!

    ## [1] "-------------------------Iteration in process--------------------------"
    ## [1] "Total iterations: 1"
    ## [1] "Compressing and Genome screening..."
    ## [1] "-------Mixed model with Kinship-----------------------------"
    ## [1] "Genotype file: 1, SNP: 1000 "
    ## [1] "1 of 1 -- Vg= 0 VE= 1.8509 -2LL= 642.91   Clustering= average   Group number= 1   Group kinship= Mean"
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "--------------------Final results presentations------------------------"
    ## [1] "Generating summary"
    ## [1] "Genomic Breeding Values (GBV) ..."
    ## [1] "Writing GBV and Acc..."
    ## [1] "GBV and accuracy distribution..."
    ## [1] "Compression portfolios..."
    ## [1] "Compression Visualization done"
    ##      Type   Cluster   Group REML               VA  VE                
    ## [1,] "Mean" "average" "1"   "642.908865892244" "0" "1.85092288386483"
    ## [1] "p3d objects transfered"
    ## [1] "Merge BLUP and BLUE"
    ## [1] "GAPIT before BLUP and BLUE"
    ## [1] "GAPIT after BLUP and BLUE"
    ## [1] "Exporting BLUP and Pred"
    ## [1] "GLM.V3 has been analyzed successfully!"
    ## [1] "The results are saved in the directory of  C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation"
    ## [1] "before ending GAPIT.Main"

    ## Loading required package: lme4

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 47.89860 46.07586
    ## [1] "GAPIT.ID in process..."
    ## [1] "GAPIT.Compression.Visualization"
    ## [1] "Pie chart"

    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 0 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 0 candidate significont markers in 3 chromosome "

    ## [1] "select 0 candidate significont markers in 4 chromosome "

    ## [1] "select 2 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 2 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 2 model in all."
    ## [1] "SUPER"
    ## [1] "Processing trait: V3"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "Zhang"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  189  individuals in phenotype data."
    ## [1] "There are  1526  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 189 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 189   4
    ## [1] "GAPIT.IC accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT.SS in process..."
    ## [1] "GAPIT will be into GWAS approach..."
    ## [1] "Zhang"
    ## [1] "The GAPIT would go into Main..."
    ## [1] "------------Examining data (QC)------------------------------------------"
    ## [1] "Try to group from and to were set to 1"
    ## [1] "group from and to were set to 1"
    ## [1] "------------Examining data (QC) done-------------------------------------"
    ## [1] "-------------------Sandwich top bun-----------------------------------"
    ## [1] "--------------------- Welcome to GAPIT ----------------------------"
    ## [1] "MLM"
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 1 model in all."
    ## [1] "MLM"
    ## [1] "GAPIT.DP in process..."
    ## [1] "GAPIT will filter marker with MAF setting !!"
    ## [1] "The markers will be filtered by SNP.MAF: 0"
    ## maf_index
    ## TRUE 
    ## 1526 
    ## NULL
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V3"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "Zhang"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  189  individuals in phenotype data."
    ## [1] "There are  1526  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There are 189 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 189   4
    ## [1] "GAPIT.IC accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT.SS in process..."
    ## [1] "GAPIT will be into GWAS approach..."
    ## [1] "Zhang"
    ## [1] "The GAPIT would go into Main..."
    ## [1] "------------Examining data (QC)------------------------------------------"
    ## [1] "Try to group from and to were set to 1"
    ## [1] "------------Examining data (QC) done-------------------------------------"
    ## [1] "-------------------Sandwich burger and dressing------------------------"

    ## The upper bound of groups is too high. It was set to the size of kinship!

    ## The lower bound of groups is too high. It was set to the size of kinship!

    ## [1] "-------------------------Iteration in process--------------------------"
    ## [1] "Total iterations: 1"
    ## [1] "Compressing and Genome screening..."
    ## [1] "-------Mixed model with Kinship-----------------------------"
    ## [1] "Genotype file: 1, SNP: 1000 "
    ## [1] "1 of 1 -- Vg= 0.1459 VE= 1.7748 -2LL= 638.44   Clustering= average   Group number= 189   Group kinship= Mean"
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "--------------------Final results presentations------------------------"
    ## [1] "Generating summary"
    ## [1] "Genomic Breeding Values (GBV) ..."
    ## [1] "GAPIT.GS accomplished successfully!"
    ## [1] "Writing GBV and Acc..."
    ## [1] "GBV and accuracy distribution..."
    ## [1] "Compression portfolios..."
    ## [1] "Compression Visualization done"
    ##      Type   Cluster   Group REML               VA                 
    ## [1,] "Mean" "average" "189" "638.444916489454" "0.145895556478027"
    ##      VE                
    ## [1,] "1.77484423566496"
    ## [1] "p3d objects transfered"
    ## [1] "Merge BLUP and BLUE"
    ## [1] "GAPIT before BLUP and BLUE"
    ## [1] "GAPIT after BLUP and BLUE"
    ## [1] "MLM.V3 has been analyzed successfully!"
    ## [1] "The results are saved in the directory of  C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation"
    ## [1] "before ending GAPIT.Main"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "GAPIT has done all analysis!!!"
    ## [1] "-------------------Sagnwich top bun: done-----------------------------"
    ## [1] "-------------------Sandwich burger and dressing------------------------"

    ## The upper bound of groups (group.to) is not sufficient. both boundries were set to a and GLM is performed!

    ## The upper bound of groups (group.to) is not sufficient. both boundries were set to their minimum and GLM is performed!

    ## [1] "-------------------------Iteration in process--------------------------"
    ## [1] "Total iterations: 1"
    ## [1] "Compressing and Genome screening..."
    ## [1] "-------The burger is SNP-----------------------------------"
    ## [1] "bin---10000---inc---10"
    ##  [1]  14  71  74 396 490 492 531 560 748 816
    ##    Type Cluster   Group    REML      VA      VE 
    ##       1   10000      10      NA      NA      NA 
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "---------------Sandwich bottom with grilled burger---------------------"
    ## [1] "---------------Sandwich bottom: reload bins ---------------------------"
    ## [1] "SNP: 500 "
    ## [1] "SNP: 1000 "
    ## [1] "SNP: 1500 "
    ## [1] "SUPER saving results..."
    ## [1] "--------------------Final results presentations------------------------"
    ## [1] "Generating summary"
    ## [1] "Genomic Breeding Values (GBV) ..."
    ## [1] "Merge BLUP and BLUE"
    ## [1] "SUPER.V3 has been analyzed successfully!"
    ## [1] "The results are saved in the directory of  C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation"
    ## [1] "before ending GAPIT.Main"
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ##  [1] 1.276044e-08 1.755053e+01 1.822702e+01 5.603647e-01 4.915592e-01
    ##  [6] 9.991344e+00 1.038430e+00 2.348105e+00 4.083590e-01 1.016975e+00
    ## [11] 6.317494e-01 1.620364e+00 2.102128e-01 1.570417e+00 4.568719e-01
    ## [16] 6.153972e-09 2.339437e-01 3.386209e-02 2.347657e-01 1.055445e+00
    ## [21] 3.021636e-01 6.591997e-01 9.385469e-01 1.287386e+00 1.621729e+01
    ## [26] 0.000000e+00 1.781716e+01 3.367626e-02 1.674063e+00 7.870699e-03
    ## [31] 4.884394e-01 9.439332e-02 4.313629e-01 4.707788e-02 8.612272e-01
    ## [36] 1.157240e-01 1.790288e-01
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 4 candidate significont markers in 1 chromosome "

    ## [1] "select 3 candidate significont markers in 2 chromosome "

    ## [1] "select 2 candidate significont markers in 3 chromosome "

    ## [1] "select 8 candidate significont markers in 4 chromosome "

    ## [1] "select 3 candidate significont markers in 5 chromosome "

    ## [1] "select 7 candidate significont markers in 6 chromosome "

    ## [1] "select 1 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 37  8
    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 3 model in all."
    ## [1] "FarmCPU"
    ## [1] "Processing trait: V3"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "FarmCPU"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  189  individuals in phenotype data."
    ## [1] "There are  1526  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 189 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 189   4
    ## [1] "GAPIT.IC accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT.SS in process..."
    ## [1] "GAPIT will be into GWAS approach..."
    ## [1] "FarmCPU"
    ## [1] "The GAPIT would go into Bus..."
    ## [1] "--------------------- Welcome to FarmCPU ----------------------------"
    ## [1] "FarmCPU Started..."
    ## [1] "Current loop: 1 out of maximum of 10"
    ## [1] "seqQTN"
    ## NULL
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 3
    ## [1] "Current loop: 2 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ## [1] 748  74 277 396
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 7
    ## [1] "Current loop: 3 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ## [1]  748  816   74  396  277   71   41 1339
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 11
    ## [1] "Current loop: 4 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1]  748  816   74  396   71 1339  277   63 1441   41
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 13
    ## [1] "Current loop: 5 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1]  748  816   74  396   71  277 1339   63 1441   41
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 13
    ## [1] "**********FarmCPU ACCOMPLISHED SUCCESSFULLY**********"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 17.2804361 17.9378160 11.4213165 16.7089601 17.4941438 17.6049289  0.3821849
    ## [1] "Creating marker p-value, MAF, estimated effect, PVE 3 plot..."

    ## [1] "FarmCPU has been done succeedly!!"
    ## [1]   71   74  277  396  748  816 1339
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 2 candidate significont markers in 1 chromosome "

    ## [1] "select 1 candidate significont markers in 2 chromosome "

    ## [1] "select 1 candidate significont markers in 3 chromosome "

    ## [1] "select 0 candidate significont markers in 4 chromosome "

    ## [1] "select 2 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 1 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 7 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 4 model in all."
    ## [1] "BLINK"
    ## [1] "Processing trait: V3"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "BLINK"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  189  individuals in phenotype data."
    ## [1] "There are  1526  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 189 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 189   4
    ## [1] "GAPIT.IC accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT.SS in process..."
    ## [1] "GAPIT will be into GWAS approach..."
    ## [1] "BLINK"
    ## [1] "The GAPIT would go into Bus..."
    ## [1] "----------------------Welcome to Blink----------------------"
    ## [1] "----------------------Iteration: 1 ----------------------"
    ## [1] "seqQTN:"
    ## NULL
    ## [1] "----------------------Iteration: 2 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 2
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 2
    ## [1] "seqQTN:"
    ## [1] 748 816
    ## [1] "----------------------Iteration: 3 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 33
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 7
    ## [1] "seqQTN:"
    ## [1]  748  816   74 1339  116   71  534
    ## [1] "----------------------Iteration: 4 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 8
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 8
    ## [1] "seqQTN:"
    ## [1]  816  748   74   71  116 1339  534  101
    ## [1] "----------------------Iteration: 5 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 9
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 9
    ## [1] "seqQTN:"
    ## [1]  816  748   74   71  116 1339  534  101   63
    ## [1] "----------------------Iteration: 6 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 9
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 9
    ## [1] "seqQTN:"
    ## [1]  816  748   74   71  116 1339  534  101   63
    ## [1] "LD.time(sec):"
    ## [1] 0 0 0 0 0 0
    ## [1] "BIC.time(sec):"
    ## [1] 0.00 0.00 0.00 0.00 0.02 0.00
    ## [1] "GLM.time(sec):"
    ## [1] 1.01 1.03 1.05 0.97 1.11 1.16
    ## [1] "-------------Blink finished successfully in 9.69 seconds!-----------------"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 17.0518051 17.2544997  3.2074873 16.0163687 11.3485214 16.5420573 17.1887545
    ## [8]  0.3357551
    ## [1] "Creating marker p-value, MAF, estimated effect, PVE 3 plot..."

    ## [1] "BLINK R is done !!!!!"
    ## [1]   71   74  101  116  534  748  816 1339
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 4 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 0 candidate significont markers in 3 chromosome "

    ## [1] "select 1 candidate significont markers in 4 chromosome "

    ## [1] "select 2 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 1 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 8 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "Reading GWAS result with GLM.V3"
    ## [1] "Reading GWAS result with SUPER.V3"
    ## [1] "Reading GWAS result with FarmCPU.V3"
    ## [1] "Reading GWAS result with BLINK.V3"

    ## [1] "GAPIT.Association.Manhattans has done !!!"
    ## [1] "GAPIT has output Multiple Manhattan figure with Symphysic type!!!"
    ## [1] "GLM"     "SUPER"   "FarmCPU" "BLINK"  
    ## [1] "Reading GWAS result with GLM.V3"
    ## [1] "Reading GWAS result with SUPER.V3"
    ## [1] "Reading GWAS result with FarmCPU.V3"
    ## [1] "Reading GWAS result with BLINK.V3"

    ## [1] "GAPIT.Association.Manhattans has done !!!"
    ## [1] "GAPIT has output Multiple Manhattan figures with Wide and High types!!!"

    ## [1] "Multracks_QQ Plotting GLM.V3..."

    ## [1] "Multracks_QQ Plotting SUPER.V3..."

    ## [1] "Multracks_QQ Plotting FarmCPU.V3..."

    ## [1] "Multracks_QQ Plotting BLINK.V3..."

    ## [1] "Multiple QQ plot has been finished!"
    ## [1] "GAPIT has output Multiple Manhattan and QQ figures with Circle types!!!"
    ## [1] "GAPIT has done all analysis!!!"
    ## [1] "Please find your all results in :"
    ## [1] "C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation"
