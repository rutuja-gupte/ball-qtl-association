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

*Important note:* GAPIT does support missing values

``` r
run_gapit <- function(vcf.name, pheno.name, f1 = "gapit1.csv", f2 = "gapit2.csv", f3 = "gapit3.csv"){
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
      # "SUPER", 
      "FarmCPU", 
      "BLINK"),
    kinship.algorithm = "VanRaden",
    kinship.cluster = "average",
    Multiple_analysis=TRUE,)
}
```

``` r
run_gapit("sample/processed_vcf.vcf.gz", "sample/traits.txt", "sample/gapit1.txt", "sample/gapit2.txt", "sample/gapit3.txt" )
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
    ## [1] "--------------------- Welcome to GAPIT ----------------------------"
    ## [1] "GLM"     "FarmCPU" "BLINK"  
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 1 model in all."
    ## [1] "GLM"
    ## [1] "GAPIT.DP in process..."
    ## [1] "GAPIT will filter marker with MAF setting !!"
    ## [1] "The markers will be filtered by SNP.MAF: 0"
    ## maf_index
    ##  TRUE 
    ## 14643 
    ## [1] "Calculating kinship..."
    ## [1] "Number of individuals and SNPs are  16  and  14643"
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
    ## [1] 16  4
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V3"
    ## [1] "GAPIT.Phenotype.View in press..."

    ## [1] "GAPIT.Phenotype.View output pdf has been generated successfully!"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "VanRaden"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  16  individuals in phenotype data."
    ## [1] "There are  14643  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 16 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 16  4
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
    ## [1] "Genotype file: 1, SNP: 2000 "
    ## [1] "Genotype file: 1, SNP: 3000 "
    ## [1] "Genotype file: 1, SNP: 4000 "
    ## [1] "Genotype file: 1, SNP: 5000 "
    ## [1] "Genotype file: 1, SNP: 6000 "
    ## [1] "Genotype file: 1, SNP: 7000 "
    ## [1] "Genotype file: 1, SNP: 8000 "
    ## [1] "Genotype file: 1, SNP: 9000 "
    ## [1] "Genotype file: 1, SNP: 10000 "
    ## [1] "Genotype file: 1, SNP: 11000 "
    ## [1] "Genotype file: 1, SNP: 12000 "
    ## [1] "Genotype file: 1, SNP: 13000 "
    ## [1] "Genotype file: 1, SNP: 14000 "
    ## [1] "1 of 1 -- Vg= 0 VE= 39.0065 -2LL= 82.02   Clustering= average   Group number= 1   Group kinship= Mean"
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "--------------------Final results presentations------------------------"
    ## [1] "Generating summary"
    ## [1] "Genomic Breeding Values (GBV) ..."
    ## [1] "Writing GBV and Acc..."
    ## [1] "GBV and accuracy distribution..."
    ## [1] "Compression portfolios..."
    ## [1] "Compression Visualization done"
    ##      Type   Cluster   Group REML               VA  VE                
    ## [1,] "Mean" "average" "1"   "82.0192785911711" "0" "39.0065461816148"
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
    ## [1] "There is no significant marker for VE !!"
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

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 0 8
    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 2 model in all."
    ## [1] "FarmCPU"
    ## [1] "Processing trait: V3"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "FarmCPU"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  16  individuals in phenotype data."
    ## [1] "There are  14643  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 16 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 16  4
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
    ## [1] "Top snps have little effect, set seqQTN to NULL!"
    ## [1] "**********FarmCPU ACCOMPLISHED SUCCESSFULLY**********"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] "There is no significant marker for VE !!"
    ## [1] "FarmCPU has been done succeedly!!"
    ## integer(0)
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 0 candidate significont markers in 1 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 0 8
    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 3 model in all."
    ## [1] "BLINK"
    ## [1] "Processing trait: V3"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "BLINK"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  16  individuals in phenotype data."
    ## [1] "There are  14643  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 16 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 16  4
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
    ## [1] "Top snps have little effect, set seqQTN to NULL!"
    ## [1] "seqQTN is:,stop here"
    ## [1] "LD.time(sec):"
    ## [1] 0 0
    ## [1] "BIC.time(sec):"
    ## [1] 0 0
    ## [1] "GLM.time(sec):"
    ## [1] 1.37 0.00
    ## [1] "-------------Blink finished successfully in 2.63 seconds!-----------------"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] "There is no significant marker for VE !!"
    ## [1] "BLINK R is done !!!!!"
    ## integer(0)
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 0 candidate significont markers in 1 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 0 8
    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "Reading GWAS result with GLM.V3"
    ## [1] "Reading GWAS result with FarmCPU.V3"
    ## [1] "Reading GWAS result with BLINK.V3"

    ## [1] "GAPIT.Association.Manhattans has done !!!"
    ## [1] "GAPIT has output Multiple Manhattan figure with Symphysic type!!!"
    ## [1] "GLM"     "FarmCPU" "BLINK"  
    ## [1] "Reading GWAS result with GLM.V3"
    ## [1] "Reading GWAS result with FarmCPU.V3"
    ## [1] "Reading GWAS result with BLINK.V3"

    ## [1] "GAPIT.Association.Manhattans has done !!!"
    ## [1] "GAPIT has output Multiple Manhattan figures with Wide and High types!!!"

    ## [1] "Multracks_QQ Plotting GLM.V3..."

    ## [1] "Multracks_QQ Plotting FarmCPU.V3..."

    ## [1] "Multracks_QQ Plotting BLINK.V3..."

    ## [1] "Multiple QQ plot has been finished!"
    ## [1] "GAPIT has output Multiple Manhattan and QQ figures with Circle types!!!"
    ## [1] "GAPIT has done all analysis!!!"
    ## [1] "Please find your all results in :"
    ## [1] "C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation"
