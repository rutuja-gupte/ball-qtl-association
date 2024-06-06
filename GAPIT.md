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
  write.table(t1, f1, sep="\t", row.names=FALSE)
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
  write.table(t2, f2, sep="\t", row.names=FALSE)
  myGD <- read.table(f2, head = TRUE)
  
  # Genetic map
  t3 <- vcf@fix[,c("POS", "CHROM")]
  t3 <- data.frame(t3)
  t3$Markers <- rownames(gt)
  t3$CHROM <- as.numeric(t3$CHROM)
  t3$POS <- as.numeric(t3$POS)
  t3 <- t3 %>% relocate(Markers, CHROM, POS)
  write.table(t3, f3, sep="\t", row.names=FALSE)
  myGM <- read.table(f3, head = TRUE)
  
  myGAPIT=GAPIT(
    Y=myY[,], #first column is ID
    GD=myGD,
    GM=myGM,
    PCA.total=3,
    model=c("GLM", "SUPER", "FarmCPU", "BLINK"),
    kinship.algorithm = "VanRaden",
    kinship.cluster = "average",
    Multiple_analysis=TRUE,)
}
```

``` r
run_gapit("TASSEL_samples/processed_vcf.vcf", "TASSEL_samples/mdp_traits2.txt", "TASSEL_samples/gapit1.csv", "TASSEL_samples/gapit2.csv", "TASSEL_samples/gapit3.csv" )
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 10
    ##   header_line: 11
    ##   variant count: 510
    ##   column count: 268
    ## Meta line 10 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 510
    ##   Character matrix gt cols: 268
    ##   skip: 0
    ##   nrows: 510
    ##   row_num: 0
    ## Processed variant: 510
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
    ##  510 
    ## [1] "Calculating kinship..."
    ## [1] "Number of individuals and SNPs are  259  and  510"
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
    ## [1] 259   4
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V3"
    ## [1] "GAPIT.Phenotype.View in press..."

    ## [1] "GAPIT.Phenotype.View output pdf has been generated successfully!"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "VanRaden"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ## [1] "1 of 1 -- Vg= 0 VE= 310.4195 -2LL= 1933.56   Clustering= average   Group number= 1   Group kinship= Mean"
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "--------------------Final results presentations------------------------"
    ## [1] "Generating summary"
    ## [1] "Genomic Breeding Values (GBV) ..."
    ## [1] "Writing GBV and Acc..."
    ## [1] "GBV and accuracy distribution..."
    ## [1] "Compression portfolios..."
    ## [1] "Compression Visualization done"
    ##      Type   Cluster   Group REML               VA  VE               
    ## [1,] "Mean" "average" "1"   "1933.55538679059" "0" "310.41951370786"
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
    ## [1] 6.845510e+00 4.461989e-06 4.121230e+01 1.106401e+01
    ## [1] "GAPIT.ID in process..."
    ## [1] "GAPIT.Compression.Visualization"
    ## [1] "Pie chart"

    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 1 candidate significont markers in 1 chromosome "

    ## [1] "select 1 candidate significont markers in 2 chromosome "

    ## [1] "select 1 candidate significont markers in 3 chromosome "

    ## [1] "select 0 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 1 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 4 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V4"
    ## [1] "GAPIT.Phenotype.View in press..."

    ## [1] "GAPIT.Phenotype.View output pdf has been generated successfully!"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "VanRaden"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ## [1] "1 of 1 -- Vg= 0 VE= 19.3587 -2LL= 1309.23   Clustering= average   Group number= 1   Group kinship= Mean"
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "--------------------Final results presentations------------------------"
    ## [1] "Generating summary"
    ## [1] "Genomic Breeding Values (GBV) ..."
    ## [1] "Writing GBV and Acc..."
    ## [1] "GBV and accuracy distribution..."
    ## [1] "Compression portfolios..."
    ## [1] "Compression Visualization done"
    ##      Type   Cluster   Group REML               VA  VE                
    ## [1,] "Mean" "average" "1"   "1309.22958169097" "0" "19.3587264406748"
    ## [1] "p3d objects transfered"
    ## [1] "Merge BLUP and BLUE"
    ## [1] "GAPIT before BLUP and BLUE"
    ## [1] "GAPIT after BLUP and BLUE"
    ## [1] "Exporting BLUP and Pred"
    ## [1] "GLM.V4 has been analyzed successfully!"
    ## [1] "The results are saved in the directory of  C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation"
    ## [1] "before ending GAPIT.Main"
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 1.672436e+01 2.372402e+00 5.879391e-07 2.492373e+00 1.580459e+01
    ## [6] 3.103793e+01 5.193979e+00
    ## [1] "Creating marker p-value, MAF, estimated effect, PVE 3 plot..."

    ## [1] "GAPIT.ID in process..."
    ## [1] "GAPIT.Compression.Visualization"
    ## [1] "Pie chart"

    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 3 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 1 candidate significont markers in 3 chromosome "

    ## [1] "select 1 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 1 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 7 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V5"
    ## [1] "GAPIT.Phenotype.View in press..."

    ## [1] "GAPIT.Phenotype.View output pdf has been generated successfully!"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "VanRaden"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ## [1] "1 of 1 -- Vg= 0 VE= 15.8975 -2LL= 1264.91   Clustering= average   Group number= 1   Group kinship= Mean"
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "--------------------Final results presentations------------------------"
    ## [1] "Generating summary"
    ## [1] "Genomic Breeding Values (GBV) ..."
    ## [1] "Writing GBV and Acc..."
    ## [1] "GBV and accuracy distribution..."
    ## [1] "Compression portfolios..."
    ## [1] "Compression Visualization done"
    ##      Type   Cluster   Group REML               VA  VE                
    ## [1,] "Mean" "average" "1"   "1264.90923535469" "0" "15.8975336304348"
    ## [1] "p3d objects transfered"
    ## [1] "Merge BLUP and BLUE"
    ## [1] "GAPIT before BLUP and BLUE"
    ## [1] "GAPIT after BLUP and BLUE"
    ## [1] "Exporting BLUP and Pred"
    ## [1] "GLM.V5 has been analyzed successfully!"
    ## [1] "The results are saved in the directory of  C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation"
    ## [1] "before ending GAPIT.Main"
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 3.678950e-09 5.603609e+00 6.767971e+01
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

    ## [1] "select 1 candidate significont markers in 2 chromosome "

    ## [1] "select 0 candidate significont markers in 3 chromosome "

    ## [1] "select 2 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 3 8

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
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ##  510 
    ## NULL
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V3"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "Zhang"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ## [1] "1 of 1 -- Vg= 122.5599 VE= 112.6871 -2LL= 1893.45   Clustering= average   Group number= 259   Group kinship= Mean"
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
    ## [1,] "Mean" "average" "259" "1893.44921464192" "122.559864499093"
    ##      VE                
    ## [1,] "112.687056620943"
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
    ##  [1]  79  84 187 195 196 263 288 429 442 448
    ##    Type Cluster   Group    REML      VA      VE 
    ##       1   10000      10      NA      NA      NA 
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "---------------Sandwich bottom with grilled burger---------------------"
    ## [1] "---------------Sandwich bottom: reload bins ---------------------------"
    ## [1] "SNP: 500 "
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
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 21.58092 30.98765
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 0 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 1 candidate significont markers in 3 chromosome "

    ## [1] "select 0 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 1 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 2 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V4"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "Zhang"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ##  510 
    ## NULL
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V4"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "Zhang"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ## [1] "1 of 1 -- Vg= 6.9558 VE= 7.5504 -2LL= 1267.19   Clustering= average   Group number= 259   Group kinship= Mean"
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
    ## [1,] "Mean" "average" "259" "1267.18980146902" "6.95579738756804"
    ##      VE                
    ## [1,] "7.55036318047671"
    ## [1] "p3d objects transfered"
    ## [1] "Merge BLUP and BLUE"
    ## [1] "GAPIT before BLUP and BLUE"
    ## [1] "GAPIT after BLUP and BLUE"
    ## [1] "MLM.V4 has been analyzed successfully!"
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
    ##  [1]  12  38  73  84 196 199 221 255 429 448
    ##    Type Cluster   Group    REML      VA      VE 
    ##       1   10000      10      NA      NA      NA 
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "---------------Sandwich bottom with grilled burger---------------------"
    ## [1] "---------------Sandwich bottom: reload bins ---------------------------"
    ## [1] "SNP: 500 "
    ## [1] "SUPER saving results..."
    ## [1] "--------------------Final results presentations------------------------"
    ## [1] "Generating summary"
    ## [1] "Genomic Breeding Values (GBV) ..."
    ## [1] "Merge BLUP and BLUE"
    ## [1] "SUPER.V4 has been analyzed successfully!"
    ## [1] "The results are saved in the directory of  C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation"
    ## [1] "before ending GAPIT.Main"
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 4.500120e-12 3.562644e+01 2.850609e+01 6.018517e+00
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 1 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 1 candidate significont markers in 3 chromosome "

    ## [1] "select 1 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 1 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 4 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V5"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "Zhang"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ##  510 
    ## NULL
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V5"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "Zhang"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ## [1] "1 of 1 -- Vg= 5.964 VE= 5.6119 -2LL= 1215.92   Clustering= average   Group number= 259   Group kinship= Mean"
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
    ## [1,] "Mean" "average" "259" "1215.92226612903" "5.96398136934955"
    ##      VE                
    ## [1,] "5.61188330619594"
    ## [1] "p3d objects transfered"
    ## [1] "Merge BLUP and BLUE"
    ## [1] "GAPIT before BLUP and BLUE"
    ## [1] "GAPIT after BLUP and BLUE"
    ## [1] "MLM.V5 has been analyzed successfully!"
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
    ##  [1]  57 122 125 136 252 273 305 401 425 474
    ##    Type Cluster   Group    REML      VA      VE 
    ##       1   10000      10      NA      NA      NA 
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "---------------Sandwich bottom with grilled burger---------------------"
    ## [1] "---------------Sandwich bottom: reload bins ---------------------------"
    ## [1] "SNP: 500 "
    ## [1] "SUPER saving results..."
    ## [1] "--------------------Final results presentations------------------------"
    ## [1] "Generating summary"
    ## [1] "Genomic Breeding Values (GBV) ..."
    ## [1] "Merge BLUP and BLUE"
    ## [1] "SUPER.V5 has been analyzed successfully!"
    ## [1] "The results are saved in the directory of  C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation"
    ## [1] "before ending GAPIT.Main"
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1]  5.60397 67.67918
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 0 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 0 candidate significont markers in 3 chromosome "

    ## [1] "select 2 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

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
    ## [1] "The 3 model in all."
    ## [1] "FarmCPU"
    ## [1] "Processing trait: V3"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "FarmCPU"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ## [1] 200 442  84 319 501  48 288 429 259
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 12
    ## [1] "Current loop: 3 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 442 319 190  79 142 448  84 261 195 395 215 288 259 429 200  48 501
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 20
    ## [1] "Current loop: 4 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 442 190 448  84 429 353 259 215 319 142  79 195 261 288  48 395 200 501
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 21
    ## [1] "Current loop: 5 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 442 190  84 429 448 353  96 259  79 142 215 195 319 261 288  48 395 501 200
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 22
    ## [1] "Current loop: 6 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 190 442 448 429 353  84 259  96 142 215  79 319 195  48 261 288 501 395 200
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 22
    ## [1] "**********FarmCPU ACCOMPLISHED SUCCESSFULLY**********"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 18.97583  6.59696 20.35775 21.35212
    ## [1] "FarmCPU has been done succeedly!!"
    ## [1] 190 429 442 448
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 0 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 1 candidate significont markers in 3 chromosome "

    ## [1] "select 0 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 2 candidate significont markers in 8 chromosome "

    ## [1] "select 1 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 4 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V4"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "FarmCPU"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ##  [1] 200 429  43  84 221  12  38 343 208 319
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 13
    ## [1] "Current loop: 3 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 221 215 343  12 208  79  73 200 319  38 429  84  43
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 16
    ## [1] "Current loop: 4 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 221 215 343  12 429  84  79 208 319 222  38  73 200  43
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 17
    ## [1] "Current loop: 5 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 221 215 343  84 429  12 319  79 222 199 208  38 200  73  43
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 18
    ## [1] "Current loop: 6 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 221  84 215 429 222 319 199  79 343  12  38  73 208 200  43
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 18
    ## [1] "**********FarmCPU ACCOMPLISHED SUCCESSFULLY**********"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1]  8.748984 19.537591 44.060470
    ## [1] "FarmCPU has been done succeedly!!"
    ## [1]  84 215 221
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 1 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 1 candidate significont markers in 3 chromosome "

    ## [1] "select 1 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 3 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V5"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "FarmCPU"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ##  [1] 252 143 429 440 493 401 471  57 200 122
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 13
    ## [1] "Current loop: 3 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 473 122 136 252 487 425 199 493 401 440 471 200  57 429 143
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 18
    ## [1] "Current loop: 4 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 122 252 473 487 136 411 200 493 498 390 440 425  57 401 199 471 429 143
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 21
    ## [1] "Current loop: 5 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 122 252 473 136 487 493 425 440 411 422 390 498  57 200 471 199 401 429 143
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 22
    ## [1] "Current loop: 6 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 122 473 252 136 487 425 422 440 411 493 498 390 200  57 471 401 429 199 143
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 22
    ## [1] "**********FarmCPU ACCOMPLISHED SUCCESSFULLY**********"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1]  5.911176 26.130627 53.640432  1.057617
    ## [1] "FarmCPU has been done succeedly!!"
    ## [1] 122 136 252 473
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 0 candidate significont markers in 1 chromosome "

    ## [1] "select 2 candidate significont markers in 2 chromosome "

    ## [1] "select 0 candidate significont markers in 3 chromosome "

    ## [1] "select 1 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "select 1 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 4 8

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
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ## [1] 4
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 3
    ## [1] "seqQTN:"
    ## [1] 200 442
    ## [1] "----------------------Iteration: 3 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 3
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 3
    ## [1] "seqQTN:"
    ## [1] 200 442
    ## [1] "LD.time(sec):"
    ## [1] 0 0 0
    ## [1] "BIC.time(sec):"
    ## [1] 0.00 0.01 0.00
    ## [1] "GLM.time(sec):"
    ## [1] 0.80 0.75 0.86
    ## [1] "-------------Blink finished successfully in 3.93 seconds!-----------------"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 51.16137 12.84662
    ## [1] "BLINK R is done !!!!!"
    ## [1] 200 442
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 0 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 1 candidate significont markers in 3 chromosome "

    ## [1] "select 0 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 1 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 2 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V4"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "BLINK"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ## [1] 7
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 7
    ## [1] "seqQTN:"
    ## [1] 200 429  43  84 221  12 199
    ## [1] "----------------------Iteration: 3 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 10
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 10
    ## [1] "seqQTN:"
    ## [1] 221 199 200 429  43
    ## [1] "----------------------Iteration: 4 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 10
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 10
    ## [1] "seqQTN:"
    ## [1] 221 199  12 429  79  84  43 196 210
    ## [1] "----------------------Iteration: 5 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 6
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 6
    ## [1] "seqQTN:"
    ## [1] 429 221 199  12  79  84  43 196 210
    ## [1] "LD.time(sec):"
    ## [1] 0.00 0.02 0.00 0.00 0.00
    ## [1] "BIC.time(sec):"
    ## [1] 0.00 0.00 0.00 0.02 0.00
    ## [1] "GLM.time(sec):"
    ## [1] 0.91 0.83 0.81 0.85 0.77
    ## [1] "-------------Blink finished successfully in 6.47 seconds!-----------------"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 46.79476 12.89910
    ## [1] "BLINK R is done !!!!!"
    ## [1] 221 429
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 0 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 0 candidate significont markers in 3 chromosome "

    ## [1] "select 1 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 1 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 2 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V5"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "BLINK"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  229  individuals in phenotype data."
    ## [1] "There are  510  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 229 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 229   4
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
    ## [1] 5
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 5
    ## [1] "seqQTN:"
    ## [1] 252 248 143 429 440
    ## [1] "----------------------Iteration: 3 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 3
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 3
    ## [1] "seqQTN:"
    ## [1] 440 136 401 252
    ## [1] "----------------------Iteration: 4 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 5
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 5
    ## [1] "seqQTN:"
    ## [1] 252 440 136 401 473
    ## [1] "----------------------Iteration: 5 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 8
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 8
    ## [1] "seqQTN:"
    ## [1] 252 440 473 401 135 122 136
    ## [1] "----------------------Iteration: 6 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 7
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 7
    ## [1] "seqQTN:"
    ## [1] 252 473 135 122 440 487 136 401
    ## [1] "----------------------Iteration: 7 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 7
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 7
    ## [1] "seqQTN:"
    ## [1] 252 473 122 487 135 440 136 401
    ## [1] "LD.time(sec):"
    ## [1] 0 0 0 0 0 0 0
    ## [1] "BIC.time(sec):"
    ## [1] 0.00 0.00 0.01 0.00 0.00 0.00 0.00
    ## [1] "GLM.time(sec):"
    ## [1] 0.81 0.78 0.73 0.75 0.78 0.75 0.75
    ## [1] "-------------Blink finished successfully in 8.34 seconds!-----------------"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1]  7.370994 72.658530  1.791658
    ## [1] "BLINK R is done !!!!!"
    ## [1] 122 252 473
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 0 candidate significont markers in 1 chromosome "

    ## [1] "select 1 candidate significont markers in 2 chromosome "

    ## [1] "select 0 candidate significont markers in 3 chromosome "

    ## [1] "select 1 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "select 1 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] FALSE
    ## [1] 3 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "Reading GWAS result with GLM.V3"
    ## [1] "Reading GWAS result with GLM.V4"
    ## [1] "Reading GWAS result with GLM.V5"
    ## [1] "Reading GWAS result with SUPER.V3"
    ## [1] "Reading GWAS result with SUPER.V4"
    ## [1] "Reading GWAS result with SUPER.V5"
    ## [1] "Reading GWAS result with FarmCPU.V3"
    ## [1] "Reading GWAS result with FarmCPU.V4"
    ## [1] "Reading GWAS result with FarmCPU.V5"
    ## [1] "Reading GWAS result with BLINK.V3"
    ## [1] "Reading GWAS result with BLINK.V4"
    ## [1] "Reading GWAS result with BLINK.V5"

    ## [1] "GAPIT.Association.Manhattans has done !!!"
    ## [1] "GAPIT has output Multiple Manhattan figure with Symphysic type!!!"
    ## [1] "GAPIT has done all analysis!!!"
    ## [1] "Please find your all results in :"
    ## [1] "C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation"
