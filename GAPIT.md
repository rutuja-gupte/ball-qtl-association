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
to the requirements of this function. The main assumption being made
here is that chromosomes are numbers and can be converted to integers

*Important note:* GAPIT does support missing values

``` r
gapit_files <- function(vcf.name, pheno.name, f1 = "gapit1.txt", f2 = "gapit2.hmp.txt"){
  vcf <- read.vcfR(vcf.name)
  traits <- read.table(pheno.name)
  gt <- extract.gt(vcf, element="GT")
  
  # Phenotype data
  t1 <- data.frame(var1 = colnames(gt))
  t1 <- left_join(t1, traits, by = c("var1"="V1"))
  t1 <- drop_na(t1)
  t1 <- data.frame(t1) %>% mutate(
    across(everything(), ~replace_na(.x, 0))
  )
  write.table(t1, f1, sep="\t", row.names=FALSE, quote=FALSE)
  
  # Genotype data
  hapmap <- vcfR2hapmap(vcf)
  write.table(hapmap, 
              file = f2,
              sep = "\t", 
              row.names = FALSE,
              col.names = FALSE)
}
```

``` r
run_gapit <- function(wd, f1 = "gapit1.txt", f2 = "gapit2.hmp.txt"){
  myY <- read.table(f1, head = TRUE)
  myG <- read.delim(f2, head = FALSE)
  
  setwd(wd)
  
  setwd("gapitGLM/")
  myGAPITGLM=GAPIT(
    Y=myY[,],
    G=myG,
    model="GLM",
    group.from = 1,
    group.to = 1,
    )
  
  setwd(wd)
  
  setwd("gapitMLMM/")
  myGAPITGLM=GAPIT(
    Y=myY[,], 
    G=myG,
    model="MLMM",
  )
  
  setwd(wd)

  setwd("gapitFarmCPU/")
  myGAPITGLM=GAPIT(
    Y=myY[,], 
    G=myG,
    model="FarmCPU",
  )
  
  setwd(wd)
  
  setwd("gapitBLINK/")
  
  myGAPITGLM=GAPIT(
    Y=myY[,], 
    G=myG,
    model="BLINK",
  )
  
  setwd(wd)
}
```

``` r
wd <- getwd()
gapit_files("processed.vcf.gz", "traits.txt")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 10
    ##   header_line: 11
    ##   variant count: 2828
    ##   column count: 280
    ## Meta line 10 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 2828
    ##   Character matrix gt cols: 280
    ##   skip: 0
    ##   nrows: 2828
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant: 2828
    ## All variants processed

``` r
run_gapit(wd)
```

    ## [1] "--------------------- Welcome to GAPIT ----------------------------"
    ## [1] "GLM"
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 1 model in all."
    ## [1] "GLM"
    ## [1] "GAPIT.DP in process..."
    ## [1] "Converting genotype..."
    ## [1] "Converting HapMap format to numerical under model of Middle"
    ## [1] "Perform numericalization"
    ## [1] "Succesfuly finished converting HapMap which has bits of 2"
    ## [1] "Converting genotype done."
    ## [1] "GAPIT will filter marker with MAF setting !!"
    ## [1] "The markers will be filtered by SNP.MAF: 0"
    ## maf_index
    ## TRUE 
    ## 2828 
    ## [1] "Calculating kinship..."
    ## [1] "Number of individuals and SNPs are  271  and  2828"
    ## [1] "Calculating ZHANG relationship defined by Zhiwu Zhang..."
    ## [1] "substracting mean..."
    ## [1] "Getting X'X..."
    ## [1] "Adjusting..."
    ## [1] "Adjustment by the minimum diagonal"
    ## [1] "Calculating kinship with Zhang method: done"
    ## [1] "kinship calculated"
    ## [1] "Creating heat map for kinship..."

    ## [1] "Kinship heat map created"
    ## [1] "Adding IDs to kinship..."
    ## [1] "Writing kinship to file..."
    ## [1] "Kinship save as file"
    ## [1] "Kinship created!"
    ## [1] "Filting marker for GAPIT.Genotype.View function ..."

    ## [1] "GAPIT.Genotype.View . pdfs generate.successfully!"
    ## NULL
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V2"
    ## [1] "GAPIT.Phenotype.View in press..."

    ## [1] "GAPIT.Phenotype.View output pdf has been generated successfully!"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "Zhang"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  271  individuals in phenotype data."
    ## [1] "There are  2828  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 271 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 271   2
    ## [1] "GAPIT.IC accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT.SS in process..."
    ## [1] "GAPIT will be into GWAS approach..."
    ## [1] "Zhang"
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
    ## [1] "1 of 1 -- Vg= 0 VE= 394.7591 -2LL= 2381.36   Clustering= average   Group number= 1   Group kinship= Mean"
    ## [1] "---------------------Sandwich bottom bun-------------------------------"
    ## [1] "--------------------Final results presentations------------------------"
    ## [1] "Generating summary"
    ## [1] "Genomic Breeding Values (GBV) ..."
    ## [1] "Writing GBV and Acc..."
    ## [1] "GBV and accuracy distribution..."
    ## [1] "Compression portfolios..."
    ## [1] "Compression Visualization done"
    ##      Type   Cluster   Group REML               VA  VE                
    ## [1,] "Mean" "average" "1"   "2381.36122847456" "0" "394.759070870578"
    ## [1] "p3d objects transfered"
    ## [1] "Merge BLUP and BLUE"
    ## [1] "GAPIT before BLUP and BLUE"
    ## [1] "GAPIT after BLUP and BLUE"
    ## [1] "Exporting BLUP and Pred"
    ## [1] "GLM.V2 has been analyzed successfully!"
    ## [1] "The results are saved in the directory of  C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitGLM"
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

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ##  [1]  2.262314  6.975510  1.573541  6.493246  0.000000 13.913383  0.000000
    ##  [8]  3.628309  1.368232  6.431126
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
    ## [1] "select 5 candidate significont markers in 1 chromosome "

    ## [1] "select 1 candidate significont markers in 2 chromosome "

    ## [1] "select 1 candidate significont markers in 3 chromosome "

    ## [1] "select 0 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 1 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 1 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] TRUE
    ## [1] 10  8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "Reading GWAS result with GLM.V2"

    ## [1] "GAPIT.Association.Manhattans has done !!!"
    ## [1] "GAPIT has output Multiple Manhattan figure with Symphysic type!!!"
    ## [1] "GAPIT has done all analysis!!!"
    ## [1] "Please find your all results in :"
    ## [1] "C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitGLM"
    ## [1] "--------------------- Welcome to GAPIT ----------------------------"
    ## [1] "MLMM"
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 1 model in all."
    ## [1] "MLMM"
    ## [1] "GAPIT.DP in process..."
    ## [1] "Converting genotype..."
    ## [1] "Converting HapMap format to numerical under model of Middle"
    ## [1] "Perform numericalization"
    ## [1] "Succesfuly finished converting HapMap which has bits of 2"
    ## [1] "Converting genotype done."
    ## [1] "GAPIT will filter marker with MAF setting !!"
    ## [1] "The markers will be filtered by SNP.MAF: 0"
    ## maf_index
    ## TRUE 
    ## 2828 
    ## [1] "Filting marker for GAPIT.Genotype.View function ..."

    ## [1] "GAPIT.Genotype.View . pdfs generate.successfully!"
    ## NULL
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V2"
    ## [1] "GAPIT.Phenotype.View in press..."

    ## [1] "GAPIT.Phenotype.View output pdf has been generated successfully!"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "MLMM"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  271  individuals in phenotype data."
    ## [1] "There are  2828  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 271 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 271   2
    ## [1] "GAPIT.IC accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT.SS in process..."
    ## [1] "GAPIT will be into GWAS approach..."
    ## [1] "MLMM"
    ## [1] "The GAPIT would go into Bus..."
    ## [1] "GWAS by MLMM method !!"
    ## [1] "Calculating kinship with VanRaden method..."
    ## [1] "substracting P..."
    ## [1] "Getting X'X..."
    ## [1] "Adjusting..."
    ## [1] "Calculating kinship with VanRaden method: done"
    ## null model done! pseudo-h= 0.668 
    ## step 1 done! pseudo-h= 0.757 
    ## step  2  done! pseudo-h= 0.763 
    ## step  3  done! pseudo-h= 0.777 
    ## step  4  done! pseudo-h= 0.808 
    ## step  5  done! pseudo-h= 0.757 
    ## step  6  done! pseudo-h= 0.735 
    ## step  7  done! pseudo-h= 0.759 
    ## step  8  done! pseudo-h= 0.745 
    ## step  9  done! pseudo-h= 0.801 
    ## backward analysis 
    ## creating output 
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] TRUE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 76.44275
    ## [1] 2137
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

    ## [1] "select 0 candidate significont markers in 4 chromosome "

    ## [1] "select 0 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 1 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 0 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] TRUE
    ## [1] 1 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "Reading GWAS result with MLMM.V2"

    ## [1] "GAPIT.Association.Manhattans has done !!!"
    ## [1] "GAPIT has output Multiple Manhattan figure with Symphysic type!!!"
    ## [1] "GAPIT has done all analysis!!!"
    ## [1] "Please find your all results in :"
    ## [1] "C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitMLMM"
    ## [1] "--------------------- Welcome to GAPIT ----------------------------"
    ## [1] "FarmCPU"
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 1 model in all."
    ## [1] "FarmCPU"
    ## [1] "GAPIT.DP in process..."
    ## [1] "Converting genotype..."
    ## [1] "Converting HapMap format to numerical under model of Middle"
    ## [1] "Perform numericalization"
    ## [1] "Succesfuly finished converting HapMap which has bits of 2"
    ## [1] "Converting genotype done."
    ## [1] "GAPIT will filter marker with MAF setting !!"
    ## [1] "The markers will be filtered by SNP.MAF: 0"
    ## maf_index
    ## TRUE 
    ## 2828 
    ## [1] "Filting marker for GAPIT.Genotype.View function ..."

    ## [1] "GAPIT.Genotype.View . pdfs generate.successfully!"
    ## NULL
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V2"
    ## [1] "GAPIT.Phenotype.View in press..."

    ## [1] "GAPIT.Phenotype.View output pdf has been generated successfully!"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "FarmCPU"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  271  individuals in phenotype data."
    ## [1] "There are  2828  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 271 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 271   2
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
    ## [1] 0
    ## [1] "Current loop: 2 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1]  214  411  265 1089 2451 1949  808  180 1654 2142 1845
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 11
    ## [1] "Current loop: 3 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1]  411  327 1394 1643  214  529  974  436  242 2726  495 1089  808 1654  180
    ## [16] 2142 1949 2451  265 1845
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 20
    ## [1] "Current loop: 4 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1]  411 1654 1949 2609 1089  220  737 1625 1725 1119 1643 2726  974  180  495
    ## [16]  808  214  327 1394 1845  529  436  242  265 2451 2142
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 26
    ## [1] "Current loop: 5 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1]  411 1643 2726 1654 1949 1089 1725  974  436  214 1050 1625 1394  180  495
    ## [16]  242  737  327  808  529 2609  220 1845  265 1119 2142 2451
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 27
    ## [1] "Current loop: 6 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1] 1643  411 2726 1949 1725 1089 1654 1050  974 1625  436  242 1394  214  495
    ## [16]  327  180  808  737 2609 1845  529  220  265 1119 2142 2451
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 27
    ## [1] "**********FarmCPU ACCOMPLISHED SUCCESSFULLY**********"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 0 0
    ## [1] "FarmCPU has been done succeedly!!"
    ## [1]  411 1643
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 1 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 0 candidate significont markers in 3 chromosome "

    ## [1] "select 0 candidate significont markers in 4 chromosome "

    ## [1] "select 1 candidate significont markers in 5 chromosome "

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
    ## [1] TRUE
    ## [1] 2 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "Reading GWAS result with FarmCPU.V2"

    ## [1] "GAPIT.Association.Manhattans has done !!!"
    ## [1] "GAPIT has output Multiple Manhattan figure with Symphysic type!!!"
    ## [1] "GAPIT has done all analysis!!!"
    ## [1] "Please find your all results in :"
    ## [1] "C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitFarmCPU"
    ## [1] "--------------------- Welcome to GAPIT ----------------------------"
    ## [1] "BLINK"
    ## [1] "--------------------Processing traits----------------------------------"
    ## [1] "Phenotype provided!"
    ## [1] "The 1 model in all."
    ## [1] "BLINK"
    ## [1] "GAPIT.DP in process..."
    ## [1] "Converting genotype..."
    ## [1] "Converting HapMap format to numerical under model of Middle"
    ## [1] "Perform numericalization"
    ## [1] "Succesfuly finished converting HapMap which has bits of 2"
    ## [1] "Converting genotype done."
    ## [1] "GAPIT will filter marker with MAF setting !!"
    ## [1] "The markers will be filtered by SNP.MAF: 0"
    ## maf_index
    ## TRUE 
    ## 2828 
    ## [1] "Calculating kinship..."
    ## [1] "Number of individuals and SNPs are  271  and  2828"
    ## [1] "Kinship created!"
    ## [1] "Filting marker for GAPIT.Genotype.View function ..."

    ## [1] "GAPIT.Genotype.View . pdfs generate.successfully!"
    ## NULL
    ## [1] "GAPIT.DP accomplished successfully for multiple traits. Results are saved"
    ## [1] "Processing trait: V2"
    ## [1] "GAPIT.Phenotype.View in press..."

    ## [1] "GAPIT.Phenotype.View output pdf has been generated successfully!"
    ## [1] "--------------------Phenotype and Genotype ----------------------------------"
    ## [1] "BLINK"
    ## [1] TRUE
    ## [1] "There are  1  traits in phenotype data."
    ## [1] "There are  271  individuals in phenotype data."
    ## [1] "There are  2828  markers in genotype data."
    ## [1] "Phenotype and Genotype are test OK !!"
    ## [1] "--------------------GAPIT Logical Done----------------------------------"
    ## [1] "GAPIT.IC in process..."
    ## [1] "There is 0 Covarinces."
    ## [1] "There are 271 common individuals in genotype , phenotype and CV files."
    ## [1] "The dimension of total CV is "
    ## [1] 271   2
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
    ## [1] 19
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 16
    ## [1] "seqQTN:"
    ## [1]  214  411  265 1089 2451 1949  808
    ## [1] "----------------------Iteration: 3 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 4
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 4
    ## [1] "seqQTN:"
    ## [1]  411 1643  214 1089  265 2451 1949  808
    ## [1] "----------------------Iteration: 4 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 4
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 4
    ## [1] "seqQTN:"
    ## [1]  411 1643  214 2744 1089  265 2451 1949  808
    ## [1] "----------------------Iteration: 5 ----------------------"
    ## [1] "LD remove is working...."
    ## [1] "Number SNPs for LD remove:"
    ## [1] 5
    ## [1] "Model selection based on BIC is working...."
    ## [1] "Number of SNPs for BIC selection:"
    ## [1] 5
    ## [1] "seqQTN:"
    ## [1]  411 1643 2744  214 1089  265 2451 1949  808
    ## [1] "LD.time(sec):"
    ## [1] 0 0 0 0 0
    ## [1] "BIC.time(sec):"
    ## [1] 0 0 0 0 0
    ## [1] "GLM.time(sec):"
    ## [1] 0.83 1.28 1.04 1.17 0.98
    ## [1] "-------------Blink finished successfully in 8.22 seconds!-----------------"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 13.257563  7.909256 12.929304
    ## [1] "BLINK R is done !!!!!"
    ## [1]  411 1643 2744
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 1 candidate significont markers in 1 chromosome "

    ## [1] "select 0 candidate significont markers in 2 chromosome "

    ## [1] "select 0 candidate significont markers in 3 chromosome "

    ## [1] "select 0 candidate significont markers in 4 chromosome "

    ## [1] "select 1 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 0 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "select 0 candidate significont markers in 9 chromosome "

    ## [1] "select 1 candidate significont markers in 10 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] TRUE
    ## [1] 3 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "Reading GWAS result with BLINK.V2"

    ## [1] "GAPIT.Association.Manhattans has done !!!"
    ## [1] "GAPIT has output Multiple Manhattan figure with Symphysic type!!!"
    ## [1] "GAPIT has done all analysis!!!"
    ## [1] "Please find your all results in :"
    ## [1] "C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitBLINK"
