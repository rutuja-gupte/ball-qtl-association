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
run_gapit <- function(f1 = "gapit1.txt", f2 = "gapit2.hmp.txt"){
  myY <- read.table(f1, head = TRUE)
  myG <- read.delim(f2, head = FALSE) 
  
  # setwd("C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitGLM/")
  # 
  # myGAPITGLM=GAPIT(
  #   Y=myY[,], #first column is ID
  #   G=myG,
  #   PCA.total=5,
  #   model="GLM",
  #   group.from = 1,
  #   group.to = 1,
  #   )
  
  # setwd("C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitMLM/")
  # 
  # myGAPITGLM=GAPIT(
  #   Y=myY[,], #first column is ID
  #   G=myG,
  #   PCA.total=5,
  #   model="MLM",
  # )
  
  # setwd("C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitCMLM/")
  # 
  # myGAPITGLM=GAPIT(
  #   Y=myY[,], #first column is ID
  #   G=myG,
  #   PCA.total=5,
  #   model="CMLM",
  # )
  
  setwd("C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitMLMM/")

  myGAPITGLM=GAPIT(
    Y=myY[,], #first column is ID
    G=myG,
    PCA.total=3,
    model="MLMM",
  )
  
  setwd("C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitFarmCPU/")

  myGAPITGLM=GAPIT(
    Y=myY[,], #first column is ID
    G=myG,
    PCA.total=3,
    model="FarmCPU",
  )
  
  setwd("C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitBLINK/")
  
  myGAPITGLM=GAPIT(
    Y=myY[,], #first column is ID
    G=myG,
    PCA.total=3,
    model="BLINK",
  )
}

gapit_files <- function(vcf.name, pheno.name, f1 = "gapit1.txt", f2 = "gapit2.hmp.txt"){
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
setwd("C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation")
gapit_files("gt5382/gt5382_processed.vcf.gz", "gt5382/gt5382_traits.txt", "gt5382/gapit1.txt", "gt5382/gapit2.hmp.txt")
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

``` r
setwd("C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation")
run_gapit("gt5382/gapit1.txt", "gt5382/gapit2.hmp.txt")
```

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
    ## 1526 
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
    ## [1] "MLMM"
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
    ## [1] "MLMM"
    ## [1] "The GAPIT would go into Bus..."
    ## [1] "GWAS by MLMM method !!"
    ## [1] "Calculating kinship with VanRaden method..."
    ## [1] "substracting P..."
    ## [1] "Getting X'X..."
    ## [1] "Adjusting..."
    ## [1] "Calculating kinship with VanRaden method: done"
    ## null model done! pseudo-h= 0.062 
    ## step 1 done! pseudo-h= 0.048 
    ## step  2  done! pseudo-h= 0 
    ## step  3  done! pseudo-h= 0 
    ## backward analysis 
    ## creating output

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
    ## [1] 33.4898930 32.1701906 30.7096044  0.2318068  0.3995438
    ## [1] "Creating marker p-value, MAF, estimated effect, PVE 3 plot..."

    ## [1]   74  748  816 1338 1339
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

    ## [1] "select 2 candidate significont markers in 5 chromosome "

    ## [1] "select 0 candidate significont markers in 6 chromosome "

    ## [1] "select 1 candidate significont markers in 7 chromosome "

    ## [1] "select 0 candidate significont markers in 8 chromosome "

    ## [1] "manhattan plot on chromosome finished"
    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Association table..."
    ## [1] "Joining tvalue and stderr"
    ## [1] "GAPIT Phenotype distribution with significant markers in process..."
    ## [1] TRUE
    ## [1] 5 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "Reading GWAS result with MLMM.V3"

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
    ## 1526 
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
    ## [1] 748  74 277 116
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 7
    ## [1] "Current loop: 3 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ## [1]  748  816   74  116  277   71   41 1339
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 11
    ## [1] "Current loop: 4 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1]  748  816   74  116   71 1339  277   63 1441   41
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 13
    ## [1] "Current loop: 5 out of maximum of 10"
    ## [1] "optimizing possible QTNs..."
    ## [1] "seqQTN"
    ##  [1]  748  816   74  116   71  277 1339   63 1441   41
    ## [1] "scanning..."
    ## [1] "number of covariates in current loop is:"
    ## [1] 13
    ## [1] "**********FarmCPU ACCOMPLISHED SUCCESSFULLY**********"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 0.000000 0.000000 2.008537 3.982478
    ## [1] "FarmCPU has been done succeedly!!"
    ## [1]   71   74  116  277  748  816 1339
    ## [1] "GAPIT.ID in process..."
    ## [1] "Filtering SNPs with MAF..."
    ## [1] "Calculating FDR..."
    ## [1] "QQ plot..."

    ## [1] "Manhattan plot (Genomewise)..."

    ## [1] "GAPIT.Manhattan accomplished successfully!zw"
    ## [1] "Manhattan plot (Chromosomewise)..."
    ## [1] "select 3 candidate significont markers in 1 chromosome "

    ## [1] "select 1 candidate significont markers in 2 chromosome "

    ## [1] "select 0 candidate significont markers in 3 chromosome "

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
    ## [1] TRUE
    ## [1] 7 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "Reading GWAS result with FarmCPU.V3"

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
    ## 1526 
    ## [1] "Calculating kinship..."
    ## [1] "Number of individuals and SNPs are  189  and  1526"
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
    ## [1] 0.00 0.01 0.00 0.00 0.00 0.00
    ## [1] "GLM.time(sec):"
    ## [1] 1.03 1.00 0.97 0.81 0.83 0.80
    ## [1] "-------------Blink finished successfully in 8.49 seconds!-----------------"
    ## [1] "Calculating Original GWAS result..."
    ## [1] "GAPIT.RandomModel beginning..."
    ## [1] FALSE
    ## [1] "Candidate Genes could Phenotype_Variance_Explained(%) :"
    ## [1] 17.0669833 17.2345861  3.1999939 16.0224033 11.3681144 16.5201575 17.1971080
    ## [8]  0.3356582
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
    ## [1] TRUE
    ## [1] 8 8

    ## [1] "GAPIT.ID accomplished successfully for multiple traits. Results are saved"
    ## [1] "GAPIT accomplished successfully for multiple traits. Result are saved"
    ## [1] "Reading GWAS result with BLINK.V3"

    ## [1] "GAPIT.Association.Manhattans has done !!!"
    ## [1] "GAPIT has output Multiple Manhattan figure with Symphysic type!!!"
    ## [1] "GAPIT has done all analysis!!!"
    ## [1] "Please find your all results in :"
    ## [1] "C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation/gapitBLINK"

``` r
setwd("C:/Users/RGupte/OneDrive - Ball Horticultural Company/Association_Documentation")
```
