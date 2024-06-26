TASSEL
================
Rutuja Gupte

## Installation

Starting with dependencies

    sudo apt install openjdk-19-jre-headless

Now the actual repository

    git clone https://bitbucket.org/tasseladmin/tassel-5-standalone.git 

cd into the new directory and check if it works

    # Meant to throw an error
    ./run_pipeline.pl

    # Download their tutorial data if you havent already. The repository may not contain the vcf file in the tutorial data.

## File processing

I am assuming that the vcf file has already been preprocessed according
to the requirements of this function. The main assumptions being made
here are:  
1. Chromosomes are numbers and can be converted to integers  
2. There are no missing values  
3. The phenotype file has 2 columns of names and no headers

``` r
traits <- read.table("gt5382/gt5382_traits.txt")
traits <- traits[,-1]
cnames <- colnames(traits)
cnames[1] <- "<Trait>"
colnames(traits) <- cnames
write.table(traits, "gt5382/gt5382_traits_final.txt", quote=FALSE, sep="\t", row.names=FALSE)
```

For TASSEL, I still havent figured out how to properly get around the
covariate thing. I am currently passing the phenotype file as covariate
but that shouldnt be the case. This is the current way I am running it.
I think it should be fine because I have the exclude last trait flag so
I am hoping it just ignores all covariates.

    ./run_pipeline.pl -fork1 -importGuess gt5382_processed.vcf.gz -fork2 -importGuess gt5382_traits_final.txt -fork3 -importGuess gt5382_traits_final.txt -excludeLastTrait -combine5 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endPlugin -export glm_output

*Hopefully better ways to do this*

The following two give the same output

The first one might be the better option

    ./run_pipeline.pl -fork1 -importGuess gt5382_processed.vcf.gz -fork2 -importGuess gt5382_traits_final.txt -combine3 -input1 -input2 -intersect -FixedEffectLMPlugin -endPlugin -export gt5382_outputtest1

    ./run_pipeline.pl -fork1 -importGuess gt5382_processed.vcf.gz -fork2 -importGuess gt5382_traits_final.txt -fork3 -importGuess gt5382_traits_final.txt -excludeLastTrait -combine5 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endPlugin -export gt5382_outputtest2
