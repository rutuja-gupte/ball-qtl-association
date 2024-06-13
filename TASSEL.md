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

For TASSEL, I still havent figured out how to properly get around the
covariate thing. I am currently passing the phenotype file as covariate
but that shouldnt be the case. This is the current way I am running it.
I think it should be fine because I have the exclude last trait flag so
I am hoping it just ignores all covariates.

    ./run_pipeline.pl -fork1 -importGuess processed.vcf.gz -fork2 -importGuess traits_final.txt -fork3 -importGuess traits_final.txt -excludeLastTrait -combine5 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endPlugin -export glm_output
