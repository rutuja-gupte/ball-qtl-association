PLINK
================
Rutuja Gupte

This is my favorite one because it is entirely command-line based and
does not need any extra processing other than the common processing we
have done so far. Even within the common pre-processing, if you do not
want to do some of that, it is totally fine because that can be
accounted for using flags.

*Important note:* PLINK does not allow floating point phenotypes. They
have to be integers.

``` r
# temp <- read.table("TASSEL_samples/mdp_traits2.txt")
# temp <- temp %>% mutate(across(3:ncol(temp), round))
# temp <- temp %>% mutate(
#     across(everything(), ~replace_na(.x, 0))
#   )
# write.table(temp, "TASSEL_samples/traits_final.txt", row.names=FALSE, col.names=FALSE, sep = "\t", quote=FALSE)
```

The main command that I would use to run PLINK is:

    ./plink --vcf gt5382_processed.vcf.gz --double-id --pheno gt5382_traits.txt --mpheno 1 --assoc --out gt5382 --allow-no-sex

For multiple phenotypes:  
(Putting this in here because I am terrible at guessing what the change
would be - spoiler alert: it is not mpheno that changes, I need to throw
in another flag)

    ./plink --vcf TASSEL_samples_processed.vcf.gz --double-id --pheno TASSEL_samples_traits.txt --mpheno 1 --assoc --out TASSEL_samples --allow-no-sex --all-pheno

The beauty of it all is that it works in both PowerShell and Linux. Make
sure all the appropriate files are in the directory containing the
executable PLINK file.

Move the .qassoc file into the working directory to visualize the data
and make Manhattan plots. We will require the R package qqman.

``` r
plink <- read.table("gt5382/gt5382.qassoc", header = TRUE)
plink <- drop_na(plink)
manhattan(plink)
```

![](PLINK_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
