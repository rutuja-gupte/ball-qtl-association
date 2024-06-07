Model Comparison
================
Rutuja Gupte

# Data

``` r
gapit <- read.csv("GC1_round1/GAPIT.Association.Filter_GWAS_results.csv")
gmodel <- read.csv("GC1_round1/out.csv", skip=15)
plink <- read.table("GC1_round1/GC_analysis.assoc", header = TRUE)
tassel <- read.table("GC1_round1/tassel_results.txt", header = TRUE)

gmodel <- drop_na(gmodel)
gmodel$CHR <- gmodel$Chrom
gmodel$BP <- as.numeric(str_extract(gmodel$Marker, "_(.*)", group = 1))
gmodel$P <- as.numeric(gmodel$p.value)
gmodel$P <- replace_na(gmodel$P, 1e-7)
gmodel$SNP <- gmodel$Marker


gapit$CHR <- gapit$Chr
gapit$BP <- gapit$Pos
gapit$P <- gapit$P.value

plink <- drop_na(plink)


gapit$method <- gapit$traits
plink$method <- "PLINK"
gmodel$method <- "GModel"

tassel$CHR <- tassel$Chr
tassel$BP <- tassel$Pos
tassel$SNP <- tassel$Marker
tassel$P <- tassel$p
tassel$method <- "TASSEL"
```

# Some plots first

``` r
manhattan(gapit)
```

![](Comparisons_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
manhattan(plink)
```

![](Comparisons_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
manhattan(gmodel)
```

![](Comparisons_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
manhattan(tassel)
```

![](Comparisons_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

# Comparisons

``` r
data.long <- bind_rows(gapit[,c("CHR", "BP", "P", "method")], 
                  plink[,c("CHR", "BP", "P", "method")],
                  gmodel[,c("CHR", "BP", "P", "method")],
                  tassel[,c("CHR", "BP", "P", "method")]
                  )
data.wide <- pivot_wider(distinct(data.long), names_from = method, values_from = P)
```

``` r
data.long %>% filter(P < 1e-5)
```

    ##     CHR        BP            P     method
    ## 1     3 151318433 9.685196e-11     GLM.V3
    ## 2     3 151520667 1.604024e-13     GLM.V3
    ## 3     3 151520905 2.119468e-13     GLM.V3
    ## 4     3 151847862 2.248703e-12     GLM.V3
    ## 5     3 151847869 8.194591e-13     GLM.V3
    ## 6     3 152101422 4.339464e-14     GLM.V3
    ## 7     3 152170961 1.451882e-14     GLM.V3
    ## 8     3 152171110 6.541081e-16     GLM.V3
    ## 9     3 152171122 6.541081e-16     GLM.V3
    ## 10    3 157186920 2.064995e-08     GLM.V3
    ## 11    3 157186980 2.064995e-08     GLM.V3
    ## 12    3 157282394 3.345208e-09     GLM.V3
    ## 13    3 157545359 1.428748e-07     GLM.V3
    ## 14    3 157545410 1.428748e-07     GLM.V3
    ## 15    3 158122283 1.878646e-06     GLM.V3
    ## 16    3 158122291 1.878646e-06     GLM.V3
    ## 17    3 158615242 7.255676e-06     GLM.V3
    ## 18    3 158672450 5.542389e-06     GLM.V3
    ## 19    3 158672498 5.542389e-06     GLM.V3
    ## 20    3 158844731 3.940975e-06     GLM.V3
    ## 21    1   6622274 2.027071e-06   SUPER.V3
    ## 22    1  65808071 4.556467e-07   SUPER.V3
    ## 23    3  44686136 2.027071e-06   SUPER.V3
    ## 24    3 151318433 1.404514e-38   SUPER.V3
    ## 25    3 151520667 4.522109e-45   SUPER.V3
    ## 26    3 151520905 6.651088e-45   SUPER.V3
    ## 27    3 151847862 1.683714e-42   SUPER.V3
    ## 28    3 151847869 1.829033e-43   SUPER.V3
    ## 29    3 152101422 5.053939e-46   SUPER.V3
    ## 30    3 152170961 4.339043e-47   SUPER.V3
    ## 31    3 152171110 6.597248e-50   SUPER.V3
    ## 32    3 152171122 6.597248e-50   SUPER.V3
    ## 33    3 157186920 2.826744e-32   SUPER.V3
    ## 34    3 157186980 2.826744e-32   SUPER.V3
    ## 35    3 157282394 1.062527e-34   SUPER.V3
    ## 36    3 157545359 6.209433e-30   SUPER.V3
    ## 37    3 157545410 6.209433e-30   SUPER.V3
    ## 38    3 157887387 2.535741e-20   SUPER.V3
    ## 39    3 158073671 9.996984e-22   SUPER.V3
    ## 40    3 158122283 1.894514e-26   SUPER.V3
    ## 41    3 158122291 1.894514e-26   SUPER.V3
    ## 42    3 158615242 1.262343e-24   SUPER.V3
    ## 43    3 158672399 9.256817e-08   SUPER.V3
    ## 44    3 158672450 6.431244e-25   SUPER.V3
    ## 45    3 158672498 6.431244e-25   SUPER.V3
    ## 46    3 158844731 7.960017e-10   SUPER.V3
    ## 47    4  44571820 1.911027e-06   SUPER.V3
    ## 48    4  64566768 1.911027e-06   SUPER.V3
    ## 49    4 114870976 1.640691e-06   SUPER.V3
    ## 50    4 114870991 1.640691e-06   SUPER.V3
    ## 51    4 124758701 4.112962e-06   SUPER.V3
    ## 52    4 130423081 1.640691e-06   SUPER.V3
    ## 53    4 156327879 1.640691e-06   SUPER.V3
    ## 54    4 161183825 1.911027e-06   SUPER.V3
    ## 55    5 120677425 2.027071e-06   SUPER.V3
    ## 56    6   6626093 2.027071e-06   SUPER.V3
    ## 57    7 155867495 4.112962e-06   SUPER.V3
    ## 58    7 252306803 1.911027e-06   SUPER.V3
    ## 59    1   6622274 7.254936e-18 FarmCPU.V3
    ## 60    1  10151897 7.705992e-17 FarmCPU.V3
    ## 61    1  54946246 8.656241e-19 FarmCPU.V3
    ## 62    1  65808071 3.714959e-18 FarmCPU.V3
    ## 63    2   3758575 9.165700e-17 FarmCPU.V3
    ## 64    3 152171110 2.789937e-73 FarmCPU.V3
    ## 65    4  44571820 4.105301e-15 FarmCPU.V3
    ## 66    4  64566768 4.105301e-15 FarmCPU.V3
    ## 67    4 124758701 2.714172e-18 FarmCPU.V3
    ## 68    4 130423081 2.619430e-18 FarmCPU.V3
    ## 69    4 161183825 4.105301e-15 FarmCPU.V3
    ## 70    7 252306803 4.105301e-15 FarmCPU.V3
    ## 71    1   6622274 7.062690e-22   BLINK.V3
    ## 72    1   7569772 1.232398e-22   BLINK.V3
    ## 73    1  10151897 1.699887e-21   BLINK.V3
    ## 74    1  65808071 1.341552e-20   BLINK.V3
    ## 75    2   3758575 1.736894e-19   BLINK.V3
    ## 76    3 152171110 4.506745e-86   BLINK.V3
    ## 77    4  44571820 6.693506e-21   BLINK.V3
    ## 78    4 124758701 1.454847e-22   BLINK.V3
    ## 79    4 130423081 7.027732e-21   BLINK.V3
    ## 80    3 151318433 3.989000e-13      PLINK
    ## 81    3 151520667 4.030000e-14      PLINK
    ## 82    3 151520905 4.030000e-14      PLINK
    ## 83    3 151847862 8.291000e-15      PLINK
    ## 84    3 151847869 1.699000e-15      PLINK
    ## 85    3 152101422 1.174000e-17      PLINK
    ## 86    3 152170961 7.322000e-16      PLINK
    ## 87    3 152171110 1.438000e-16      PLINK
    ## 88    3 152171122 1.438000e-16      PLINK
    ## 89    3 157186920 2.999000e-08      PLINK
    ## 90    3 157186980 2.999000e-08      PLINK
    ## 91    3 157282394 7.589000e-08      PLINK
    ## 92    3 157545359 7.020000e-07      PLINK
    ## 93    3 157545410 7.020000e-07      PLINK
    ## 94    3 158122283 3.151000e-06      PLINK
    ## 95    3 158122291 3.151000e-06      PLINK
    ## 96    1   7950814 5.500000e-06     GModel
    ## 97    2    643492 1.000000e-07     GModel
    ## 98    3    123715 1.000000e-07     GModel
    ## 99    3 152171110 1.000000e-07     GModel
    ## 100   4 135104193 6.600000e-06     GModel
    ## 101   4 151940350 4.400000e-06     GModel
    ## 102   5    321231 1.000000e-07     GModel
    ## 103   6    175345 1.000000e-07     GModel
    ## 104   7 207225414 1.500000e-06     GModel
    ## 105   3 151318433 1.349000e-18     TASSEL
    ## 106   3 151520667 6.899200e-26     TASSEL
    ## 107   3 151520905 6.899200e-26     TASSEL
    ## 108   3 151847862 4.805500e-24     TASSEL
    ## 109   3 151847869 6.677800e-26     TASSEL
    ## 110   3 152101422 8.661700e-34     TASSEL
    ## 111   3 152170961 1.190400e-31     TASSEL
    ## 112   3 152171110 6.804600e-39     TASSEL
    ## 113   3 152171122 6.804600e-39     TASSEL
    ## 114   3 157186920 6.557700e-11     TASSEL
    ## 115   3 157186980 6.557700e-11     TASSEL
    ## 116   3 157282394 5.911900e-11     TASSEL
    ## 117   3 157545359 1.555900e-08     TASSEL
    ## 118   3 157545410 1.555900e-08     TASSEL
    ## 119   3 158122283 4.699200e-07     TASSEL
    ## 120   3 158122291 4.699200e-07     TASSEL
    ## 121   3 158615242 6.865800e-06     TASSEL
    ## 122   3 158672450 4.532400e-06     TASSEL
    ## 123   3 158672498 4.532400e-06     TASSEL
    ## 124   3 158844731 4.145400e-06     TASSEL

``` r
data.wide %>% filter(if_any(3:ncol(data.wide), ~ . < 1e-5)) %>% print(n=Inf, width = Inf)
```

    ## # A tibble: 50 × 9
    ##      CHR        BP    GLM.V3  SUPER.V3 FarmCPU.V3  BLINK.V3     PLINK     GModel
    ##    <int>     <dbl>     <dbl>     <dbl>      <dbl>     <dbl>     <dbl>      <dbl>
    ##  1     3 151318433  9.69e-11  1.40e-38  NA        NA         3.99e-13 NA        
    ##  2     3 151520667  1.60e-13  4.52e-45  NA        NA         4.03e-14 NA        
    ##  3     3 151520905  2.12e-13  6.65e-45  NA        NA         4.03e-14 NA        
    ##  4     3 151847862  2.25e-12  1.68e-42  NA        NA         8.29e-15 NA        
    ##  5     3 151847869  8.19e-13  1.83e-43  NA        NA         1.70e-15 NA        
    ##  6     3 152101422  4.34e-14  5.05e-46  NA        NA         1.17e-17 NA        
    ##  7     3 152170961  1.45e-14  4.34e-47  NA        NA         7.32e-16 NA        
    ##  8     3 152171110  6.54e-16  6.60e-50   2.79e-73  4.51e-86  1.44e-16  0.0000001
    ##  9     3 152171122  6.54e-16  6.60e-50  NA        NA         1.44e-16 NA        
    ## 10     3 157186920  2.06e- 8  2.83e-32  NA        NA         3.00e- 8 NA        
    ## 11     3 157186980  2.06e- 8  2.83e-32  NA        NA         3.00e- 8 NA        
    ## 12     3 157282394  3.35e- 9  1.06e-34  NA        NA         7.59e- 8 NA        
    ## 13     3 157545359  1.43e- 7  6.21e-30  NA        NA         7.02e- 7 NA        
    ## 14     3 157545410  1.43e- 7  6.21e-30  NA        NA         7.02e- 7 NA        
    ## 15     3 158122283  1.88e- 6  1.89e-26  NA        NA         3.15e- 6 NA        
    ## 16     3 158122291  1.88e- 6  1.89e-26  NA        NA         3.15e- 6 NA        
    ## 17     3 158615242  7.26e- 6  1.26e-24  NA        NA         9.08e- 5 NA        
    ## 18     3 158672450  5.54e- 6  6.43e-25  NA        NA         3.21e- 5 NA        
    ## 19     3 158672498  5.54e- 6  6.43e-25  NA        NA         3.21e- 5 NA        
    ## 20     3 158844731  3.94e- 6  7.96e-10  NA        NA         1.67e- 5 NA        
    ## 21     1   6622274 NA         2.03e- 6   7.25e-18  7.06e-22 NA        NA        
    ## 22     1  65808071 NA         4.56e- 7   3.71e-18  1.34e-20  3.09e- 1 NA        
    ## 23     3  44686136 NA         2.03e- 6  NA        NA        NA        NA        
    ## 24     3 157887387 NA         2.54e-20  NA        NA         3.72e- 4 NA        
    ## 25     3 158073671 NA         1.00e-21  NA        NA         5.11e- 5 NA        
    ## 26     3 158672399 NA         9.26e- 8  NA        NA         9.08e- 5 NA        
    ## 27     4  44571820 NA         1.91e- 6   4.11e-15  6.69e-21 NA        NA        
    ## 28     4  64566768 NA         1.91e- 6   4.11e-15 NA        NA        NA        
    ## 29     4 114870976 NA         1.64e- 6  NA        NA        NA        NA        
    ## 30     4 114870991 NA         1.64e- 6  NA        NA        NA        NA        
    ## 31     4 124758701 NA         4.11e- 6   2.71e-18  1.45e-22 NA        NA        
    ## 32     4 130423081 NA         1.64e- 6   2.62e-18  7.03e-21 NA        NA        
    ## 33     4 156327879 NA         1.64e- 6  NA        NA        NA        NA        
    ## 34     4 161183825 NA         1.91e- 6   4.11e-15 NA        NA        NA        
    ## 35     5 120677425 NA         2.03e- 6  NA        NA        NA        NA        
    ## 36     6   6626093 NA         2.03e- 6  NA        NA        NA        NA        
    ## 37     7 155867495 NA         4.11e- 6  NA        NA        NA        NA        
    ## 38     7 252306803 NA         1.91e- 6   4.11e-15 NA        NA        NA        
    ## 39     1  10151897 NA        NA          7.71e-17  1.70e-21 NA        NA        
    ## 40     1  54946246 NA        NA          8.66e-19 NA        NA        NA        
    ## 41     2   3758575 NA        NA          9.17e-17  1.74e-19 NA        NA        
    ## 42     1   7569772 NA        NA         NA         1.23e-22 NA        NA        
    ## 43     1   7950814 NA        NA         NA        NA        NA         0.0000055
    ## 44     2    643492 NA        NA         NA        NA        NA         0.0000001
    ## 45     3    123715 NA        NA         NA        NA        NA         0.0000001
    ## 46     4 135104193 NA        NA         NA        NA        NA         0.0000066
    ## 47     4 151940350 NA        NA         NA        NA        NA         0.0000044
    ## 48     5    321231 NA        NA         NA        NA        NA         0.0000001
    ## 49     6    175345 NA        NA         NA        NA        NA         0.0000001
    ## 50     7 207225414 NA        NA         NA        NA        NA         0.0000015
    ##       TASSEL
    ##        <dbl>
    ##  1  1.35e-18
    ##  2  6.90e-26
    ##  3  6.90e-26
    ##  4  4.81e-24
    ##  5  6.68e-26
    ##  6  8.66e-34
    ##  7  1.19e-31
    ##  8  6.80e-39
    ##  9  6.80e-39
    ## 10  6.56e-11
    ## 11  6.56e-11
    ## 12  5.91e-11
    ## 13  1.56e- 8
    ## 14  1.56e- 8
    ## 15  4.70e- 7
    ## 16  4.70e- 7
    ## 17  6.87e- 6
    ## 18  4.53e- 6
    ## 19  4.53e- 6
    ## 20  4.15e- 6
    ## 21 NA       
    ## 22 NA       
    ## 23 NA       
    ## 24  1.04e- 5
    ## 25  2.37e- 5
    ## 26  2.06e- 5
    ## 27 NA       
    ## 28 NA       
    ## 29 NA       
    ## 30 NA       
    ## 31 NA       
    ## 32 NA       
    ## 33 NA       
    ## 34 NA       
    ## 35 NA       
    ## 36 NA       
    ## 37 NA       
    ## 38 NA       
    ## 39 NA       
    ## 40 NA       
    ## 41 NA       
    ## 42 NA       
    ## 43 NA       
    ## 44 NA       
    ## 45 NA       
    ## 46 NA       
    ## 47 NA       
    ## 48 NA       
    ## 49 NA       
    ## 50 NA

``` r
data.wide %>% filter(if_all(3:ncol(data.wide), ~ . < 1e-5)) %>% print(n=Inf, width = Inf)
```

    ## # A tibble: 1 × 9
    ##     CHR        BP   GLM.V3 SUPER.V3 FarmCPU.V3 BLINK.V3    PLINK    GModel
    ##   <int>     <dbl>    <dbl>    <dbl>      <dbl>    <dbl>    <dbl>     <dbl>
    ## 1     3 152171110 6.54e-16 6.60e-50   2.79e-73 4.51e-86 1.44e-16 0.0000001
    ##     TASSEL
    ##      <dbl>
    ## 1 6.80e-39
