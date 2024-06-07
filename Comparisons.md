Model Comparison
================
Rutuja Gupte

# Data

``` r
gapit <- read.csv("TASSEL_samples/GAPIT.Association.Filter_GWAS_results.csv")
gmodel <- read.csv("TASSEL_samples/out.csv", skip=15)
plink <- read.table("TASSEL_samples/tassel_sample.qassoc", header = TRUE)

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
# gmodel$method <- "GModel"
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
# manhattan(gmodel)
```

# Comparisons

``` r
data.long <- bind_rows(gapit[,c("CHR", "BP", "P", "method")], 
                  plink[,c("CHR", "BP", "P", "method")],
                  # gmodel[,c("CHR", "BP", "P", "method")]
                  )
data.wide <- pivot_wider(distinct(data.long), names_from = method, values_from = P)
```

``` r
data.long %>% filter(P < 1e-5)
```

    ##    CHR        BP            P     method
    ## 1    3 180532745 2.522444e-06     GLM.V3
    ## 2    8 160536300 4.620555e-06     GLM.V3
    ## 3    3 180532745 9.299295e-07     GLM.V4
    ## 4    8 134813437 4.577439e-06     GLM.V4
    ## 5    4 167650309 1.511361e-07     GLM.V5
    ## 6    4 167650309 5.452788e-09   SUPER.V5
    ## 7    3 168689563 5.659851e-06 FarmCPU.V3
    ## 8    8 160536300 8.287069e-06 FarmCPU.V3
    ## 9    4  29209539 5.245487e-06 FarmCPU.V4
    ## 10   2  44606596 3.026087e-06 FarmCPU.V5
    ## 11   3 180532745 2.128248e-06   BLINK.V3
    ## 12   8 160536300 3.990245e-06   BLINK.V3
    ## 13   8 134813437 6.605887e-06   BLINK.V4
    ## 14   4 167650309 3.432674e-12   BLINK.V5
    ## 15   9 124067174 4.178212e-06   BLINK.V5

``` r
data.wide %>% filter(if_any(3:ncol(data.wide), ~ . < 1e-5)) %>% print(n=Inf, width = Inf)
```

    ## # A tibble: 8 × 15
    ##     CHR        BP      GLM.V3       GLM.V4       GLM.V5   SUPER.V3   SUPER.V4
    ##   <int>     <int>       <dbl>        <dbl>        <dbl>      <dbl>      <dbl>
    ## 1     3 180532745  0.00000252  0.000000930 NA           NA          0.0000113
    ## 2     8 160536300  0.00000462 NA           NA            0.0000145 NA        
    ## 3     4  29209539 NA           0.0000538   NA           NA          0.0000292
    ## 4     8 134813437 NA           0.00000458  NA           NA          0.0000655
    ## 5     4 167650309 NA          NA            0.000000151 NA         NA        
    ## 6     3 168689563 NA          NA           NA            0.0000292 NA        
    ## 7     2  44606596 NA          NA           NA           NA         NA        
    ## 8     9 124067174 NA          NA           NA           NA         NA        
    ##   SUPER.V5  FarmCPU.V3  FarmCPU.V4  FarmCPU.V5    BLINK.V3    BLINK.V4  BLINK.V5
    ##      <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>     <dbl>
    ## 1 NA       NA          NA          NA           0.00000213 NA          NA       
    ## 2 NA        0.00000829 NA          NA           0.00000399 NA          NA       
    ## 3 NA       NA           0.00000525 NA          NA           0.0000108  NA       
    ## 4 NA        0.0000836  NA          NA          NA           0.00000661 NA       
    ## 5  5.45e-9 NA          NA           0.0000190  NA          NA           3.43e-12
    ## 6 NA        0.00000566 NA          NA          NA          NA          NA       
    ## 7 NA       NA          NA           0.00000303 NA          NA           7.16e- 5
    ## 8 NA       NA          NA           0.0000122  NA          NA           4.18e- 6
    ##    PLINK
    ##    <dbl>
    ## 1 0.617 
    ## 2 0.116 
    ## 3 0.541 
    ## 4 0.791 
    ## 5 0.241 
    ## 6 0.246 
    ## 7 0.0213
    ## 8 0.700

``` r
data.wide %>% filter(if_all(3:ncol(data.wide), ~ . < 1e-5)) %>% print(n=Inf, width = Inf)
```

    ## # A tibble: 0 × 15
    ## # ℹ 15 variables: CHR <int>, BP <int>, GLM.V3 <dbl>, GLM.V4 <dbl>, GLM.V5 <dbl>, SUPER.V3 <dbl>, SUPER.V4 <dbl>, SUPER.V5 <dbl>, FarmCPU.V3 <dbl>, FarmCPU.V4 <dbl>, FarmCPU.V5 <dbl>, BLINK.V3 <dbl>, BLINK.V4 <dbl>, BLINK.V5 <dbl>, PLINK <dbl>
