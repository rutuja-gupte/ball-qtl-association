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
gmodel$BP <- as.numeric(str_extract(gmodel$Marker, "_(.*)_", group = 1))
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
data.long <- bind_rows(gapit[,c("SNP", "CHR", "BP", "P", "method")], 
                  plink[,c("SNP", "CHR", "BP", "P", "method")],
                  # gmodel[,c("SNP", "CHR", "BP", "P", "method")]
                  )
data.long$SNP <- as.character(data.long$SNP)
data.wide <- pivot_wider(distinct(data.long), names_from = method, values_from = P)
```

``` r
data.long %>% filter(P < 1e-5)
```

    ##            SNP CHR        BP            P     method
    ## 1   PZA03734.2   3 180532745 2.522444e-06     GLM.V3
    ## 2   PZB00811.1   8 160536300 4.620555e-06     GLM.V3
    ## 3   PZA03734.2   3 180532745 9.299295e-07     GLM.V4
    ## 4   PZA03591.1   8 134813437 4.577439e-06     GLM.V4
    ## 5   PZA03379.2   4 167650309 1.511361e-07     GLM.V5
    ## 6   PZA03379.2   4 167650309 5.452788e-09   SUPER.V5
    ## 7   PZD00040.3   3 168689563 5.659851e-06 FarmCPU.V3
    ## 8   PZB00811.1   8 160536300 8.287069e-06 FarmCPU.V3
    ## 9   PZA03693.1   4  29209539 5.245487e-06 FarmCPU.V4
    ## 10 PZA02808.12   2  44606596 3.026087e-06 FarmCPU.V5
    ## 11  PZA03734.2   3 180532745 2.128248e-06   BLINK.V3
    ## 12  PZB00811.1   8 160536300 3.990245e-06   BLINK.V3
    ## 13  PZA03591.1   8 134813437 6.605887e-06   BLINK.V4
    ## 14  PZA03379.2   4 167650309 3.432674e-12   BLINK.V5
    ## 15  PZB01432.1   9 124067174 4.178212e-06   BLINK.V5

``` r
data.wide %>% filter(if_any(3:ncol(data.wide), ~ . < 1e-5)) %>% print(n=Inf)
```

    ## # A tibble: 8 × 16
    ##   SNP           CHR     BP   GLM.V3   GLM.V4   GLM.V5 SUPER.V3 SUPER.V4 SUPER.V5
    ##   <chr>       <int>  <int>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
    ## 1 PZA03734.2      3 1.81e8  2.52e-6  9.30e-7 NA       NA        1.13e-5 NA      
    ## 2 PZB00811.1      8 1.61e8  4.62e-6 NA       NA        1.45e-5 NA       NA      
    ## 3 PZA03693.1      4 2.92e7 NA        5.38e-5 NA       NA        2.92e-5 NA      
    ## 4 PZA03591.1      8 1.35e8 NA        4.58e-6 NA       NA        6.55e-5 NA      
    ## 5 PZA03379.2      4 1.68e8 NA       NA        1.51e-7 NA       NA        5.45e-9
    ## 6 PZD00040.3      3 1.69e8 NA       NA       NA        2.92e-5 NA       NA      
    ## 7 PZA02808.12     2 4.46e7 NA       NA       NA       NA       NA       NA      
    ## 8 PZB01432.1      9 1.24e8 NA       NA       NA       NA       NA       NA      
    ## # ℹ 7 more variables: FarmCPU.V3 <dbl>, FarmCPU.V4 <dbl>, FarmCPU.V5 <dbl>,
    ## #   BLINK.V3 <dbl>, BLINK.V4 <dbl>, BLINK.V5 <dbl>, PLINK <dbl>
