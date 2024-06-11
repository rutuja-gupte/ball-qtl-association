Model Comparison
================
Rutuja Gupte

# Data

``` r
gapit <- read.csv("gt5382/GAPIT.Association.Filter_GWAS_results.csv")
gmodel <- read.csv("gt5382/out.csv", skip=15)
plink <- read.table("gt5382/gt5382.qassoc", header = TRUE)
tassel <- read.table("gt5382/5382_out1.txt", header = TRUE)

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
tassel <- drop_na(tassel)

gapit$method <- gapit$traits
plink$method <- "PLINK"
# gmodel$method <- "GModel"

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
# manhattan(gmodel)
manhattan(tassel)
```

![](Comparisons_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

# Comparisons

``` r
data.long <- bind_rows(gapit[,c("CHR", "BP", "P", "method")], 
                  plink[,c("CHR", "BP", "P", "method")],
                  # gmodel[,c("CHR", "BP", "P", "method")],
                  tassel[,c("CHR", "BP", "P", "method")]
                  )

data.long <- data.long %>%
  group_by(CHR, BP, method) %>%
  summarise(P = mean(P))
```

    ## `summarise()` has grouped output by 'CHR', 'BP'. You can override using the
    ## `.groups` argument.

``` r
data.wide <- pivot_wider(distinct(data.long), names_from = method, values_from = P)
```

``` r
data.long %>% filter(P < 1e-5)
```

    ## # A tibble: 29 × 4
    ## # Groups:   CHR, BP [14]
    ##      CHR       BP method            P
    ##    <int>    <int> <chr>         <dbl>
    ##  1     1  1157973 SUPER.V3   9.05e- 6
    ##  2     1 39613580 BLINK.V3   2.74e- 9
    ##  3     1 39613580 FarmCPU.V3 4.51e- 8
    ##  4     1 42721991 BLINK.V3   2.39e- 9
    ##  5     1 42721991 FarmCPU.V3 5.23e-10
    ##  6     1 42721991 SUPER.V3   4.80e- 6
    ##  7     1 70326408 BLINK.V3   6.25e- 9
    ##  8     2 69565451 FarmCPU.V3 3.65e- 7
    ##  9     2 69565451 SUPER.V3   8.34e- 6
    ## 10     3  6555509 FarmCPU.V3 1.76e- 9
    ## # ℹ 19 more rows

``` r
data.wide %>% filter(if_any(1:(ncol(data.wide)-2), ~ . < 1e-5)) %>% print(n=Inf, width = Inf)
```

    ## # A tibble: 14 × 8
    ## # Groups:   CHR, BP [14]
    ##      CHR       BP       PLINK      TASSEL  SUPER.V3  BLINK.V3 FarmCPU.V3
    ##    <int>    <int>       <dbl>       <dbl>     <dbl>     <dbl>      <dbl>
    ##  1     1  1157973 0.000568    NA           9.05e- 6 NA         NA       
    ##  2     1 39613580 0.000568    NA           1.01e- 5  2.74e- 9   4.51e- 8
    ##  3     1 42721991 0.000568    NA           4.80e- 6  2.39e- 9   5.23e-10
    ##  4     1 70326408 0.000568    NA           1.09e- 5  6.25e- 9  NA       
    ##  5     2 69565451 0.000568    NA           8.34e- 6 NA          3.65e- 7
    ##  6     3  6555509 0.000568    NA           1.09e- 5 NA          1.76e- 9
    ##  7     3 69978891 0.577        0.00000301 NA        NA         NA       
    ##  8     4  4488295 0.000568    NA           8.34e- 6  1.12e- 6  NA       
    ##  9     4 10274863 0.000568    NA           9.05e- 6 NA         NA       
    ## 10     5  1504132 0.000000668 NA           4.98e-10  1.93e-14   1.69e-15
    ## 11     5 11206077 0.000568    NA           5.93e- 6 NA         NA       
    ## 12     5 19864625 0.000000668 NA           5.14e-10  4.72e-15   1.67e-14
    ## 13     7 34313233 0.00260      0.00971     7.28e- 7 NA         NA       
    ## 14     7 34313291 0.00202      0.00759     5.21e- 7  2.25e- 7   5.59e- 7
    ##         GLM.V3
    ##          <dbl>
    ##  1 NA         
    ##  2 NA         
    ##  3 NA         
    ##  4 NA         
    ##  5 NA         
    ##  6 NA         
    ##  7 NA         
    ##  8 NA         
    ##  9 NA         
    ## 10  0.00000242
    ## 11 NA         
    ## 12  0.00000305
    ## 13 NA         
    ## 14 NA

``` r
data.wide %>% filter(if_all(1:(ncol(data.wide)-2), ~ . < 1e-5)) %>% print(n=Inf, width = Inf)
```

    ## # A tibble: 0 × 8
    ## # Groups:   CHR, BP [0]
    ## # ℹ 8 variables: CHR <int>, BP <int>, PLINK <dbl>, TASSEL <dbl>, SUPER.V3 <dbl>, BLINK.V3 <dbl>, FarmCPU.V3 <dbl>, GLM.V3 <dbl>
