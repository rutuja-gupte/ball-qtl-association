Model Comparisons
================
Rutuja Gupte

# Data

``` r
l <- grep("gapit", list.dirs(), value = TRUE)
dfs <- lapply(l, function(r){
d <- read.csv(paste(r, "GAPIT.Association.Filter_GWAS_results.csv", sep = "/"), header=TRUE)
return(d)
})
gapit <- do.call("rbind", dfs)

plink <- read.table("gt5382/gt5382.PHENO1.glm.linear")
colnames(plink) <- c("CHROM",   "POS",  "ID",   "REF",  "ALT",  "PROVISIONAL_REF",  "A1",   "OMITTED",  "A1_FREQ",  "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P",    "ERRCODE")
tassel <- read.table("gt5382/gt5382_1.txt", header = TRUE)

gapit$CHR <- gapit$Chr
gapit$BP <- gapit$Pos
gapit$P <- gapit$P.value

plink <- drop_na(plink)
plink$CHR <- plink$CHROM
plink$BP <- plink$POS
plink$SNP <- paste(plink$CHROM, "_", plink$POS, sep="")

tassel <- drop_na(tassel)

gapit$method <- str_split_i(paste("GAPIT.", gapit$traits), "\\.V3", 1)
plink$method <- "PLINK"

tassel$CHR <- tassel$Chr
tassel$BP <- tassel$Pos
tassel$SNP <- tassel$Marker
tassel$P <- tassel$p
tassel$method <- "TASSEL"

gmodel <- read.csv("gt5382/out.csv", skip=15)
gmodel <- drop_na(gmodel)
gmodel$CHR <- gmodel$Chrom
gmodel$BP <- as.numeric(str_extract(gmodel$Marker, "_(.*)_", group = 1))
gmodel$P <- as.numeric(gmodel$p.value)
gmodel$P <- replace_na(gmodel$P, 1e-7)
gmodel$SNP <- gmodel$Marker
gmodel$method <- "GModel"
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
data.long <- bind_rows(
                  gapit[,c("CHR", "BP", "P", "method")], 
                  plink[,c("CHR", "BP", "P", "method")],
                  gmodel[,c("CHR", "BP", "P", "method")],
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

More plots

``` r
data.long %>% filter(P < 1e-5) %>% ggplot() +
  geom_point(aes(x=BP, y=-log10(P))) +
  geom_hline(yintercept=5, color='red') +
  facet_grid(cols=vars(CHR), rows=vars(method)) +
  theme(axis.text.x = element_text(angle = 90))
```

![](Comparisons_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
data.long %>% filter(P < 1e-5)
```

    ## # A tibble: 33 × 4
    ## # Groups:   CHR, BP [20]
    ##      CHR       BP method                P
    ##    <int>    <dbl> <chr>             <dbl>
    ##  1     1   452827 GModel         1   e- 7
    ##  2     1 39613580 GAPIT. BLINK   5.54e- 9
    ##  3     1 39613580 GAPIT. FarmCPU 2.92e- 7
    ##  4     1 42721991 GAPIT. BLINK   5.21e- 9
    ##  5     1 42721991 GAPIT. FarmCPU 7.88e-10
    ##  6     1 70326408 GAPIT. BLINK   1.92e- 8
    ##  7     1 70326408 GAPIT. FarmCPU 5.07e- 8
    ##  8     1 78113465 GModel         6   e- 7
    ##  9     2  2072516 GModel         1   e- 7
    ## 10     2 69565451 GAPIT. BLINK   8.54e- 7
    ## # ℹ 23 more rows

``` r
data.wide %>% filter(if_any(1:(ncol(data.wide)-2), ~ . < 1e-5)) %>% print(n=Inf, width = Inf)
```

    ## # A tibble: 20 × 8
    ## # Groups:   CHR, BP [20]
    ##      CHR       BP        PLINK      TASSEL     GModel `GAPIT. BLINK`
    ##    <int>    <dbl>        <dbl>       <dbl>      <dbl>          <dbl>
    ##  1     1   452827 NA           NA           0.0000001      NA       
    ##  2     1 39613580  0.000568    NA          NA               5.54e- 9
    ##  3     1 42721991  0.000568    NA          NA               5.21e- 9
    ##  4     1 70326408  0.000568    NA          NA               1.92e- 8
    ##  5     1 78113465  0.0405       0.0563      0.0000006      NA       
    ##  6     2  2072516 NA           NA           0.0000001      NA       
    ##  7     2 69565451  0.000568    NA          NA               8.54e- 7
    ##  8     2 76205547  0.0803       0.183       0.0000001      NA       
    ##  9     3   229097 NA           NA           0.0000001      NA       
    ## 10     3 69978891  0.577        0.00000301 NA              NA       
    ## 11     4   637778 NA           NA           0.0000001      NA       
    ## 12     4 60695692  0.00953      0.0259      0.0000015      NA       
    ## 13     5   222932 NA           NA           0.0000001      NA       
    ## 14     5  1504132  0.000000668 NA          NA               5.64e-14
    ## 15     5 19864625  0.000000668 NA          NA               3.09e-14
    ## 16     6   112235 NA           NA           0.0000001      NA       
    ## 17     7    71273 NA           NA           0.0000001      NA       
    ## 18     7 34313233  0.00260      0.00971    NA              NA       
    ## 19     7 34313291  0.00202      0.00759     0.0000001       5.41e- 8
    ## 20     8    46336 NA           NA           0.0000001      NA       
    ##    `GAPIT. FarmCPU` `GAPIT. MLMM`
    ##               <dbl>         <dbl>
    ##  1        NA        NA           
    ##  2         2.92e- 7  0.0000291   
    ##  3         7.88e-10  0.0000249   
    ##  4         5.07e- 8 NA           
    ##  5        NA        NA           
    ##  6        NA        NA           
    ##  7         2.04e- 7  0.0000324   
    ##  8        NA        NA           
    ##  9        NA        NA           
    ## 10        NA        NA           
    ## 11        NA        NA           
    ## 12        NA        NA           
    ## 13        NA        NA           
    ## 14         2.27e-13  0.0000000165
    ## 15         7.59e-14  0.0000000352
    ## 16        NA        NA           
    ## 17        NA        NA           
    ## 18         1.86e- 7  0.00000427  
    ## 19        NA         0.00000322  
    ## 20        NA        NA

``` r
data.wide %>% filter(if_all(1:(ncol(data.wide)-2), ~ . < 1e-5)) %>% print(n=Inf, width = Inf)
```

    ## # A tibble: 0 × 8
    ## # Groups:   CHR, BP [0]
    ## # ℹ 8 variables: CHR <int>, BP <dbl>, PLINK <dbl>, TASSEL <dbl>, GModel <dbl>, GAPIT. BLINK <dbl>, GAPIT. FarmCPU <dbl>, GAPIT. MLMM <dbl>

Analyzing the significant sites further

``` r
hits <- data.wide %>% filter(if_any(1:(ncol(data.wide)-2), ~ . < 1e-5))
vcf <- read.vcfR("gt5382/gt5382_processed.vcf.gz")
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
hits$snp <- paste(hits$CHR, "_", hits$BP)
vcf.snps <- paste(vcf@fix[,"CHROM"], "_", vcf@fix[,"POS"])
bool_ser <- vcf.snps %in% hits$snp
vcf@gt <- vcf@gt[bool_ser,]
vcf@fix <- vcf@fix[bool_ser,]
gt <- extract.gt(vcf, element="GT")
gt <- data.frame(gt)
gt$CHR <- vcf@fix[,"CHROM"]
gt$BP <- vcf@fix[,"POS"]
traits <- read.table("gt5382/gt5382_traits.txt")
colnames(traits) <- c("V1", "sample", "value")
really_long_hits <- pivot_longer(gt, names_to="sample", values_to="gt", cols = c(-CHR, -BP))
really_long_hits <- left_join(really_long_hits, traits[,c("sample", "value")], by="sample")
really_long_hits$snp <- paste(really_long_hits$CHR, "_", really_long_hits$BP, sep="")
```

Now trying to plot them all horizontally

``` r
really_long_hits %>%
  ggplot() +
  geom_jitter(aes(x=snp, y=value, color=gt), width=0.2, height=0.5) +
  theme(axis.text.x = element_text(angle = 90))
```

![](Comparisons_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
