Processing of tuber bruising phenotype
================
Olivia Angelin-Bonnet
June 29, 2022

``` r
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             legend.position = "bottom",
             text = element_text(size = 16))
```

## Data path

``` r
samples_info_file = here("data/samples_info.csv")
bruising_file = here("data/raw_bruising_scores.csv")
```

``` r
samples_info = read_csv(samples_info_file, col_types = cols())
```

## Reading the bruising phenotype

We’ll read the bruising score dataset and get an average bruising score
for each genotype (= average over 6 tubers, 3 from each biological
replicate).

``` r
bruising_df = read_csv(bruising_file, show_col_types = FALSE) %>% 
  select(Sample, Rep, starts_with("Bruise")) %>% 
  pivot_longer(cols = starts_with("Bruise"),
               names_to = "Bruise",
               values_to = "bruising_score") %>% 
  group_by(Sample) %>% 
  summarise(bruising_score_mean = mean(bruising_score, na.rm = TRUE))
```

We can assess whether each phenotype is normally distributed using the
Shapiro-Wilk test.

``` r
shapiro.test(bruising_df$bruising_score_mean)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  bruising_df$bruising_score_mean
    ## W = 0.98149, p-value = 0.03066

The test shows evidence of non-normality.

## Normalisation of phenotypes using `bestNormalize`

The package `bestNormalize` tests different normalisation methods on a
given dataset and automatically selects the best performing one. It can
be used as follows:

``` r
set.seed(1)
norm_vals <- bestNormalize(bruising_df$bruising_score_mean)
norm_vals
```

    ## Best Normalizing transformation with 160 Observations
    ##  Estimated Normality Statistics (Pearson P / df, lower => more normal):
    ##  - arcsinh(x): 1.1975
    ##  - Center+scale: 1.2019
    ##  - Exp(x): 3.1225
    ##  - Log_b(x+a): 2.1275
    ##  - orderNorm (ORQ): 1.1713
    ##  - sqrt(x + a): 1.2806
    ##  - Yeo-Johnson: 1.1538
    ## Estimation method: Out-of-sample via CV with 10 folds and 5 repeats
    ##  
    ## Based off these, bestNormalize chose:
    ## Standardized Yeo-Johnson Transformation with 160 nonmissing obs.:
    ##  Estimated statistics:
    ##  - lambda = 0.4944419 
    ##  - mean (before standardization) = 1.144867 
    ##  - sd (before standardization) = 0.5031638

``` r
tibble(Sample = bruising_df$Sample,
       bruising_score_mean = norm_vals$x.t) %>% 
  write_csv(here("genomics_gwas/processed_data/pheno_df_trans.csv"))
```

``` r
shapiro.test(norm_vals$x.t)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  norm_vals$x.t
    ## W = 0.98881, p-value = 0.2338
