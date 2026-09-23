---
title: "Untitled"
output: 
  html_document: 
    keep_md: true
date: "2026-04-18"
---



## TO DO

- Predict on all test data and filter afterwards 

## Load Functions and Data

### Functions


``` r
source(here::here("Functions",
                    "mv_imspe_funs.R"))

source(here::here("Functions",
                    "kernel_funs.R"))

source(here::here("Functions",
                    "coords_tranform_funs.R"))

source(here::here("Functions",
                    "prediction_sampling_funs.R"))

source(here::here("Functions",
                    "orth_basis_funs.R"))

source(here::here("Functions",
                    "rm_mvgps_mat_mu_pred_functions.R"))
```


### Training Data

Column names and properties 


``` r
metals <- c("Mo","Nb","Ta","V","W")
properties <- c("log_max_strain",
                "bc_1","bc_2","bc_3",
                "bc_4","bc_5","bc_6")

props <- c("ultimate_stress","ultimate_strain",
           "toughness","young_mod")

prop_ltx <- unname(TeX(c("$\\sigma_u$","$\\epsilon_u$",
                         "$U_T$","$E$")))

coefs_ltx <- unname(TeX(c("$\\ln(\\epsilon_u)$",
                                  "$\\beta_1$","$\\beta_2$","$\\beta_3$",
                               "$\\beta_4$", "$\\beta_5$", "$\\beta_6$")))


n_prop <- length(properties)
n_dim <- length(metals) - 1
```

Data to compute transformation bounds (ILR)


``` r
mixture_dat <- gen.mixture(11,5)

eps_ilr <- 0.02

xx <- ilr_bounded_minmax_additive(as.matrix(mixture_dat),
                                  b = 0,
                                  epsilon = eps_ilr)
```

Initial data set


``` r
dat_ss_c <- rio::import(here::here("Results",
                                 "ss_all_ort_coefs_c_uns.csv"))

dat_ss_c <- as_tibble(dat_ss_c)%>%
              mutate(alloy_cnt = alloy_cnt.x)%>%
              select(-alloy_cnt.x, - alloy_cnt.y)

dat_ss_c <- dat_ss_c%>%
            filter(!(alloy %in% metals))%>%
            mutate(id_row = row_number())

dat_ss_t <- rio::import(here::here("Results",
                                 "ss_all_ort_coefs_t_uns.csv"))

dat_ss_t <- as_tibble(dat_ss_t)%>%
              mutate(alloy_cnt = alloy_cnt.x)%>%
              select(-alloy_cnt.x, - alloy_cnt.y)

dat_ss_t <- dat_ss_t%>%
            filter(!(alloy %in% metals))%>%
            mutate(id_row = row_number())

dat_ss_all <- dat_ss_t%>%
                filter(!(alloy %in% metals))%>%
                mutate(type = "t")%>%
                bind_rows(dat_ss_c%>%
                            mutate(type ="c"))%>%
                unite("alloy_id",
                      Mo:W,remove = FALSE)%>%
                arrange(alloy_id, rep, type)#%>%
                  #filter(! alloy %in% metals)

dat_ss_fil <- dat_ss_all%>%
              ungroup()%>%
              group_by(alloy_id)%>%
              arrange(desc(rep))%>%
              slice_head(n=1)%>%
  mutate(elements= (Mo>0) + (Nb>0) + (Ta>0) + (V>0) + (W>0))
```

Quaternaries data


``` r
dat_ss_ter_c <- rio::import(here::here("Results",
                                 "ss_ter_ort_coefs_c_uns.csv"))

dat_ss_ter_t <- rio::import(here::here("Results",
                                 "ss_ter_ort_coefs_t_uns.csv"))

dat_ss_ter_c <- as_tibble(dat_ss_ter_c)%>%
  mutate(elements= (Mo>0) + (Nb>0) + (Ta>0) + (V>0) + (W>0))

dat_ss_ter_t <- as_tibble(dat_ss_ter_t)%>%
  mutate(elements= (Mo>0) + (Nb>0) + (Ta>0) + (V>0) + (W>0))

max_ac <- max(dat_ss_fil$alloy_cnt)
dat_ss_t_all <- dat_ss_t%>%
  bind_rows(dat_ss_ter_t%>%
              mutate(alloy_cnt = alloy_cnt + max_ac))%>%
  arrange(alloy_cnt)

dat_ss_c_all <- dat_ss_c%>%
  bind_rows(dat_ss_ter_c%>%
              mutate(alloy_cnt = alloy_cnt + max_ac))%>%
  arrange(alloy_cnt)

dat_ss_fil_all <- dat_ss_fil%>%
  bind_rows(dat_ss_ter_c%>%
              mutate(alloy_cnt = alloy_cnt + max_ac)%>%
              group_by(alloy_cnt)%>%
              arrange(desc(rep))%>%
              slice_head(n=1)%>%
              ungroup())%>%
  arrange(alloy_cnt)


dat_ss_ter_c <- as_tibble(dat_ss_ter_c)%>%
                  pivot_longer(cols = all_of(properties),
                               names_to = "property",
                               values_to = "val_true")%>%
                  arrange(alloy_cnt)

dat_ss_ter_t <- as_tibble(dat_ss_ter_t)%>%
                  pivot_longer(cols = all_of(properties),
                               names_to = "property",
                               values_to = "val_true")%>%
                  arrange(alloy_cnt)

X_ter <- dat_ss_ter_c%>%
  group_by(alloy_cnt)%>%
  slice_head(n=1)
X_ter <- as.matrix(X_ter[,metals])

X_ter_B <- dat_ss_ter_c%>%
  group_by(alloy_cnt, rep)%>%
  slice_head(n=1)

X_ter_B <- as.matrix(X_ter_B[,metals])
```

## Test Heterosedaskticity


``` r
gauss_test <- function(x){
  nortest::ad.test(x)$p.value
}

bar_test <- function(x, g){
  g <- factor(g)
  df = data.frame(x = x,
                  g = g)
 res <- bartlett.test(x ~ g,data = df)
 res$p.value
}
```



``` r
dat_ss_c_val_err <- dat_ss_c_all%>%
  mutate(elements= (Mo>0) + (Nb>0) + (Ta>0) + (V>0) + (W>0))%>%
  filter(elements > 1)%>%
  group_by(alloy_cnt)%>%
  add_tally()%>%
  filter(n > 1)%>%
  pivot_longer(cols = all_of(properties),
                               names_to = "property",
                               values_to = "val_true")%>%
  group_by(alloy_cnt, property)%>%
  mutate(val_err = val_true - mean(val_true))%>%
  ungroup()


dat_ss_c_val_err%>%
  ggplot(aes(sample = val_err))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~property,
             scales = "free")
```

![](Predictions_19-05---Copy_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

``` r
dat_ss_c_val_err%>%
  group_by(property)%>%
  summarise(p_val = gauss_test(val_err))
```

```
## # A tibble: 7 × 2
##   property          p_val
##   <chr>             <dbl>
## 1 bc_1           1.78e- 7
## 2 bc_2           6.01e-11
## 3 bc_3           3.7 e-24
## 4 bc_4           3.7 e-24
## 5 bc_5           3.7 e-24
## 6 bc_6           3.7 e-24
## 7 log_max_strain 5.80e- 9
```

``` r
dat_ss_c_val_err%>%
  filter(elements < 3)%>%
  group_by(property, alloy)%>%
  summarise(p_val = gauss_test(val_err),
            p_val_bt = bar_test(val_err,alloy_cnt))%>%
  mutate(ad_pass = p_val > 0.05,
         bt_pass = p_val_bt > 0.05)
```

```
## `summarise()` has grouped output by 'property'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 70 × 6
## # Groups:   property [7]
##    property alloy   p_val   p_val_bt ad_pass bt_pass
##    <chr>    <chr>   <dbl>      <dbl> <lgl>   <lgl>  
##  1 bc_1     MoNb  0.0351  0.0136     FALSE   FALSE  
##  2 bc_1     MoTa  0.425   0.00912    TRUE    FALSE  
##  3 bc_1     MoV   0.170   0.0428     TRUE    FALSE  
##  4 bc_1     MoW   0.0374  0.113      FALSE   TRUE   
##  5 bc_1     NbTa  0.331   0.129      TRUE    TRUE   
##  6 bc_1     NbV   0.870   0.714      TRUE    TRUE   
##  7 bc_1     NbW   0.505   0.171      TRUE    TRUE   
##  8 bc_1     TaV   0.0692  0.0238     TRUE    FALSE  
##  9 bc_1     TaW   0.00567 0.00000500 FALSE   FALSE  
## 10 bc_1     VW    0.410   0.135      TRUE    TRUE   
## # ℹ 60 more rows
```




``` r
dat_ss_t_val_err <- dat_ss_t_all%>%
  mutate(elements= (Mo>0) + (Nb>0) + (Ta>0) + (V>0) + (W>0))%>%
  filter(elements > 1)%>%
  group_by(alloy_cnt)%>%
  add_tally()%>%
  filter(n > 1)%>%
  pivot_longer(cols = all_of(properties),
                               names_to = "property",
                               values_to = "val_true")%>%
  group_by(alloy_cnt, property)%>%
  mutate(val_err = val_true - mean(val_true))%>%
  ungroup()


dat_ss_t_val_err%>%
  ggplot(aes(sample = val_err))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~property,
             scales = "free")
```

![](Predictions_19-05---Copy_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

``` r
dat_ss_t_val_err%>%
  group_by(property)%>%
  summarise(p_val = gauss_test(val_err))
```

```
## # A tibble: 7 × 2
##   property          p_val
##   <chr>             <dbl>
## 1 bc_1           1.56e- 2
## 2 bc_2           3.7 e-24
## 3 bc_3           2.69e-14
## 4 bc_4           3.7 e-24
## 5 bc_5           3.7 e-24
## 6 bc_6           3.7 e-24
## 7 log_max_strain 1.68e-16
```

``` r
dat_ss_t_val_err%>%
  filter(elements < 3)%>%
  group_by(property, alloy)%>%
  summarise(p_val = gauss_test(val_err),
            p_val_bt = bar_test(val_err,alloy_cnt))%>%
  mutate(ad_pass = p_val > 0.05,
         bt_pass = p_val_bt > 0.05)
```

```
## `summarise()` has grouped output by 'property'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 70 × 6
## # Groups:   property [7]
##    property alloy   p_val p_val_bt ad_pass bt_pass
##    <chr>    <chr>   <dbl>    <dbl> <lgl>   <lgl>  
##  1 bc_1     MoNb  0.550   0.723    TRUE    TRUE   
##  2 bc_1     MoTa  0.373   0.192    TRUE    TRUE   
##  3 bc_1     MoV   0.761   0.833    TRUE    TRUE   
##  4 bc_1     MoW   0.156   0.00324  TRUE    FALSE  
##  5 bc_1     NbTa  0.0494  0.000414 FALSE   FALSE  
##  6 bc_1     NbV   0.988   0.420    TRUE    TRUE   
##  7 bc_1     NbW   0.00720 0.912    FALSE   TRUE   
##  8 bc_1     TaV   0.967   0.291    TRUE    TRUE   
##  9 bc_1     TaW   0.445   0.00230  TRUE    FALSE  
## 10 bc_1     VW    0.987   0.0830   TRUE    TRUE   
## # ℹ 60 more rows
```




``` r
dat_ss_c_higher_val_err <- dat_ss_c_val_err %>%
  filter(elements > 2)

dat_ss_c_higher_val_err %>%
  ggplot(aes(sample = val_err)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~property, scales = "free")
```

![](Predictions_19-05---Copy_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

``` r
dat_ss_c_higher_val_err %>%
  group_by(property) %>%
  summarise(p_val = gauss_test(val_err))
```

```
## # A tibble: 7 × 2
##   property          p_val
##   <chr>             <dbl>
## 1 bc_1           4.51e- 2
## 2 bc_2           6.94e- 2
## 3 bc_3           5.34e- 4
## 4 bc_4           1.51e-18
## 5 bc_5           3.7 e-24
## 6 bc_6           3.7 e-24
## 7 log_max_strain 2.24e- 1
```

``` r
dat_ss_c_higher_val_err %>%
  group_by(property, elements) %>%
  summarise(
    p_val = gauss_test(val_err),
    p_val_bt = bar_test(val_err, alloy_cnt)
  ) %>%
  mutate(
    ad_pass = p_val > 0.05,
    bt_pass = p_val_bt > 0.05
  )
```

```
## `summarise()` has grouped output by 'property'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 14 × 6
## # Groups:   property [7]
##    property       elements    p_val p_val_bt ad_pass bt_pass
##    <chr>             <int>    <dbl>    <dbl> <lgl>   <lgl>  
##  1 bc_1                  4 7.82e- 1 2.72e- 1 TRUE    TRUE   
##  2 bc_1                  5 8.72e- 2 1.15e- 1 TRUE    TRUE   
##  3 bc_2                  4 6.56e- 1 3.17e- 1 TRUE    TRUE   
##  4 bc_2                  5 7.88e- 2 9.08e- 2 TRUE    TRUE   
##  5 bc_3                  4 2.86e- 1 5.79e- 2 TRUE    TRUE   
##  6 bc_3                  5 3.94e- 3 1.55e- 2 FALSE   FALSE  
##  7 bc_4                  4 7.60e- 4 2.11e- 3 FALSE   FALSE  
##  8 bc_4                  5 9.11e-15 3.03e- 9 FALSE   FALSE  
##  9 bc_5                  4 2.22e- 5 5.28e- 4 FALSE   FALSE  
## 10 bc_5                  5 3.7 e-24 9.87e-18 FALSE   FALSE  
## 11 bc_6                  4 3.29e- 6 2.12e- 6 FALSE   FALSE  
## 12 bc_6                  5 3.7 e-24 4.72e-17 FALSE   FALSE  
## 13 log_max_strain        4 9.04e- 1 4.80e- 1 TRUE    TRUE   
## 14 log_max_strain        5 2.03e- 1 1.41e- 1 TRUE    TRUE
```



``` r
dat_ss_t_higher_val_err <- dat_ss_t_val_err %>%
  filter(elements > 2)

dat_ss_t_higher_val_err %>%
  ggplot(aes(sample = val_err)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~property, scales = "free")
```

![](Predictions_19-05---Copy_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

``` r
dat_ss_t_higher_val_err %>%
  group_by(property) %>%
  summarise(p_val = gauss_test(val_err))
```

```
## # A tibble: 7 × 2
##   property          p_val
##   <chr>             <dbl>
## 1 bc_1           2.02e- 1
## 2 bc_2           2.83e-19
## 3 bc_3           2.25e- 2
## 4 bc_4           2.36e-14
## 5 bc_5           1.36e-22
## 6 bc_6           3.57e-12
## 7 log_max_strain 1.29e-12
```

``` r
dat_ss_t_higher_val_err %>%
  group_by(property, elements) %>%
  summarise(
    p_val = gauss_test(val_err),
    p_val_bt = bar_test(val_err, alloy_cnt)
  ) %>%
  mutate(
    ad_pass = p_val > 0.05,
    bt_pass = p_val_bt > 0.05
  )
```

```
## `summarise()` has grouped output by 'property'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 14 × 6
## # Groups:   property [7]
##    property       elements    p_val p_val_bt ad_pass bt_pass
##    <chr>             <int>    <dbl>    <dbl> <lgl>   <lgl>  
##  1 bc_1                  4 7.63e- 1 8.48e- 1 TRUE    TRUE   
##  2 bc_1                  5 1.00e- 1 1.12e- 3 TRUE    FALSE  
##  3 bc_2                  4 8.73e- 1 6.94e- 1 TRUE    TRUE   
##  4 bc_2                  5 1.16e-15 1.85e-12 FALSE   FALSE  
##  5 bc_3                  4 6.80e- 1 2.55e- 1 TRUE    TRUE   
##  6 bc_3                  5 2.74e- 2 4.11e- 3 FALSE   FALSE  
##  7 bc_4                  4 8.93e- 1 2.62e- 1 TRUE    TRUE   
##  8 bc_4                  5 1.40e-11 4.75e-10 FALSE   FALSE  
##  9 bc_5                  4 2.18e- 1 6.81e- 2 TRUE    TRUE   
## 10 bc_5                  5 5.25e-16 1.30e-10 FALSE   FALSE  
## 11 bc_6                  4 2.79e- 2 1.61e- 2 FALSE   FALSE  
## 12 bc_6                  5 1.17e- 9 1.35e- 8 FALSE   FALSE  
## 13 log_max_strain        4 9.45e- 1 7.80e- 1 TRUE    TRUE   
## 14 log_max_strain        5 6.51e-12 9.52e-61 FALSE   FALSE
```



##fitted nugget  


``` r
load(here::here("Results", "Predictions_het.RData"))
load(here::here("Results", "Predictions_hom.RData"))

pred_ids <- dat_pm %>%
  distinct(id = id_obs, alloy_cnt, Mo, Nb, Ta, V, W) %>%
  filter((Mo > 0) + (Nb > 0) + (Ta > 0) + (V > 0) + (W > 0) == 2) %>%
  select(id, alloy_cnt)

stopifnot(nrow(pred_ids) == 40)


dat_het <- dat_nug %>%
  filter(
    iter == 0, Kernel == "m32",
    (Test == "Compression" & mean_type == "pol") |
      (Test == "Tension" & mean_type == "int")
  ) %>%
  inner_join(pred_ids, by = "id") %>%
  pivot_longer(
    cols = all_of(properties),
    names_to = "property",
    values_to = "nugget_het"
  ) %>%
  select(Test, alloy_cnt, property, nugget_het)

dat_hom <- dat_nug_hom %>%
  filter(
    iter == 0, Kernel == "m32",
    (Test == "Compression-Hom" & mean_type == "pol") |
      (Test == "Tension-Hom" & mean_type == "int")
  ) %>%
  mutate(Test = sub("-Hom$", "", Test)) %>%
  inner_join(pred_ids, by = "id") %>%
  pivot_longer(
    cols = all_of(properties),
    names_to = "property",
    values_to = "nugget_hom"
  ) %>%
  select(Test, alloy_cnt, property, nugget_hom)

dat_observed <- bind_rows(
  dat_ss_c_all %>% mutate(Test = "Compression"),
  dat_ss_t_all %>% mutate(Test = "Tension")
) %>%
  filter(alloy_cnt %in% pred_ids$alloy_cnt)%>%
  pivot_longer(
    cols = all_of(properties),
    names_to = "property",
    values_to = "val_true"
  ) %>%
  group_by(Test, alloy_cnt, property) %>%
  summarise(
    replications = n(),
    observed_var = var(val_true),
    sum_sq = sum((val_true - mean(val_true))^2),
    .groups = "drop"
  )

dat_compare <- dat_observed %>%
  inner_join(dat_het, by = c("Test", "alloy_cnt", "property")) %>%
  inner_join(dat_hom, by = c("Test", "alloy_cnt", "property")) %>%
  mutate(
    log_score_gain = -0.5 * (
      (replications - 1) * log(nugget_het / nugget_hom) +
        sum_sq * (1 / nugget_het - 1 / nugget_hom)
    )
  )
```




``` r
lambda_results <- dat_compare %>%
  group_by(Test, property) %>%
  summarise(
    compositions = n(),
    nugget_max_min = max(nugget_het) / min(nugget_het),
    score_gain_vs_fitted_hom = sum(log_score_gain),
    shape_gain_vs_best_constant =
      n() * log(mean(observed_var) /
                mean(observed_var / nugget_het)) -
      sum(log(nugget_het)),
    .groups = "drop"
  )

lambda_results %>%
  filter(Test == "Compression") %>%
  print(n = Inf, width = Inf)
```

```
## # A tibble: 7 × 6
##   Test        property       compositions nugget_max_min
##   <chr>       <chr>                 <int>          <dbl>
## 1 Compression bc_1                     40           1.00
## 2 Compression bc_2                     40           1.00
## 3 Compression bc_3                     40          19.9 
## 4 Compression bc_4                     40           7.53
## 5 Compression bc_5                     40           1.00
## 6 Compression bc_6                     40           1.00
## 7 Compression log_max_strain           40           1.00
##   score_gain_vs_fitted_hom shape_gain_vs_best_constant
##                      <dbl>                       <dbl>
## 1                   0.0813                   -0.00402 
## 2                   0.0763                    0.00181 
## 3                  15.3                       6.01    
## 4                   1.95                     -1.23    
## 5                   0.309                    -0.000755
## 6                  -2.71                     -0.00293 
## 7                 195.                        0.00229
```

``` r
lambda_results %>%
  filter(Test == "Tension") %>%
  print(n = Inf, width = Inf)
```

```
## # A tibble: 7 × 6
##   Test    property       compositions nugget_max_min score_gain_vs_fitted_hom
##   <chr>   <chr>                 <int>          <dbl>                    <dbl>
## 1 Tension bc_1                     40           1.00                 -0.569  
## 2 Tension bc_2                     40           1.00                  0.00457
## 3 Tension bc_3                     40           1.00                235.     
## 4 Tension bc_4                     40           1.00                 20.1    
## 5 Tension bc_5                     40           1.00                -21.6    
## 6 Tension bc_6                     40          14.4                  17.4    
## 7 Tension log_max_strain           40         575.                  238.     
##   shape_gain_vs_best_constant
##                         <dbl>
## 1                 -0.00419   
## 2                 -0.0000336 
## 3                  0.000144  
## 4                 -0.00126   
## 5                  0.00000157
## 6                 12.4       
## 7                -16.1
```









## Training

Training set


``` r
temp <- dat_ss_fil%>%
  ungroup()%>%
  filter(elements == 2,
         round(pmax(Mo,Nb,Ta,V,W),1) %in% c(0.1,0.3,0.5,0.7,0.9) 
         )%>%
  arrange(alloy_cnt)

#temp <- temp[-which(temp$alloy_cnt == 59),] #59 is the quanary with x =0.2 for all elements


id_train <- temp$alloy_cnt

dat_ss_train <- dat_ss_fil%>%
      filter((alloy_cnt %in% id_train))%>%
  arrange(alloy_cnt)


dat_bin_test <- dat_ss_fil%>%
            filter(! alloy_cnt %in% id_train,
                   elements == 2)

dat_train <- dat_ss_all%>%
      filter(alloy_cnt %in% id_train)%>%
      arrange(alloy_cnt)

dat_train_c <- dat_train%>%
      filter(type == "c")
    
dat_train_t <- dat_train%>%
      filter(type == "t")

out_nn_c <- compute_mean_var_fun(Y = as.matrix(dat_train_c[,properties]),
                              a_mult = rep(3,nrow(temp)))

out_nn_t <- compute_mean_var_fun(Y = as.matrix(dat_train_t[,properties]),
                              a_mult = rep(3,nrow(temp)))

dat_train <- dat_train%>%
                  pivot_longer(cols = all_of(properties),
                               names_to = "property",
                               values_to = "val_true")


dat_train_B <- dat_train%>%
    group_by(alloy_cnt, rep)%>%
  slice_head(n=1)
    
dat_train_c <- dat_train%>%
      filter(type == "c")
    
dat_train_t <- dat_train%>%
      filter(type == "t")
```

Candidates for active learning


``` r
dat_cand <- dat_ss_fil_all%>%
  filter(elements > 2,
         rep > 1)

ids_cand <- dat_cand$alloy_cnt

dat_ss_c_cand <- dat_ss_c_all%>%
                filter(alloy_cnt %in% ids_cand)

dat_ss_t_cand <- dat_ss_t_all%>%
                filter(alloy_cnt %in% ids_cand)

dat_rep_1 <- dat_ss_fil_all%>%
  filter(elements > 2,
         rep < 2)
```


``` r
X_ini <- as.matrix(dat_ss_train[,metals])
X_ini_B <- as.matrix(dat_train_B[,metals])

X_cand_ini <- as.matrix(dat_cand[,metals])


X_cand_ini_B <- as.matrix(dat_ss_c_cand[,metals])

dat_ss_c_cand <- dat_ss_c_cand%>%
                  pivot_longer(cols = all_of(properties),
                               names_to = "property",
                               values_to = "val_true")

dat_ss_t_cand <- dat_ss_t_cand%>%
                  pivot_longer(cols = all_of(properties),
                               names_to = "property",
                               values_to = "val_true")
```


## Predictions for Overall Performance

### Predict Models


``` r
kernels_v <- c("m32", "m52")
mean_type <- c("pol", "int")
```


``` r
np_id_h_all <- numeric(0)
for(mt in mean_type){
  for(kt in kernels_v){
    fn <- paste0("ALV2_26_prior_gam_np_h_", mt, "_", kt, ".RData")
    fp <- here::here("Results","Fixed Data",
                       fn)
    
    if(file.exists(fp)){
      load(fp)
      np_id_h_all <- c(np_id_h_all, 
                       np_id_h)
                       #np_id_h[-length(np_id_h)])
    }
    
  }
}
np_id_h_all <- sort(unique(np_id_h_all))
length(np_id_h_all)
```

```
## [1] 11
```

``` r
for(mt in mean_type){
  for(kt in kernels_v){
    fn <- paste0("ALV2_26_prior_gam_np_h_gen_sig_", mt, "_", kt, ".RData")
    fp <- here::here("Results","Fixed Data",
                       fn)
    
    if(file.exists(fp)){
      load(fp)
       np_id_h_all <- c(np_id_h_all, 
                        np_id_h)
                       #np_id_h[-length(np_id_h)])
    }
    
  }
}

np_id_h_all <- sort(unique(np_id_h_all))
length(np_id_h_all)
```

```
## [1] 18
```

``` r
id_rl <- dat_cand$alloy_cnt[np_id_h_all]

dat_pred <- rbind(dat_cand,
                  dat_bin_test,
                  dat_rep_1)%>%
            arrange(alloy_cnt)

X_pred <- dat_pred[,metals] |> as.matrix()
X_pred_irl <- ilr_bounded_minmax_additive(X_pred,
                                             epsilon = eps_ilr,
                                             z_min = xx$z_min,
                                             z_max = xx$z_max)$u

Basis_X_pred <- create_poly_basis(X_pred,
                                   n_out = n_prop,
                                   prop_names = properties,
                                   sparse_mat = T,
                               degree = 1,
                               omit_last_var = T,
                                   var_names = metals)
```



``` r
st <- Sys.time()
fn <- here::here("Results",
                 "Predictions_het.RData")

pred_type = "median" # mean or median

if(file.exists(fn)){
  load(fn)
}else{
  
dat_mu <- dat_lb <- dat_ub <- dat_var <- dat_nug <- dat_pm <- vector(mode = "list", 
                                                                     length = 2)

mods <- c("Tension", "Compression")


k <- 0
impse_dat <- vector(mode = "list")

for(mt in mean_type){
  BX <- NULL
  if(mt == "pol") BX <- Basis_X_pred
  message(Sys.time())
  for(kt in kernels_v){
    for(i in 0:5){
      mod_name <- paste(mt, kt, sep="_")
      file_name <- paste0("train_bin_rl_prior_gam_",
                          mt, "_", kt,
                          "_2026_05_19_RLV2_iter_",i,".RData")
      
      fp <- here::here("Results",
                       "Fixed Data",
                       file_name)
      
      if(file.exists(fp)){
        for(j in mods){
          k <- k + 1
           results <- import_predict_mod_gp_sampling_fun(newX = X_pred_irl,
                                                 Xdesign = dat_pred[,c(metals,"alloy_cnt")],
                                                 newX_pred = BX,
                                                 file_path = fp,
                                                 step_id = mod_name,
                                                 mod_id = j,
                                                 rseed = 1323,
                                                 n_funs = length(properties) - 1,
                                                 n_samples = 1e4,
                                                 prop_names = properties,
                                                 pred_type = pred_type)
           
           #pareto_mean <- pareto_scores_prop_fun(as.matrix(results$mp_preds$mean))
        
            gc()
           
           #results$mp_mean$pareto_score = pareto_mean$pareto_score
           dat_mu[[k]] <- results$gp_mu %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_var[[k]] <- results$gp_var %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_nug[[k]] <- results$gp_nug %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_pm[[k]] <- results$mp_mean %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_lb[[k]] <- results$mp_lb %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_ub[[k]] <- results$mp_ub %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           impse_dat[[k]] <- tibble(impse = results$imspe,
                                    iter = i,
                                    Model = mod_name,
                                    mean_type = mt,
                                    Kernel = kt,
                                    Type = j)
        }
      }
      
    }
  }
}


dat_mu <- bind_rows(dat_mu)
dat_var <- bind_rows(dat_var)
dat_nug <- bind_rows(dat_nug)
dat_pm <- bind_rows(dat_pm)
dat_lb <- bind_rows(dat_lb)
dat_ub <- bind_rows(dat_ub)
impse_dat <- bind_rows(impse_dat)

save(dat_mu,
     dat_var,
     dat_nug,
     dat_pm,
     dat_lb,
     dat_ub,
     impse_dat,
     file = fn)
ft <- Sys.time()
print(ft - st)
}
```


``` r
st <- Sys.time()

fn <- here::here("Results",
                 "Predictions_hom.RData")

if(file.exists(fn)){
  load(fn)
}else{
  

dat_mu_hom <- dat_lb_hom <- dat_ub_hom <- dat_var_hom <- dat_nug_hom <- dat_pm_hom <- vector(mode = "list", 
                                                                     length = 2)

mods <- c("Tension-Hom", "Compression-Hom")

k <- 0
impse_dat_hom <- vector(mode = "list")

for(mt in mean_type){
  BX <- NULL
  if(mt == "pol") BX <- Basis_X_pred
  message(Sys.time())
  for(kt in kernels_v){
    for(i in 0:5){
      mod_name <- paste(mt, kt, sep="_")
      file_name <- paste0("train_bin_rl_prior_gam_",
                          mt, "_", kt,
                          "_2026_05_19_RLV2_iter_",i,".RData")
      
      fp <- here::here("Results",
                       "Fixed Data",
                       file_name)
      
      if(file.exists(fp)){
        for(j in mods){
          k <- k + 1
           results <- import_predict_mod_gp_sampling_fun(newX = X_pred_irl,
                                                 Xdesign = dat_pred[,c(metals,"alloy_cnt")],
                                                 newX_pred = BX,
                                                 file_path = fp,
                                                 step_id = mod_name,
                                                 mod_id = j,
                                                 gp_type = "Hom",
                                                 rseed = 1323,
                                                 n_funs = length(properties) - 1,
                                                 n_samples = 1e4,
                                                 prop_names = properties,
                                                 pred_type = pred_type)
           
           #pareto_mean <- pareto_scores_prop_fun(as.matrix(results$mp_preds$mean))
        
            gc()
           
           #results$mp_mean$pareto_score = pareto_mean$pareto_score
           dat_mu_hom[[k]] <- results$gp_mu %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_var_hom[[k]] <- results$gp_var %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_nug_hom[[k]] <- results$gp_nug %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_pm_hom[[k]] <- results$mp_mean %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_lb_hom[[k]] <- results$mp_lb %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_ub_hom[[k]] <- results$mp_ub %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           impse_dat_hom[[k]] <- tibble(impse = results$imspe,
                                    iter = i,
                                    Model = mod_name,
                                    mean_type = mt,
                                    Kernel = kt,
                                    Type = j)
        }
      }
      
    }
  }
}


dat_mu_hom <- bind_rows(dat_mu_hom)
dat_var_hom <- bind_rows(dat_var_hom)
dat_nug_hom <- bind_rows(dat_nug_hom)
dat_pm_hom <- bind_rows(dat_pm_hom)
dat_lb_hom <- bind_rows(dat_lb_hom)
dat_ub_hom <- bind_rows(dat_ub_hom)
impse_dat_hom <- bind_rows(impse_dat_hom)


save(dat_mu_hom,
     dat_var_hom,
     dat_nug_hom,
     dat_pm_hom,
     dat_lb_hom,
     dat_ub_hom,
     impse_dat_hom,
     file = fn)
}


ft <- Sys.time()
print(ft - st)
```

```
## Time difference of 0.02070212 secs
```


``` r
st <- Sys.time()
fn <- here::here("Results",
                 "Predictions_lr_het.RData")

dat_mu_gen <- dat_lb_gen <- dat_ub_gen <- dat_var_gen <-
  dat_nug_gen <- dat_pm_gen <- vector(mode = "list",  length = 2)

mods <- c("Tension", "Compression")

pred_type = "median" # mean or median

k <- 0
impse_dat_gen <- vector(mode = "list")

if(file.exists(fn)){
  load(fn)
}else{

for(mt in mean_type){
  BX <- NULL
  if(mt == "pol") BX <- Basis_X_pred
  message(Sys.time())
  for(kt in kernels_v){
    for(i in 0:5){
      mod_name <- paste(mt, kt, sep="_")
      file_name <- paste0("train_bin_rl_low_rank_prior_gam_",
                          mt, "_", kt,
                          "_2026_09_21_RLV2_iter_",i,".RData")

      fp <- here::here("Results",
                       "Fixed Data",
                       file_name)

      if(file.exists(fp)){
        for(j in mods){
          k <- k + 1
           results <- import_predict_mod_gp_sampling_fun(newX = X_pred_irl,
                                                 Xdesign = dat_pred[,c(metals,"alloy_cnt")],
                                                 newX_pred = BX,
                                                 file_path = fp,
                                                 step_id = mod_name,
                                                 mod_id = j,
                                                 imspe_accuracy = 0.01,
                                                 rseed = 1323,
                                                 n_funs = length(properties) - 1,
                                                 n_samples = 1e4,
                                                 prop_names = properties,
                                                 pred_type = pred_type)

           #pareto_mean <- pareto_scores_prop_fun(as.matrix(results$mp_preds$mean))

            gc()

           #results$mp_mean$pareto_score = pareto_mean$pareto_score
           dat_mu_gen[[k]] <- results$gp_mu %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_var_gen[[k]] <- results$gp_var %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_nug_gen[[k]] <- results$gp_nug %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_pm_gen[[k]] <- results$mp_mean %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_lb_gen[[k]] <- results$mp_lb %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_ub_gen[[k]] <- results$mp_ub %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           impse_dat_gen[[k]] <- tibble(impse = results$imspe,
                                    iter = i,
                                    Model = mod_name,
                                    mean_type = mt,
                                    Kernel = kt,
                                    Type = j)
        }
      }

    }
  }
}


dat_mu_gen <- bind_rows(dat_mu_gen)
dat_var_gen <- bind_rows(dat_var_gen)
dat_nug_gen <- bind_rows(dat_nug_gen)
dat_pm_gen <- bind_rows(dat_pm_gen)
dat_lb_gen <- bind_rows(dat_lb_gen)
dat_ub_gen <- bind_rows(dat_ub_gen)
impse_dat_gen <- bind_rows(impse_dat_gen)

save(dat_mu_gen,
     dat_var_gen,
     dat_nug_gen,
     dat_pm_gen,
     dat_lb_gen,
     dat_ub_gen,
     impse_dat_gen,
     file = fn)
}
ft <- Sys.time()
print(ft - st)
```

```
## Time difference of 0.02436399 secs
```


``` r
st <- Sys.time()

dat_mu_hom_gen <- dat_lb_hom_gen <- dat_ub_hom_gen <- dat_var_hom_gen <-
  dat_nug_hom_gen <- dat_pm_hom_gen <- vector(mode = "list", length = 2)

mods <- c("Tension-Hom", "Compression-Hom")

k <- 0
impse_dat_hom_gen <- vector(mode = "list")

fn <- here::here("Results",
                 "Predictions_lr_hom.RData")


if(file.exists(fn)){
  load(fn)
}else{
  
for(mt in mean_type){
  BX <- NULL
  if(mt == "pol") BX <- Basis_X_pred
  message(Sys.time())
  for(kt in kernels_v){
    for(i in 0:5){
      mod_name <- paste(mt, kt, sep="_")
      file_name <- paste0("train_bin_rl_low_rank_prior_gam_",
                          mt, "_", kt,
                          "_2026_09_21_RLV2_iter_",i,".RData")

      fp <- here::here("Results",
                       "Fixed Data",
                       file_name)

      if(file.exists(fp)){
        for(j in mods){
          k <- k + 1
           results <- import_predict_mod_gp_sampling_fun(newX = X_pred_irl,
                                                 Xdesign = dat_pred[,c(metals,"alloy_cnt")],
                                                 newX_pred = BX,
                                                 file_path = fp,
                                                 step_id = mod_name,
                                                 mod_id = j,
                                                 gp_type = "Hom",
                                                 rseed = 1323,
                                                 n_funs = length(properties) - 1,
                                                 n_samples = 1e4,
                                                 imspe_accuracy = 0.01,
                                                 prop_names = properties,
                                                 pred_type = pred_type)

           #pareto_mean <- pareto_scores_prop_fun(as.matrix(results$mp_preds$mean))

            gc()

           #results$mp_mean$pareto_score = pareto_mean$pareto_score
           dat_mu_hom_gen[[k]] <- results$gp_mu %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_var_hom_gen[[k]] <- results$gp_var %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_nug_hom_gen[[k]] <- results$gp_nug %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_pm_hom_gen[[k]] <- results$mp_mean %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_lb_hom_gen[[k]] <- results$mp_lb %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           dat_ub_hom_gen[[k]] <- results$mp_ub %>% mutate(iter = i,
                                                   Kernel = kt,
                                                   mean_type = mt)
           impse_dat_hom_gen[[k]] <- tibble(impse = results$imspe,
                                    iter = i,
                                    Model = mod_name,
                                    mean_type = mt,
                                    Kernel = kt,
                                    Type = j)
        }
      }

    }
  }
}


dat_mu_hom_gen <- bind_rows(dat_mu_hom_gen)
dat_var_hom_gen <- bind_rows(dat_var_hom_gen)
dat_nug_hom_gen <- bind_rows(dat_nug_hom_gen)
dat_pm_hom_gen <- bind_rows(dat_pm_hom_gen)
dat_lb_hom_gen <- bind_rows(dat_lb_hom_gen)
dat_ub_hom_gen <- bind_rows(dat_ub_hom_gen)
impse_dat_hom_gen <- bind_rows(impse_dat_hom_gen)

ft <- Sys.time()
print(ft - st)

save(dat_mu_hom_gen,
     dat_var_hom_gen,
     dat_nug_hom_gen,
     dat_pm_hom_gen,
     dat_lb_hom_gen,
     dat_ub_hom_gen,
     impse_dat_hom_gen,
     file = fn)

}
```

### Comparison MAPE


``` r
dat_sim <- bind_rows(dat_ss_c_all%>%
                           mutate(Test = "Compression"),
                         dat_ss_t_all%>%
                           mutate(Test = "Tension"))  %>%
  filter(alloy_cnt %in% dat_pred$alloy_cnt)%>%
  mutate(ultimate_stress = exp(log_max_stress),
         ultimate_strain = exp(log_max_strain))%>%
  pivot_longer(cols= all_of(props),
               values_to = "val",
               names_to = "property")
```



``` r
dat_comp <- dat_pm_hom%>%
            mutate(Test = ifelse(Test == "Tension-Hom", 
                                 "Tension", "Compression" ),
                   Model = "Hom",
                   Sig = "Diag")%>%
            bind_rows(dat_pm%>%
                        mutate(Model = "Het",
                               Sig = "Diag"),
                      dat_pm_gen%>%
                        mutate(Model = "Het",
                               Sig = "Low-Rank"),
                      dat_pm_hom_gen%>%
                          mutate(Test = ifelse(Test == "Tension-Hom", 
                                               "Tension", "Compression" ),
                                 Model = "Hom",
                                 Sig = "Low-Rank")
                      )%>%
           pivot_longer(cols= all_of(props),
               values_to = "pred",
               names_to = "property")%>%
            inner_join(dat_sim%>%
                         select(alloy_cnt,Test,property,val)%>%
                         group_by(alloy_cnt,Test,property)%>%
                         summarise(val = mean(val)),
                       by = join_by(alloy_cnt,Test,property))%>%
  mutate(elements= (Mo>0) + (Nb>0) + (Ta>0) + (V>0) + (W>0))
```

```
## `summarise()` has grouped output by 'alloy_cnt', 'Test'. You can override using
## the `.groups` argument.
```



``` r
mod_names_nn <- c("int_m32",
                  "pol_m32",
                  "int_m52",
               "pol_m52",
               "DL")
mod_labels_nn <- c("M3/2",
                  "Lin M3/2",
                  "M5/2",
               "Lin M5/2",
               "Experts")

mod_labels_ltx <- unname(TeX(c("M$3/2$","Lin M$3/2$",
                         "M$5/2$","Lin M$5/2$","Experts")))

dat_comp%>%
  filter(elements < 3,
         iter < 1,
         Step %in% mod_names_nn,
         Model == "Het",
         Test == "Compression")%>%
  mutate(Test = factor(Test,
                       levels = c("Tension","Compression")),
         Step = factor(Step,
                       levels = mod_names_nn,
                       labels = mod_labels_ltx),
         Kernel = factor(Kernel,
                         levels = kernels_v,
                         labels = c("M3/2", "M5/2")),
         mean_type = factor(mean_type,
                            levels = c("int","pol"),
                            labels = c("Constant", "Linear")))%>%
  filter(!(Test == "Compression" & property =="toughness"))%>%
  mutate(rer = (pred-val)/val,
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  ggplot(aes(y= 100*rer,
             colour = Sig))+
  geom_boxplot(aes(x = Sig))+
  scale_color_brewer("Model", palette = "Dark2")+
  scale_x_discrete(labels = rep(" ", 5))+
  ggh4x::facet_grid2(Step~property,
             scales = "free",
             independent = "y",
             labeller = label_parsed
             )+
  labs(y = "Relative Error (%)",
       x = "Model",
       title = "Accuracy on Binary Composition Alloys-Compression")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))
```

![](Predictions_19-05---Copy_files/figure-html/acc_bin_comp-1.png)<!-- -->


``` r
mod_names_nn <- c("int_m32",
                  "pol_m32",
                  "int_m52",
               "pol_m52",
               "DL")
mod_labels_nn <- c("M3/2",
                  "Lin M3/2",
                  "M5/2",
               "Lin M5/2",
               "Experts")

mod_labels_ltx <- unname(TeX(c("M$3/2$","Lin M$3/2$",
                         "M$5/2$","Lin M$5/2$","Experts")))

dat_comp%>%
  filter(elements < 3,
         iter < 1,
         Step %in% mod_names_nn,
         Model == "Het",
         Test == "Tension")%>%
  mutate(Test = factor(Test,
                       levels = c("Tension","Compression")),
         Step = factor(Step,
                       levels = mod_names_nn,
                       labels = mod_labels_ltx),
         Kernel = factor(Kernel,
                         levels = kernels_v,
                         labels = c("M3/2", "M5/2")),
         mean_type = factor(mean_type,
                            levels = c("int","pol"),
                            labels = c("Constant", "Linear")))%>%
  filter(!(Test == "Compression" & property =="toughness"))%>%
  mutate(rer = (pred-val)/val,
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  ggplot(aes(y= 100*rer,
             colour = Sig))+
  geom_boxplot(aes(x = Sig))+
  scale_color_brewer("Model", palette = "Dark2")+
  scale_x_discrete(labels = rep(" ", 5))+
  ggh4x::facet_grid2(Step~property,
             scales = "free",
             independent = "y",
             labeller = label_parsed
             )+
  labs(y = "Relative Error (%)",
       x = "Model",
       title = "Accuracy on Binary Composition Alloys-Tension")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))
```

![](Predictions_19-05---Copy_files/figure-html/acc_bin_tension-1.png)<!-- -->


``` r
mod_names_nn <- c("int_m32",
                  "pol_m32",
                  "int_m52",
               "pol_m52",
               "DL")
mod_labels_nn <- c("M3/2",
                  "Lin M3/2",
                  "M5/2",
               "Lin M5/2",
               "Experts")

mod_labels_ltx <- unname(TeX(c("M$3/2$","Lin M$3/2$",
                         "M$5/2$","Lin M$5/2$","Experts")))

dat_comp%>%
  filter(elements < 3,
         iter < 1,
         Step %in% mod_names_nn,
         Sig == "Diag",
         Test == "Tension")%>%
  mutate(Test = factor(Test,
                       levels = c("Tension","Compression")),
         Step = factor(Step,
                       levels = mod_names_nn,
                       labels = mod_labels_ltx),
         Kernel = factor(Kernel,
                         levels = kernels_v,
                         labels = c("M3/2", "M5/2")),
         mean_type = factor(mean_type,
                            levels = c("int","pol"),
                            labels = c("Constant", "Linear")))%>%
  filter(!(Test == "Compression" & property =="toughness"))%>%
  mutate(rer = (pred-val)/val,
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  ggplot(aes(y= 100*rer,
             colour = Model))+
  geom_boxplot(aes(x = Model))+
  scale_color_brewer("Model", palette = "Dark2")+
  scale_x_discrete(labels = rep(" ", 5))+
  ggh4x::facet_grid2(Step~property,
             scales = "free",
             independent = "y",
             labeller = label_parsed
             )+
  labs(y = "Relative Error (%)",
       x = "Model",
       title = "Accuracy on Binary Composition Alloys-Tension")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))
```

![](Predictions_19-05---Copy_files/figure-html/acc_bin_hom_het_ten-1.png)<!-- -->


``` r
mod_names_nn <- c("int_m32",
                  "pol_m32",
                  "int_m52",
               "pol_m52",
               "DL")
mod_labels_nn <- c("M3/2",
                  "Lin M3/2",
                  "M5/2",
               "Lin M5/2",
               "Experts")

mod_labels_ltx <- unname(TeX(c("M$3/2$","Lin M$3/2$",
                         "M$5/2$","Lin M$5/2$","Experts")))

dat_comp%>%
  filter(elements < 3,
         iter < 1,
         Step %in% mod_names_nn,
         Sig == "Diag",
         Test != "Tension")%>%
  mutate(Test = factor(Test,
                       levels = c("Tension","Compression")),
         Step = factor(Step,
                       levels = mod_names_nn,
                       labels = mod_labels_ltx),
         Kernel = factor(Kernel,
                         levels = kernels_v,
                         labels = c("M3/2", "M5/2")),
         mean_type = factor(mean_type,
                            levels = c("int","pol"),
                            labels = c("Constant", "Linear")))%>%
  filter(!(Test == "Compression" & property =="toughness"))%>%
  mutate(rer = (pred-val)/val,
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  ggplot(aes(y= 100*rer,
             colour = Model))+
  geom_boxplot(aes(x = Model))+
  scale_color_brewer("Model", palette = "Dark2")+
  scale_x_discrete(labels = rep(" ", 5))+
  ggh4x::facet_grid2(Step~property,
             scales = "free",
             independent = "y",
             labeller = label_parsed
             )+
  labs(y = "Relative Error (%)",
       x = "Model",
       title = "Accuracy on Binary Composition Alloys-Compression")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))
```

![](Predictions_19-05---Copy_files/figure-html/acc_bin_hom_het_comp-1.png)<!-- -->




``` r
dat_comp%>%
  filter(elements > 2,
         ! alloy_cnt %in% id_rl,
         Model == "Het",
         Test == "Tension")%>%
  mutate(Test = factor(Test,
                       levels = c("Tension","Compression")),
         Kernel = factor(Kernel,
                         levels = kernels_v,
                         labels = c("M3/2", "M5/2")),
         mean_type = factor(mean_type,
                            levels = c("int","pol"),
                            labels = c("Constant", "Linear")))%>%
  filter(!(Test == "Compression" & property =="toughness"))%>%
  mutate(rer = abs(pred-val)/val,
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  group_by(mean_type,Kernel, Sig, iter,Test,property)%>%
  summarise(rer= mean(rer))%>%
  ggplot(aes(x= iter,
             y= 100*rer,
             colour = Kernel,
             linetype = Sig))+
  geom_line(linewidth = 1.25)+
  # scale_color_brewer("Property",
  #                    labels = prop_ltx,
  #                    palette = "Dark2")+
  scale_color_brewer("Kernel", palette = "Dark2")+
  ggh4x::facet_grid2(mean_type~property,
             scales = "free",
             independent = "y",
             labeller = label_parsed)+
  #geom_hline(yintercept = 0, colour = "red",linetype = "dashed")+
  labs(y = "MAPE (%)",
       x = "Iteration",
       linetype = "Covariance",
       title="Accuracy of MVGPR on General Composition Alloys (Tension)")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))+ 
    guides(colour = guide_legend(nrow = 1))
```

```
## `summarise()` has grouped output by 'mean_type', 'Kernel', 'Sig', 'iter',
## 'Test'. You can override using the `.groups` argument.
```

![](Predictions_19-05---Copy_files/figure-html/mape_al_tension-1.png)<!-- -->



``` r
dat_comp%>%
  filter(elements > 2,
         ! alloy_cnt %in% id_rl,
         Model == "Hom",
         Test == "Compression")%>%
  mutate(Test = factor(Test,
                       levels = c("Tension","Compression")),
         Kernel = factor(Kernel,
                         levels = kernels_v,
                         labels = c("M3/2", "M5/2")),
         mean_type = factor(mean_type,
                            levels = c("int","pol"),
                            labels = c("Constant", "Linear")))%>%
  filter(!(Test == "Compression" & property =="toughness"))%>%
  mutate(rer = abs(pred-val)/val,
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  group_by(mean_type,Kernel, Sig, iter,Test,property)%>%
  summarise(rer= mean(rer))%>%
  ggplot(aes(x= iter,
             y= 100*rer,
             colour = Kernel,
             linetype = Sig))+
  geom_line(linewidth = 1.25)+
  # scale_color_brewer("Property",
  #                    labels = prop_ltx,
  #                    palette = "Dark2")+
  scale_color_brewer("Kernel", palette = "Dark2")+
  facet_grid(mean_type~property,
             scales = "free",
             labeller = label_parsed)+
  geom_hline(yintercept = 0, colour = "red",
             linetype = "dashed")+
  labs(y = "MAPE (%)",
       x = "Iteration",
       linetype = "Covariance",
       title="Accuracy of MVGPR on General Composition Alloys (Compression)")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))+ 
    guides(colour = guide_legend(nrow = 1))
```

```
## `summarise()` has grouped output by 'mean_type', 'Kernel', 'Sig', 'iter',
## 'Test'. You can override using the `.groups` argument.
```

![](Predictions_19-05---Copy_files/figure-html/mape_al_hom-1.png)<!-- -->


``` r
dat_comp%>%
  filter(elements > 2,
         ! alloy_cnt %in% id_rl,
         Model == "Hom",
         Test == "Tension")%>%
  mutate(Test = factor(Test,
                       levels = c("Tension","Compression")),
         Kernel = factor(Kernel,
                         levels = kernels_v,
                         labels = c("M3/2", "M5/2")),
         mean_type = factor(mean_type,
                            levels = c("int","pol"),
                            labels = c("Constant", "Linear")))%>%
  filter(!(Test == "Compression" & property =="toughness"))%>%
  mutate(rer = abs(pred-val)/val,
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  group_by(mean_type,Kernel, Sig, iter,Test,property)%>%
  summarise(rer= mean(rer))%>%
  ggplot(aes(x= iter,
             y= 100*rer,
             colour = Kernel,
             linetype = Sig))+
  geom_line(linewidth = 1.25)+
  # scale_color_brewer("Property",
  #                    labels = prop_ltx,
  #                    palette = "Dark2")+
  scale_color_brewer("Kernel", palette = "Dark2")+
  facet_grid(mean_type~property,
             scales = "free",
             labeller = label_parsed)+
  geom_hline(yintercept = 0, colour = "red",
             linetype = "dashed")+
  labs(y = "MAPE (%)",
       x = "Iteration",
       linetype = "Covariance",
       title="Accuracy of MVGPR on General Composition Alloys (Compression)")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))+ 
    guides(colour = guide_legend(nrow = 1))
```

```
## `summarise()` has grouped output by 'mean_type', 'Kernel', 'Sig', 'iter',
## 'Test'. You can override using the `.groups` argument.
```

![](Predictions_19-05---Copy_files/figure-html/mape_al_hom_ten-1.png)<!-- -->



## IMSPE


``` r
impse_dat_all <- impse_dat%>%
  mutate(Sig = "Diag")%>%
  bind_rows(impse_dat_hom%>%
              mutate(Sig = "Diag"),
            impse_dat_gen%>%
              mutate(Sig = "Gen"),
            impse_dat_hom%>%
              mutate(Sig = "Gen"))%>%
  #filter(impse < 1e4)%>%
  mutate(Model_Type = factor(Type,
                       levels  = c("Tension","Compression",
                                   "Tension-Hom","Compression-Hom"),
                       labels = c("Het","Het",
                                  "Hom","Hom")),
         Test_Type = factor(Type,
                       levels  = c("Tension","Compression",
                                   "Tension-Hom","Compression-Hom"),
                       labels = c("Tension","Compression",
                                  "Tension","Compression")),
         Kernel = factor(Kernel,
                         levels = kernels_v,
                         labels = c("M3/2", "M5/2")),
         mean_type = factor(mean_type,
                            levels = c("int","pol"),
                            labels = c("Constant", "Linear")))
```




``` r
impse_dat_all%>%
  filter(Sig == "Diag")%>%
  ggplot(aes(x= iter,
             y= impse,
             colour = Model_Type,
             linetype = Kernel))+
  geom_line(linewidth = 1.25)+
  # scale_color_brewer("Property",
  #                    labels = prop_ltx,
  #                    palette = "Dark2")+
  scale_color_brewer("Error", palette = "Dark2")+
  ggh4x::facet_grid2(Test_Type~mean_type,
             scales = "free",
             independent = "y")+
  labs(y = "IMSPE",
       x = "Iteration",
       linetype = "Kernel",
       title="IMSPE During the AL Algorithm")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))+ 
    guides(colour = guide_legend(nrow = 1))
```

![](Predictions_19-05---Copy_files/figure-html/imspe_al-1.png)<!-- -->


``` r
impse_dat_all%>%
  filter(Sig == "Gen")%>%
  ggplot(aes(x= iter,
             y= impse,
             colour = Model_Type,
             linetype = Kernel))+
  geom_line(linewidth = 1.25)+
  # scale_color_brewer("Property",
  #                    labels = prop_ltx,
  #                    palette = "Dark2")+
  scale_color_brewer("Error", palette = "Dark2")+
  ggh4x::facet_grid2(Test_Type~mean_type,
             scales = "free",
             independent = "y")+
  labs(y = "IMSPE",
       x = "Iteration",
       linetype = "Kernel",
       title="IMSPE During the AL Algorithm")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))+ 
    guides(colour = guide_legend(nrow = 1))
```

![](Predictions_19-05---Copy_files/figure-html/imspe_al_gen-1.png)<!-- -->



``` r
p1 <- impse_dat%>%
  bind_rows(impse_dat_hom)%>%
  mutate(Model_Type = factor(Type,
                       levels  = c("Tension","Compression",
                                   "Tension-Hom","Compression-Hom"),
                       labels = c("Het","Het",
                                  "Hom","Hom")),
         Test_Type = factor(Type,
                       levels  = c("Tension","Compression",
                                   "Tension-Hom","Compression-Hom"),
                       labels = c("Tension","Compression",
                                  "Tension","Compression")),
         Kernel = factor(Kernel,
                         levels = kernels_v,
                         labels = c("M3/2", "M5/2")),
         mean_type = factor(mean_type,
                            levels = c("int","pol"),
                            labels = c("Constant", "Linear")))%>%
  group_by(Model_Type, Kernel, mean_type, iter)%>%
  summarise(impse = sum(impse))%>%
  ggplot(aes(x= iter,
             y= log(impse),
             colour = Kernel))+
  geom_line(linewidth = 1.25)+
  #scale_y_log10() +
  # scale_color_brewer("Property",
  #                    labels = prop_ltx,
  #                    palette = "Dark2")+
  scale_color_brewer("Kernel", palette = "Dark2")+
  ggh4x::facet_grid2(Model_Type~mean_type,
             scales = "free",
             independent = "y")+
  labs(y = "Log(IMSPE)",
       x = "Iteration",
       title="IMSPE During the AL Algorithm")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))+ 
    guides(colour = guide_legend(nrow = 1))
```

```
## `summarise()` has grouped output by 'Model_Type', 'Kernel', 'mean_type'. You
## can override using the `.groups` argument.
```

``` r
p1 
```

![](Predictions_19-05---Copy_files/figure-html/imspe_al_sum-1.png)<!-- -->



### Predictions NN





``` r
nn_bin_t <- read.csv(here::here("Results", "NN_Experts",
                                 "dat_test_bin_even_t_predictions.csv"))
nn_bin_c <- read.csv(here::here("Results", "NN_Experts",
                                 "dat_test_bin_even_c_predictions.csv"))
nn_qq_t <- read.csv(here::here("Results", "NN_Experts",
                                 "dat_test_qq_t_predictions.csv"))
nn_qq_c <- read.csv(here::here("Results", "NN_Experts",
                                 "dat_test_qq_c_predictions.csv"))

mu_ids <- 7 + 1:7
var_ids <- max(mu_ids) + 1:7
alloy_cnt_col <- 6

gp_mu_nn_t <- bind_rows(nn_bin_t[,c(alloy_cnt_col,mu_ids)],
                        nn_qq_t[,c(alloy_cnt_col,mu_ids)])%>%
              arrange(alloy_cnt)
gp_mu_nn_t <- as.matrix(gp_mu_nn_t[,-1])
gp_mu_nn_t <- (gp_mu_nn_t%r*% out_nn_t$rs_mult) %r+% out_nn_t$rs_mean


gp_var_nn_t<- bind_rows(nn_bin_t[,c(alloy_cnt_col,var_ids)],
                        nn_qq_t[,c(alloy_cnt_col,var_ids)])%>%
              arrange(alloy_cnt)
gp_var_nn_t <- exp(as.matrix(gp_var_nn_t[,-1])) %r*%
                    out_nn_t$rs_mult^2


gp_mu_nn_c <- bind_rows(nn_bin_c[,c(alloy_cnt_col,mu_ids)],
                        nn_qq_c[,c(alloy_cnt_col,mu_ids)])%>%
              arrange(alloy_cnt)
gp_mu_nn_c <- as.matrix(gp_mu_nn_c[,-1])
gp_mu_nn_c <- (gp_mu_nn_c%r*% out_nn_c$rs_mult) %r+% out_nn_c$rs_mean

gp_var_nn_c<- bind_rows(nn_bin_c[,c(alloy_cnt_col,var_ids)],
                        nn_qq_c[,c(alloy_cnt_col,var_ids)])%>%
              arrange(alloy_cnt)
gp_var_nn_c <- exp(as.matrix(gp_var_nn_c[,-1]))  %r*%
                    out_nn_c$rs_mult^2
```



``` r
res_nn_t <- predict_mod_gp_sampling_fun(Xdesign = dat_pred[,c(metals,"alloy_cnt")],
                                            gp_mu = gp_mu_nn_t,
                                           gp_var = gp_var_nn_t,
                                           step_id = "NN",
                                           test_type = "Tension",
                                           rseed = 1323,
                                           n_funs = length(properties) - 1,
                                           n_samples = 1e4,
                                           prop_names = properties,
                                                 pred_type = pred_type)

res_nn_c <- predict_mod_gp_sampling_fun(Xdesign = dat_pred[,c(metals,"alloy_cnt")],
                                            gp_mu = gp_mu_nn_c,
                                           gp_var = gp_var_nn_c,
                                           step_id = "NN",
                                           test_type = "Compression",
                                           rseed = 1323,
                                           n_samples = 1e4,
                                           n_funs = length(properties) - 1,
                                           prop_names = properties,
                                                 pred_type = pred_type)
```




``` r
dl_name <- "CL-REMoE"

dat_comp_nn <- bind_rows(res_nn_c$mp_mean ,
                          res_nn_t$mp_mean )%>% 
              mutate(iter = 0,
                     Step = dl_name,
                     Kernel = dl_name,
                     mean_type = dl_name,
                     Model = dl_name,
                     elements= (Mo>0) + (Nb>0) + (Ta>0) + (V>0) + (W>0))%>%
           pivot_longer(cols= all_of(props),
               values_to = "pred",
               names_to = "property")%>%
            inner_join(dat_sim%>%
                         select(alloy_cnt,Test,property,val)%>%
                         group_by(alloy_cnt,Test,property)%>%
                         summarise(val = mean(val)),
                       by = join_by(alloy_cnt,Test,property))
```

```
## `summarise()` has grouped output by 'alloy_cnt', 'Test'. You can override using
## the `.groups` argument.
```

``` r
dat_comp_nn <- dat_comp_nn %>%
                bind_rows(dat_comp)
```


``` r
prop_ltx <- unname(TeX(c("$\\sigma_u$","$\\epsilon_f$",
                         "$U_T$","$E$")))
mod_names_nn <- c("int_m32",
                  "pol_m32",
                  "int_m52",
               "pol_m52",
               dl_name)
mod_labels_nn <- c("M3/2",
                  "Lin M3/2",
                  "M5/2",
               "Lin M5/2",
               dl_name)

p1 <- dat_comp_nn%>%
  filter(Step %in% mod_names_nn,
         Model %in% c("Het", dl_name),
         iter < 1,
         elements < 3)%>%
  filter(!(Test == "Compression" & property =="toughness"))%>%
  mutate(rer = (pred-val)/val,
         Step = factor(Step,
                       levels = mod_names_nn,
                       labels = mod_labels_nn),
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  ggplot(aes(y= 100*rer,
             colour = Step))+
  geom_boxplot(aes(x = Step))+
  scale_color_brewer("Model", palette = "Dark2")+
  scale_x_discrete(labels = rep(" ", 5))+
  facet_grid(Test~property,
             scales = "free",
             labeller = label_parsed)+
  labs(y = "Relative Error (%)",
       x = "Model",
       title = "Accuracy on Even Binary Alloys")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))

p1
```

![](Predictions_19-05---Copy_files/figure-html/acc_bin_mods_nn-1.png)<!-- -->


``` r
p1 <- dat_comp_nn%>%
  filter(Step %in% mod_names_nn,
         Model %in% c("Het", dl_name),
         iter < 1,
         elements > 2)%>%
  filter(!(Test == "Compression" & property =="toughness"))%>%
  mutate(rer = (pred-val)/val,
         Step = factor(Step,
                       levels = mod_names_nn,
                       labels = mod_labels_nn),
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  ggplot(aes(y= 100*rer,
             colour = Step))+
  geom_boxplot(aes(x = Step))+
  scale_color_brewer("Model", palette = "Dark2")+
  scale_x_discrete(labels = rep(" ", 5))+
  facet_grid(Test~property,
             scales = "free",
             labeller = label_parsed)+
  labs(y = "Relative Error (%)",
       x = "Model",
       title = "Accuracy on General Composition Alloys")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10))

p1
```

![](Predictions_19-05---Copy_files/figure-html/acc_qq_mods_nn-1.png)<!-- -->



``` r
dat_comp_nn%>%
  filter(Step %in% mod_names_nn,
         Model %in% c("Het", dl_name),
         iter < 1,
         #! alloy_cnt %in% id_rl,
         elements > 2)%>%
  mutate(rer = abs(pred-val)/val,
         Step = factor(Step,
                       levels = mod_names_nn,
                       labels = mod_labels_nn),
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  group_by(Model, Test, Step,property)%>%
  summarise(rer= 100*mean(rer))%>%
  pivot_wider(names_from = "property",
              values_from = "rer")%>%
  arrange(Test)%>%
  kableExtra::kbl("latex",booktabs = T,digits = 2)
```

```
## `summarise()` has grouped output by 'Model', 'Test', 'Step'. You can override
## using the `.groups` argument.
```


\begin{tabular}[t]{lllrrrr}
\toprule
Model & Test & Step & sigma[u] & epsilon[f] & U[T] & E\\
\midrule
CL-REMoE & Compression & CL-REMoE & 5.33 & 6.60 & 10.99 & 8.77\\
Het & Compression & M3/2 & 18.74 & 14.29 & 30.87 & 13.71\\
Het & Compression & Lin M3/2 & 13.96 & 12.86 & 26.59 & 9.93\\
Het & Compression & M5/2 & 22.11 & 27.60 & 50.51 & 12.48\\
Het & Compression & Lin M5/2 & 18.85 & 20.27 & 39.47 & 10.70\\
\addlinespace
CL-REMoE & Tension & CL-REMoE & 6.51 & 9.09 & 8.76 & 8.47\\
Het & Tension & M3/2 & 12.49 & 17.35 & 16.71 & 15.87\\
Het & Tension & Lin M3/2 & 23.77 & 13.67 & 24.81 & 17.32\\
Het & Tension & M5/2 & 16.52 & 22.51 & 31.08 & 15.71\\
Het & Tension & Lin M5/2 & 19.18 & 15.90 & 27.85 & 37.63\\
\bottomrule
\end{tabular}



``` r
dat_comp_nn%>%
  filter(Step %in% mod_names_nn,
         Model %in% c("Het", dl_name),
         iter < 1,
         elements <3)%>%
  mutate(rer = abs(pred-val)/val,
         Step = factor(Step,
                       levels = mod_names_nn,
                       labels = mod_labels_nn),
         property = factor(property,
                           levels = props,
                           labels = prop_ltx))%>%
  group_by(Model, Test, Step,property)%>%
  summarise(rer= 100*mean(rer))%>%
  pivot_wider(names_from = "property",
              values_from = "rer")%>%
  arrange(Test)%>%
  kableExtra::kbl("latex",booktabs = T,digits = 2)
```

```
## `summarise()` has grouped output by 'Model', 'Test', 'Step'. You can override
## using the `.groups` argument.
```


\begin{tabular}[t]{lllrrrr}
\toprule
Model & Test & Step & sigma[u] & epsilon[f] & U[T] & E\\
\midrule
CL-REMoE & Compression & CL-REMoE & 4.33 & 4.22 & 6.88 & 5.55\\
Het & Compression & M3/2 & 1.22 & 6.79 & 7.46 & 6.73\\
Het & Compression & Lin M3/2 & 1.48 & 1.09 & 2.18 & 1.41\\
Het & Compression & M5/2 & 8.37 & 4.51 & 16.76 & 4.65\\
Het & Compression & Lin M5/2 & 1.11 & 2.73 & 3.14 & 2.80\\
\addlinespace
CL-REMoE & Tension & CL-REMoE & 2.68 & 2.98 & 3.96 & 3.35\\
Het & Tension & M3/2 & 8.82 & 1.19 & 9.39 & 7.41\\
Het & Tension & Lin M3/2 & 1.85 & 1.27 & 2.37 & 1.86\\
Het & Tension & M5/2 & 1.43 & 1.32 & 2.75 & 1.50\\
Het & Tension & Lin M5/2 & 1.92 & 5.62 & 6.08 & 6.17\\
\bottomrule
\end{tabular}


## Binary Plots

### MoNb


``` r
in_pr <- c(0.2,0.4,0.6,0.8)
alloys_in <- c("Nb","Mo", "MoNb")
dat_sb <- dat_ss_t%>%
  mutate(Training = ifelse(Nb %in% in_pr,"Test","Train"))%>%
  filter(alloy %in% alloys_in)

x <-seq(0,1, by = 0.025)
X_design <- tibble(Mo = x,
                   Nb = 1-x,
                   Ta = 0,
                   V = 0,
                   W = 0)
X_pred <-as.matrix(X_design)

X_pred_irl <- ilr_bounded_minmax_additive(X_pred,
                                             epsilon = eps_ilr,
                                             z_min = xx$z_min,
                                             z_max = xx$z_max)$u

BX <- create_poly_basis(X_pred,
                                   n_out = n_prop,
                                   prop_names = properties,
                                   sparse_mat = T,
                               degree = 1,
                               omit_last_var = T,
                                   var_names = metals)

X_design <- X_design%>%
  unite("alloy_id",
        Mo:W,
        remove = FALSE)

steps_in <- c(1,6)
n_files <- length(steps_in)
```


``` r
dat_mu <- dat_lb <- dat_ub <- dat_var <- dat_nug <- dat_pm <- vector(mode = "list", length = n_files)
imspe_v <- rep(0,n_files)
for(i in steps_in){
  fn <- paste0("train_bin_rl_prior_gam_pol_m52_2026_05_19_RLV2_iter_",i-1,".RData")
  fp <- here::here("Results","Fixed Data",fn)
  
  results <- import_predict_mod_gp_sampling_fun(newX = X_pred_irl,
                                                 Xdesign = X_design,
                                                 newX_pred = BX,
                                           file_path = fp,
                                           n_funs = length(properties) - 1,
                                           step_id = i-1,
                                           mod_id = 1,
                                           rseed = 1323,
                                           n_samples = 1e4,
                                           prop_names = properties,
                                           pred_type = pred_type)
  
 dat_mu[[i]] <- results$gp_mu
 dat_var[[i]] <- results$gp_var
 dat_nug[[i]] <- results$gp_nug
 dat_pm[[i]] <- results$mp_mean
 dat_lb[[i]] <- results$mp_lb
 dat_ub[[i]] <- results$mp_ub
 imspe_v[i] <- results$imspe
}


dat_mu <- bind_rows(dat_mu)
dat_var <- bind_rows(dat_var)
dat_nug <- bind_rows(dat_nug)
dat_pm <- bind_rows(dat_pm)
dat_lb <- bind_rows(dat_lb)
dat_ub <- bind_rows(dat_ub)
```


``` r
dat_lb2 <- dat_lb%>%
  pivot_longer(all_of(props),
               values_to = "val",
               names_to = "property")%>%
  mutate(property = factor(property, levels = props,
                           labels = prop_ltx))

dat_ub2 <- dat_ub%>%
  pivot_longer(all_of(props),
               values_to = "val",
               names_to = "property")%>%
  mutate(property = factor(property, levels = props,
                           labels = prop_ltx))

dat_pm2 <- dat_pm%>%
  pivot_longer(all_of(props),
               values_to = "val",
               names_to = "property")%>%
  mutate(Step = factor(Step),
         ub = dat_ub2$val,
         lb = dat_lb2$val)%>%
  mutate(property = factor(property, levels = props,
                           labels = prop_ltx))

dat_sb2 <- dat_sb%>%
  mutate(ultimate_stress = exp(log_max_stress),
         ultimate_strain = exp(log_max_strain))%>%
  pivot_longer(all_of(props),
               values_to = "val",
               names_to = "property")%>%
  mutate(property = factor(property, levels = props,
                           labels = prop_ltx))

dat_pm2%>%
  ggplot(aes(x= Nb,
             y = val))+
  geom_line(aes(color = Step),
            linewidth = 1.25)+
  scale_color_brewer("Iteration", palette = "Dark2")+
  geom_ribbon(aes(ymin = lb,
                  ymax = ub,
                  fill = Step),
              alpha =0.5)+
  scale_fill_brewer("Iteration", palette = "Dark2")+
  geom_point(size = 1.5,
             data = dat_sb2%>%
                    filter(Training == "Train"),
             color = "blue")+
   geom_point(size = 1.5,
             data = dat_sb2%>%
                    filter(Training == "Test"),
             color = "red")+
  labs(y = "Value",
       title = "Mechanical Properties of MoNb Binaries Under Tension")+
  facet_wrap(~property,labeller = label_parsed,
             scales = "free", nrow = 1)+
  theme(legend.position = "bottom")
```

![](Predictions_19-05---Copy_files/figure-html/pred_prop_monb-1.png)<!-- -->


``` r
dat_pm2%>%
  filter(property %in% c("epsilon[f]", "U[T]"),
         Step == 0)%>%
  ggplot(aes(x= Nb,
             y = val))+
  geom_line(aes(color = Step),
            linewidth = 1.25)+
  scale_color_brewer("Iteration", palette = "Dark2")+
  geom_ribbon(aes(ymin = lb,
                  ymax = ub,
                  fill = Step),
              alpha =0.5)+
  scale_fill_brewer("Iteration", palette = "Dark2")+
  geom_point(size = 1.5,
             data = dat_sb2%>%
                    filter(property %in% c("epsilon[f]", "U[T]")),
             color = "blue")+
  labs(y = "Value",
       title = "d) Properties Prediction of MoNb Alloys")+
  facet_wrap(~property,labeller = label_parsed,
             scales = "free", nrow = 1)+
  theme(legend.position = "none")
```

![](Predictions_19-05---Copy_files/figure-html/pred_prop_monb_small-1.png)<!-- -->


### VW Binary


``` r
in_pr <- c(0.2,0.4,0.6,0.8)
alloys_in <- c("V","W", "VW")

x <-seq(0.01,0.99, by = 0.01)
X_vw <- tibble(Mo = 0,
                   Nb = 0,
                   Ta = 0,
                   V = x,
                   W = 1-x)
X_vw_pred <- as.matrix(X_vw)

X_vw <- X_vw%>%
  unite("alloy_id",
        Mo:W,
        remove = FALSE)
```

### VW


``` r
in_pr <- c(0.2,0.4,0.6,0.8)
alloys_in <- c("V","W", "VW")

dat_sb <- dat_ss_c%>%
  mutate(Training = ifelse(V %in% in_pr,"Test","Train"))%>%
  filter(alloy %in% alloys_in)

x <-seq(0.01,0.99, by = 0.01)
X_vw <- tibble(Mo = 0,
                   Nb = 0,
                   Ta = 0,
                   V = x,
                   W = 1-x)

X_pred <-as.matrix(X_vw)

X_pred_irl <- ilr_bounded_minmax_additive(X_pred,
                                             epsilon = eps_ilr,
                                             z_min = xx$z_min,
                                             z_max = xx$z_max)$u

BX <- create_poly_basis(X_pred,
                                   n_out = n_prop,
                                   prop_names = properties,
                                   sparse_mat = T,
                               degree = 1,
                               omit_last_var = T,
                                   var_names = metals)

X_design <- X_vw%>%
  unite("alloy_id",
        Mo:W,
        remove = FALSE)

steps_in <- c(1,6)
n_files <- length(steps_in)
```


``` r
dat_mu <- dat_lb <- dat_ub <- dat_var <- dat_nug <- dat_pm <- vector(mode = "list", length = n_files)
imspe_v <- rep(0,n_files)
for(i in steps_in){
  fn <- paste0("train_bin_rl_prior_gam_pol_m52_2026_05_19_RLV2_iter_",i-1,".RData")
  fp <- here::here("Results","Fixed Data",fn)
  
  results <- import_predict_mod_gp_sampling_fun(newX = X_pred_irl,
                                                 Xdesign = X_design,
                                                 newX_pred = BX,
                                           file_path = fp,
                                           n_funs = length(properties) - 1,
                                           step_id = i-1,
                                           mod_id = 2,
                                           rseed = 1323,
                                           n_samples = 1e4,
                                           prop_names = properties,
                                           pred_type = pred_type)
  
 dat_mu[[i]] <- results$gp_mu
 dat_var[[i]] <- results$gp_var
 dat_nug[[i]] <- results$gp_nug
 dat_pm[[i]] <- results$mp_mean
 dat_lb[[i]] <- results$mp_lb
 dat_ub[[i]] <- results$mp_ub
 imspe_v[i] <- results$imspe
}


dat_mu <- bind_rows(dat_mu)
dat_var <- bind_rows(dat_var)
dat_nug <- bind_rows(dat_nug)
dat_pm <- bind_rows(dat_pm)
dat_lb <- bind_rows(dat_lb)
dat_ub <- bind_rows(dat_ub)
```


``` r
dat_lb2 <- dat_lb%>%
  pivot_longer(all_of(props),
               values_to = "val",
               names_to = "property")%>%
  mutate(property = factor(property, levels = props,
                           labels = prop_ltx))

dat_ub2 <- dat_ub%>%
  pivot_longer(all_of(props),
               values_to = "val",
               names_to = "property")%>%
  mutate(property = factor(property, levels = props,
                           labels = prop_ltx))

dat_pm2 <- dat_pm%>%
  pivot_longer(all_of(props),
               values_to = "val",
               names_to = "property")%>%
  mutate(Step = factor(Step),
         ub = dat_ub2$val,
         lb = dat_lb2$val)%>%
  mutate(property = factor(property, levels = props,
                           labels = prop_ltx))

dat_sb2 <- dat_sb%>%
  mutate(ultimate_stress = exp(log_max_stress),
         ultimate_strain = exp(log_max_strain))%>%
  pivot_longer(all_of(props),
               values_to = "val",
               names_to = "property")%>%
  mutate(property = factor(property, levels = props,
                           labels = prop_ltx))

dat_pm2%>%
  ggplot(aes(x= V,
             y = val))+
  geom_line(aes(color = Step),
            linewidth = 1.25)+
  scale_color_brewer("Iteration", palette = "Dark2")+
  geom_ribbon(aes(ymin = lb,
                  ymax = ub,
                  fill = Step),
              alpha =0.25)+
  scale_fill_brewer("Iteration", palette = "Dark2")+
  geom_point(size = 1.5,
             data = dat_sb2%>%
                    filter(Training == "Train"),
             color = "blue")+
   geom_point(size = 1.5,
             data = dat_sb2%>%
                    filter(Training == "Test"),
             color = "red")+
  labs(y = "Value",
       title = "Mechanical Properties of VW Binaries Under Compression")+
  facet_wrap(~property,labeller = label_parsed,
             scales = "free", nrow = 1)+
  theme(legend.position = "bottom")
```

![](Predictions_19-05---Copy_files/figure-html/pred_prop_vw-1.png)<!-- -->



``` r
dat_pm2%>%
  filter(property %in% c("epsilon[f]", "U[T]"),
         Step == 0)%>%
  ggplot(aes(x= V,
             y = val))+
  geom_line(linewidth = 1.25,
            colour = "blue")+
  geom_ribbon(aes(ymin = lb,
                  ymax = ub),
              fill = "steelblue",
              alpha =0.5)+
  geom_point(size = 1.5,
             data = dat_sb2%>%
                    filter(property %in% c("epsilon[f]", "U[T]")),
             color = "red")+
  labs(y = "Value",
       title = "d) Properties Prediction of VW Alloys")+
  facet_wrap(~property,labeller = label_parsed,
             scales = "free", nrow = 1)+
  theme(legend.position = "none",
        text=element_text(size = 12),
                  legend.box.margin = margin(-10,-10,-10,-10),
                  legend.text = element_text(size=10))
```

![](Predictions_19-05---Copy_files/figure-html/pred_prop_vw_small-1.png)<!-- -->


### Beta


``` r
X_vw <- X_vw%>%
  mutate(id = row_number())

var_ms_long <- dat_var%>%
            pivot_longer(all_of(properties),
                         values_to = "var",
                         names_to = "property")

pred_ms_c <- dat_mu%>%
            pivot_longer(all_of(properties),
                         values_to = "pm",
                         names_to = "property")%>%
              mutate(var_ms = var_ms_long$var)%>%
  inner_join(X_vw%>%
               select(V,W,id),
             by = c("id"),
             relationship = "many-to-one")%>%
        mutate(property = factor(property,
                                 levels = properties,
                                 labels = coefs_ltx))
dat_sb_3_c <- dat_ss_all%>%
  filter(alloy %in% alloys_in,
         type == "c")%>%
  mutate(Training = ifelse(round(V,2) %in% in_pr,"Test","Train"))%>%
            pivot_longer(all_of(properties),
                         values_to = "pm",
                         names_to = "property")%>%
        mutate(property = factor(property,
                                 levels = properties,
                                 labels = coefs_ltx))

dat_sb_3_t <- dat_ss_all%>%
  filter(alloy %in% alloys_in,
         type == "t")%>%
  mutate(Training = ifelse(round(V,2) %in% in_pr,"Test","Train"))%>%
            pivot_longer(all_of(properties),
                         values_to = "pm",
                         names_to = "property")%>%
        mutate(property = factor(property,
                                 levels = properties,
                                 labels = coefs_ltx))
```



``` r
dat_mu <- dat_mu_c <- vector(mode = "list", length = 2)
for(i in steps_in){
  fn <- paste0("train_bin_rl_prior_gam_pol_m52_2026_05_19_RLV2_iter_",i-1,".RData")
  fp <- here::here("Results","Fixed Data",fn)
   results <- import_predict_mod_gp_sampling_fun(newX = X_pred_irl,
                                                 Xdesign = X_design,
                                                 newX_pred = BX,
                                           file_path = fp,
                                           n_funs = length(properties) - 1,
                                           step_id = i-1,
                                           mod_id = 2,
                                           rseed = 1323,
                                           n_samples = 1e4,
                                           prop_names = properties,
                                           pred_type = pred_type)
  
  dat_mu_c[[i]] <- results$gp_mu
 
  results <- import_predict_mod_gp_sampling_fun(newX = X_pred_irl,
                                                 Xdesign = X_design,
                                                 newX_pred = BX,
                                           file_path = fp,
                                           n_funs = length(properties) - 1,
                                           step_id = i-1,
                                           mod_id = 1,
                                           rseed = 1323,
                                           n_samples = 1e4,
                                           prop_names = properties)
  
 dat_mu[[i]] <- results$gp_mu
 

}


dat_mu <- bind_rows(dat_mu)

pred_ms_t <- dat_mu%>%
            pivot_longer(all_of(properties),
                         values_to = "pm",
                         names_to = "property")%>%
              mutate(var_ms = var_ms_long$var)%>%
  inner_join(X_vw%>%
               select(V,W,id),
             by = c("id"),
             relationship = "many-to-one")%>%
        mutate(property = factor(property,
                                 levels = properties,
                                 labels = coefs_ltx))

pred_ms_c <- bind_rows(dat_mu_c)%>%
            pivot_longer(all_of(properties),
                         values_to = "pm",
                         names_to = "property")%>%
              mutate(var_ms = var_ms_long$var)%>%
  inner_join(X_vw%>%
               select(V,W,id),
             by = c("id"),
             relationship = "many-to-one")%>%
        mutate(property = factor(property,
                                 levels = properties,
                                 labels = coefs_ltx))
```





``` r
qz <- qnorm(0.05/2)

pc <- pred_ms_c%>%
  mutate(lb = pm - qz*sqrt(var_ms),
        ub = pm + qz*sqrt(var_ms))%>%
  mutate(Step = factor(Step))%>%
  ggplot(aes(x= V,
             y = pm))+
  geom_line(aes(color = Step),
            linewidth = 1.25)+
  scale_color_brewer("Iteration", palette = "Dark2")+
  geom_ribbon(aes(ymin = lb,
                  ymax = ub,
                  fill = Step),
              alpha =0.25)+
  scale_fill_brewer("Iteration", palette = "Dark2")+
  geom_point(size = 1.5,
             data = dat_sb_3_c%>%
                    filter(Training == "Train"),
             color = "blue")+
   geom_point(size = 1.5,
             data = dat_sb_3_c%>%
                    filter(Training == "Test"),
             color = "red")+
  facet_wrap(~property , scales = "free",
             nrow = 1,
             labeller = label_parsed)+
  labs(y = "Value",
       title = "Compression")+
  theme(legend.position = "right",
        text=element_text(size = 12),
                  legend.box.margin = margin(-10,-10,-10,-10),
                  legend.text = element_text(size=10))

pt <- pred_ms_t%>%
  mutate(lb = pm - qz*sqrt(var_ms),
        ub = pm + qz*sqrt(var_ms))%>%
  mutate(Step = factor(Step))%>%
  ggplot(aes(x= V,
             y = pm))+
  geom_line(aes(color = Step),
            linewidth = 1.25)+
  scale_color_brewer("Iteration", palette = "Dark2")+
  geom_ribbon(aes(ymin = lb,
                  ymax = ub,
                  fill = Step),
              alpha =0.25)+
  scale_fill_brewer("Iteration", palette = "Dark2")+
  geom_point(size = 1.5,
             data = dat_sb_3_t%>%
                    filter(Training == "Train"),
             color = "blue")+
   geom_point(size = 1.5,
             data = dat_sb_3_t%>%
                    filter(Training == "Test"),
             color = "red")+
  facet_wrap(~property , scales = "free",
             nrow = 1,
             labeller = label_parsed)+
  labs(y = "Value",
       title = "Tension")+
  theme(legend.position = "right",
        text=element_text(size = 12),
                  legend.box.margin = margin(-10,-10,-10,-10),
                  legend.text = element_text(size=10))

gridExtra::grid.arrange(pc,pt,ncol = 1)
```

![](Predictions_19-05---Copy_files/figure-html/pred_beta_v2-1.png)<!-- -->


``` r
qz <- qnorm(0.05/2)
props_in <- c("ln(epsilon[f])", "beta[1]")
pc <- pred_ms_c%>%
  filter(Step == 0,
         property %in% props_in)%>%
  mutate(lb = pm - qz*sqrt(var_ms),
        ub = pm + qz*sqrt(var_ms))%>%
  ggplot(aes(x= V,
             y = pm))+
  geom_line(colour = "blue",
            linewidth = 1.25)+
  geom_ribbon(aes(ymin = lb,
                  ymax = ub),
              fill = "steelblue",
              alpha =0.25)+
   geom_point(size = 1.5,
             data = dat_sb_3_c%>%
                    filter(property %in% props_in),
             color = "red")+
  facet_wrap(~property , scales = "free",
             nrow = 1,
             labeller = label_parsed)+
  labs(y = "Value",
       title = "c) Basis Coefficients")+
  theme(legend.position = "right",
        text=element_text(size = 12),
                  legend.box.margin = margin(-10,-10,-10,-10),
                  legend.text = element_text(size=10))
pc
```

![](Predictions_19-05---Copy_files/figure-html/pred_beta_small-1.png)<!-- -->



## Predictions for SS curves

### Plot SS


``` r
load(file = here::here("ss_data_all.RData"))

ss_dat_all <- ss_dat_all%>%
  mutate(test_type = ifelse(test_type=="Tension", 
                            "Tension",
                            "Compression"))%>%
  unite("alloy_id",
        Mo:W,
        remove = FALSE)
```



``` r
alloys_in <- c("V","W", "VW")
dat_sb <- ss_dat_all%>%
  mutate(W = round(W,1))%>%
  filter(alloy %in% alloys_in,
         W %in% c(0.2,0.5,0.8),
         test_type =="Tension")

ids_in <- X_vw %>%
          mutate(W = round(W,2))%>%
            filter(W %in% c(0.2,0.5,0.8))

gp_mu <- results$gp_mu

gp_mu <- gp_mu%>%
            filter(id %in% ids_in$id)

var_gp <- results$gp_var
var_gp <- var_gp%>%
  filter(id %in% ids_in$id)
```


``` r
test <- predict_ss_curve_sampling_fun(mu_gp =
                                        as.matrix(gp_mu[,1:5]),
                             var_gp = as.matrix(var_gp[,1:5]),
                             n_log_s = 1e3,
                             n_cond_s = 1e3,
                             max_strain = 0.08,
                             pred_type = "mean")

gc()
```

```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  3247051 173.5    6060757  323.7   6060757  323.7
## Vcells 19702618 150.4  197653732 1508.0 247051233 1884.9
```




``` r
dat_ss_ci <- tibble(strain = rep(test$strain, ncol(test$mean)),
                    V = rep(ids_in$V, each = nrow(test$mean)),
                    W = rep(ids_in$W, each = nrow(test$mean)),
                    pred = as.vector(test$mean),
                    ub = as.vector(test$ub),
                    lb = as.vector(test$lb),
                    )

dat_h <- tibble(W = ids_in$W,
                med = exp(as.vector(gp_mu$log_max_strain) + 
                            qnorm(0.05) * sqrt(as.vector(var_gp$log_max_strain))))
```





``` r
dat_ss_ci%>%
  ggplot(aes(x = strain))+
  geom_ribbon(aes(ymax = ub,
                  ymin = lb),
              fill="blue",
              alpha = 0.3)+
  geom_line(aes(y = pred),
            color = "blue",
            linewidth = 1.25)+
  # geom_vline(aes(xintercept = med),
  #            data = dat_h,
  #             color = "red",
  #            linewidth = 1.25,
  #           linetype = "dashed")+
  geom_line(aes(x=strain,
                y = stress,
                group = rep),
            data = dat_sb%>%
              filter(strain < 0.08),
            color = "black")+
  facet_wrap(~V,
             scales = "free")+
  labs(y = TeX("Stress $\\sigma$ (GPa)"),
       x = TeX("Strain $\\epsilon$"),
       title = "Predicted SS Curve for VW Binary Alloys under Tension")+
  theme_bw()
```

```
## Warning: Removed 37 rows containing missing values or values outside the scale range
## (`geom_ribbon()`).
```

![](Predictions_19-05---Copy_files/figure-html/ss_ci-1.png)<!-- -->



