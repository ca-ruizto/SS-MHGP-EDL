library(Matrix)
library(collapse)
library(tidyverse)
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


n_dim <- 4
n_locs <- 50
n_props <- 7

mixture_dat <- AlgDesign::gen.mixture(21,5)

eps_ilr <-0.02

xx <- ilr_bounded_minmax_additive(as.matrix(mixture_dat),
                                  b = 0,
                                  epsilon = eps_ilr,
                                  method = "basic")

metals <- c("Mo","Nb","Ta","V","W")
dat_ss_c <- rio::import(here::here("Results",
                                   "ss_all_ort_coefs_c_uns.csv"))

dat_ss_fil <- dat_ss_c%>%
  unite("alloy_id",
        Mo:W,remove = FALSE)%>%
  ungroup()%>%
  group_by(alloy_id)%>%
  arrange(desc(rep))%>%
  slice_head(n=1)%>%
  mutate(elements= (Mo>0) + (Nb>0) + (Ta>0) + (V>0) + (W>0))

temp <- dat_ss_fil%>%
  ungroup()%>%
  filter(elements == 2,
         round(pmax(Mo,Nb,Ta,V,W),1) %in% c(0.1,0.3,0.5,0.7,0.9) 
  )%>%
  arrange(alloy_cnt)

X <- ilr_bounded_minmax_additive(as.matrix(temp[,metals]),
                                 epsilon = eps_ilr,
                                 method = "basic",
                                 z_min = xx$z_min,
                                 z_max = xx$z_max)$u


X <- matrix(runif(n_dim * n_locs), nrow = n_locs,
            ncol = n_dim)
th <- runif(n_dim, min = 0.5, max = 3)

W_hgp <- Wij(mu1 = as.matrix(X[,1], ncol = 1),
             theta = th[1], type = "Matern3_2")

W_int <- Wij_int(X = X[,1], 
                 theta = th[1], cov_fun = matern_5_2_fun )

all.equal(W_hgp, W_int)

D0 <- mv_dist_fun(X)

cor_fun_test <- function(theta, D_mat) general_cor_fun(theta = theta,
                                                       D_mat = D_mat,
                                                       kernel_fun = matern_5_2_fun)

C <- cor_fun_test(D_mat = D0, 
                theta = log(th))
Sigtest <- diag(c(1,runif(n_props-1)))
CS <- kronecker(C,Sigtest)
CS_lam <- CS + Diagonal(x = runif(n_locs * n_props) + 1e-10)
Ki <- chol2inv(chol(CS_lam))


pred_obj <- list(Ki = Ki,
                 theta = log(th),
                 X0 = X,
                 cor_fun_gp = cor_fun_test,
                 Sig = Sigtest,
                 tau = 2.1,
                 cov_type = "Matern5_2",
                 n_properties = n_props)

IMSPE_MV(pred_obj)

IMSPE_MV_NI(pred_obj)


```{r}
file_name <- paste0("train_bin_rl_prior_gam_",
                    "pol", "_", "m32",
                    "_2026_05_19_RLV2_iter_",0,".RData")

fp <- here::here("Results",
                 "Fixed Data",
                 file_name)

load(fp)
```

Test one case with negative imspe

```{r}
n_pol <- c("Tension-Hom","Compression-Hom", "Tension","Compression")

imspe_samp_dat <- vector(mode = "list")
dat_ls <- vector(mode = "list")
val_int <- val_ni <- numeric(4)
k <- 0
np_vals <- c(3,5,10)
for(p in n_pol){
  message(paste("Starting", p, "at", Sys.time()))
  k <- k + 1
  val <- numeric(0)
  val_int[k] <- IMSPE_MV(pol[[p]])
  val_ni[k] <- IMSPE_MV_NI(pol[[p]], max_eval = 10^4)
  
}

impse_analytical <- tibble(imspe = val_int,
                           Model = n_pol)
impse_ni <- tibble(imspe = val_ni,
                   Model = n_pol)

dat_ls <- bind_rows(dat_ls)
```

```{r fig.height=3, fig.width=6}
plot_bad <- imspe_samp_dat%>%
  ggplot()+
  geom_line(aes(x = n_samples,
                y = imspe))+
  geom_hline(aes(yintercept = imspe ),
             data = impse_analytical,
             color = "red")+
  geom_hline(aes(yintercept = imspe ),
             data = impse_ni,
             color = "blue")+
  facet_wrap(~Model,scales = "free")

plot_bad
```

```{r}
file_name <- paste0("train_bin_rl_prior_gam_",
                    "pol", "_", "m32",
                    "_2026_05_19_RLV2_iter_",1,".RData")

fp <- here::here("Results",
                 "Fixed Data",
                 file_name)

load(fp)
```




```{r}
np <- 3

a <- IMSPE_MV( pol[[np]])

b <- IMSPE_MV_NI(pol[[np]], max_eval = 1e2, return_object = T)
a
b

abs(a - b$integral) < b$error

err_v <- rep(0,4)
err_v[1] <- b$error

for(i in 2:4){
  d <- IMSPE_MV_NI(pol[[np]], max_eval = 10^(i + 1), return_object = T)
  err_v[i] <- d$error
}

plot(err_v)
#IMSPE_MV_NI(pol[[np]], max_eval = 3e4, return_object = F)


```
