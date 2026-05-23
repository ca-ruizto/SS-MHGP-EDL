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


