# ###############################################################
# Description: companion software to replicate the real
#              data application of the article
#              «A goodness-of-fit test for functional time series
#              with applications to diffusion processes»
# Authors: A. López-Pérez and J. Álvarez-Liébana
# ###############################################################

# ###########################
# PACKAGES
# ###########################
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
repos <- "http://cran.us.r-project.org"
if(!require(ggplot2)) install.packages("ggplot2", repos = repos)
if(!require(fda.usc)) install.packages("fda.usc", repos = repos)
if(!require(tidyverse)) install.packages("tidyverse", repos = repos)
if(!require(lubridate)) install.packages("lubridate", repos = repos)
if(!require(glue)) install.packages("glue", repos = repos)
if(!require(reshape2)) install.packages("reshape2", repos = repos)

# Companion software
source("./OU_estimation_test.R")

# ###########################
# DATA
# ###########################

# URL of data freely available at GitHub repository
# https://github.com/dadosdelaplace/gof-test-arh-ou-process
repo <- paste0("https://raw.githubusercontent.com/dadosdelaplace/",
               "gof-test-arh-ou-process/main/data")
url <- c(glue("{repo}/EURGBP_2019_5min.csv"),
         glue("{repo}/EURUSD_2019_5min.csv"),
         glue("{repo}/GBPUSD_2019_5min.csv"))

# Data was extracted from histdata.com/download-free-forex-data/
EURGBP <- read_delim(file = url[1], delim = ";")
EURUSD <- read_delim(file = url[2], delim = ";")
GBPUSD <- read_delim(file = url[3], delim = ";")

# Preparing dates
dates <- as.Date(as.character(EURGBP$date), format = "%Y%m%d")
hours <- floor(as.numeric(EURGBP$hour) / 1e4)
minutes <- floor((as.numeric(EURGBP$hour) - hours * 1e4) / 1e2)
hourtime <- hms::as_hms(hours * 60 * 60 + minutes*60)

# Build data
data <- data.frame("dates" = dates, "hourtime" = hourtime,
                   "EURGBP" = EURGBP[, 3], "EURUSD" = EURUSD[, 3],
                   "GBPUSD" = GBPUSD[, 3])
names(data) <- c("dates", "hourtime", "EURGBP", "EURUSD", "GBPUSD")

# Settings of working space: 288 daily curves recorded each 5 minutes
t <- seq(0, 1, by = 1/(60*24/5 - 1))
fda_data <- list()
fda_data$EURGBP <-
  fdata(t(matrix(data$EURGBP, nrow = length(t))), argvals = t)
fda_data$EURUSD <-
  fdata(t(matrix(data$EURUSD, nrow = length(t))), argvals = t)
fda_data$GBPUSD <-
  fdata(t(matrix(data$GBPUSD, nrow = length(t))), argvals = t)


# ###########################
# PLOTS
# ###########################


# # plot SDE
# sde_datos <- data.frame(sde_proc, len1 = 1:length(sde_proc)/288)
# ggplot(sde_datos, aes(x = len1, y = sde_proc)) +
#   geom_line() + 
#   labs(x = "t", y = "exchange rate") 
# 
# # plot FD
# data      <- as.data.frame(t(fstoch_proc$data)) 
# data$id   <- t
# plot_data <- melt(data, id.var = "id")
# ind       <- sort(rep(1:255, nrow(data)))
# ggplot(plot_data, aes(x = id, y = value, group = variable, colour = ind)) +
#   geom_line()  + 
#   labs(x = "t", y = "exchange rate", color = "day") 
# 

# ###########################
# TESTING EURGBP DATA
# ###########################

# Parameters
set.seed(123)
z <- 1 # order to test
hyp_simp <- FALSE # we test composite hypothesis
hyp_comp <- TRUE
est_method  <- "fpcr_l1s" # FLMFR estimation method
thre_p <- thre_q <- 0.995 # Initial (p, q) for capturing 99.5% of variance
cv_1se <- TRUE # optimal lambda be the lambda.1se, as returned by cv.glmnet?
B <- 1000 # Bootstrap replicates

# STAGE 1: Is X_n an ARH(1) process (that is, X_n vs X_{n-1} FLMFR)?
X <- list("fstoch_proc" = fda_data$EURGBP)
testing_ARHz <- ARHz_test(X, z = z, B = B, hyp_simp = hyp_simp,
                          hyp_comp = hyp_comp, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = cv_1se)
testing_ARHz$testing_comp$p.value

# STAGE 2: F-test OU vs ARH(1) (unrestrictedd)
testing_F_OU <- F_stat_OU(testing_ARHz$X_flmfr, X, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = FALSE)
F_stat <- testing_F_OU$F_stat
F_stat_ast <- theta_est_ast <- rep(0, B) 

# Parameter estimates
beta_est <- testing_F_OU$pred_OU$theta_est
sigma_est <- abs(optimize(mce, r = as.vector(data$EURGBP),
                          Delta = 1/(length(t) - 1),
                          interval = c(0, 2))$minimum)
par.sde_ast <- list("alpha" = 0, "beta" = beta_est,
                    "sigma" = sigma_est) # centered process, alpha = mu = 0

# Parametric boostrap
pb <- txtProgressBar(max = B, style = 3)
for (b in 1:B) {
  
  # Simulating bootstrap replicates of X
  X_ast <- r_stoch_proc(dim(fda_data$EURGBP$data)[1], t = t,
                        par.sde = par.sde_ast, mu = 0, X0 = 0,
                        warm_up = -1, type = "CKLS", model = "OU",
                        verbose = FALSE, plot = FALSE)
  
  # Converting X_ast into a FLMFR X_n vs X_{n-1}
  X_flmfr_ast <- ARH_to_FLMFR(X_ast[["fstoch_proc"]], z)
  
  # Computing bootstrap F-statistics
  testing_F_OU_ast <- F_stat_OU(X_flmfr_ast, X_ast, est_method = est_method,
                                thre_p = thre_p, thre_q = thre_q,
                                cv_1se = FALSE, verbose = FALSE)
  theta_est_ast[b] <- testing_F_OU_ast$pred_OU$theta_est
  F_stat_ast[b]    <- testing_F_OU_ast$F_stat
  
  setTxtProgressBar(pb, b)
}

# Approximation of the p-value by MC
p_value_F_test <- mean(F_stat < F_stat_ast)

# p-values
c(testing_ARHz$testing_comp$p.value, p_value_F_test)


# ###########################
# TESTING EURUSD DATA
# ###########################

# Parameters
set.seed(123)
z <- 1 # order to test
hyp_simp <- FALSE # we test composite hypothesis
hyp_comp <- TRUE
est_method  <- "fpcr_l1s" # FLMFR estimation method
thre_p <- thre_q <- 0.995 # Initial (p, q) for capturing 99.5% of variance
cv_1se <- TRUE # optimal lambda be the lambda.1se, as returned by cv.glmnet?
B <- 1000 # Bootstrap replicates

# STAGE 1: Is X_n an ARH(1) process (that is, X_n vs X_{n-1} FLMFR)?
X <- list("fstoch_proc" = fda_data$EURUSD)
testing_ARHz <- ARHz_test(X, z = z, B = B, hyp_simp = hyp_simp,
                          hyp_comp = hyp_comp, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = cv_1se)
testing_ARHz$testing_comp$p.value

# STAGE 2: F-test OU vs ARH(1) (unrestrictedd)
testing_F_OU <- F_stat_OU(testing_ARHz$X_flmfr, X, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = FALSE)
F_stat <- testing_F_OU$F_stat
F_stat_ast <- theta_est_ast <- rep(0, B) 

# Parameter estimates
beta_est <- testing_F_OU$pred_OU$theta_est
sigma_est <- abs(optimize(mce, r = as.vector(data$EURUSD),
                          Delta = 1/(length(t) - 1),
                          interval = c(0, 2))$minimum)
par.sde_ast <- list("alpha" = 0, "beta" = beta_est,
                    "sigma" = sigma_est) # centered process, alpha = mu = 0

# Parametric boostrap
pb <- txtProgressBar(max = B, style = 3)
for (b in 1:B) {
  
  # Simulating bootstrap replicates of X
  X_ast <- r_stoch_proc(dim(fda_data$EURUSD$data)[1], t = t,
                        par.sde = par.sde_ast, mu = 0, X0 = 0,
                        warm_up = -1, type = "CKLS", model = "OU",
                        verbose = FALSE, plot = FALSE)
  
  # Converting X_ast into a FLMFR X_n vs X_{n-1}
  X_flmfr_ast <- ARH_to_FLMFR(X_ast[["fstoch_proc"]], z)
  
  # Computing bootstrap F-statistics
  testing_F_OU_ast <- F_stat_OU(X_flmfr_ast, X_ast, est_method = est_method,
                                thre_p = thre_p, thre_q = thre_q,
                                cv_1se = FALSE, verbose = FALSE)
  theta_est_ast[b] <- testing_F_OU_ast$pred_OU$theta_est
  F_stat_ast[b]    <- testing_F_OU_ast$F_stat
  
  setTxtProgressBar(pb, b)
}

# Approximation of the p-value by MC
p_value_F_test <- mean(F_stat < F_stat_ast)

# p-values
c(testing_ARHz$testing_comp$p.value, p_value_F_test)


# ###########################
# TESTING GBPUSD DATA
# ###########################

# Parameters
set.seed(123)
z <- 1 # order to test
hyp_simp <- FALSE # we test composite hypothesis
hyp_comp <- TRUE
est_method  <- "fpcr_l1s" # FLMFR estimation method
thre_p <- thre_q <- 0.995 # Initial (p, q) for capturing 99.5% of variance
cv_1se <- TRUE # optimal lambda be the lambda.1se, as returned by cv.glmnet?
B <- 1000 # Bootstrap replicates

# STAGE 1: Is X_n an ARH(1) process (that is, X_n vs X_{n-1} FLMFR)?
X <- list("fstoch_proc" = fda_data$GBPUSD)
testing_ARHz <- ARHz_test(X, z = z, B = B, hyp_simp = hyp_simp,
                          hyp_comp = hyp_comp, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = cv_1se)
testing_ARHz$testing_comp$p.value

# STAGE 2: F-test OU vs ARH(1) (unrestrictedd)
testing_F_OU <- F_stat_OU(testing_ARHz$X_flmfr, X, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, cv_1se = FALSE)
F_stat <- testing_F_OU$F_stat
F_stat_ast <- theta_est_ast <- rep(0, B) 

# Parameter estimates
beta_est <- testing_F_OU$pred_OU$theta_est
sigma_est <- abs(optimize(mce, r = as.vector(data$GBPUSD),
                          Delta = 1/(length(t) - 1),
                          interval = c(0, 2))$minimum)
par.sde_ast <- list("alpha" = 0, "beta" = beta_est,
                    "sigma" = sigma_est) # centered process, alpha = mu = 0

# Parametric boostrap
pb <- txtProgressBar(max = B, style = 3)
for (b in 1:B) {
  
  # Simulating bootstrap replicates of X
  X_ast <- r_stoch_proc(dim(fda_data$GBPUSD$data)[1], t = t,
                        par.sde = par.sde_ast, mu = 0, X0 = 0,
                        warm_up = -1, type = "CKLS", model = "OU",
                        verbose = FALSE, plot = FALSE)
  
  # Converting X_ast into a FLMFR X_n vs X_{n-1}
  X_flmfr_ast <- ARH_to_FLMFR(X_ast[["fstoch_proc"]], z)
  
  # Computing bootstrap F-statistics
  testing_F_OU_ast <- F_stat_OU(X_flmfr_ast, X_ast, est_method = est_method,
                                thre_p = thre_p, thre_q = thre_q,
                                cv_1se = FALSE, verbose = FALSE)
  theta_est_ast[b] <- testing_F_OU_ast$pred_OU$theta_est
  F_stat_ast[b]    <- testing_F_OU_ast$F_stat
  
  setTxtProgressBar(pb, b)
}

# Approximation of the p-value by MC
p_value_F_test <- mean(F_stat < F_stat_ast)

# p-values
c(testing_ARHz$testing_comp$p.value, p_value_F_test)
