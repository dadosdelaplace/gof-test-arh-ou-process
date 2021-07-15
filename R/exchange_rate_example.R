library(fda.usc)
library(ggplot2)
library(reshape2)

### Read data --------------------------------------------------------------------------------------

datos <- read.csv2("EURGBP_2019_5min.csv", header = TRUE, colClasses = c(rep("character",2), rep("numeric",1)), dec = ".")
# datos <- read.csv2("EURUSD_2019_5min.csv", header = TRUE, colClasses = c(rep("character",2), rep("numeric",1)), dec = ".")
# datos <- read.csv2("GBPUSD_2019_5min.csv", header = TRUE, colClasses = c(rep("character",2), rep("numeric",1)), dec = ".")

t           <- seq(0, 1, by = 1/287)      # 288 "5 minutes" per day
sde_proc    <- datos[, 3]
fstoch_proc <- fdata(t(matrix(sde_proc, nrow = length(t))), argvals = t)


### Plot data --------------------------------------------------------------------------------------

# plot SDE
sde_datos <- data.frame(sde_proc, len1 = 1:length(sde_proc)/288)
ggplot(sde_datos, aes(x = len1, y = sde_proc)) +
  geom_line() + 
  labs(x = "t", y = "exchange rate") 

# plot FD
data      <- as.data.frame(t(fstoch_proc$data)) 
data$id   <- t
plot_data <- melt(data, id.var = "id")
ind       <- sort(rep(1:255, nrow(data)))
ggplot(plot_data, aes(x = id, y = value, group = variable, colour = ind)) +
  geom_line()  + 
  labs(x = "t", y = "exchange rate", color = "day") 


### Test data --------------------------------------------------------------------------------------

z           <- 1
hyp_simp    <- FALSE
hyp_comp    <- TRUE
est_method  <- "fpcr_l1s"
thre_p      <- 0.995
thre_q      <- 0.995
lambda      <- NULL
boot_scores <- TRUE
cv_1se      <- TRUE
B           <- 1000

t           <- seq(0, 1, by = 1/287)      # 288 "5 minutes" per day
sde_proc    <- datos[, 3]
sde_proc    <- sde_proc - mean(sde_proc)  # centered process

## SDE as a set of functional trajectories (assumed to be centered)
fstoch_proc <- fdata(t(matrix(sde_proc, nrow = length(t))), argvals = t)
X           <- list("fstoch_proc" = fstoch_proc)


## STAGE 1: Is X_n an ARH(1) process (that is, X_n vs X_{n-1} FLMFR)?
testing_ARHp <- ARHp_test(X, z = z, B = B, hyp_simp = hyp_simp,
                          hyp_comp = hyp_comp, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, lambda = lambda,
                          boot_scores = boot_scores, plot_dens = FALSE,
                          plot_proc = FALSE, save_fit_flm = FALSE,
                          save_boot_stats = FALSE, cv_1se = cv_1se)

testing_ARHp$testing_comp$p.value


## STAGE 2: F-test OU vs ARH(1) (unrestrictedd)
testing_F_OU <- F_stat_OU(testing_ARHp$X_flmfr, X, est_method = est_method,
                          thre_p = thre_p, thre_q = thre_q, lambda = lambda,
                          cv_1se = FALSE)
F_stat       <- testing_F_OU$F_stat
F_stat_ast   <- theta_est_ast <- rep(0, B)

# Parameter estimates
beta_est     <- testing_F_OU$pred_OU$theta_est
sigma_est    <- abs(optimize(mce, r = as.vector(sde_proc), Delta = 1/(length(t) - 1), interval = c(0,2))$minimum)
par.sde_ast  <- list("alpha" = 0, "beta" = beta_est, "sigma" = sigma_est) # centered process, alpha = mu = 0

for (b in 1:B) {
  # Simulating bootstrap replicates of X
  X_ast <- r_stoch_proc(dim(fstoch_proc$data)[1], t = t, par.sde = par.sde_ast, mu = 0, X0 = 0,
                        warm_up = -1, type = "CKLS", model = "OU", verbose = FALSE, plot = FALSE)
  
  # Converting X_ast into a FLMFR X_n vs X_{n-1}
  X_flmfr_ast <- ARH_to_FLMFR(X_ast[["fstoch_proc"]], z)
  
  # Computing bootstrap F-statistics
  testing_F_OU_ast <- F_stat_OU(X_flmfr_ast, X_ast, est_method = est_method,
                                thre_p = thre_p, thre_q = thre_q,
                                lambda = lambda, cv_1se = FALSE)
  theta_est_ast[b] <- testing_F_OU_ast$pred_OU$theta_est
  F_stat_ast[b]    <- testing_F_OU_ast$F_stat
}

# Approximation of the p-value by MC
p_value_F_test <- mean(F_stat < F_stat_ast)

## p-values
c(testing_ARHp$testing_comp$p.value, p_value_F_test)

