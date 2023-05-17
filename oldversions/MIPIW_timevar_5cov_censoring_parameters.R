# The purpose of this sensitivity analysis is to evaluate the 
# robustness of the non-informative censoring assumption.


# G: Treatment assignment indicator
# W_1(t): Treatment confounder
# W_2: Treatment confounder
# Z_1: Observation times confounder
# Z_2(t): Observation time variable (not a confounder)
# Y_i(t): Time-varying outcome
# t: observation time
# phi: random effect
# epsilon: random error, function of time

# Data generation model: Y_i(t)  = (2-t) + beta1G_i + beta2[W_i1(t) - E(w_i(t)|G_i)] + 
#                                     beta3[W_i2 - E(W_i2|G_i)] + beta4[Z_i1 - E(Z_i1|G_i)] + 
#                                     beta5[Z_i2(t) - E(Z_i2(t) | G_i)] + phi_i + epsilon_i(t)
#
# Marginal model: E(Y_i(t) | G_i)  = (2-t) + beta1G_i 
#
# Observation intensity: lambda_i(t) = eta_i*sqrt(t)/2*exp(gamma1G_i + gamma2W_i1(t) + gamma3W_i2 + gamma3Z_i1 + gamma4Z_i2(t))
# Treatment assignment: P(G_i = 1) = expit(alpha0 + alpha1W_i1 + alpha2W_i2) 
#
# Censoring hazard: h(t) = h_0*exp(kappa1D_i + kappa2W_i + kappa3Z_i)
#
# tau: maximum follow up time, independent censoring.
# N: number of simulations
#

library(knitr)
library(latex2exp)
library(kableExtra)

setwd("C:\\Users\\grace\\Documents\\WaterlooPHD\\Research\\IIW_Simulation\\MIPIW\\Simulation1")
source("MIPIW_timevar_5cov_censoring_functions.R")
# 
# Violating independence between all three processes
set.seed(999)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
beta4 = 1
beta5 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
gamma4 = 0.2
gamma5 = 0.3
alpha0 = -1
alpha1 = 3
alpha2 = 2
kappa1 = -0.5
kappa2 = 0.1
kappa3 = -0.2
tau = 7
N = 1000

# #
# 
# test1 <- simulateResultsMIPI_timevar_SENS(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                          kappa1, kappa2, kappa3,
#                                          alpha0, alpha1, tau)
# saveRDS(test1, "test1.rds")

test1 <- readRDS("test1.rds")


set.seed(333)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
kappa1vec = c(0,  -0.5)
kappa2vec = c(0,  -0.5)
kappa3vec = c(0,  -0.5)
alpha0 = -1
alpha1 = 3
tau = 7
N = 1000

# sensanalysis2 <- simulateALLMIPI_timevar_SENS(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                          kappa1vec, kappa2vec, kappa3vec,
#                                          alpha0, alpha1, tau)
# 
# saveRDS(sensanalysis2, "sensanalysis2.rds")
sensanalysis2 <- readRDS("sensanalysis2.rds")

kable(sensanalysis2$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))




set.seed(333)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
kappa1vec = c(0,  0.5)
kappa2vec = c(0,  0.5)
kappa3vec = c(0,  0.5)
alpha0 = -1
alpha1 = 3
tau = 7
N = 1000

sensanalysis3 <- simulateALLMIPI_timevar_SENS(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
                                         kappa1vec, kappa2vec, kappa3vec,
                                         alpha0, alpha1, tau)

saveRDS(sensanalysis3, "sensanalysis3.rds")
sensanalysis3 <- readRDS("sensanalysis3.rds")

kable(sensanalysis3$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))



set.seed(39203)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
kappa1vec = c(0,  -0.2)
kappa2vec = c(0,  -0.2)
kappa3vec = c(0,  -1)
alpha0 = -1
alpha1 = 3
tau = 7
N = 1000

# sensanalysis4 <- simulateALLMIPI_timevar_SENS(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                               kappa1vec, kappa2vec, kappa3vec,
#                                               alpha0, alpha1, tau)
# 
# saveRDS(sensanalysis4, "sensanalysis4.rds")
sensanalysis4 <- readRDS("sensanalysis4.rds")

kable(sensanalysis4$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))



set.seed(6578643)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
kappa1vec = c(0,  -0.2)
kappa2vec = c(0,  0.2)
kappa3vec = c(0, 1)
alpha0 = -1
alpha1 = 3
tau = 7
N = 1000
# 
# sensanalysis5 <- simulateALLMIPI_timevar_SENS(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                               kappa1vec, kappa2vec, kappa3vec,
#                                               alpha0, alpha1, tau)
# 
# saveRDS(sensanalysis5, "sensanalysis5.rds")
sensanalysis5 <- readRDS("sensanalysis5.rds")

kable(sensanalysis5$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))




set.seed(465123)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
kappa1vec = c(0,  -0.75)
kappa2vec = c(0,  0.80)
kappa3vec = c(0, 1.5)
alpha0 = -1
alpha1 = 3
tau = 7
N = 1000

# sensanalysis6 <- simulateALLMIPI_timevar_SENS(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                               kappa1vec, kappa2vec, kappa3vec,
#                                               alpha0, alpha1, tau)
# 
# saveRDS(sensanalysis6, "sensanalysis6.rds")
sensanalysis6 <- readRDS("sensanalysis6.rds")

kable(sensanalysis6$resultsmat, booktabs = T, digits = 3,format = "latex")%>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))

# ### comparison
# boxplot(sensanalysis6$fulldatalist[[1]], main = TeX(r'(Simulation Results: $\kappa_1 = -0.75$, $\kappa_2 = 0.80$, $\kappa_3 = 1.5$)', bold=TRUE),
#         ylim = c(-0.5, 1.6))
# abline(h = 0.5, lty = 3, col = "red")
# 
# boxplot(resultsMIPI_all_timevar2$fulldatalist[[12]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0.3$, $\gamma_3 = 0.6$, $\alpha_1 = 3$)', bold=TRUE),
#         ylim = c(-0.5, 1.6))
# abline(h = 0.5, col = "red", lty = 3)
# #####




set.seed(28724)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0
gamma3 = 0.6
kappa1vec = c(0,  -0.75)
kappa2vec = c(0,  0.80)
kappa3vec = c(0, 1.5)
alpha0 = -1
alpha1 = 3
tau = 7
N = 1000

sensanalysis7 <- simulateALLMIPI_timevar_SENS(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
                                              kappa1vec, kappa2vec, kappa3vec,
                                              alpha0, alpha1, tau)

saveRDS(sensanalysis7, "sensanalysis7.rds")
sensanalysis7 <- readRDS("sensanalysis7.rds")

kable(sensanalysis7$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))

### comparison

library(latex2exp)
boxplot(sensanalysis6$fulldatalist[[1]], main = TeX(r'(Simulation Results: $\kappa_1 = -0.75$, $\kappa_2 = 0.80$, $\kappa_3 = 1.5$)', bold=TRUE),
        ylim = c(-0.5, 1.6))
abline(h = 0.5, lty = 3, col = "red")

boxplot(resultsMIPI_all_timevar2$fulldatalist[[12]],
        main = TeX(r'(Simulation Results: $\gamma_2 = 0.3$, $\gamma_3 = 0.6$, $\alpha_1 = 3$)', bold=TRUE),
        ylim = c(-0.5, 1.6))
abline(h = 0.5, col = "red", lty = 3)
#####



