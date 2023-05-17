
# G: Treatment assignment indicator
# W: Treatment confounder
# Z: Observation times confounder
# Y: Time-varying outcome
# t: observation time
# phi: random effect
# epsilon: random error, function of time

# Data generation model: Y_i(t)  = (2-t) + beta1G_i + beta2[W_i - E(w_i|G_i)] + 
#                                     beta3[Z_i - E(Z_i|G_i)] + phi_i + epsilon_i(t)
#
# Marginal model: E(Y_i(t) | G_i) = E{(2-t) + beta1G_i + beta2[W_i - E(w_i|G_i)] + beta3(Z_i - E(Z_i|G_i)) | G_i]
#                                 = (2-t) + beta1G_i + betaE{2[W_i - W(W_i|G_i) | G_i] + beta3E{Z_i - E(Z_i|G_i)|G_i}
#                                 = (2-t) + beta1G_i 
#
# Observation intensity: lambda_i(t) = eta_i*exp(gamma1G_i + gamma2W_i + gamma3Z_i)
# Treatment assignment: P(G_i = 1) = expit(alpha0 + alpha1W_i + alpha2Z_i) 
#
#
# tau: maximum follow up time, independent censoring.
# N: number of simulations
#

setwd("C:\\Users\\grace\\Documents\\WaterlooPHD\\Research\\IIW_Simulation\\FIPTIW\\Simulation1")
source("FIPTIW_timeinvar_functions.R")


set.seed(500)
n = 100
beta1 = 0.5
beta2 = 1
beta3 = 1
gamma1 = 0.5
gamma2 = 0
gamma3 = 0.3
alpha0 = 0.5
alpha1 = 0.4
tau = 7
N = 1000
# 
# # 
# resultsMIPI <- simulateResultsMIPI(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau)
# resultsMIPI
  


# set.seed(500)
# n = 100
# beta1 = 0.5
# beta2 = 1
# beta3 = 1
# gamma1 = 0.5
# gamma2vec = c(0, 0.5)
# gamma3vec = c(0, 0.6)
# alpha0 = 0.5
# alpha1vec = c(0, 0.4)
# tau = 7
# N = 1000
# 
# resultsMIPI_all <- simulateALLMIPI(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, alpha0, alpha1vec, tau)
# saveRDS(resultsMIPI_all, "resultsMIPI1.rds")
resultsMIPI_all <- readRDS("resultsMIPI1.rds")
kable(resultsMIPI_all$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))

kable(resultsMIPI_all$resultsmat_stab, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))






