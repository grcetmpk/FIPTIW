
#~~~~~~~~~~~~~~~~~~ SIMULATiONS I and II ~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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
source("FIPTIW_timevarW_functions.R")

# parallelization stuff
ncores <- detectCores() - 1
nclusters <- makeCluster(ncores)

require(knitr)
require(kableExtra)


#~~~~~~~~~~ FIPTIW simulation, assumptions satisfied (not shown in paper).

#### n = 100 ###

set.seed(500)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3)
gamma3vec = c(0, 0.6)
alpha0 = -1
alpha1vec = c(0, 1, 3)
tau = 7
N = 1000

# resultsFIPTIW_n100 <- simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
#                                                         alpha0, alpha1vec, tau, outputfulldatalist = FALSE)
# saveRDS(resultsFIPTIW_n100, "resultsFIPTIW_n100.rds")
resultsFIPTIW_n100 <- readRDS("resultsFIPTIW_n100.rds")

kable(resultsFIPTIW_n100$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "FIPTIW" = 3))




### n = 50 ###

set.seed(6234)
n = 50
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3)
gamma3vec = c(0, 0.6)
alpha0 = -1
alpha1vec = c(0, 1, 3)
tau = 7
N = 1000

# resultsFIPTIW_n50 <- simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
#                                                        alpha0, alpha1vec, tau, outputfulldatalist = FALSE)
#saveRDS(resultsFIPTIW_n50, "resultsFIPTIW_n50.rds")
resultsFIPTIW_n50 <- readRDS("resultsFIPTIW_n50.rds")

kable(resultsFIPTIW_n50$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "FIPTIW" = 3))



### n = 500 ###


set.seed(23424)
n = 500
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3)
gamma3vec = c(0, 0.6)
alpha0 = -1
alpha1vec = c(0, 1, 3)
tau = 7
N = 1000

# resultsFIPTIW_n500 <- simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
#                                                         alpha0, alpha1vec, tau, outputfulldatalist = FALSE)
# saveRDS(resultsFIPTIW_n500, "resultsFIPTIW_n500.rds")
resultsFIPTIW_n500 <- readRDS("resultsFIPTIW_n500.rds")

kable(resultsFIPTIW_n500$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "FIPTIW" = 3))





#### ~~~~~~~~~~~~ Variable selection simulation (Sim II) ~~~~~~~~~~~~~~~~~~~ #

set.seed(54641230)
n = 100
beta1 = 0.5
beta2vec = c(0, 2)
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3)
gamma3 = 0.6
alpha0 = 0
alpha1 = 0
tau = 7
N = 1000

# varsel_n100 <- simulateALLVarSel(N, n, beta1, beta2vec, beta3, gamma1, gamma2vec, gamma3, alpha0, alpha1, tau)
# saveRDS(varsel_n100, "varsel_n100.rds")
varsel_n100 <- readRDS("varsel_n100.rds")

kable(varsel_n100$resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
  add_header_above(c(" " = 3, "Variables used to estimate intensity" = 8))



set.seed(6412)
n = 50
beta1 = 0.5
beta2vec = c(0, 2)
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3)
gamma3 = 0.6
alpha0 = 0
alpha1 = 0
tau = 7
N = 1000

# varsel_n50 <- simulateALLVarSel(N, n, beta1, beta2vec, beta3, gamma1, gamma2vec, gamma3, alpha0, alpha1, tau)
#saveRDS(varsel_n50, "varsel_n50.rds")
varsel_n50 <- readRDS("varsel_n50.rds")
kable(varsel_n50$resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
  add_header_above(c(" " = 3, "Variables used to estimate intensity" = 8))



set.seed(5465)
n = 500
beta1 = 0.5
beta2vec = c(0, 2)
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3)
gamma3 = 0.6
alpha0 = 0
alpha1 = 0
tau = 7
N = 1000

varsel_timevar_n500 <- simulateALLVarSel(N, n, beta1, beta2vec, beta3, 
                                         gamma1, gamma2vec, gamma3, alpha0, alpha1, tau,
                                         inParallel = T, nclusters = nclusters)
saveRDS(varsel_n500, "varsel_timevar_n1000.rds")

kable(varsel_n500$resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
  add_header_above(c(" " = 3, "Variables used to estimate intensity" = 8))








#~~~~~~~~~~~~~~ Violations of Noninform Censoring Sim (Sim I)~~~~~~~~~#
# 
# set.seed(3402385)
# n = 100
# beta1 = 0.5
# beta2 = 2
# beta3 = 1
# gamma1 = 0.5
# gamma2 = 0.3
# gamma3 = 0.6
# alpha0 = -1
# alpha1 = 1
# tau = 7
# N = 1000
# censinform = F 
# eta1vec = c(0, 0.1, 0.4)
# eta2vec = c(0, 0.2, 0.5)
# eta3vec = c(0, 0.4, 0.6)
# 
# results_censoring_n100 <- simulateALLFIPTIW_CENS(N, n, beta1, beta2, beta3, gamma1,
#             gamma2vec, gamma3vec, alpha0, alpha1, tau, outputfulldatalist = F,
#             censinform = T, eta1, eta2, eta3)
# saveRDS(results_censoring_n100, "results_censoring_n100.rds")
# results_censoring_n100 <- readRDS("results_censoring_n100.rds")
# 
# kable(results_censoring_n100$resultsmat[, c(1,2,3, 4,6,7,9,10,12, 13, 15)], booktabs = T, 
#       digits = 3,format = "latex") %>%
#       add_header_above(c(" " = 3, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2))
# 
# 
# 
# 
# 
# 
# set.seed(6875)
# n = 50
# beta1 = 0.5
# beta2 = 2
# beta3 = 1
# gamma1 = 0.5
# gamma2 = 0.3
# gamma3 = 0.6
# alpha0 = -1
# alpha1 = 1
# eta1 = 0.1
# eta2 = 0.2
# eta3 = 0.4
# tau = 7
# N = 1000
# censinform = F
# eta1vec = c(0, 0.1, 0.4)
# eta2vec = c(0, 0.2, 0.5)
# eta3vec = c(0, 0.4, 0.6)
# 
# # results_censoring_n50 <- simulateALLFIPTIW_CENS(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, alpha0, alpha1, tau, outputfulldatalist = F, censinform = T, eta1, eta2, eta3)
# # saveRDS(results_censoring_n50, "results_censoring_n50.rds")
# 
# results_censoring_n50 <- readRDS("results_censoring_n50.rds")
# 
# kable(results_censoring_n50$resultsmat[, c(1,2,3, 4,6,7,9,10,12, 13, 15)], booktabs = T, digits = 3,format = "latex") %>%
#   add_header_above(c(" " = 3, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2))
# 
# 
# # n = 500
# 
# set.seed(984512)
# n = 500
# beta1 = 0.5
# beta2 = 2
# beta3 = 1
# gamma1 = 0.5
# gamma2 = 0.3
# gamma3 = 0.6
# alpha0 = -1
# alpha1 = 1
# eta1 = -0.1
# eta2 = -0.2
# eta3 = -0.4
# tau = 7
# N = 1000
# censinform = T
# eta1vec = c(0, 0.1, 0.4)
# eta2vec = c(0, 0.2, 0.5)
# eta3vec = c(0, 0.4, 0.6)
# # 
# 
# results_censoring_IPCW_n500 <- simulateALLFIPTICW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, alpha0, alpha1, tau, outputfulldatalist = F, censinform = T, eta1, eta2, eta3)
# # saveRDS(results_censoring_n500, "results_censoring_n500.rds")
# 
# results_censoring_n500 <- readRDS("results_censoring_n500.rds")
# 
# kable(results_censoring_n500$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
#   add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "FIPTIW" = 3))
# 


#~~~~~~~~~~~~~~ Violations of Noninform Censoring Sim with IPCW~~~~~~~~~#

#### n = 100

# set.seed(234324)
# n = 100
# beta1 = 0.5
# beta2 = 2
# beta3 = 1
# gamma1 = 0.5
# gamma2 = 0.3
# gamma3 = 0.6
# alpha0 = -1
# alpha1 = 1
# tau = 7
# N = 1000
# censinform = F 
# eta1vec = c(0, 0.1, 0.4)
# eta2vec = c(0, 0.2, 0.5)
# eta3vec = c(0, 0.4, 0.6)
# 
# results_censoring_n100_ipcw <- simulateALLFIPTICW_CENS(N, n, beta1, beta2, beta3, gamma1,
#                                                  gamma2vec, gamma3vec, alpha0, alpha1, tau, outputfulldatalist = F,
#                                                  censinform = T, eta1, eta2, eta3)
saveRDS(results_censoring_n100_ipcw, "results_censoring_n100_ipcw.rds")
results_censoring_n100_ipcw <- readRDS("results_censoring_n100_ipcw.rds")

kable(results_censoring_n100_ipcw$resultsmat[, c(1,2,3, 4,6,7,9,10,12, 13, 15, 16, 18)], booktabs = T, 
      digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2, "FIPTICW" = 2))



##### n = 50

set.seed(3028)
n = 50
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = -1
alpha1 = 1
tau = 7
N = 1000
censinform = F
eta1vec = c(0, 0.1, 0.4)
eta2vec = c(0, 0.2, 0.5)
eta3vec = c(0, 0.4, 0.6)

results_censoring_n50_ipcw <- simulateALLFIPTICW_CENS(N, n, beta1, beta2, beta3, gamma1,
                                                 gamma2vec, gamma3vec, alpha0, alpha1, tau, outputfulldatalist = F,
                                                 censinform = T, eta1, eta2, eta3)
saveRDS(results_censoring_n50_ipcw, "results_censoring_n50_ipcw.rds")
results_censoring_n50_ipcw <- readRDS("results_censoring_n50_ipcw.rds")

kable(results_censoring_n50_ipcw$resultsmat[, c(1,2,3, 4,6,7,9,10,12, 13, 15, 16, 18)], booktabs = T, 
      digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2, "FIPTICW" = 2))



##### n = 500

set.seed(234234)
n = 500
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = -1
alpha1 = 1
tau = 7
N = 1000
censinform = F
eta1vec = c(0, 0.1, 0.4)
eta2vec = c(0, 0.2, 0.5)
eta3vec = c(0, 0.4, 0.6)

results_censoring_n500_ipcw <- simulateALLFIPTICW_CENS(N, n, beta1, beta2, beta3, gamma1,
                                                      gamma2vec, gamma3vec, alpha0, alpha1, tau, outputfulldatalist = F,
                                                      censinform = T, eta1, eta2, eta3)
saveRDS(results_censoring_n500_ipcw, "results_censoring_n500_ipcw.rds")
results_censoring_n500_ipcw <- readRDS("results_censoring_n50 fdswjhmv c0_ipcw.rds")

kable(results_censoring_n50_ipcw$resultsmat[, c(1,2,3, 4,6,7,9,10,12, 13, 15, 16, 18)], booktabs = T, 
      digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2, "FIPTICW" = 2))


