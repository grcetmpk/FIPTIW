
#~~~~~~~~~~~~~~~~~~ SIMULATiONS I and II ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
######~~~~~~~~~~~~~~CONTINUOUS OUTCOME~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

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
source("FIPTIW_timevarW_functions.R")

# parallelization stuff
ncores <- detectCores() - 3
nclusters <- makeCluster(ncores)

require(knitr)
require(kableExtra)
require(ggplot2)
require(ggpubr)


#~~~~~~~~~~ FIPTIW simulation, assumptions satisfied (not shown in paper).

#### n = 100 ###

set.seed(2345)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3)
gamma3vec = c(0, 0.6)
alpha0 = -1
alpha1vec = c(0, 1)
tau = 7
N = 50
ncutpts = 2

resultsFIPTIW_n100 <- simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec,
                                                        alpha0, alpha1vec, tau, outputfulldatalist = FALSE, ncutpts = 2,
                                        inParallel = F,  outcome = "continuous")
saveRDS(resultsFIPTIW_n100, "resultsFIPTIW_n100.rds")

kable(resultsFIPTIW_n100$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 5, "IIW" = 5, "IPTW" = 5, "FIPTIW" = 5))





# 
# 
# ### n = 50 ###
# 
# set.seed(6234)
# n = 50
# beta1 = 0.5
# beta2 = 2
# beta3 = 1
# gamma1 = 0.5
# gamma2vec = c(0, 0.3)
# gamma3vec = c(0, 0.6)
# alpha0 = -1
# alpha1vec = c(0, 1, 3)
# tau = 7
# N = 1000
# 
# # resultsFIPTIW_n50 <- simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
# #                                                        alpha0, alpha1vec, tau, outputfulldatalist = FALSE)
# #saveRDS(resultsFIPTIW_n50, "resultsFIPTIW_n50.rds")
# resultsFIPTIW_n50 <- readRDS("resultsFIPTIW_n50.rds")
# 
# 
# ### n = 500 ###
# 
# 
# set.seed(23424)
# n = 500
# beta1 = 0.5
# beta2 = 2
# beta3 = 1
# gamma1 = 0.5
# gamma2vec = c(0, 0.3)
# gamma3vec = c(0, 0.6)
# alpha0 = -1
# alpha1vec = c(0, 1, 3)
# tau = 7
# N = 1000
# 
# # resultsFIPTIW_n500 <- simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
# #                                                         alpha0, alpha1vec, tau, outputfulldatalist = FALSE)
# # saveRDS(resultsFIPTIW_n500, "resultsFIPTIW_n500.rds")
# resultsFIPTIW_n500 <- readRDS("resultsFIPTIW_n500.rds")
# 
# kable(resultsFIPTIW_n500$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
#   add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "FIPTIW" = 3))
# 
# 



#### ~~~~~~~~~~~~ Variable selection simulation (Sim II) ~~~~~~~~~~~~~~~~~~~ #

set.seed(54641230)
n = 100
beta1 = 0.5
beta2vec = c(0, 2)
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3, 0.6)
gamma3 = 0.6
alpha0 = 0
alpha1 = 0
tau = 7
N = 1000
# 
# varsel_n100 <- simulateALLVarSel(N, n, beta1, beta2vec, beta3, gamma1, gamma2vec, 
#                                  gamma3, alpha0, alpha1, tau, ncutpts = 2, 
#                                  inParallel = T, nclusters = nclusters)
saveRDS(varsel_n100, "varsel_n100.rds")
varsel_n100 <- readRDS("varsel_n100.rds")

kable(varsel_n100$resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
  add_header_above(c(" " = 3, "Variables used to estimate intensity" = 8))



set.seed(39487)
n = 50
beta1 = 0.5
beta2vec = c(0, 2)
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3, 0.6)
gamma3 = 0.6
alpha0 = 0
alpha1 = 0
tau = 7
N = 1000
ncutpts = 2
#
# varsel_n50 <- simulateALLVarSel(N, n, beta1, beta2vec, beta3, gamma1, gamma2vec, gamma3, 
#                                 alpha0, alpha1, tau, ncutpts = 1, inParallel = T, nclusters = nclusters)
# saveRDS(varsel_n50, "varsel_n50.rds")
varsel_n50 <- readRDS("varsel_n50.rds")
kable(varsel_n50$resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
  add_header_above(c(" " = 3, "Variables used to estimate intensity" = 8))



set.seed(5465)
n = 500
beta1 = 0.5
beta2vec = c(0, 2)
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3, 0.6)
gamma3 = 0.6
alpha0 = 0
alpha1 = 0
tau = 7
N = 1000
ncutpts = 2

# varsel_timevar_n500 <- simulateALLVarSel(N, n, beta1, beta2vec,
#                                          beta3, gamma1, gamma2vec, gamma3, alpha0, alpha1, tau,
#                                          inParallel = T, nclusters, ncutpts = 2)
# saveRDS(varsel_timevar_n500, "varsel_timevar_n1000.rds")

kable(varsel_timevar_n500$resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
  add_header_above(c(" " = 3, "Variables used to estimate intensity" = 8))




#~~~~~~~~~~~~~~ Violations of Noninform Censoring Sim with IPCW~~~~~~~~~#

#### n = 100

set.seed(234)
n = 100
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
censinform = T
eta1vec = c(0)
eta2vec = c(0, 0.2, 0.5)
eta3vec = c(0, 0.4, 0.6)
ncutpts = 2

results_censoring_n100 <- simulateALLFIPTICW_CENS(N, n, beta1, beta2, beta3, gamma1,
                                                 gamma2vec, gamma3vec, alpha0, alpha1, tau, outputfulldatalist = T,
                                                 censinform = T, eta1vec, eta2vec, eta3vec, ncutpts)
saveRDS(results_censoring_n100, "results_censoring_n100.rds")
results_censoring_n100 <- readRDS("results_censoring_n100.rds")

plot_censoring_n100 <- makeBoxPlotsATE(results_censoring_n100$fulldatalist, eta2vec)
annotate_figure(plot_censoring_n100, top = text_grob("Combined Results for n = 100", face = "bold", size = 20))


coveragemat_n100 <- results_censoring_n100$resultsmat[, c(2,3,7,8,12,13,17,18,22,23, 27,28)]
kable(coveragemat_n100, booktabs = T, 
      digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2, "FIPTICW" = 2))




##### n = 50

set.seed(2432)
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
eta1vec = c(0)
eta2vec = c(0, 0.2, 0.5)
eta3vec = c(0, 0.4, 0.6)


results_censoring_n50 <- simulateALLFIPTICW_CENS(N, n, beta1, beta2, beta3, gamma1,
                                                 gamma2vec, gamma3vec, alpha0, alpha1, tau, outputfulldatalist = T,
                                                 censinform = T, eta1vec, eta2vec, eta3vec, ncutpts = 1)
saveRDS(results_censoring_n50, "results_censoring_n50.rds")
results_censoring_n50 <- readRDS("results_censoring_n50.rds")

plot_censoring_n50 <- makeBoxPlotsATE(results_censoring_n50$fulldatalist, eta2vec)
annotate_figure(plot_censoring_n50, top = text_grob("Combined Results for n = 50", face = "bold", size = 20))



coveragemat_n50 <- results_censoring_n50$resultsmat[, c(2,3,7,8,12,13,17,18,22,23, 27,28)]
kable(coveragemat_n50, booktabs = T, 
      digits = 3,format = "latex") %>%
  add_header_above(c(" " = 2, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2, "FIPTICW" = 2))



##### n = 500

set.seed(444)
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
eta1vec = c(0)
eta2vec = c(0, 0.2, 0.5)
eta3vec = c(0, 0.4, 0.6)
ncutpts = 2

results_censoring_n500 <- simulateALLFIPTICW_CENS(N, n, beta1, beta2, beta3, gamma1,
                                                      gamma2vec, gamma3vec, alpha0, alpha1, tau, outputfulldatalist = T,
                                                      censinform = T, eta1vec, eta2vec, eta3vec, ncutpts)
saveRDS(results_censoring_n500, "results_censoring_n500.rds")
results_censoring_n500 <- readRDS("results_censoring_n500.rds")

plot_censoring_n500 <- makeBoxPlotsATE(results_censoring_n500$fulldatalist, eta2vec)
annotate_figure(plot_censoring_n500, top = text_grob("Combined Results for n = 500", face = "bold", size = 20))


coveragemat <- results_censoring_n500$resultsmat[, c(2,3,7,8,12,13,17,18,22,23, 27,28)]
kable(coveragemat, booktabs = T, 
      digits = 3,format = "latex") %>%
  add_header_above(c(" " = 2, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2, "FIPTICW" = 2))




######################3 BINARY RESULTS

source("FIPTIW_timevarW_functions.R")

# parallelization stuff
ncores <- detectCores() - 3
nclusters <- makeCluster(ncores)

require(knitr)
require(kableExtra)
require(ggplot2)
require(ggpubr)


#~~~~~~~~~~ FIPTIW simulation, assumptions satisfied (not shown in paper).

#### n = 100 ###

set.seed(2345)
n = 100
beta1 = 0.5
beta2 = 0.5 
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0)
gamma3vec = c(0)
alpha0 = -1
alpha1vec = c(0)
tau = 7
N = 100
ncutpts = 2

resultsFIPTIW_n100_binary <- simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec,
                                        alpha0, alpha1vec, tau, outputfulldatalist = FALSE, ncutpts = 2,
                                        inParallel = F,  outcome = "binary")


