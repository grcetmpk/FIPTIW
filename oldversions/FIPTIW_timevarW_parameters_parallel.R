
#~~~~~~~~~~~~~~~~~~ MAIN SIMULATION STUDY ~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

source("MIPIW_timevarW_functions_parallel1.R")
library(knitr)
library(kableExtra)

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
# 
# resultsMIPI_all_timevar_n100_parallel <- simulateALLMIPI_timevar(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
#                                                         alpha0, alpha1vec, tau, outputfulldatalist = FALSE)
# saveRDS(resultsMIPI_all_timevar_n100_parallel, "resultsMIPI_timevar_n100_parallel.rds")
resultsMIPI_all_timevar_n100_parallel <- readRDS("resultsMIPI_timevar_n100_parallel.rds")


# library(latex2exp)
# names(resultsMIPI_all_timevar_n100_parallel$fulldatalist)
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0_gamma3_0_alpha1_0.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[1]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0$, $\gamma_3 = 0$,$\alpha_1 = 0$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0_gamma3_0_alpha1_1.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[2]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0$, $\gamma_3 = 0$,$\alpha_1 = 1$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0_gamma3_0_alpha1_3.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[3]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0$, $\gamma_3 = 0$,$\alpha_1 = 3$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0_gamma3_0.6_alpha1_0.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[4]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0$, $\gamma_3 = 0.6$,$\alpha_1 = 0$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0_gamma3_0.6_alpha1_1.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[5]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0$, $\gamma_3 = 0.6$,$\alpha_1 = 1$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0_gamma3_0.6_alpha1_3.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[6]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0$, $\gamma_3 = 0.6$,$\alpha_1 = 3$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0.3_gamma3_0_alpha1_0.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[7]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0.3$, $\gamma_3 = 0$,$\alpha_1 = 0$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0.3_gamma3_0_alpha1_1.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[8]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0.3$, $\gamma_3 = 0$,$\alpha_1 = 1$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0.3_gamma3_0_alpha1_3.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[9]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0.3$, $\gamma_3 = 0$,$\alpha_1 = 3$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0.3_gamma3_0.6_alpha1_0.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[10]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0.3$, $\gamma_3 = 0.6$, $\alpha_1 = 0$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# 
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0.3_gamma3_0.6_alpha1_1.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[11]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0.3$, $\gamma_3 = 0.6$, $\alpha_1 = 1$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# png(file = ".//simulationresultsplot//fulldata_gamma2_0.3_gamma3_0.6_alpha1_3.png")
# boxplot(resultsMIPI_all_timevar_n100_parallel$fulldatalist[[12]],
#         main = TeX(r'(Simulation Results: $\gamma_2 = 0.3$, $\gamma_3 = 0.6$, $\alpha_1 = 3$)', bold=TRUE))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()



kable(resultsMIPI_all_timevar_n100_parallel$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))




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
# 
# resultsMIPI_all_timevar_n50_parallel <- simulateALLMIPI_timevar(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
#                                                        alpha0, alpha1vec, tau, outputfulldatalist = FALSE)
# saveRDS(resultsMIPI_all_timevar_n50_parallel, "resultsMIPI_timevar_n50_parallel.rds")
resultsMIPI_all_timevar_n50_parallel <- readRDS("resultsMIPI_timevar_n50_parallel.rds")

kable(resultsMIPI_all_timevar_n50_parallel$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))



### n = 500 ###

source("MIPIW_timevarW_functions.R") #using non-parallelized code as its crashing
set.seed(41563)
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

# resultsMIPI_all_timevar_n500 <- simulateALLMIPI_timevar(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
#                                                          alpha0, alpha1vec, tau, outputfulldatalist = FALSE)
# saveRDS(resultsMIPI_all_timevar_n500, "resultsMIPI_timevar_n500.rds")
resultsMIPI_all_timevar_n500 <- readRDS("resultsMIPI_timevar_n500.rds")

kable(resultsMIPI_all_timevar_n500$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))







# ~~~~~~~~~~~~ Variable selection simulation (Appendix) ~~~~~~~~~~~~~~~~~~~ #
# 
# set.seed(500)
# n = 100
# beta1 = 0.5
# beta2 = 2
# beta3 = 1
# gamma1 = 0.5
# gamma2 = 0
# gamma3vec = c(0,0.6)
# alpha0 = 0.5
# alpha1 = 0
# tau = 7
# N = 1000
# 
# varsel_timevar <- simulateALLVarSel_timevar(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau)
# kable(varsel_timevar$resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
#   add_header_above(c("Non-stabilized " = 1, "Variables used to estimate intensity" = 9))


