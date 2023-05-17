#### Sensitivity Analysis for Informative Censoring

setwd("C:\\Users\\grace\\Documents\\WaterlooPHD\\Research\\IIW_Simulation\\MIPIW\\Simulation1")
source("MIPIW_timevarW_informativecensoring_functions.R")

library(knitr)
library(kableExtra)
require(reshape2)
require(ggplot2)
library(latex2exp)
#### n = 50 ###


set.seed(5465)
n = 50
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = -1
alpha1 = 1
var_phi_vec <- c(0, 0.1, 0.25, 0.5)
mu_c_vec <- c(0, 0.5, 1)
eta1 = 0
eta2 = 0
eta3 = 0
tau = 7
N = 1000

# results_censoring_n50_nocov <- simulateALLMIPI_timevar(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                         alpha0, alpha1, eta1, eta2, eta3, var_phi_vec, mu_c_vec, 
#                                                  tau, outputfulldatalist = T)

# saveRDS(results_censoring_n50_nocov, "results_censoring_n50_nocov.rds")
results_censoring_n50_nocov <- readRDS("results_censoring_n50_nocov.rds")
kable(results_censoring_n50_nocov$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 2, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))


set.seed(1235)
eta1 = 0.2
eta2 = 0.1
eta3 = 0.1

# results_censoring_n50_cov <- simulateALLMIPI_timevar(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                        alpha0, alpha1, eta1, eta2, eta3, var_phi_vec, mu_c_vec, 
#                                                        tau, outputfulldatalist = T)
# saveRDS(results_censoring_n50_cov, "results_censoring_n50_cov.rds")

results_censoring_n50_cov <- readRDS("results_censoring_n50_cov.rds")
kable(results_censoring_n50_cov$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 2, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))





#### n = 100 ###

set.seed(654321)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = -1
alpha1 = 1
var_phi_vec <- c(0, 0.1, 0.25, 0.5)
mu_c_vec <- c(0, 0.5, 1)
eta1 = 0
eta2 = 0
eta3 = 0
tau = 7
N = 1000

# results_censoring_n100_nocov <- simulateALLMIPI_timevar(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                         alpha0, alpha1, eta1, eta2, eta3, var_phi_vec, mu_c_vec,
#                                                  tau, outputfulldatalist = T)
# 
# saveRDS(results_censoring_n100_nocov, "results_censoring_n100_nocov.rds")
results_censoring_n100_nocov <- readRDS("results_censoring_n100_nocov.rds")
kable(results_censoring_n100_nocov$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 2, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))




set.seed(654891)
eta1 = 0.2
eta2 = 0.1
eta3 = 0.1

# results_censoring_n100_cov <- simulateALLMIPI_timevar(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                        alpha0, alpha1, eta1, eta2, eta3, var_phi_vec, mu_c_vec,
#                                                        tau, outputfulldatalist = T)
# saveRDS(results_censoring_n100_cov, "results_censoring_n100_cov.rds")

results_censoring_n100_cov <- readRDS("results_censoring_n100_cov.rds")
kable(results_censoring_n100_cov$resultsmat, booktabs = T, digits = 3,format = "latex") %>%
  add_header_above(c(" " = 2, "Unweighted" = 3, "IIW" = 3, "IPTW" = 3, "MIPIW" = 3))




# graphs for ENAR poster
# focus only on n = 100 and mu_c = 0

results_no_covariates <- results_censoring_n100_nocov$fulldatalist 
results_no_covariates_muc0 <- data.frame(rbind(results_no_covariates[[1]], results_no_covariates[[2]], 
                                               results_no_covariates[[3]], results_no_covariates[[4]]))
results_no_covariates_muc0$varphi <- as.factor(c(rep(0, 1000), rep(0.1, 1000), rep(0.25, 1000), rep(0.5, 1000)))

df_nocovariates <- melt(results_no_covariates_muc0, id.vars = "varphi")
colnames(df_nocovariates) <- c("varphi", "Method", "ATE")


ggplot(df_nocovariates, aes(x = Method, y = ATE, fill = varphi)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 1) + 
  scale_fill_discrete(name = TeX(r'($\sigma^2_{\epsilon}$)', bold=TRUE)) +
  ggtitle("Results: No Covariates Related to Censoring Hazard") +
  theme(text = element_text(size=20)) +
  ylim(c(-2, 3))
  




results_covariates <- results_censoring_n100_cov$fulldatalist 
results_covariates_muc0 <- data.frame(rbind(results_covariates[[1]], results_covariates[[2]], 
                                               results_covariates[[3]], results_covariates[[4]]))
results_covariates_muc0$varphi <- as.factor(c(rep(0, 1000), rep(0.1, 1000), rep(0.25, 1000), rep(0.5, 1000)))

df_covariates <- melt(results_covariates_muc0, id.vars = "varphi")
colnames(df_covariates) <- c("varphi", "Method", "ATE")

ggplot(df_covariates, aes(x = Method, y = ATE, fill = varphi)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 1) + 
  scale_fill_discrete(name = TeX(r'($\sigma^2_{\epsilon}$)', bold=TRUE)) +
  ggtitle("Results: Covariates Related to Censoring Hazard") +
  theme(text = element_text(size=20)) +
  ylim(c(-2, 3))





# 
# 
# # make graphs for presentation
# 
# 
# p1 <- results_censoring_n100$fulldatalist[[1]]
# colnames(p1) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# 
# png(file = ".//simulationresultsplot//fulldata_muc_0_varphi_0.png")
# boxplot(p1,
#         main = TeX(r'(Simulation Results: $\mu_c = 0$, $\sigma^2_{\phi} = 0$)', bold=TRUE),
#         ylim = c(-2, 2.5),
#         ylab = "ATE")
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# p2 <- results_censoring_n100$fulldatalist[[2]]
# colnames(p2) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# png(file = ".//simulationresultsplot//fulldata_muc_0_varphi_0.1.png")
# boxplot(p2,
#         main = TeX(r'(Simulation Results: $\mu_c = 0$, $\sigma^2_{\phi} = 0.10$)', bold=TRUE),
#         ylim = c(-2, 2.5),
#         ylab = "ATE")
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# 
# 
# p3 <- results_censoring_n100$fulldatalist[[3]]
# colnames(p3) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# png(file = ".//simulationresultsplot//fulldata_muc_0_varphi_0.25.png")
# boxplot(p3,
#         main = TeX(r'(Simulation Results: $\mu_c = 0$, $\sigma^2_{\phi} = 0.25$)', bold=TRUE),
#                    ylim = c(-2, 2.5),
#                    ylab = "ATE")
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# p4 <- results_censoring_n100$fulldatalist[[4]]
# colnames(p4) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# png(file = ".//simulationresultsplot//fulldata_muc_0_varphi_0.50.png")
# boxplot(p4,
#         main = TeX(r'(Simulation Results: $\mu_c = 0$, $\sigma^2_{\phi} = 0.50$)', bold=TRUE),
#         ylim = c(-2, 2.5),
#         ylab = "ATE")
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 
# 
# 
# 
# # 
# # p5 <- results_censoring_n100$fulldatalist[[5]]
# # colnames(p5) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# # png(file = ".//simulationresultsplot//fulldata_muc_0.5_varphi_0.png")
# # boxplot(p5,
# #         main = TeX(r'(Simulation Results: $\mu_c = 0.5$, $\sigma^2_{\phi} = 0$)', bold=TRUE),
# #         ylim = c(-2, 2.5),
# #         ylab = "ATE")
# # abline(h = 0.5, col = "red", lty = 3)
# # dev.off()
# # 
# # 
# # p6 <- results_censoring_n100$fulldatalist[[6]]
# # colnames(p6) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# # png(file = ".//simulationresultsplot//fulldata_muc_0.5_varphi_0.1.png")
# # boxplot(p6,
# #         main = TeX(r'(Simulation Results: $\mu_c = 0.5$, $\sigma^2_{\phi} = 0.1$)', bold=TRUE),
# #         ylim = c(-2, 2.5),
# #         ylab = "ATE")
# # abline(h = 0.5, col = "red", lty = 3)
# # dev.off()
# # 
# # p7 <- results_censoring_n100$fulldatalist[[7]]
# # colnames(p7) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# # png(file = ".//simulationresultsplot//fulldata_muc_0.5_varphi_0.25.png")
# # boxplot(p7,
# #         main = TeX(r'(Simulation Results: $\mu_c = 0.5$, $\sigma^2_{\phi} = 0.25$)', bold=TRUE),
# #         ylim = c(-2, 2.5),
# #         ylab = "ATE")
# # abline(h = 0.5, col = "red", lty = 3)
# # dev.off()
# # 
# # 
# # p8 <- results_censoring_n100$fulldatalist[[8]]
# # colnames(p8) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# # png(file = ".//simulationresultsplot//fulldata_muc_0.5_varphi_0.5.png")
# # boxplot(p8,
# #         main = TeX(r'(Simulation Results: $\mu_c = 0.5$, $\sigma^2_{\phi} = 0.5$)', bold=TRUE),
# #         ylim = c(-2, 2.5),
# #         ylab = "ATE")
# # abline(h = 0.5, col = "red", lty = 3)
# # dev.off()
# # 
# # 
# # p9 <- results_censoring_n100$fulldatalist[[9]]
# # colnames(p9) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# # png(file = ".//simulationresultsplot//fulldata_muc_1_varphi_0.png")
# # boxplot(p9,
# #         main = TeX(r'(Simulation Results: $\mu_c = 1$, $\sigma^2_{\phi} = 0$)', bold=TRUE),
# #         ylim = c(-2, 2.5),
# #         ylab = "ATE")
# # abline(h = 0.5, col = "red", lty = 3)
# # dev.off()
# # 
# # p10 <- results_censoring_n100$fulldatalist[[10]]
# # colnames(p10) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# # png(file = ".//simulationresultsplot//fulldata_muc_1_varphi_0.1.png")
# # boxplot(p10,
# #         main = TeX(r'(Simulation Results: $\mu_c = 1$, $\sigma^2_{\phi} = 0.1$)', bold=TRUE),
# #         ylim = c(-2, 2.5),
# #         ylab = "ATE")
# # abline(h = 0.5, col = "red", lty = 3)
# # dev.off()
# # 
# # 
# # p11 <- results_censoring_n100$fulldatalist[[11]]
# # colnames(p11) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# # png(file = ".//simulationresultsplot//fulldata_muc_1_varphi_0.25.png")
# # boxplot(p11,
# #         main = TeX(r'(Simulation Results: $\mu_c = 1$, $\sigma^2_{\phi} = 0.25$)', bold=TRUE),
# #         ylim = c(-2, 2.5),
# #         ylab = "ATE")
# # abline(h = 0.5, col = "red", lty = 3)
# # dev.off()
# # 
# # p12 <- results_censoring_n100$fulldatalist[[12]]
# # colnames(p12) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# # png(file = ".//simulationresultsplot//fulldata_muc_1_varphi_0.5.png")
# # boxplot(p12,
# #         main = TeX(r'(Simulation Results: $\mu_c = 1$, $\sigma^2_{\phi} = 0.5$)', bold=TRUE),
# #         ylim = c(-2, 2.5),
# #         ylab = "ATE")
# # abline(h = 0.5, col = "red", lty = 3)
# # dev.off()
# 
# 
# # example with assumptions satisfied
# resultsMIPI_all_timevar_n100 <- readRDS("resultsMIPI_timevar_n100.rds")
# test <- resultsMIPI_all_timevar_n100$fulldatalist[[5]]
# colnames(test) <- c("Naive", "IIW", "IPTW", "FIPTIW")
# png(file = ".//simulationresultsplot//tester.png")
# boxplot(test,
#         main = TeX(r'(Simulation Results: Noninformative Censoring)', bold=TRUE),
#         ylab = "ATE", ylim = c(-2, 2.5))
# abline(h = 0.5, col = "red", lty = 3)
# dev.off()
# 







