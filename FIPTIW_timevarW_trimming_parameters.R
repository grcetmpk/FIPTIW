#~~~~~~~~~~~~~ WEIGHT TRIMMING ~~~~~~~~~~~~~~~~~~#

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

require(knitr)
require(kableExtra)
require(reshape2)
require(ggplot2)
require(latex2exp)
require(patchwork)
require(grid)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ low IPTW weights, low IIW weights
#### n = 100 ###

set.seed(13546)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 0.5
tau = 7
N = 1000
baseline = F


# weighttrimresults_lowIPTW_lowIIW<- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                             alpha0, alpha1, tau, baseline)
# saveRDS(weighttrimresults_lowIPTW_lowIIW, "weighttrimresults_lowIPTW_lowIIW.rds")
weighttrimresults_lowIPTW_lowIIW <- readRDS("weighttrimresults_lowIPTW_lowIIW.rds")


# view average weight distribution
weighttrimresults_lowIPTW_lowIIW$othermat_iptw
weighttrimresults_lowIPTW_lowIIW$othermat_iiw
weighttrimresults_lowIPTW_lowIIW$othermat_fiptiw




### create dataset for data vis
biasmat <- weighttrimresults_lowIPTW_lowIIW$biasmat
colnames(biasmat) <- c("Trimmed Before", "Trimmed After")
biasmat <- melt(biasmat)
biasmat$value <- biasmat$value/beta1
varmat <- weighttrimresults_lowIPTW_lowIIW$varmat
colnames(varmat) <- c("Trimmed Before", "Trimmed After")
varmat <- melt(varmat)
msemat <- weighttrimresults_lowIPTW_lowIIW$msemat
colnames(msemat) <- c("Trimmed Before", "Trimmed After")
msemat <- melt(msemat)
resultsmat_lowIPTW_lowIIW <- cbind(biasmat, varmat$value, msemat$value)
colnames(resultsmat_lowIPTW_lowIIW) <- c("Percentile", "Group", "Bias", "Variance", "MSE")


# bias plot
p1_lowIPTW_lowIIW <- ggplot(resultsmat_lowIPTW_lowIIW, aes(x = Percentile, y = Bias, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+ 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)") +
  ylab("Rel. Bias")+
  ylim(-0.3, 0.9) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))

### Variance plots

p2_lowIPTW_lowIIW<- ggplot(resultsmat_lowIPTW_lowIIW, aes(x = Percentile, y = Variance, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) + 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)")+
  ylim(0.015, 0.03)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))

### MSE plots


p3_lowIPTW_lowIIW <- ggplot(resultsmat_lowIPTW_lowIIW, aes(x = Percentile, y = MSE, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8)+ 
  theme(legend.title=element_blank()) +
  xlab("Threshold (Percentile)")+
  ylim(0, 0.2)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))


combined_lowIPTW_lowIIW <- p1_lowIPTW_lowIIW + p2_lowIPTW_lowIIW + p3_lowIPTW_lowIIW & theme(legend.position = "bottom")
# combined_lowIPTW_lowIIW + 
#   plot_layout(guides = "collect") +
#   plot_annotation(title = "Results: Low IPTW, Low IIW")
# 

### determine minimum bias, var, mse, etc

#minimum bias
resultsmat_lowIPTW_lowIIW[resultsmat_lowIPTW_lowIIW$Bias == min(abs(resultsmat_lowIPTW_lowIIW$Bias)), ]
resultsmat_lowIPTW_lowIIW[resultsmat_lowIPTW_lowIIW$Variance == min(abs(resultsmat_lowIPTW_lowIIW$Variance)), ]
resultsmat_lowIPTW_lowIIW[resultsmat_lowIPTW_lowIIW$MSE == min(abs(resultsmat_lowIPTW_lowIIW$MSE)), ]





####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Moderate IPTW, low IIW

#### n = 100 ###

set.seed(442)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 3.5
tau = 7
N = 1000

weightrimresults_moderateIPTW_lowIIW <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
                                           alpha0, alpha1, tau, baseline = F)
saveRDS(weightrimresults_moderateIPTW_lowIIW, "weightrimresults_moderateIPTW_lowIIW.rds")
weightrimresults_moderateIPTW_lowIIW <- readRDS("weightrimresults_moderateIPTW_lowIIW.rds")


# view average weight distribution

weightrimresults_moderateIPTW_lowIIW$othermat_iptw
weightrimresults_moderateIPTW_lowIIW$othermat_iiw
weightrimresults_moderateIPTW_lowIIW$othermat_fiptiw


### create df for plots

biasmat_moderateIPTW_lowIIW <- weightrimresults_moderateIPTW_lowIIW$biasmat
colnames(biasmat_moderateIPTW_lowIIW) <- c("Trimmed Before", "Trimmed After")
biasmat_moderateIPTW_lowIIW <- melt(biasmat_moderateIPTW_lowIIW)
biasmat_moderateIPTW_lowIIW$value <- biasmat_moderateIPTW_lowIIW$value/beta1

varmat_moderateIPTW_lowIIW <- weightrimresults_moderateIPTW_lowIIW$varmat
colnames(varmat_moderateIPTW_lowIIW) <- c("Trimmed Before", "Trimmed After")
varmat_moderateIPTW_lowIIW <- melt(varmat_moderateIPTW_lowIIW)

msemat_moderateIPTW_lowIIW <- weightrimresults_moderateIPTW_lowIIW$msemat
colnames(msemat_moderateIPTW_lowIIW) <- c("Trimmed Before", "Trimmed After")
msemat_moderateIPTW_lowIIW <- melt(msemat_moderateIPTW_lowIIW)

resultsmat_moderateIPTW_lowIIW <- cbind(biasmat_moderateIPTW_lowIIW, varmat_moderateIPTW_lowIIW$value, msemat_moderateIPTW_lowIIW$value)
colnames(resultsmat_moderateIPTW_lowIIW) <- c("Percentile", "Group", "Bias", "Variance", "MSE")

# bias plot
p1_moderateIPTW_lowIIW<- ggplot(resultsmat_moderateIPTW_lowIIW, aes(x = Percentile, y = Bias, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+ 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)") +
  ylab("Rel. Bias")+
  ylim(-0.3, 0.9) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))


### Variance plots
p2_moderateIPTW_lowIIW<- ggplot(resultsmat_moderateIPTW_lowIIW, aes(x = Percentile, y = Variance, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) + 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)")+
  ylim(0.015, 0.03)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))



### MSE plots
p3_moderateIPTW_lowIIW <- ggplot(resultsmat_moderateIPTW_lowIIW, aes(x = Percentile, y = MSE, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8)+ 
  theme(legend.title=element_blank()) +
  xlab("Threshold (Percentile)")+
  ylim(0, 0.2)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))


combined_moderateIPTW_lowIIW <- p1_moderateIPTW_lowIIW + p2_moderateIPTW_lowIIW + p3_moderateIPTW_lowIIW & theme(legend.position = "bottom")
# combined_moderateIPTW_lowIIW + 
#   plot_layout(guides = "collect") +
#   plot_annotation(title = "Results: Moderate IPTW, Low IIW")

### determine minimum bias, var, mse, etc
resultsmat_moderateIPTW_lowIIW_first <- resultsmat_moderateIPTW_lowIIW[resultsmat_moderateIPTW_lowIIW$Group == "Trimmed Before", ]
resultsmat_moderateIPTW_lowIIW_after <- resultsmat_moderateIPTW_lowIIW[resultsmat_moderateIPTW_lowIIW$Group == "Trimmed After", ]

#trimming first
resultsmat_moderateIPTW_lowIIW_first[resultsmat_moderateIPTW_lowIIW_first$Bias == min(abs(resultsmat_moderateIPTW_lowIIW_first$Bias)), ]
resultsmat_moderateIPTW_lowIIW_first[resultsmat_moderateIPTW_lowIIW_first$Variance == min(abs(resultsmat_moderateIPTW_lowIIW_first$Variance)), ]
resultsmat_moderateIPTW_lowIIW_first[resultsmat_moderateIPTW_lowIIW_first$MSE == min(abs(resultsmat_moderateIPTW_lowIIW_first$MSE)), ]


#trimming after
resultsmat_moderateIPTW_lowIIW_after[resultsmat_moderateIPTW_lowIIW_after$Bias == min(abs(resultsmat_moderateIPTW_lowIIW_after$Bias)), ]
resultsmat_moderateIPTW_lowIIW_after[resultsmat_moderateIPTW_lowIIW_after$Variance == min(abs(resultsmat_moderateIPTW_lowIIW_after$Variance)), ]
resultsmat_moderateIPTW_lowIIW_after[resultsmat_moderateIPTW_lowIIW_after$MSE == min(abs(resultsmat_moderateIPTW_lowIIW_after$MSE)), ]









######## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~high IPTW, low IIW
#### n = 100 ###

set.seed(6984)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 5.5
tau = 7
N = 1000

weightrimresults_highIPTW_lowIIW <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
                                          alpha0, alpha1, tau, baseline = F)
saveRDS(weightrimresults_highIPTW_lowIIW, "weightrimresults_highIPTW_lowIIW.rds")
weightrimresults_highIPTW_lowIIW <- readRDS("weightrimresults_highIPTW_lowIIW.rds")

weightrimresults_highIPTW_lowIIW$othermat_iptw
weightrimresults_highIPTW_lowIIW$othermat_iiw
weightrimresults_highIPTW_lowIIW$othermat_fiptiw




### create df for plots

biasmat_highIPTW_lowIIW <- weightrimresults_highIPTW_lowIIW$biasmat
colnames(biasmat_highIPTW_lowIIW) <- c("Trimmed Before", "Trimmed After")
biasmat_highIPTW_lowIIW<- melt(biasmat_highIPTW_lowIIW)
biasmat_highIPTW_lowIIW$value <- biasmat_highIPTW_lowIIW$value/beta1


varmat_highIPTW_lowIIW<- weightrimresults_highIPTW_lowIIW$varmat
colnames(varmat_highIPTW_lowIIW) <- c("Trimmed Before", "Trimmed After")
varmat_highIPTW_lowIIW <- melt(varmat_highIPTW_lowIIW)

msemat_highIPTW_lowIIW <- weightrimresults_highIPTW_lowIIW$msemat
colnames(msemat_highIPTW_lowIIW) <- c("Trimmed Before", "Trimmed After")
msemat_highIPTW_lowIIW <- melt(msemat_highIPTW_lowIIW)

resultsmat_highIPTW_lowIIW <- cbind(biasmat_highIPTW_lowIIW, varmat_highIPTW_lowIIW$value, msemat_highIPTW_lowIIW$value)
colnames(resultsmat_highIPTW_lowIIW) <- c("Percentile", "Group", "Bias", "Variance", "MSE")

# bias plot
p1_highIPTW_lowIIW <- ggplot(resultsmat_highIPTW_lowIIW, aes(x = Percentile, y = Bias, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+ 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)") +
  ylab("Rel. Bias") +
  ylim(-0.3, 0.9) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))

### Variance plots
p2_highIPTW_lowIIW <- ggplot(resultsmat_highIPTW_lowIIW, aes(x = Percentile, y = Variance, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) + 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)")+
  ylim(0.015, 0.03)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))



### MSE plots
p3_highIPTW_lowIIW <- ggplot(resultsmat_highIPTW_lowIIW, aes(x = Percentile, y = MSE, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8)+ 
  theme(legend.title=element_blank()) +
  xlab("Threshold (Percentile)")+
  ylim(0, 0.2)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))


combined_highIPTW_lowIIW <- p1_highIPTW_lowIIW + p2_highIPTW_lowIIW + p3_highIPTW_lowIIW & theme(legend.position = "bottom")
# combined_highIPTW_lowIIW + 
#   plot_layout(guides = "collect") +
#   plot_annotation(title = "Results: High IPTW, Low IIW")


### determine minimum bias, var, mse, etc
resultsmat_highIPTW_lowIIW_first <- resultsmat_highIPTW_lowIIW[resultsmat_highIPTW_lowIIW$Group == "Trimmed Before", ]
resultsmat_highIPTW_lowIIW_after <- resultsmat_highIPTW_lowIIW[resultsmat_highIPTW_lowIIW$Group == "Trimmed After", ]


#trimming first
resultsmat_highIPTW_lowIIW_first[resultsmat_highIPTW_lowIIW_first$Bias == min(abs(resultsmat_highIPTW_lowIIW_first$Bias)), ]
resultsmat_highIPTW_lowIIW_first[resultsmat_highIPTW_lowIIW_first$Variance == min(abs(resultsmat_highIPTW_lowIIW_first$Variance)), ]
resultsmat_highIPTW_lowIIW_first[resultsmat_highIPTW_lowIIW_first$MSE == min(abs(resultsmat_highIPTW_lowIIW_first$MSE)), ]

#trimming after
resultsmat_highIPTW_lowIIW_after[resultsmat_highIPTW_lowIIW_after$Bias == min(abs(resultsmat_highIPTW_lowIIW_after$Bias)), ]
resultsmat_highIPTW_lowIIW_after[resultsmat_highIPTW_lowIIW_after$Variance == min(abs(resultsmat_highIPTW_lowIIW_after$Variance)), ]
resultsmat_highIPTW_lowIIW_after[resultsmat_highIPTW_lowIIW_after$MSE == min(abs(resultsmat_highIPTW_lowIIW_after$MSE)), ]




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~low IPTW and moderate IIW
#### n = 100 ###

set.seed(462)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = -.75
alpha0 = 0
alpha1 = 0.5
tau = 7
N = 1000

weighttrimresults_lowIPTW_moderateIIW<- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
                                          alpha0, alpha1, tau, baseline = F)
saveRDS(weighttrimresults_lowIPTW_moderateIIW, "weighttrimresults_lowIPTW_moderateIIW.rds")
weighttrimresults_lowIPTW_moderateIIW <- readRDS("weighttrimresults_lowIPTW_moderateIIW.rds")


# view average weight distribution
weighttrimresults_lowIPTW_moderateIIW$othermat_iptw
weighttrimresults_lowIPTW_moderateIIW$othermat_iiw
weighttrimresults_lowIPTW_moderateIIW$othermat_fiptiw


biasmat_lowIPTW_moderateIIW <- weighttrimresults_lowIPTW_moderateIIW$biasmat
colnames(biasmat_lowIPTW_moderateIIW) <- c("Trimmed Before", "Trimmed After")
biasmat_lowIPTW_moderateIIW<- melt(biasmat_lowIPTW_moderateIIW)

varmat_lowIPTW_moderateIIW<- weighttrimresults_lowIPTW_moderateIIW$varmat
colnames(varmat_lowIPTW_moderateIIW) <- c("Trimmed Before", "Trimmed After")
varmat_lowIPTW_moderateIIW <- melt(varmat_lowIPTW_moderateIIW)

msemat_lowIPTW_moderateIIW <- weighttrimresults_lowIPTW_moderateIIW$msemat
colnames(msemat_lowIPTW_moderateIIW) <- c("Trimmed Before", "Trimmed After")
msemat_lowIPTW_moderateIIW <- melt(msemat_lowIPTW_moderateIIW)

resultsmat_lowIPTW_moderateIIW <- cbind(biasmat_lowIPTW_moderateIIW, varmat_lowIPTW_moderateIIW$value, msemat_lowIPTW_moderateIIW$value)
colnames(resultsmat_lowIPTW_moderateIIW) <- c("Percentile", "Group", "Bias", "Variance", "MSE")

# bias plot
p1_lowIPTW_moderateIIW <- ggplot(resultsmat_lowIPTW_moderateIIW, aes(x = Percentile, y = Bias, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+ 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)")+
  ylab("Rel. Bias")+
  ylim(-0.3, 0.9)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))

### Variance plots
p2_lowIPTW_moderateIIW <- ggplot(resultsmat_lowIPTW_moderateIIW, aes(x = Percentile, y = Variance, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) + 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)")+
  ylim(0.015, 0.03)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))


### MSE plots
p3_lowIPTW_moderateIIW <- ggplot(resultsmat_lowIPTW_moderateIIW, aes(x = Percentile, y = MSE, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8)+ 
  theme(legend.title=element_blank()) +
  xlab("Threshold (Percentile)")+
  ylim(0, 0.2)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))


combined_lowIPTW_moderateIIW <- p1_lowIPTW_moderateIIW + p2_lowIPTW_moderateIIW + p3_lowIPTW_moderateIIW & theme(legend.position = "bottom")
# combined_lowIPTW_moderateIIW + 
#   plot_layout(guides = "collect") +
#   plot_annotation(title = "Results: Moderate IPTW, Moderate IIW")


### determine minimum bias, var, mse, etc
resultsmat_lowIPTW_moderateIIW_first <- resultsmat_lowIPTW_moderateIIW[resultsmat_lowIPTW_moderateIIW$Group == "Trimmed Before", ]
resultsmat_lowIPTW_moderateIIW_after <- resultsmat_lowIPTW_moderateIIW[resultsmat_lowIPTW_moderateIIW$Group == "Trimmed After", ]



#trimming first
resultsmat_lowIPTW_moderateIIW_first[resultsmat_lowIPTW_moderateIIW_first$Bias == -min(abs(resultsmat_lowIPTW_moderateIIW_first$Bias)), ]
resultsmat_lowIPTW_moderateIIW_first[resultsmat_lowIPTW_moderateIIW_first$Variance == min(abs(resultsmat_lowIPTW_moderateIIW_first$Variance)), ]
resultsmat_lowIPTW_moderateIIW_first[resultsmat_lowIPTW_moderateIIW_first$MSE == min(abs(resultsmat_lowIPTW_moderateIIW_first$MSE)), ]

#trimming after
resultsmat_lowIPTW_moderateIIW_after[resultsmat_lowIPTW_moderateIIW_after$Bias == -min(abs(resultsmat_lowIPTW_moderateIIW_after$Bias)), ]
resultsmat_lowIPTW_moderateIIW_after[resultsmat_lowIPTW_moderateIIW_after$Variance == min(abs(resultsmat_lowIPTW_moderateIIW_after$Variance)), ]
resultsmat_lowIPTW_moderateIIW_after[resultsmat_lowIPTW_moderateIIW_after$MSE == min(abs(resultsmat_lowIPTW_moderateIIW_after$MSE)), ]




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ low IPTW and high IIW
#### n = 100 ###

set.seed(5321)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = -1
alpha0 = 0
alpha1 = 0.5
tau = 7
N = 1000
# 
weighttrimresults_lowIPTW_highIIW<- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
                                                                      alpha0, alpha1, tau, baseline = F)
saveRDS(weighttrimresults_lowIPTW_highIIW, "weighttrimresults_lowIPTW_highIIW.rds")
weighttrimresults_lowIPTW_highIIW <- readRDS("weighttrimresults_lowIPTW_highIIW.rds")


# view average weight distribution
weighttrimresults_lowIPTW_highIIW$othermat_iptw
weighttrimresults_lowIPTW_highIIW$othermat_iiw
weighttrimresults_lowIPTW_highIIW$othermat_fiptiw


biasmat_lowIPTW_highIIW <- weighttrimresults_lowIPTW_highIIW$biasmat
colnames(biasmat_lowIPTW_highIIW) <- c("Trimmed Before", "Trimmed After")
biasmat_lowIPTW_highIIW<- melt(biasmat_lowIPTW_highIIW)

varmat_lowIPTW_highIIW<- weighttrimresults_lowIPTW_highIIW$varmat
colnames(varmat_lowIPTW_highIIW) <- c("Trimmed Before", "Trimmed After")
varmat_lowIPTW_highIIW <- melt(varmat_lowIPTW_highIIW)

msemat_lowIPTW_highIIW <- weighttrimresults_lowIPTW_highIIW$msemat
colnames(msemat_lowIPTW_highIIW) <- c("Trimmed Before", "Trimmed After")
msemat_lowIPTW_highIIW <- melt(msemat_lowIPTW_highIIW)

resultsmat_lowIPTW_highIIW <- cbind(biasmat_lowIPTW_highIIW, varmat_lowIPTW_highIIW$value, msemat_lowIPTW_highIIW$value)
colnames(resultsmat_lowIPTW_highIIW) <- c("Percentile", "Group", "Bias", "Variance", "MSE")

# bias plot
p1_lowIPTW_highIIW <- ggplot(resultsmat_lowIPTW_highIIW, aes(x = Percentile, y = Bias, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+ 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)")+
  ylab("Rel. Bias")+
  ylim(-0.3, 0.9)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))

### Variance plots
p2_lowIPTW_highIIW <- ggplot(resultsmat_lowIPTW_highIIW, aes(x = Percentile, y = Variance, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) + 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)")+
  ylim(0.015, 0.03)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))

### MSE plots
p3_lowIPTW_highIIW <- ggplot(resultsmat_lowIPTW_highIIW, aes(x = Percentile, y = MSE, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8)+ 
  theme(legend.title=element_blank()) +
  xlab("Threshold (Percentile)")+
  ylim(0, 0.2)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))


combined_lowIPTW_highIIW <- p1_lowIPTW_highIIW + p2_lowIPTW_highIIW + p3_lowIPTW_highIIW & theme(legend.position = "bottom")
# combined_lowIPTW_highIIW + 
#   plot_layout(guides = "collect") +
#   plot_annotation(title = "Results: Moderate IPTW, Moderate IIW")


### determine minimum bias, var, mse, etc
resultsmat_lowIPTW_highIIW_first <- resultsmat_lowIPTW_highIIW[resultsmat_lowIPTW_highIIW$Group == "Trimmed Before", ]
resultsmat_lowIPTW_highIIW_after <- resultsmat_lowIPTW_highIIW[resultsmat_lowIPTW_highIIW$Group == "Trimmed After", ]




#trimming first
resultsmat_lowIPTW_highIIW_first[resultsmat_lowIPTW_highIIW_first$Bias == -min(abs(resultsmat_lowIPTW_highIIW_first$Bias)), ]
resultsmat_lowIPTW_highIIW_first[resultsmat_lowIPTW_highIIW_first$Variance == min(abs(resultsmat_lowIPTW_highIIW_first$Variance)), ]
resultsmat_lowIPTW_highIIW_first[resultsmat_lowIPTW_highIIW_first$MSE == min(abs(resultsmat_lowIPTW_highIIW_first$MSE)), ]

#trimming after
resultsmat_lowIPTW_highIIW_after[resultsmat_lowIPTW_highIIW_after$Bias == -min(abs(resultsmat_lowIPTW_highIIW_after$Bias)), ]
resultsmat_lowIPTW_highIIW_after[resultsmat_lowIPTW_highIIW_after$Variance == min(abs(resultsmat_lowIPTW_highIIW_after$Variance)), ]
resultsmat_lowIPTW_highIIW_after[resultsmat_lowIPTW_highIIW_after$MSE == min(abs(resultsmat_lowIPTW_highIIW_after$MSE)), ]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ moderate IPTW and moderate IIW
#### n = 100 ###

set.seed(984654)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = -1
alpha0 = 0
alpha1 = 3.5
tau = 7
N = 1000

weighttrimresults_moderateIPTW_moderateIIW<- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
                                                                  alpha0, alpha1, tau, baseline = F)
saveRDS(weighttrimresults_moderateIPTW_moderateIIW, "weighttrimresults_moderateIPTW_moderateIIW.rds")
weighttrimresults_moderateIPTW_moderateIIW <- readRDS("weighttrimresults_moderateIPTW_moderateIIW.rds")


# view average weight distribution
weighttrimresults_moderateIPTW_moderateIIW$othermat_iptw
weighttrimresults_moderateIPTW_moderateIIW$othermat_iiw
weighttrimresults_moderateIPTW_moderateIIW$othermat_fiptiw


biasmat_moderateIPTW_moderateIIW <- weighttrimresults_moderateIPTW_moderateIIW$biasmat
colnames(biasmat_moderateIPTW_moderateIIW) <- c("Trimmed Before", "Trimmed After")
biasmat_moderateIPTW_moderateIIW<- melt(biasmat_moderateIPTW_moderateIIW)

varmat_moderateIPTW_moderateIIW<- weighttrimresults_moderateIPTW_moderateIIW$varmat
colnames(varmat_moderateIPTW_moderateIIW) <- c("Trimmed Before", "Trimmed After")
varmat_moderateIPTW_moderateIIW <- melt(varmat_moderateIPTW_moderateIIW)

msemat_moderateIPTW_moderateIIW <- weighttrimresults_moderateIPTW_moderateIIW$msemat
colnames(msemat_moderateIPTW_moderateIIW) <- c("Trimmed Before", "Trimmed After")
msemat_moderateIPTW_moderateIIW <- melt(msemat_moderateIPTW_moderateIIW)

resultsmat_moderateIPTW_moderateIIW <- cbind(biasmat_moderateIPTW_moderateIIW, varmat_moderateIPTW_moderateIIW$value, msemat_moderateIPTW_moderateIIW$value)
colnames(resultsmat_moderateIPTW_moderateIIW) <- c("Percentile", "Group", "Bias", "Variance", "MSE")

# bias plot
p1_moderateIPTW_moderateIIW <- ggplot(resultsmat_moderateIPTW_moderateIIW, aes(x = Percentile, y = Bias, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+ 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)")+
  ylab("Rel. Bias")+
  ylim(-0.3, 0.9)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))

### Variance plots
p2_moderateIPTW_moderateIIW <- ggplot(resultsmat_moderateIPTW_moderateIIW, aes(x = Percentile, y = Variance, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8) + 
  theme(legend.title=element_blank())+
  xlab("Threshold (Percentile)") +
  ylim(0.015, 0.03)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))

### MSE plots
p3_moderateIPTW_moderateIIW <- ggplot(resultsmat_moderateIPTW_moderateIIW, aes(x = Percentile, y = MSE, color = Group, shape = Group)) + 
  geom_line() + 
  geom_point(size = 0.8)+ 
  theme(legend.title=element_blank()) +
  xlab("Threshold (Percentile)")+
  ylim(0, 0.2)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))


combined_moderateIPTW_moderateIIW <- p1_moderateIPTW_moderateIIW + p2_moderateIPTW_moderateIIW + p3_moderateIPTW_moderateIIW & theme(legend.position = "bottom")
# combined_moderateIPTW_moderateIIW + 
#   plot_layout(guides = "collect") +
#   plot_annotation(title = "Results: Moderate IPTW, Moderate IIW")


### determine minimum bias, var, mse, etc
resultsmat_moderateIPTW_moderateIIW_first <- resultsmat_moderateIPTW_moderateIIW[resultsmat_moderateIPTW_moderateIIW$Group == "Trimmed Before", ]
resultsmat_moderateIPTW_moderateIIW_after <- resultsmat_moderateIPTW_moderateIIW[resultsmat_moderateIPTW_moderateIIW$Group == "Trimmed After", ]


#trimming first
resultsmat_moderateIPTW_moderateIIW_first[resultsmat_moderateIPTW_moderateIIW_first$Bias == min(abs(resultsmat_moderateIPTW_moderateIIW_first$Bias)), ]
resultsmat_moderateIPTW_moderateIIW_first[resultsmat_moderateIPTW_moderateIIW_first$Variance == min(abs(resultsmat_moderateIPTW_moderateIIW_first$Variance)), ]
resultsmat_moderateIPTW_moderateIIW_first[resultsmat_moderateIPTW_moderateIIW_first$MSE == min(abs(resultsmat_moderateIPTW_moderateIIW_first$MSE)), ]

#trimming after
resultsmat_moderateIPTW_moderateIIW_after[resultsmat_moderateIPTW_moderateIIW_after$Bias == min(abs(resultsmat_moderateIPTW_moderateIIW_after$Bias)), ]
resultsmat_moderateIPTW_moderateIIW_after[resultsmat_moderateIPTW_moderateIIW_after$Variance == min(abs(resultsmat_moderateIPTW_moderateIIW_after$Variance)), ]
resultsmat_moderateIPTW_moderateIIW_after[resultsmat_moderateIPTW_moderateIIW_after$MSE == min(abs(resultsmat_moderateIPTW_moderateIIW_after$MSE)), ]




#~~~~~~~~~~~ PLOT ~~~~~~~~~~~~~~

allplots <- wrap_elements(grid::textGrob('Low Trt \n Low Obs', gp = gpar( fontsize = 10))) + p1_lowIPTW_lowIIW + p2_lowIPTW_lowIIW  + p3_lowIPTW_lowIIW +
  wrap_elements(grid::textGrob('Moderate Trt \n Low Obs', gp = gpar( fontsize = 10))) + p1_moderateIPTW_lowIIW + p2_moderateIPTW_lowIIW + p3_moderateIPTW_lowIIW +
  wrap_elements(grid::textGrob('High Trt \n Low Obs', gp = gpar( fontsize = 10))) + p1_highIPTW_lowIIW + p2_highIPTW_lowIIW + p3_highIPTW_lowIIW +
  wrap_elements(grid::textGrob('Low Trt \n Moderate Obs', gp = gpar( fontsize = 10))) + p1_lowIPTW_moderateIIW + p2_lowIPTW_moderateIIW + p3_lowIPTW_moderateIIW +
  wrap_elements(grid::textGrob('Low Trt \n High Obs', gp = gpar( fontsize = 10))) + p1_lowIPTW_highIIW +p2_lowIPTW_highIIW +p3_lowIPTW_highIIW +
  wrap_elements(grid::textGrob('Moderate Trt \n Moderate Obs', gp = gpar( fontsize = 10))) + p1_moderateIPTW_moderateIIW + p2_moderateIPTW_moderateIIW + p3_moderateIPTW_moderateIIW & theme(legend.position = "bottom")


allplots + 
  plot_layout(guides = "collect", ncol = 4) +
  plot_annotation(title = "Simulation III Results")


