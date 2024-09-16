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
require(ggpubr)
require(dplyr)

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
ncutpts = 2


# weighttrimresults_lowIPTW_lowIIW<- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                             alpha0, alpha1, tau, ncutpts = 2)
# saveRDS(weighttrimresults_lowIPTW_lowIIW, "weighttrimresults_lowIPTW_lowIIW.rds")
weighttrimresults_lowIPTW_lowIIW <- readRDS("weighttrimresults_lowIPTW_lowIIW.rds")


# view average weight distribution
weighttrimresults_lowIPTW_lowIIW$othermat_iptw
weighttrimresults_lowIPTW_lowIIW$othermat_iiw
weighttrimresults_lowIPTW_lowIIW$othermat_fiptiw


### determine minimum bias, var, mse, etc

#bias closest to zero
bias_lowlow <- weighttrimresults_lowIPTW_lowIIW$biasmat
mse_lowlow <- weighttrimresults_lowIPTW_lowIIW$msemat
coverage_lowlow <- weighttrimresults_lowIPTW_lowIIW$coveragemat

bias_lowlow[which(abs(bias_lowlow$Bias) == min(abs(bias_lowlow$Bias))),] 
# 99th percentile, trimmed first
# bias = 0.000

bias_lowlow[bias_lowlow$Percentile==95, ]
# compared to bias = -0.086 for trimmed first, -0.068 for trimmed after


mse_lowlow[which(abs(mse_lowlow$MSE) == min(abs(mse_lowlow$MSE))),] 
#97-100th percentile (MSE 0.129) for both trimmed first and after

mse_lowlow[mse_lowlow$Percentile==95, ]
# MSE = 0.133 for trimmed first, 
# MSE = 0.133 for trimmed after

coverage_lowlow[coverage_lowlow$Percentile==95, ]
# 0.949 trimmed first,
# 0.949 trimmed after

coverage_lowlow[coverage_lowlow$Percentile==100, ]
# 0.953

lowlowPlots<- plotWeightTrimmingResults(biasmat = weighttrimresults_lowIPTW_lowIIW$biasmat,
                                                         varmat = weighttrimresults_lowIPTW_lowIIW$varmat,
                                                         msemat = weighttrimresults_lowIPTW_lowIIW$msemat,
                                                         coveragemat = weighttrimresults_lowIPTW_lowIIW$coveragemat)




####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Moderate IPTW, low IIW

#### n = 100 ###

set.seed(342)
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
ncutpts = 2

weightrimresults_moderateIPTW_lowIIW <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
                                           alpha0, alpha1, tau, ncutpts = 2)
saveRDS(weightrimresults_moderateIPTW_lowIIW, "weightrimresults_moderateIPTW_lowIIW.rds")
weightrimresults_moderateIPTW_lowIIW <- readRDS("weightrimresults_moderateIPTW_lowIIW.rds")


# view average weight distribution

weightrimresults_moderateIPTW_lowIIW$othermat_iptw
weightrimresults_moderateIPTW_lowIIW$othermat_iiw
weightrimresults_moderateIPTW_lowIIW$othermat_fiptiw


### determine minimum bias, var, mse, etc

#bias closest to zero
bias_modlow <- weightrimresults_moderateIPTW_lowIIW$biasmat
mse_modlow <- weightrimresults_moderateIPTW_lowIIW$msemat
coverage_modlow <- weightrimresults_moderateIPTW_lowIIW$coveragemat

bias_modlow[which(abs(bias_modlow$Bias) == min(abs(bias_modlow$Bias))),] 
# 99 percentile, trimmed after
# bias = 0.098

bias_modlow[bias_modlow$Percentile==95, ]
# compared to bias =  for trimmed first,  for trimmed after
# Trimmedfirst = 0.215
# Trimmed after = 0.165

mse_modlow[which(abs(mse_modlow$MSE) == min(abs(mse_modlow$MSE))),] 
#87th percentile (MSE = 0.305) for trimmed after

mse_modlow[mse_modlow$Percentile==95, ]
# MSE = 0.325 for trimmed first, 
# MSE = 0.329 for trimmed after

coverage_modlow[coverage_modlow$Percentile==95, ]
#  0.875 trimmed first,
#  0.882 trimmed after

coverage_modlow[coverage_modlow$Percentile==100, ]
# 0.874

modlowplots <- plotWeightTrimmingResults(biasmat = weightrimresults_moderateIPTW_lowIIW$biasmat,
                                                         varmat = weightrimresults_moderateIPTW_lowIIW$varmat,
                                                         msemat = weightrimresults_moderateIPTW_lowIIW$msemat,
                                                         coveragemat = weightrimresults_moderateIPTW_lowIIW$coveragemat)









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
                                          alpha0, alpha1, tau, ncutpts = 2)
saveRDS(weightrimresults_highIPTW_lowIIW, "weightrimresults_highIPTW_lowIIW.rds")
weightrimresults_highIPTW_lowIIW <- readRDS("weightrimresults_highIPTW_lowIIW.rds")

weightrimresults_highIPTW_lowIIW$othermat_iptw
weightrimresults_highIPTW_lowIIW$othermat_iiw
weightrimresults_highIPTW_lowIIW$othermat_fiptiw



### determine minimum bias, var, mse, etc

#bias closest to zero
bias_highlow <- weightrimresults_highIPTW_lowIIW$biasmat
mse_highlow <- weightrimresults_highIPTW_lowIIW$msemat
coverage_highlow <- weightrimresults_highIPTW_lowIIW$coveragemat

bias_highlow[which(abs(bias_highlow$Bias) == min(abs(bias_highlow$Bias))),] 
# 98th percentile, trimmed after
# bias = 0.083


bias_highlow[bias_highlow$Percentile==95, ]
# compared to bias =  0.127 for trimmed first, 0.103 for trimmed after


mse_highlow[which(abs(mse_highlow$MSE) == min(abs(mse_highlow$MSE))),] 
#57-66th percentile (MSE = 0.428) for both trimmed after

mse_highlow[mse_highlow$Percentile==95, ]
# MSE =  0.500 for trimmed first, 
# MSE =  for 0.524 trimmed after

coverage_highlow[coverage_highlow$Percentile==95, ]
# 0.869trimmed first,
# 0.853 trimmed after

coverage_highlow[coverage_highlow$Percentile==100, ]
#0.799

highlowplots <- plotWeightTrimmingResults(biasmat = weightrimresults_highIPTW_lowIIW$biasmat,
                                                         varmat = weightrimresults_highIPTW_lowIIW$varmat,
                                                         msemat = weightrimresults_highIPTW_lowIIW$msemat,
                                                         coveragemat = weightrimresults_highIPTW_lowIIW$coveragemat)





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
                                          alpha0, alpha1, tau, ncutpts = 2)
saveRDS(weighttrimresults_lowIPTW_moderateIIW, "weighttrimresults_lowIPTW_moderateIIW.rds")
weighttrimresults_lowIPTW_moderateIIW <- readRDS("weighttrimresults_lowIPTW_moderateIIW.rds")


# view average weight distribution
weighttrimresults_lowIPTW_moderateIIW$othermat_iptw
weighttrimresults_lowIPTW_moderateIIW$othermat_iiw
weighttrimresults_lowIPTW_moderateIIW$othermat_fiptiw




### determine minimum bias, var, mse, etc

#bias closest to zero
bias_lowmod <- weighttrimresults_lowIPTW_moderateIIW$biasmat
mse_lowmod <- weighttrimresults_lowIPTW_moderateIIW$msemat
coverage_lowmod <- weighttrimresults_lowIPTW_moderateIIW$coveragemat

bias_lowmod[which(abs(bias_lowmod$Bias) == min(abs(bias_lowmod$Bias))),] 
# 100th percentile, trimmed first
# bias = 0.075

bias_lowmod[bias_lowmod$Percentile==95, ]
# compared to bias = 0.141 for trimmed first, 0.141 for trimmed after


mse_lowmod[which(abs(mse_lowmod$MSE) == min(abs(mse_lowmod$MSE))),] 
#98th percentile (MSE 0.455) for both trimmed first and after

mse_lowmod[mse_lowmod$Percentile==95, ]
# MSE = 0.280 for trimmed first, 
# MSE = 0.362 for trimmed after

coverage_lowmod[coverage_lowmod$Percentile==95, ]
# 0.908 trimmed first,
# 0.866 trimmed after

coverage_lowmod[coverage_lowmod$Percentile==100, ]
# 0.925 trimmed first,


lowmodplots <- plotWeightTrimmingResults(biasmat = weighttrimresults_lowIPTW_moderateIIW$biasmat,
                                                         varmat = weighttrimresults_lowIPTW_moderateIIW$varmat,
                                                         msemat = weighttrimresults_lowIPTW_moderateIIW$msemat,
                                                         coveragemat = weighttrimresults_lowIPTW_moderateIIW$coveragemat)




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

# weighttrimresults_lowIPTW_highIIW<- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                       alpha0, alpha1, tau, ncutpts = 2)
# saveRDS(weighttrimresults_lowIPTW_highIIW, "weighttrimresults_lowIPTW_highIIW.rds")
weighttrimresults_lowIPTW_highIIW <- readRDS("weighttrimresults_lowIPTW_highIIW.rds")


# view average weight distribution
weighttrimresults_lowIPTW_highIIW$othermat_iptw
weighttrimresults_lowIPTW_highIIW$othermat_iiw
weighttrimresults_lowIPTW_highIIW$othermat_fiptiw




### determine minimum bias, var, mse, etc

#bias closest to zero
bias_lowhigh <- weighttrimresults_lowIPTW_highIIW$biasmat
mse_lowhigh <- weighttrimresults_lowIPTW_highIIW$msemat
coverage_lowhigh <- weighttrimresults_lowIPTW_highIIW$coveragemat

bias_lowhigh[which(abs(bias_lowhigh$Bias) == min(abs(bias_lowhigh$Bias))),] 
# 100th percentile, trimmed first
# bias = 0.095

bias_lowhigh[bias_lowhigh$Percentile==95, ]
# compared to bias = 0.233 for trimmed first, 0.401 for trimmed after


mse_lowhigh[which(abs(mse_lowhigh$MSE) == min(abs(mse_lowhigh$MSE))),] 
#96-98th percentile (MSE =0.385) for trimmed first

mse_lowhigh[mse_lowhigh$Percentile==95, ]
# MSE = 0.388 for trimmed first, 
# MSE = 0.461 for trimmed after

coverage_lowhigh[coverage_lowhigh$Percentile==95, ]
#  0.884 trimmed first,
# 0.835 trimmed after

lowhighplots <- plotWeightTrimmingResults(biasmat = weighttrimresults_lowIPTW_highIIW$biasmat,
                                                         varmat = weighttrimresults_lowIPTW_highIIW$varmat,
                                                         msemat = weighttrimresults_lowIPTW_highIIW$msemat,
                                                         coveragemat = weighttrimresults_lowIPTW_highIIW$coveragemat)










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
# 
# weighttrimresults_moderateIPTW_moderateIIW<- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, ncutpts = 2)
# saveRDS(weighttrimresults_moderateIPTW_moderateIIW, "weighttrimresults_moderateIPTW_moderateIIW.rds")
weighttrimresults_moderateIPTW_moderateIIW <- readRDS("weighttrimresults_moderateIPTW_moderateIIW.rds")


# view average weight distribution
weighttrimresults_moderateIPTW_moderateIIW$othermat_iptw
weighttrimresults_moderateIPTW_moderateIIW$othermat_iiw
weighttrimresults_moderateIPTW_moderateIIW$othermat_fiptiw




### determine minimum bias, var, mse, etc

#bias closest to zero
bias_modmod <- weighttrimresults_moderateIPTW_moderateIIW$biasmat
mse_modmod <- weighttrimresults_moderateIPTW_moderateIIW$msemat
coverage_modmod <- weighttrimresults_moderateIPTW_moderateIIW$coveragemat

bias_modmod[which(abs(bias_modmod$Bias) == min(abs(bias_modmod$Bias))),] 
# 99th percentile, trimmed first
# bias = 

bias_modmod[bias_modmod$Percentile==95, ]
# compared to bias =  for trimmed first,  for trimmed after


mse_modmod[which(abs(mse_modmod$MSE) == min(abs(mse_modmod$MSE))),] 
#98th percentile (MSE ) for both trimmed first and after

mse_modmod[mse_modmod$Percentile==95, ]
# MSE =  for trimmed first, 
# MSE =  for trimmed after

coverage_modmod[coverage_modmod$Percentile==95, ]
#  trimmed first,
#  trimmed after

modmodplots <- plotWeightTrimmingResults(biasmat = weighttrimresults_modIPTW_modIIW$biasmat,
                                                         varmat = weighttrimresults_modIPTW_modIIW$varmat,
                                                         msemat = weighttrimresults_modIPTW_modIIW$msemat,
                                                         coveragemat = weighttrimresults_modIPTW_modIIW$coveragemat)








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


