
source("FIPTIW_timevarW_functions.R")

require(parallel)
# parallelization stuff
# ncores <- detectCores() - 3
nclusters <- makeCluster(45)

require(knitr)
require(kableExtra)
require(ggplot2)
require(ggpubr)
require(knitr)
require(kableExtra)
require(reshape2)
require(latex2exp)
require(patchwork)
require(grid)
require(dplyr)


#~~~~~~~~~~ FIPTIW simulation, assumptions satisfied comparing spline vs offset (APPENDIX).

### BINARY

n = 100
beta1 = 0.5
beta2 = 0.5 
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3)
gamma3vec = c(0, 0.6)
alpha0 = -1
alpha1vec = c(0, 1)
tau = 7
N = 500


set.seed(8659)
# resultsFIPTIW_n100_binary <- simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec,
#                                                alpha0, alpha1vec, tau, outcome = "binary",
#                                                censinform = F, eta1 = NULL, eta2 = NULL, eta3= NULL,
#                                                usesplines = F, ncutpts = NULL, outputfulldatalist = T,
#                                                inParallel = T, nclusters = nclusters)
# 
# 
# saveRDS(resultsFIPTIW_n100_binary, "resultsFIPTIW_n100_binary_satisfied.rds")
resultsFIPTIW_n100_binary <- readRDS("resultsFIPTIW_n100_binary_satisfied.rds")



set.seed(100)
# resultsFIPTIW_n100_binary_spline <- simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec,
#                                                alpha0, alpha1vec, tau, outcome = "binary",
#                                                censinform = F, eta1 = NULL, eta2 = NULL, eta3= NULL,
#                                                usesplines = T, ncutpts = 2, outputfulldatalist = T,
#                                                inParallel = T, nclusters = nclusters)
# saveRDS(resultsFIPTIW_n100_binary_spline, "resultsFIPTIW_n100_binary_satisfied_spline.rds")
resultsFIPTIW_n100_binary_spline <- readRDS("resultsFIPTIW_n100_binary_satisfied_spline.rds")


results_n100_binary_both <- cbind(resultsFIPTIW_n100_binary$resultsmat, resultsFIPTIW_n100_binary_spline$resultsmat[,5:8])

kable(results_n100_binary_both, booktabs = T, digits = 3 ,format = "latex") %>%
  add_header_above(c(" " = 4, "Parametric Offset Used" = 4, "Cubic Splines Used" = 4))


#### CONTINUOUS



set.seed(234)
# resultsFIPTIW_n100_normal <-  simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
#                                                  alpha0, alpha1vec, tau, outcome = "continuous", 
#                                                  censinform = F, eta1 = NULL, eta2 = NULL, eta3= NULL,
#                                                  usesplines = F, ncutpts = NULL, outputfulldatalist = T,
#                                                  inParallel = T, nclusters = nclusters)
# 
# 
# 
# saveRDS(resultsFIPTIW_n100_normal, "resultsFIPTIW_n100_continuous_satisfied.rds")
resultsFIPTIW_n100_normal <- readRDS("resultsFIPTIW_n100_continuous_satisfied.rds")


set.seed(2348)
# resultsFIPTIW_n100_normal_spline <-  simulateALLFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
#                                             alpha0, alpha1vec, tau, outcome = "continuous", 
#                                             censinform = F, eta1 = NULL, eta2 = NULL, eta3= NULL,
#                                             usesplines = T, ncutpts = 2, outputfulldatalist = T,
#                                             inParallel = T, nclusters = nclusters)
# saveRDS(resultsFIPTIW_n100_normal_spline, "resultsFIPTIW_n100_continuous_satisfied_spline.rds")
resultsFIPTIW_n100_normal_spline <- readRDS("resultsFIPTIW_n100_continuous_satisfied_spline.rds")


results_n100_normal_both <- cbind(resultsFIPTIW_n100_normal$resultsmat, resultsFIPTIW_n100_normal_spline$resultsmat[,5:8])

kable(results_n100_normal_both, booktabs = T, digits = 3 ,format = "latex") %>%
  add_header_above(c(" " = 3, "Parametric Offset Used" = 5, "Cubic Splines Used" = 5))










###############################################################################
################ SIMULATION I #################################################
########## NONINFORMATIVE CENSORING ###########################################
###############################################################################

#### ~~~~~~~~~~~~ Binary Outcome (n = 50, 100, 500) ~~~~~~~~~~~~~~~~~~~ #



beta1 = 0.5
beta2 = 0.5
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = -1
alpha1 = 1
tau = 7
N = 500
censinform = T
eta1vec = c(0)
eta2vec = c(0, 0.2, 0.5)
eta3vec = c(0, 0.4, 0.6)



### binary n = 100, with splines

set.seed(23487)
# results_censoring_n100_binary <- simulateALLFIPTICW_CENS(N, n = 100, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                    alpha0, alpha1, tau, outcome = "binary",
#                                    censinform = T, eta1vec, eta2vec, eta3vec,
#                                    usesplines = T, ncutpts =  2, outputfulldatalist = T,
#                                    inParallel = T, nclusters = nclusters)
# 
# 
# saveRDS(results_censoring_n100_binary, "results_censoring_n100_binary.rds")
results_censoring_n100_binary <- readRDS("results_censoring_n100_binary.rds")

plot_censoring_n100_binary <- makeBoxPlotsATE(results_censoring_n100_binary$fulldatalist, eta2vec)
annotate_figure(plot_censoring_n100_binary, top = text_grob("Combined Results for n = 100,  Binary Outcome", face = "bold", size = 20))


coveragemat_n100_binary <- results_censoring_n100_binary$resultsmat[, c(2,3,7,8,12,13,17,18,22,23, 27,28)]
kable(coveragemat_n100_binary, booktabs = T, 
      digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2, "FIPTICW" = 2))





### binary n = 50, with splines

set.seed(7354)
# results_censoring_n50_binary <- simulateALLFIPTICW_CENS(N, n = 50, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                                   alpha0, alpha1, tau, outcome = "binary",
#                                                   censinform = T, eta1vec, eta2vec, eta3vec,
#                                                   usesplines = T, ncutpts =  2, outputfulldatalist = T,
#                                                   inParallel = T, nclusters = nclusters)
# 
# 
# saveRDS(results_censoring_n50_binary, "results_censoring_n50_binary.rds")
results_censoring_n50_binary <- readRDS("results_censoring_n50_binary.rds")

plot_censoring_n50_binary <- makeBoxPlotsATE(results_censoring_n50_binary$fulldatalist, eta2vec)
annotate_figure(plot_censoring_n50_binary, top = text_grob("Combined Results for n = 50, Binary Outcome", face = "bold", size = 20))


coveragemat_n50_binary <- results_censoring_n50_binary$resultsmat[, c(2,3,7,8,12,13,17,18,22,23, 27,28)]
kable(coveragemat_n50_binary, booktabs = T, 
      digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2, "FIPTICW" = 2))



### binary n = 500, with splines

set.seed(643)
results_censoring_n500_binary <- simulateALLFIPTICW_CENS(N, n = 500, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
                                                 alpha0, alpha1, tau, outcome = "binary",
                                                 censinform = T, eta1vec, eta2vec, eta3vec,
                                                 usesplines = T, ncutpts =  2, outputfulldatalist = T,
                                                 inParallel = T, nclusters = nclusters)


# saveRDS(results_censoring_n500_binary, "results_censoring_n500_binary.rds")
# results_censoring_n500_binary <- readRDS("results_censoring_n500_binary.rds")

plot_censoring_n500_binary <- makeBoxPlotsATE(results_censoring_n500_binary$fulldatalist, eta2vec)
annotate_figure(plot_censoring_n500_binary, top = text_grob("Combined Results for n = 500, Binary Outcome", face = "bold", size = 20))


coveragemat_n500_binary <- results_censoring_n500_binary$resultsmat[, c(2,3,7,8,12,13,17,18,22,23, 27,28)]
kable(coveragemat_n500_binary, booktabs = T, 
      digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2, "FIPTICW" = 2))





#### ~~~~~~~~~~~~ Continuous Outcome (n = 100) ~~~~~~~~~~~~~~~~~~~ #



beta1 = 0.5
beta2 = 2
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = -1
alpha1 = 1
tau = 7
N = 500
censinform = T
eta1vec = c(0)
eta2vec = c(0, 0.2, 0.5)
eta3vec = c(0, 0.4, 0.6)



### continuous n = 100, with splines

set.seed(2345)
# results_censoring_n100_continuous <- simulateALLFIPTICW_CENS(N, n = 100, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                                   alpha0, alpha1, tau, outcome = "continuous",
#                                                   censinform = T, eta1vec, eta2vec, eta3vec,
#                                                   usesplines = T, ncutpts =  2, outputfulldatalist = T,
#                                                   inParallel = T, nclusters = nclusters)
# 
# 
# saveRDS(results_censoring_n100_continuous, "results_censoring_n100_continuous.rds")
results_censoring_n100_continuous <- readRDS("results_censoring_n100_continuous.rds")

plot_censoring_n100_continuous <- makeBoxPlotsATE(results_censoring_n100_continuous$fulldatalist, eta2vec)
annotate_figure(plot_censoring_n100_continuous, top = text_grob("Combined Results for n = 100, Normally Distributed Outcome", face = "bold", size = 20))


coveragemat_n100_continuous <- results_censoring_n100_continuous$resultsmat[, c(2,3,7,8,12,13,17,18,22,23, 27,28)]
kable(coveragemat_n100_continuous, booktabs = T, 
      digits = 3,format = "latex") %>%
  add_header_above(c(" " = 3, "Unweighted" = 2, "IIW" = 2, "IPTW" = 2, "FIPTIW" = 2, "FIPTICW" = 2))










###############################################################################
################ SIMULATION II ################################################
############## VARIABLE SELECTION #############################################
###############################################################################

#### ~~~~~~~~~~~~ Binary Outcome (n = 50, 100, 500) ~~~~~~~~~~~~~~~~~~~ #


beta1 = 0.5
beta2vec = c(0, 0.5)
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3, 0.6)
gamma3 = 0.6
alpha0 = 0
alpha1 = 0
tau = 7
N = 500

#### n = 100, binary, splines

set.seed(54641230)
# varsel_n100_binary <- simulateALLVarSel(N, n = 100, beta1, beta2vec, beta3, gamma1, gamma2vec, gamma3,
#                   alpha0, alpha1, tau, outcome = "binary",
#                   usesplines = T, ncutpts = 2,
#                   inParallel = T, nclusters = nclusters)
# 
# 
# saveRDS(varsel_n100_binary , "varsel_n100_binary.rds")
varsel_n100_binary<- readRDS("varsel_n100_binary.rds")

kable(varsel_n100_binary $resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
  add_header_above(c(" " = 3, "Variables used to estimate intensity" = 8))




### n = 50, binary, splines
set.seed(39487)
# varsel_n50_binary <- simulateALLVarSel(N, n = 50, beta1, beta2vec, beta3, gamma1, gamma2vec, gamma3, 
#                                         alpha0, alpha1, tau, outcome = "binary",
#                                         usesplines = T, ncutpts = 2,
#                                         inParallel = T, nclusters = nclusters)
# 
# 
# saveRDS(varsel_n50_binary , "varsel_n50_binary.rds")
varsel_n50_binary<- readRDS("varsel_n50_binary.rds")

kable(varsel_n50_binary $resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
  add_header_above(c(" " = 3, "Variables used to estimate intensity" = 8))




### n = 500, binary, splines
set.seed(3958)
# varsel_n500_binary <- simulateALLVarSel(N, n = 500, beta1, beta2vec, beta3, gamma1, gamma2vec, gamma3, 
#                                        alpha0, alpha1, tau, outcome = "binary",
#                                        usesplines = T, ncutpts = 2,
#                                        inParallel = T, nclusters = nclusters)
# 
# saveRDS(varsel_n500_binary , "varsel_n500_binary.rds")
varsel_n500_binary<- readRDS("varsel_n500_binary.rds")

kable(varsel_n500_binary $resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
  add_header_above(c(" " = 3, "Variables used to estimate intensity" = 8))


#### ~~~~~~~~~~~~ Normal Outcome (n = 100) ~~~~~~~~~~~~~~~~~~~ #


beta1 = 0.5
beta2vec = c(0, 2)
beta3 = 1
gamma1 = 0.5
gamma2vec = c(0, 0.3, 0.6)
gamma3 = 0.6
alpha0 = 0
alpha1 = 0
tau = 7
N = 500

#### n = 100, continuous, splines
set.seed(54641230)
# varsel_n100_continuous <- simulateALLVarSel(N, n = 100, beta1, beta2vec, beta3, gamma1, gamma2vec, gamma3, 
#                                         alpha0, alpha1, tau, outcome = "continuous",
#                                         usesplines = T, ncutpts = 2,
#                                         inParallel = T, nclusters = nclusters)
# 
# 
# saveRDS(varsel_n100_continuous , "varsel_n100_continuous.rds")
varsel_n100_continuous<- readRDS("varsel_n100_continuous.rds")

kable(varsel_n100_continuous$resultsmat_beta01, booktabs = T, digits = 2,format = "latex") %>%
  add_header_above(c(" " = 3, "Variables used to estimate intensity" = 8))






###############################################################################
################ SIMULATION III ###############################################
############### WEIGHT TRIMMING  ##############################################
###############################################################################

#### ~~~~~~~~~~~~ Binary Outcome w Splines (n = 100) ~~~~~~~~~~~~~~~~~~~ #

### LOW-LOW


n = 100
beta1 = 0.5
beta2 = 0.5
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 0.5
tau = 7
N = 500

set.seed(13546)

# weighttrimresults_lowIPTW_lowIIW <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                           alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_lowIPTW_lowIIW, "weighttrimresults_lowIPTW_lowIIW.rds")
weighttrimresults_lowIPTW_lowIIW <- readRDS("weighttrimresults_lowIPTW_lowIIW.rds")
compareMinimums(weighttrimresults_lowIPTW_lowIIW)

plots_lowIPTW_lowIIW <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_lowIIW)


### MODERATE-LOW


gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 3.5

set.seed(243)
# weighttrimresults_modIPTW_lowIIW <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                                                   alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_modIPTW_lowIIW, "weighttrimresults_modIPTW_lowIIW.rds")
weighttrimresults_modIPTW_lowIIW <- readRDS("weighttrimresults_modIPTW_lowIIW.rds")
compareMinimums(weighttrimresults_modIPTW_lowIIW)

plots_modIPTW_lowIIW <- plotWeightTrimmingResults(weighttrimresults_modIPTW_lowIIW)



### HIGH-LOW
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 5.5

set.seed(23489)
# weighttrimresults_highIPTW_lowIIW <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                                                   alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_highIPTW_lowIIW, "weighttrimresults_highIPTW_lowIIW.rds")
weighttrimresults_highIPTW_lowIIW <- readRDS("weighttrimresults_highIPTW_lowIIW.rds")
compareMinimums(weighttrimresults_highIPTW_lowIIW)


plots_highIPTW_lowIIW <- plotWeightTrimmingResults(weighttrimresults_highIPTW_lowIIW)


### LOW-MOD

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.40
alpha0 = 0
alpha1 = 0.5

set.seed(3874)
# weighttrimresults_lowIPTW_modIIW <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                    alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_lowIPTW_modIIW, "weighttrimresults_lowIPTW_modIIW.rds")
weighttrimresults_lowIPTW_modIIW <- readRDS("weighttrimresults_lowIPTW_modIIW.rds")
compareMinimums(weighttrimresults_lowIPTW_modIIW)

plots_lowIPTW_modIIW <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_modIIW)


### LoW-HIGH

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.6
alpha0 = 0
alpha1 = 0.5

set.seed(23487)
# weighttrimresults_lowIPTW_highIIW <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_lowIPTW_highIIW, "weighttrimresults_lowIPTW_highIIW.rds")
weighttrimresults_lowIPTW_highIIW <- readRDS("weighttrimresults_lowIPTW_highIIW.rds")
compareMinimums(weighttrimresults_lowIPTW_highIIW)

plots_lowIPTW_highIIW <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_highIIW)




### MOD-MOD

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.4
alpha0 = 0
alpha1 = 3.5

set.seed(99)
# weighttrimresults_modIPTW_modIIW <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                                                    alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_modIPTW_modIIW, "weighttrimresults_modIPTW_modIIW.rds")
weighttrimresults_modIPTW_modIIW <- readRDS("weighttrimresults_modIPTW_modIIW.rds")
compareMinimums(weighttrimresults_modIPTW_modIIW)

plots_modIPTWmodIIW <- plotWeightTrimmingResults(weighttrimresults_modIPTW_modIIW)



######### Binary Outcome No Splines ######


### LOW-LOW


n = 100
beta1 = 0.5
beta2 = 0.5
beta3 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 0.5
tau = 7
N = 500

set.seed(13546)

# weighttrimresults_lowIPTW_lowIIW_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                           alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_lowIPTW_lowIIW_nospline, "weighttrimresults_lowIPTW_lowIIW_nospline.rds")
weighttrimresults_lowIPTW_lowIIW_nospline <- readRDS("weighttrimresults_lowIPTW_lowIIW_nospline.rds")
compareMinimums(weighttrimresults_lowIPTW_lowIIW_nospline)

plots_lowIPTW_lowIIW_nospline <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_lowIIW_nospline)


### MODERATE-LOW


gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 3.5

set.seed(243)
# weighttrimresults_modIPTW_lowIIW_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_modIPTW_lowIIW_nospline, "weighttrimresults_modIPTW_lowIIW_nospline.rds")
weighttrimresults_modIPTW_lowIIW_nospline <- readRDS("weighttrimresults_modIPTW_lowIIW_nospline.rds")
compareMinimums(weighttrimresults_modIPTW_lowIIW_nospline)

plots_modIPTW_lowIIW_nospline <- plotWeightTrimmingResults(weighttrimresults_modIPTW_lowIIW_nospline)



### HIGH-LOW
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 5.5

set.seed(23489)
# weighttrimresults_highIPTW_lowIIW_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_highIPTW_lowIIW_nospline, "weighttrimresults_highIPTW_lowIIW_nospline.rds")
weighttrimresults_highIPTW_lowIIW_nospline <- readRDS("weighttrimresults_highIPTW_lowIIW_nospline.rds")
compareMinimums(weighttrimresults_highIPTW_lowIIW_nospline)


plots_highIPTW_lowIIW_nospline <- plotWeightTrimmingResults(weighttrimresults_highIPTW_lowIIW_nospline)


### LOW-MOD

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.40
alpha0 = 0
alpha1 = 0.5

set.seed(3874)
# weighttrimresults_lowIPTW_modIIW_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_lowIPTW_modIIW_nospline, "weighttrimresults_lowIPTW_modIIW_nospline.rds")
weighttrimresults_lowIPTW_modIIW_nospline <- readRDS("weighttrimresults_lowIPTW_modIIW_nospline.rds")
compareMinimums(weighttrimresults_lowIPTW_modIIW_nospline)

plots_lowIPTW_modIIW_nospline <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_modIIW_nospline)


### LoW-HIGH

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.6
alpha0 = 0
alpha1 = 0.5

set.seed(23487)
# weighttrimresults_lowIPTW_highIIW_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_lowIPTW_highIIW_nospline, "weighttrimresults_lowIPTW_highIIW_nospline.rds")
weighttrimresults_lowIPTW_highIIW_nospline <- readRDS("weighttrimresults_lowIPTW_highIIW_nospline.rds")
compareMinimums(weighttrimresults_lowIPTW_highIIW_nospline)

plots_lowIPTW_highIIW_nospline <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_highIIW_nospline)




### MOD-MOD

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.4
alpha0 = 0
alpha1 = 3.5

set.seed(99)
# weighttrimresults_modIPTW_modIIW_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                                                   alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "binary")
# saveRDS(weighttrimresults_modIPTW_modIIW_nospline, "weighttrimresults_modIPTW_modIIW_nospline.rds")
weighttrimresults_modIPTW_modIIW_nospline <- readRDS("weighttrimresults_modIPTW_modIIW_nospline.rds")
compareMinimums(weighttrimresults_modIPTW_modIIW_nospline)

plots_modIPTWmodIIW_nospline <- plotWeightTrimmingResults(weighttrimresults_modIPTW_modIIW_nospline)









#### ~~~~~~~~~~~~ Continuous Outcome with splines (n = 100) ~~~~~~~~~~~~~~~~~~~ #

### LOW-LOW


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
N = 500

set.seed(235)

# weighttrimresults_lowIPTW_lowIIW_continuous_continuous <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                           alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_lowIPTW_lowIIW_continuous_continuous, "weighttrimresults_lowIPTW_lowIIW_continuous_continuous.rds")
weighttrimresults_lowIPTW_lowIIW_continuous_continuous <- readRDS("weighttrimresults_lowIPTW_lowIIW_continuous_continuous.rds")
compareMinimums(weighttrimresults_lowIPTW_lowIIW_continuous_continuous)

plots_lowIPTW_lowIIW_continuous_continuous <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_lowIIW_continuous_continuous)


### MODERATE-LOW


gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 3.5

set.seed(764)
# weighttrimresults_modIPTW_lowIIW_continuous <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_modIPTW_lowIIW_continuous, "weighttrimresults_modIPTW_lowIIW_continuous.rds")
weighttrimresults_modIPTW_lowIIW_continuous <- readRDS("weighttrimresults_modIPTW_lowIIW_continuous.rds")
compareMinimums(weighttrimresults_modIPTW_lowIIW_continuous)

plots_modIPTW_lowIIW_continuous <- plotWeightTrimmingResults(weighttrimresults_modIPTW_lowIIW_continuous)



### HIGH-LOW
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 5.5

set.seed(1353)
# weighttrimresults_highIPTW_lowIIW_continuous <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_highIPTW_lowIIW_continuous, "weighttrimresults_highIPTW_lowIIW_continuous.rds")
weighttrimresults_highIPTW_lowIIW_continuous <- readRDS("weighttrimresults_highIPTW_lowIIW_continuous.rds")
compareMinimums(weighttrimresults_highIPTW_lowIIW_continuous)


plots_highIPTW_lowIIW_continuous <- plotWeightTrimmingResults(weighttrimresults_highIPTW_lowIIW_continuous)


### LOW-MOD

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.40
alpha0 = 0
alpha1 = 0.5

set.seed(3249)
# weighttrimresults_lowIPTW_modIIW_continuous <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_lowIPTW_modIIW_continuous, "weighttrimresults_lowIPTW_modIIW_continuous.rds")
weighttrimresults_lowIPTW_modIIW_continuous <- readRDS("weighttrimresults_lowIPTW_modIIW_continuous.rds")
compareMinimums(weighttrimresults_lowIPTW_modIIW_continuous)

plots_lowIPTW_modIIW_continuous <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_modIIW_continuous)


### LoW-HIGH

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.6
alpha0 = 0
alpha1 = 0.5

set.seed(900)
# weighttrimresults_lowIPTW_highIIW_continuous <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_lowIPTW_highIIW_continuous, "weighttrimresults_lowIPTW_highIIW_continuous.rds")
weighttrimresults_lowIPTW_highIIW_continuous <- readRDS("weighttrimresults_lowIPTW_highIIW_continuous.rds")
compareMinimums(weighttrimresults_lowIPTW_highIIW_continuous)

plots_lowIPTW_highIIW_continuous <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_highIIW_continuous)




### MOD-MOD

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.4
alpha0 = 0
alpha1 = 3.5

set.seed(188)
# weighttrimresults_modIPTW_modIIW_continuous <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                                                   alpha0, alpha1, tau, usesplines = T, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_modIPTW_modIIW_continuous, "weighttrimresults_modIPTW_modIIW_continuous.rds")
weighttrimresults_modIPTW_modIIW_continuous <- readRDS("weighttrimresults_modIPTW_modIIW_continuous.rds")
compareMinimums(weighttrimresults_modIPTW_modIIW_continuous)

plots_modIPTWmodIIW_continuous <- plotWeightTrimmingResults(weighttrimresults_modIPTW_modIIW_continuous)



#### ~~~~~~~~~~~~ Continuous Outcome with splines (n = 100) ~~~~~~~~~~~~~~~~~~~ #

### LOW-LOW


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
N = 500

set.seed(235)

# weighttrimresults_lowIPTW_lowIIW_continuous_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                           alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_lowIPTW_lowIIW_continuous_nospline, "weighttrimresults_lowIPTW_lowIIW_continuous_nospline.rds")
weighttrimresults_lowIPTW_lowIIW_continuous_nospline <- readRDS("weighttrimresults_lowIPTW_lowIIW_continuous_nospline.rds")
compareMinimums(weighttrimresults_lowIPTW_lowIIW_continuous_nospline)

plots_lowIPTW_lowIIW_continuous_nospline <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_lowIIW_continuous_nospline)


### MODERATE-LOW


gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 3.5

set.seed(764)
# weighttrimresults_modIPTW_lowIIW_continuous_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_modIPTW_lowIIW_continuous_nospline, "weighttrimresults_modIPTW_lowIIW_continuous_nospline.rds")
weighttrimresults_modIPTW_lowIIW_continuous_nospline <- readRDS("weighttrimresults_modIPTW_lowIIW_continuous_nospline.rds")
compareMinimums(weighttrimresults_modIPTW_lowIIW_continuous_nospline)

plots_modIPTW_lowIIW_continuous_nospline <- plotWeightTrimmingResults(weighttrimresults_modIPTW_lowIIW_continuous_nospline)



### HIGH-LOW
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.6
alpha0 = 0
alpha1 = 5.5

set.seed(1353)
# weighttrimresults_highIPTW_lowIIW_continuous_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_highIPTW_lowIIW_continuous_nospline, "weighttrimresults_highIPTW_lowIIW_continuous_nospline.rds")
weighttrimresults_highIPTW_lowIIW_continuous_nospline <- readRDS("weighttrimresults_highIPTW_lowIIW_continuous_nospline.rds")
compareMinimums(weighttrimresults_highIPTW_lowIIW_continuous_nospline)


plots_highIPTW_lowIIW_continuous_nospline <- plotWeightTrimmingResults(weighttrimresults_highIPTW_lowIIW_continuous_nospline)


### LOW-MOD

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.40
alpha0 = 0
alpha1 = 0.5

set.seed(3249)
# weighttrimresults_lowIPTW_modIIW_continuous_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_lowIPTW_modIIW_continuous_nospline, "weighttrimresults_lowIPTW_modIIW_continuous_nospline.rds")
weighttrimresults_lowIPTW_modIIW_continuous_nospline <- readRDS("weighttrimresults_lowIPTW_modIIW_continuous_nospline.rds")
compareMinimums(weighttrimresults_lowIPTW_modIIW_continuous_nospline)

plots_lowIPTW_modIIW_continuous_nospline <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_modIIW_continuous_nospline)


### LoW-HIGH

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.6
alpha0 = 0
alpha1 = 0.5

set.seed(900)
# weighttrimresults_lowIPTW_highIIW_continuous_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_lowIPTW_highIIW_continuous_nospline, "weighttrimresults_lowIPTW_highIIW_continuous_nospline.rds")
weighttrimresults_lowIPTW_highIIW_continuous_nospline <- readRDS("weighttrimresults_lowIPTW_highIIW_continuous_nospline.rds")
compareMinimums(weighttrimresults_lowIPTW_highIIW_continuous_nospline)

plots_lowIPTW_highIIW_continuous_nospline <- plotWeightTrimmingResults(weighttrimresults_lowIPTW_highIIW_continuous_nospline)




### MOD-MOD

gamma1 = 0.5
gamma2 = 0.3
gamma3 = -0.4
alpha0 = 0
alpha1 = 3.5

set.seed(188)
# weighttrimresults_modIPTW_modIIW_continuous_nospline <- simulateResultsWeightTrimming(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
#                                                                   alpha0, alpha1, tau, usesplines = F, ncutpts = 2, outcome = "continuous")
# saveRDS(weighttrimresults_modIPTW_modIIW_continuous_nospline, "weighttrimresults_modIPTW_modIIW_continuous_nospline.rds")
weighttrimresults_modIPTW_modIIW_continuous_nospline <- readRDS("weighttrimresults_modIPTW_modIIW_continuous_nospline.rds")
compareMinimums(weighttrimresults_modIPTW_modIIW_continuous_nospline)

plots_modIPTWmodIIW_continuous_nospline <- plotWeightTrimmingResults(weighttrimresults_modIPTW_modIIW_continuous_nospline)





