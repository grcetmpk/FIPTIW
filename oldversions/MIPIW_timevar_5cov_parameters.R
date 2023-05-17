# The purpose of this sensitivity analysis is to evaluate the 
# robustness of the conditional independence assumption


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
# Censoring: C_i ~ Unif(tau/2, tau)
# tau: maximum follow up time, independent censoring.
# N: number of simulations
#\
#
#
# We fit various weighting models:
# - Both models are correctly specified
# - Probability of treatment model is missing $W_1$
# - Probability of treatment model is missing $W_2$
# - Intensity model is missing $Z_1$
# - Intensity model is missing $Z_2(t)$
# - Both models are missing $W_1 or W_1(t)$
# - Both models are missing $W_2$
# - Probability of treatment model is missing $W_1$ and intensity model is missing $Z_1$
# - Probability of treatment model is missing $W_2$ and intensity model is missing $Z_2(t)$



library(knitr)
library(latex2exp)
library(kableExtra)

setwd("C:\\Users\\grace\\Documents\\WaterlooPHD\\Research\\IIW_Simulation\\MIPIW\\Simulation1")
source("MIPIW_timevar_5cov_functions.R")


# Violating independence between all three processes (independent censoring)
set.seed(33525)
n = 100
beta1 = 0.5
beta2 = 2
beta3 = 1
beta4 = 1
beta5 = 1
gamma1 = 0.5
gamma2 = 0.3
gamma3 = 0.3
gamma4 = 0.2
gamma5 = 0.3
alpha0 = -0.2
alpha1 = 0.3
alpha2 = 0.3
tau = 7
N = 1000


simulateResultsMIPI_timevar_5cov_full(N, n, beta1, beta2, beta3, beta4, beta5,
                                                    gamma1, gamma2, gamma3, gamma4, gamma5,
                                                    alpha0, alpha1, alpha2, tau)

