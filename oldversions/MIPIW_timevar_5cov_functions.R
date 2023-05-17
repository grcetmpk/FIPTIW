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
# Observation intensity: lambda_i(t) = eta_i*sqrt(t)/2*exp(gamma1G_i + gamma2W_i1(t) + gamma3W_i2 + gamma4Z_i1 + gamma5Z_i2(t))
# Treatment assignment: P(G_i = 1) = expit(alpha0 + alpha1W_i1 + alpha2W_i2) 
#
# Censoring: C_i ~ Unif(tau/2, tau)
#
# tau: maximum follow up time, independent censoring.
# N: number of simulations

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


require(survival)
require(geepack)
require(nleqslv)
require(knitr)
require(kableExtra)

expit <- function(x){ return(exp(x)/(1+exp(x)))}


#~~~~~~~~~~~~~ Data Generation ~~~~~~~~~~~~~~~~#

gendataMIPI_timevar_5cov_full <- function(n, beta1, beta2, beta3, beta4, beta5,
                                          gamma1, gamma2, gamma3, gamma4, gamma5,
                                          alpha0, alpha1, alpha2, tau){
  
  ## Generates a single data set given the following:
  
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
  # Observation intensity: lambda_i(t) = eta_i*sqrt(t)/2*exp(gamma1G_i + gamma2W_i1(t) + gamma3W_i2 + gamma4Z_i1 + gamma5Z_i2(t))
  # Treatment assignment: P(G_i = 1) = expit(alpha0 + alpha1W_i1 + alpha2W_i2) 
  #
  # Censoring: C_i ~ U(tau/2, tau)
  #
  # tau: maximum follow up time, independent censoring.
  # N: number of simulations
  
  mu_Z_G0 <- 2
  var_Z_G0 <- 1 
  mu_Z_G1 <- 0
  var_Z_G1 <- 0.5
  
  var_phi <- 0.25
  var_epsilon <- 1
  

  
  getObsTimesViaThinning <- function(lambda, tmax){
    # Thinning method for generating NHPP observation times
    #
    # lambda: rate function for the NHPP 
    # tmax: maximum follow up time (tau)
    
    # initialize 
    nobs <- 0
    m <- 0
    t0 <- 0
    s0 <- 0
    lambdabar <- lambda(tmax) #supremum, intensity is bounded by this HPP
    
    times <- c(t0)
    s <- s0
    
    i = m + 1 # R starts at 1 not 0, minor adjustment
    
    while(s < tmax){
      u <- runif(1)
      w <- -log(u)/lambdabar # makes w ~ exp(lambdabar)
      s <-  s + w # s_m+1 = s + w
      D <- runif(1)
      
      if(D <= lambda(s)/lambdabar){
        times <- c(times, s)
        nobs <- nobs + 1
      }
      
      i <- i + 1
      m <- m + 1
    }
    
    times <- times[-1]
    
    if(times[nobs] <= tmax){
      return(times)
    } else{
      return(times[-nobs])
    }
  }
  
  
  
  #generate subjects one by one and then combine
  id = 1
  simdatafull <- data.frame(matrix(NA, ncol = 14))
  colnames(simdatafull) <- c("id", "time", "G", "W1", "W2", "Z1", "Z2", "censortime", "censored",
                             "cexp_W1t_G", "cexp_W2_G", "cexp_Z1_G" ,"cexp_Z2t_G",  "y")
  numcensored = 0
  
  while(id <= n){
    
    ## Generate Covariates
    
    # generate treatment confounder w ~ U(0,1)
    W1 <- runif(1, 0, 1)
    
    # generate second treatment confounder
    W2 <- rnorm(1,0,1)
    
    # generate treatment assignment
    prG <- expit(alpha0 + alpha1*W1 + alpha2*W2)
    G <- rbinom(1, 1, prG)
    
    # generate observation times confounder Z ~ N( mu_Z_G0, var_Z_G0) if G = 0 and 
    # N( mu_Z_G1, var_Z_G1) if G = 1
    
    Z1 <- ifelse(G == 0, rnorm(1, mu_Z_G0, sqrt(var_Z_G0)), rnorm(1, mu_Z_G1, sqrt(var_Z_G1)))
    
    # generate independent observation times covariate
    Z2 <- runif(1, 0.5, 1.5)
    
    # generate random effect
    phi <- rnorm(1, 0, var_phi)
    
    
    ## Generate Observation Times
    
    # simulate censoring time  according to
    censortime <- runif(1, tau/2, tau)
    
    # calculate eta
    eta <- rgamma(1, shape = 100, scale = 0.01) #eta is for random effect with mean 1 and sd 0.1 
    #(obs times within subject correlated if sd !0)
    
    # Observation intensity: 
    #lambda_i(t) = eta_i*sqrt(t)/2*exp(gamma1G_i + gamma2W_i1(t) + gamma3W_i2 + gamma4Z_i1 + gamma5Z_i2(t))
    lambdafun <- function(t) eta*sqrt(t)/2*exp(gamma1*G + gamma2*W1*log(t) + gamma3*W2 + gamma4*Z1 + gamma5*Z2*log(t))
    
    # generate observation times via thinning
    obstimes <- getObsTimesViaThinning(lambdafun, tau)
    obstimes <- round(obstimes,4)
    
    
    #calculate number of uncensored events
    nevents_uncensored <- length(obstimes)
    
    # calculate the number of observed events for individual i
    nevents <- sum(obstimes < censortime)
    
    if(nevents > 0){ #those censored before getting an observation are omitted
      
      
      simdata <- data.frame(cbind(rep(id, nevents), obstimes[1:nevents], rep(G, nevents), 
                                  rep(W1, nevents), rep(W2, nevents), rep(Z1, nevents), rep(Z2, nevents),
                                  rep(censortime, nevents)), rep(0, nevents))
      colnames(simdata) <- c("id", "time", "G", "W1", "W2", "Z1", "Z2", "censortime", "censored")
      
      if(censortime < tau){
        numcensored = numcensored + 1
        simdata[nevents + 1,] <- c(id, censortime, G, W1, W2, Z1, Z2, censortime, 1) 
      }
      
      
      ## Generate outcome
      # need conditional expectation and variance of Z | X
      
      simdata$cexp_W1t_G <- 0.5*log(simdata$time) #generated from U(0,1) independently from G
      simdata$cexp_W2_G <- 0 #generated from N(0,1) independently from G
      simdata$cexp_Z1_G <- ifelse(simdata$G == 0, mu_Z_G0, mu_Z_G1)
      simdata$cexp_Z2t_G <- log(simdata$time) #generated from U(0.5,1.5) independently from G
      
      
      simdata$y <- (2 - simdata$time) + beta1*simdata$G + beta2*(simdata$W1*log(simdata$time) - simdata$cexp_W1t_G) +
        beta3*(simdata$W2 - simdata$cexp_W2_G) + beta4*(simdata$Z1 - simdata$cexp_Z1_G) +
        beta5*(simdata$Z2*log(simdata$time) - simdata$cexp_Z2t_G) +
        rep(phi, dim(simdata)[1]) + rnorm(dim(simdata)[1], 0, sqrt(var_epsilon))
      
      
      #remove duplicated times
      #if subject has two of the same times, remove one of them.
      is.dup <- duplicated(simdata[, c(1,2)])#find which are duplicates (by ID and time)
      if(sum(is.dup > 0)){
        simdata <- simdata[!is.dup,] #remove them
      }
      
      #make sure observation times >0
      simdata <- simdata[simdata$time>0, ]
      
      
      simdatafull <- rbind(simdatafull, simdata)
      
    }
    
    
    
    id = id + 1
    
  }
  
  simdata <- simdatafull[-1,] #remove empty first row
  
  
  # order dataset by ID then time
  simdata <- simdata[with(simdata, order(id, time)), ]
  
  
  #include a variable counting observation number, to be used for lagging time for Surv function
  simdata$obsnumber <- with(simdata, ave(id, id, FUN = seq_along))
  simdata[which(simdata$censored == 1), 15] <- 0 # set observation number to 0 if its indicating censoring time
  
  # #create lagged time variable
  simdata$time.lag <- simdata$time[c(nrow(simdata),1:(nrow(simdata)-1))]
  simdata$time.lag[simdata$obsnumber == 1] <- 0
  
  
  # indicator for if subject had a visit (NOT censoring/death)
  simdata$visit <- 1 - simdata$censored
  
  #check the generated data
  numevents <- summary(tapply(simdata$visit, simdata$id, sum)) 
  
  newn <- length(unique(simdata$id)) #number of people after censoring etc
  
  out <- list(simdata, numevents, newn, numcensored)
  names(out) <- c("simdata", "numevents",  "newn", "numcensored")
  
  return(out) 
  
}



#~~~~~~~~~~~~~~~~~~ Reg. Simulation Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

simulateOneMIPI_timevar_5cov_full <- function(n, beta1, beta2, beta3, beta4, beta5,
                                              gamma1, gamma2, gamma3, gamma4, gamma5,
                                              alpha0, alpha1, alpha2, tau){
  
  # Simulates one instance of the simulation, obtaining estimates for beta1 under various weighting
  # models where we have the following models:
  # - Both models are correctly specified
  # - Probability of treatment model is missing $W_1$
  # - Probability of treatment model is missing $W_2$
  # - Intensity model is missing $Z_1$
  # - Intensity model is missing $Z_2(t)$
  # - Both models are missing $W_1 or W_1(t)$
  # - Both models are missing $W_2$
  # - Probability of treatment model is missing $W_1$ and intensity model is missing $Z_1$
  # - Probability of treatment model is missing $W_2$ and intensity model is missing $Z_2(t)$
  
  
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
  # Observation intensity: lambda_i(t) = eta_i*sqrt(t)/2*exp(gamma1G_i + gamma2W_i1(t) + gamma3W_i2 + gamma4Z_i1 + gamma5Z_i2(t))
  # Treatment assignment: P(G_i = 1) = expit(alpha0 + alpha1W_i1 + alpha2W_i2) 
  #
  # Censoring: C_i ~ U(tau/2, tau)
  #
  # tau: maximum follow up time, independent censoring.
  # N: number of simulations
  

  # simulate the data set 
  singlerun <- gendataMIPI_timevar_5cov_full(n, beta1, beta2, beta3, beta4, beta5,
                                             gamma1, gamma2, gamma3, gamma4, gamma5,
                                             alpha0, alpha1, alpha2, tau)
  
  simdata <- singlerun$simdata #only observed data
  numevents <- singlerun$numevents
  newn <- dim(simdata[simdata$obsnumber == 1,])[1]
  numcensored <- singlerun$numcensored
  
  
  # NAIVE UNWEIGHTED MODEL FOR COMPARISON
  beta1_naive <- summary(geeglm(y ~  G + offset(2 - time) - 1, id = id, data = simdata))$coef[1,1] #ATE is biased
  
  
  ### MODEL 1: Correctly Specified ###
  
  # IIW Weights
  gamma.hat_mod1 <- coxph(Surv(time.lag, time, visit) ~ G + tt(W1) +W2 +Z1 + tt(Z2) - 1, 
                          tt = function(x, t, ...) x*log(t), 
                     data = simdata, timefix = FALSE)$coef
  
  delta.hat_mod1 <- coxph(Surv(time.lag, time, visit) ~ G - 1, data = simdata, timefix = FALSE)$coef
  
  iiw_mod1 <- 1/exp(cbind(simdata$G, simdata$W1*log(simdata$time), 
                          simdata$W2, simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat_mod1)
  iiw_stab_mod1 <- exp(cbind(simdata$G)%*%delta.hat_mod1)/exp(cbind(simdata$G, 
                                simdata$W1*log(simdata$time), simdata$W2, 
                                simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat_mod1)
  
  # IPW weights
  psmod_mod1 <- glm(G ~ W1 + W2, family = binomial(link = "logit"), data = simdata)
  ps_mod1 <- expit(predict(psmod_mod1))
  prGmod_mod1 <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG_mod1 <- expit(predict(prGmod_mod1))
  
  ipw_mod1 <- 1/ps_mod1*simdata$G+1/(1-ps_mod1)*(1-simdata$G)
  ipw_stab_mod1 <- prG_mod1/ps_mod1*simdata$G+(1-prG_mod1)/(1-ps_mod1)*(1-simdata$G)

  # MIPIW weights
  
  mipiw_mod1 <- ipw_mod1*iiw_mod1
  mipiw_stab_mod1 <- ipw_stab_mod1*iiw_stab_mod1
  beta1_mipiw_mod1 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata,
                                  weights = mipiw_mod1))$coef[1,1]
  beta1_mipiw_stab_mod1 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, 
                                       weights = mipiw_stab_mod1))$coef[1,1]

  
  
  ### MODEL 2: IPTW MISSING W1, IIW Correctly Specified ###
  
  # IIW Weights
  gamma.hat_mod2 <- coxph(Surv(time.lag, time, visit) ~ G + tt(W1) +W2 +Z1 + tt(Z2) - 1, 
                          tt = function(x, t, ...) x*log(t), 
                          data = simdata, timefix = FALSE)$coef
  
  delta.hat_mod2 <- coxph(Surv(time.lag, time, visit) ~ G - 1, data = simdata, timefix = FALSE)$coef
  
  iiw_mod2 <- 1/exp(cbind(simdata$G, simdata$W1*log(simdata$time), 
                          simdata$W2, simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat_mod2)
  iiw_stab_mod2 <- exp(cbind(simdata$G)%*%delta.hat_mod2)/exp(cbind(simdata$G, 
                                          simdata$W1*log(simdata$time), simdata$W2, 
                                          simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat_mod2)
  
  # IPW weights
  psmod_mod2 <- glm(G ~ W2, family = binomial(link = "logit"), data = simdata)
  ps_mod2 <- expit(predict(psmod_mod2))
  prGmod_mod2 <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG_mod2 <- expit(predict(prGmod_mod2))
  
  ipw_mod2 <- 1/ps_mod2*simdata$G+1/(1-ps_mod2)*(1-simdata$G)
  ipw_stab_mod2 <- prG_mod2/ps_mod2*simdata$G+(1-prG_mod2)/(1-ps_mod2)*(1-simdata$G)
  
  # MIPIW weights
  
  mipiw_mod2 <- ipw_mod2*iiw_mod2
  mipiw_stab_mod2 <- ipw_stab_mod2*iiw_stab_mod2
  beta1_mipiw_mod2 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata,
                                  weights = mipiw_mod2))$coef[1,1]
  beta1_mipiw_stab_mod2 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, 
                                       weights = mipiw_stab_mod2))$coef[1,1]
  
  
  ### MODEL 3: IPTW MISSING W2, IIW Correctly Specified ###
  
  # IIW Weights
  gamma.hat_mod3 <- coxph(Surv(time.lag, time, visit) ~ G + tt(W1) +W2 +Z1 + tt(Z2) - 1, 
                          tt = function(x, t, ...) x*log(t), 
                          data = simdata, timefix = FALSE)$coef
  
  delta.hat_mod3 <- coxph(Surv(time.lag, time, visit) ~ G - 1, data = simdata, timefix = FALSE)$coef
  
  iiw_mod3 <- 1/exp(cbind(simdata$G, simdata$W1*log(simdata$time), 
                          simdata$W2, simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat_mod3)
  iiw_stab_mod3 <- exp(cbind(simdata$G)%*%delta.hat_mod3)/exp(cbind(simdata$G, 
                                                                    simdata$W1*log(simdata$time), simdata$W2, 
                                                                    simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat_mod3)
  
  # IPW weights
  psmod_mod3 <- glm(G ~ W1, family = binomial(link = "logit"), data = simdata)
  ps_mod3 <- expit(predict(psmod_mod3))
  prGmod_mod3 <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG_mod3 <- expit(predict(prGmod_mod3))
  
  ipw_mod3 <- 1/ps_mod3*simdata$G+1/(1-ps_mod3)*(1-simdata$G)
  ipw_stab_mod3 <- prG_mod3/ps_mod3*simdata$G+(1-prG_mod3)/(1-ps_mod3)*(1-simdata$G)
  
  # MIPIW weights
  
  mipiw_mod3 <- ipw_mod3*iiw_mod3
  mipiw_stab_mod3 <- ipw_stab_mod3*iiw_stab_mod3
  beta1_mipiw_mod3 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata,
                                  weights = mipiw_mod3))$coef[1,1]
  beta1_mipiw_stab_mod3 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, 
                                       weights = mipiw_stab_mod3))$coef[1,1]
  
  
  
  
  ### MODEL 4: IPTW CORRECTLY SPECIFIED, IIW MISSING Z1 ###
  
  # IIW Weights
  gamma.hat_mod4 <- coxph(Surv(time.lag, time, visit) ~ G + tt(W1) + W2 + tt(Z2) - 1, 
                          tt = function(x, t, ...) x*log(t), 
                          data = simdata, timefix = FALSE)$coef
  
  delta.hat_mod4 <- coxph(Surv(time.lag, time, visit) ~ G - 1, data = simdata, timefix = FALSE)$coef
  
  iiw_mod4 <- 1/exp(cbind(simdata$G, simdata$W1*log(simdata$time), 
                          simdata$W2, simdata$Z2*log(simdata$time))%*%gamma.hat_mod4)
  
  iiw_stab_mod4 <- exp(cbind(simdata$G)%*%delta.hat_mod4)/exp(cbind(simdata$G, 
                          simdata$W1*log(simdata$time), simdata$W2, 
                            simdata$Z2*log(simdata$time))%*%gamma.hat_mod4)
  
  # IPW weights
  psmod_mod4 <- glm(G ~ W1 + W2, family = binomial(link = "logit"), data = simdata)
  ps_mod4 <- expit(predict(psmod_mod4))
  prGmod_mod4 <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG_mod4 <- expit(predict(prGmod_mod4))
  
  ipw_mod4 <- 1/ps_mod4*simdata$G+1/(1-ps_mod4)*(1-simdata$G)
  ipw_stab_mod4 <- prG_mod4/ps_mod4*simdata$G+(1-prG_mod4)/(1-ps_mod4)*(1-simdata$G)
  
  # MIPIW weights
  
  mipiw_mod4 <- ipw_mod4*iiw_mod4
  mipiw_stab_mod4 <- ipw_stab_mod4*iiw_stab_mod4
  beta1_mipiw_mod4 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata,
                                  weights = mipiw_mod4))$coef[1,1]
  beta1_mipiw_stab_mod4 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, 
                                       weights = mipiw_stab_mod4))$coef[1,1]
  
  
  
  
  ### MODEL 5: IPTW CORRECTLY SPECIFIED, IIW MISSING Z_2(t) ###
  
  # IIW Weights
  gamma.hat_mod5 <- coxph(Surv(time.lag, time, visit) ~ G + tt(W1) + W2 + Z1 - 1, 
                          tt = function(x, t, ...) x*log(t), 
                          data = simdata, timefix = FALSE)$coef
  
  delta.hat_mod5 <- coxph(Surv(time.lag, time, visit) ~ G - 1, data = simdata, timefix = FALSE)$coef
  
  iiw_mod5 <- 1/exp(cbind(simdata$G, simdata$W1*log(simdata$time), 
                          simdata$W2, simdata$Z1)%*%gamma.hat_mod5)
  
  iiw_stab_mod5 <- exp(cbind(simdata$G)%*%delta.hat_mod5)/exp(cbind(simdata$G, 
                                          simdata$W1*log(simdata$time), simdata$W2, 
                                          simdata$Z1)%*%gamma.hat_mod5)
  
  # IPW weights
  psmod_mod5 <- glm(G ~ W1 + W2, family = binomial(link = "logit"), data = simdata)
  ps_mod5 <- expit(predict(psmod_mod5))
  prGmod_mod5 <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG_mod5 <- expit(predict(prGmod_mod5))
  
  ipw_mod5 <- 1/ps_mod5*simdata$G+1/(1-ps_mod5)*(1-simdata$G)
  ipw_stab_mod5 <- prG_mod5/ps_mod5*simdata$G+(1-prG_mod5)/(1-ps_mod5)*(1-simdata$G)
  
  # MIPIW weights
  
  mipiw_mod5 <- ipw_mod5*iiw_mod5
  mipiw_stab_mod5 <- ipw_stab_mod5*iiw_stab_mod5
  beta1_mipiw_mod5 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata,
                                  weights = mipiw_mod5))$coef[1,1]
  beta1_mipiw_stab_mod5 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, 
                                       weights = mipiw_stab_mod5))$coef[1,1]
  
  
  
  
  ### MODEL 6: BOTH WEIGHTING MODELS MISSING W_1 or W_1(t)###
  
  # IIW Weights
  gamma.hat_mod6 <- coxph(Surv(time.lag, time, visit) ~ G + W2 +Z1 + tt(Z2) - 1, 
                          tt = function(x, t, ...) x*log(t), 
                          data = simdata, timefix = FALSE)$coef
  
  delta.hat_mod6 <- coxph(Surv(time.lag, time, visit) ~ G - 1, data = simdata, timefix = FALSE)$coef
  
  iiw_mod6 <- 1/exp(cbind(simdata$G, 
                          simdata$W2, simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat_mod6)
  iiw_stab_mod6 <- exp(cbind(simdata$G)%*%delta.hat_mod6)/exp(cbind(simdata$G, 
                                              simdata$W2, 
                                              simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat_mod6)
  
  # IPW weights
  psmod_mod6 <- glm(G ~  W2, family = binomial(link = "logit"), data = simdata)
  ps_mod6 <- expit(predict(psmod_mod6))
  prGmod_mod6 <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG_mod6 <- expit(predict(prGmod_mod6))
  
  ipw_mod6 <- 1/ps_mod6*simdata$G+1/(1-ps_mod6)*(1-simdata$G)
  ipw_stab_mod6 <- prG_mod6/ps_mod6*simdata$G+(1-prG_mod6)/(1-ps_mod6)*(1-simdata$G)
  
  # MIPIW weights
  
  mipiw_mod6 <- ipw_mod6*iiw_mod6
  mipiw_stab_mod6 <- ipw_stab_mod6*iiw_stab_mod6
  beta1_mipiw_mod6 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata,
                                  weights = mipiw_mod6))$coef[1,1]
  beta1_mipiw_stab_mod6 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, 
                                       weights = mipiw_stab_mod6))$coef[1,1]
  
  
  ### MODEL 6: BOTH WEIGHTING MODELS MISSING W2
  
  # IIW Weights
  gamma.hat_mod7 <- coxph(Surv(time.lag, time, visit) ~ G + tt(W1) +Z1 + tt(Z2) - 1, 
                          tt = function(x, t, ...) x*log(t), 
                          data = simdata, timefix = FALSE)$coef
  
  delta.hat_mod7 <- coxph(Surv(time.lag, time, visit) ~ G - 1, data = simdata, timefix = FALSE)$coef
  
  iiw_mod7 <- 1/exp(cbind(simdata$G, simdata$W1*log(simdata$time), 
                          simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat_mod7)
  iiw_stab_mod7 <- exp(cbind(simdata$G)%*%delta.hat_mod7)/exp(cbind(simdata$G, 
                                      simdata$W1*log(simdata$time), 
                                      simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat_mod7)
  
  # IPW weights
  psmod_mod7 <- glm(G ~ W1, family = binomial(link = "logit"), data = simdata)
  ps_mod7 <- expit(predict(psmod_mod7))
  prGmod_mod7 <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG_mod7 <- expit(predict(prGmod_mod7))
  
  ipw_mod7 <- 1/ps_mod7*simdata$G+1/(1-ps_mod7)*(1-simdata$G)
  ipw_stab_mod7 <- prG_mod7/ps_mod7*simdata$G+(1-prG_mod7)/(1-ps_mod7)*(1-simdata$G)
  
  # MIPIW weights
  
  mipiw_mod7 <- ipw_mod7*iiw_mod7
  mipiw_stab_mod7 <- ipw_stab_mod7*iiw_stab_mod7
  beta1_mipiw_mod7 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata,
                                  weights = mipiw_mod7))$coef[1,1]
  beta1_mipiw_stab_mod7 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, 
                                       weights = mipiw_stab_mod7))$coef[1,1]
  
  
  ### MODEL 8: IPTW MODEL MISSING W_1, IIW MODEL MISSING Z_1
  
  # IIW Weights
  gamma.hat_mod8 <- coxph(Surv(time.lag, time, visit) ~ G + tt(W1) + W2 + tt(Z2) - 1, 
                          tt = function(x, t, ...) x*log(t), 
                          data = simdata, timefix = FALSE)$coef
  
  delta.hat_mod8 <- coxph(Surv(time.lag, time, visit) ~ G - 1, data = simdata, timefix = FALSE)$coef
  
  iiw_mod8 <- 1/exp(cbind(simdata$G, simdata$W1*log(simdata$time), 
                          simdata$W2, simdata$Z2*log(simdata$time))%*%gamma.hat_mod8)
  iiw_stab_mod8 <- exp(cbind(simdata$G)%*%delta.hat_mod8)/exp(cbind(simdata$G, 
                                      simdata$W1*log(simdata$time), simdata$W2, 
                                      simdata$Z2*log(simdata$time))%*%gamma.hat_mod8)
  
  # IPW weights
  psmod_mod8 <- glm(G ~ W2, family = binomial(link = "logit"), data = simdata)
  ps_mod8 <- expit(predict(psmod_mod8))
  prGmod_mod8 <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG_mod8 <- expit(predict(prGmod_mod8))
  
  ipw_mod8 <- 1/ps_mod8*simdata$G+1/(1-ps_mod8)*(1-simdata$G)
  ipw_stab_mod8 <- prG_mod8/ps_mod8*simdata$G+(1-prG_mod8)/(1-ps_mod8)*(1-simdata$G)
  
  # MIPIW weights
  
  mipiw_mod8 <- ipw_mod8*iiw_mod8
  mipiw_stab_mod8 <- ipw_stab_mod8*iiw_stab_mod8
  beta1_mipiw_mod8 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata,
                                  weights = mipiw_mod8))$coef[1,1]
  beta1_mipiw_stab_mod8 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, 
                                       weights = mipiw_stab_mod8))$coef[1,1]
  
  
  ### MODEL 9: IPTW MODEL MISSING W2, IIW MODEL MISSING Z2(t)
  # IIW Weights
  gamma.hat_mod9 <- coxph(Surv(time.lag, time, visit) ~ G + tt(W1) +W2 +Z1 - 1, 
                          tt = function(x, t, ...) x*log(t), 
                          data = simdata, timefix = FALSE)$coef
  
  delta.hat_mod9 <- coxph(Surv(time.lag, time, visit) ~ G - 1, data = simdata, timefix = FALSE)$coef
  
  iiw_mod9 <- 1/exp(cbind(simdata$G, simdata$W1*log(simdata$time), 
                          simdata$W2, simdata$Z1)%*%gamma.hat_mod9)
  iiw_stab_mod9 <- exp(cbind(simdata$G)%*%delta.hat_mod9)/exp(cbind(simdata$G, 
                                      simdata$W1*log(simdata$time), simdata$W2, 
                                      simdata$Z1)%*%gamma.hat_mod9)
  
  # IPW weights
  psmod_mod9 <- glm(G ~ W1, family = binomial(link = "logit"), data = simdata)
  ps_mod9 <- expit(predict(psmod_mod9))
  prGmod_mod9 <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG_mod9 <- expit(predict(prGmod_mod9))
  
  ipw_mod9 <- 1/ps_mod9*simdata$G+1/(1-ps_mod9)*(1-simdata$G)
  ipw_stab_mod9 <- prG_mod9/ps_mod9*simdata$G+(1-prG_mod9)/(1-ps_mod9)*(1-simdata$G)
  
  # MIPIW weights
  
  mipiw_mod9 <- ipw_mod9*iiw_mod9
  mipiw_stab_mod9 <- ipw_stab_mod9*iiw_stab_mod9
  beta1_mipiw_mod9 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata,
                                  weights = mipiw_mod9))$coef[1,1]
  beta1_mipiw_stab_mod9 <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, 
                                       weights = mipiw_stab_mod9))$coef[1,1]
  
  
  
  
  out <- list(beta1,  beta1_naive,  
              beta1_mipiw_mod1, beta1_mipiw_stab_mod1,
              beta1_mipiw_mod2, beta1_mipiw_stab_mod2,
              beta1_mipiw_mod3, beta1_mipiw_stab_mod3,
              beta1_mipiw_mod4, beta1_mipiw_stab_mod4,
              beta1_mipiw_mod5, beta1_mipiw_stab_mod5,
              beta1_mipiw_mod6, beta1_mipiw_stab_mod6,
              beta1_mipiw_mod7, beta1_mipiw_stab_mod7,
              beta1_mipiw_mod8, beta1_mipiw_stab_mod8,
              beta1_mipiw_mod9, beta1_mipiw_stab_mod9)
  names(out) <- c('beta1',  'beta1_naive',  
                       'beta1_mipiw_mod1', 'beta1_mipiw_stab_mod1',
                       'beta1_mipiw_mod2', 'beta1_mipiw_stab_mod2',
                       'beta1_mipiw_mod3', 'beta1_mipiw_stab_mod3',
                       'beta1_mipiw_mod4', 'beta1_mipiw_stab_mod4',
                       'beta1_mipiw_mod5', 'beta1_mipiw_stab_mod5',
                       'beta1_mipiw_mod6', 'beta1_mipiw_stab_mod6',
                       'beta1_mipiw_mod7', 'beta1_mipiw_stab_mod7',
                       'beta1_mipiw_mod8', 'beta1_mipiw_stab_mod8',
                       'beta1_mipiw_mod9', 'beta1_mipiw_stab_mod9')
  return(out)
}



simulateResultsMIPI_timevar_5cov_full <-  function(N, n, beta1, beta2, beta3, beta4, beta5,
                                                   gamma1, gamma2, gamma3, gamma4, gamma5,
                                                   alpha0, alpha1, alpha2, tau){
  # Simulates N instances of the given scheme and produces a table that summarizes the findings 
  # (var, bias mse) 
  
  results_beta1 <- matrix(data = NA, nrow = N, ncol = 10)
  results_beta1_stab <- matrix(data = NA, nrow = N, ncol = 10)
  
  for(i in 1:N){
    if(i%%100 == 0){print(i)}
    simrun <- simulateOneMIPI_timevar_5cov_full(n, beta1, beta2, beta3,beta4, beta5,
                                                gamma1, gamma2, gamma3, gamma4, gamma5,
                                                 alpha0, alpha1, alpha2, tau)
    results_beta1[i,] <- c(simrun$beta1_naive, simrun$beta1_mipiw_mod1,
                           simrun$beta1_mipiw_mod2, simrun$beta1_mipiw_mod3,
                           simrun$beta1_mipiw_mod4, simrun$beta1_mipiw_mod5,
                           simrun$beta1_mipiw_mod6, simrun$beta1_mipiw_mod7,
                           simrun$beta1_mipiw_mod8, simrun$beta1_mipiw_mod9) 
    
    results_beta1_stab[i,] <- c(simrun$beta1_naive, simrun$beta1_mipiw_stab_mod1,
                            simrun$beta1_mipiw_stab_mod2, simrun$beta1_mipiw_stab_mod3,
                            simrun$beta1_mipiw_stab_mod4, simrun$beta1_mipiw_stab_mod5,
                            simrun$beta1_mipiw_stab_mod6, simrun$beta1_mipiw_stab_mod7, 
                            simrun$beta1_mipiw_stab_mod8, simrun$beta1_mipiw_stab_mod9) 
  }
  
  
  # calculate the bias, variance, mse for each method/weighting model
  bias_beta1 <- round(apply(results_beta1, FUN =  mean, MARGIN = 2)  - beta1 ,3)
  names(bias_beta1) <- c("Naive", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", 'Model 7', "Model 8", "Model 9")

  var_beta1 <- round(apply(results_beta1, FUN = var, MARGIN = 2), 3)
  names(var_beta1) <- c("Naive", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", 'Model 7', "Model 8", "Model 9")

  mse_beta1 <- round(apply((results_beta1 - beta1)^2, FUN = mean, MARGIN = 2), 3)
  names(mse_beta1) <- c("Naive", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", 'Model 7', "Model 8", "Model 9")

  
  bias_beta1_stab <- round(apply(results_beta1_stab, FUN =  mean, MARGIN = 2)  - beta1 ,3)
  names(bias_beta1_stab) <- c("Naive", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", 'Model 7', "Model 8", "Model 9")
  
  var_beta1_stab <- round(apply(results_beta1_stab, FUN = var, MARGIN = 2), 3)
  names(var_beta1_stab) <- c("Naive", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", 'Model 7', "Model 8", "Model 9")
  
  mse_beta1_stab <- round(apply((results_beta1_stab - beta1)^2, FUN = mean, MARGIN = 2), 3)
  names(mse_beta1_stab) <- c("Naive", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", 'Model 7', "Model 8", "Model 9")
  
  # put the results in a table
  tab <- cbind(bias_beta1, var_beta1, mse_beta1, bias_beta1_stab, var_beta1_stab, mse_beta1_stab)
  colnames(tab) <- c("Bias", "variance", "MSE", "Bias", "variance", "MSE")
  
  
  kabletab <- kable(tab, booktabs = T, digits = 3 ,format = "latex") %>%
    add_header_above(c( Model = 1, "Non-stablized" = 3, "Stabilized" = 3))
  
  out <- list(tab, kabletab)
  names(out) <- c('tab', 'kabletab')
  
  
  
  return(out)
  
}










