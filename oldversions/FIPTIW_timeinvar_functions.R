
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
# Treatment assignment: P(G_i = 1) = expit(alpha0 + alpha1W_i) 
#
#
# tau: maximum follow up time, independent censoring.
# N: number of simulations
#


require(survival)
require(geepack)
require(nleqslv)
require(knitr)
require(kableExtra)

expit <- function(x){ return(exp(x)/(1+exp(x)))}


#~~~~~~~~~~~~~ Data Generation ~~~~~~~~~~~~~~~~#

gendataMIPI <- function(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau){
  # n: number of subjects
  # tau: maximum follow up
  #
  # Data generation model: Y_i(t)  = (2-t) + beta1G_i + beta2[W_i - E(w_i|G_i)] + 
  #                                     beta3[Z_i - E(Z_i|G_i)] + phi_i + epsilon_i(t)
  #
  # Marginal model: E(Y_i(t) | G_i) = E{(2-t) + beta1G_i + beta2[W_i - E(w_i|G_i)] + beta3(Z_i - E(Z_i|G_i)) | G_i]
  #                                 = (2-t) + beta1G_i + betaE{2[W_i - W(W_i|G_i) | G_i] + beta3E{Z_i - E(Z_i|G_i)|G_i}
  #                                 = (2-t) + beta1G_i 
  #  
  # Observation intensity: lambda_i(t) = eta_i*exp(gamma1G_i + gamma2W_i + gamma3Z_i)
  # Treatment assignment: P(G_i = 1) = expit(alpha0 + alpha1W_i) 
  #
  # where
  #
  # G: Treatment assignment indicator
  # W: Treatment confounder
  # Z: Observation times confounder
  # Y: Time-varying outcome
  # t: observation time
  # phi: random effect
  # epsilon: random error, function of time
  
  
  
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
  simdatafull <- data.frame(matrix(NA, ncol = 9))
  colnames(simdatafull) <- c("id", "time", "G", "W", "Z", "censortime", "cexp_W_G", "cexp_Z_G",  "y")
  
  while(id <= n){
    
    ## Generate Covariates
    
    # generate treatment confounder w ~ N(0,1)
    W <- rnorm(1, 0, 1)
    
    # generate treatment assignment
    prG <- expit(alpha0 + alpha1*W)
    G <- rbinom(1, 1, prG)
    
    # generate observation times confounder Z ~ N( mu_Z_G0, var_Z_G0) if G = 0 and 
    # N( mu_Z_G1, var_Z_G1) if G = 1

    Z <- ifelse(G == 0, rnorm(1, mu_Z_G0, sqrt(var_Z_G0)), rnorm(1, mu_Z_G1, sqrt(var_Z_G1)))
   
    
    # generate random effect
    phi <- rnorm(1, 0, var_phi)
    
    
    ## Generate Observation Times
    
    # simulate censoring time
    censortime <- runif(1, tau/2, tau)
    
    # calculate eta
    eta <- rgamma(1, shape = 100, scale = 0.01) #eta is for random effect with mean 1 and sd 0.1 
    #(obs times within subject correlated if sd !0)
    
    
    lambdafun <- function(t) eta*sqrt(t)/2*exp(gamma1*G + gamma2*W + gamma3*Z)
    
    # generate observation times via thinning
    obstimes <- getObsTimesViaThinning(lambdafun, tau)
    obstimes <- round(obstimes,4)
    
    if(length(obstimes) > 0){
      
      
      # calculate the number of events for individual i
      nevents <- length(obstimes)
      
      simdata <- data.frame(cbind(rep(id, nevents), obstimes, rep(G, nevents), 
                                  rep(W, nevents), rep(Z, nevents), rep(censortime, nevents)))
      colnames(simdata) <- c("id", "time", "G", "W", "Z", "censortime")
      
      
      
      ## Generate outcome
      # need conditional expectation and variance of Z | X
      
      simdata$cexp_W_G <- 0 #generated from N(0,1) independently from G
      simdata$cexp_Z_G <- ifelse(simdata$G == 0, mu_Z_G0, mu_Z_G1)
      
      
      simdata$y <- (2 - simdata$time) + beta1*simdata$G + beta2*(simdata$W - simdata$cexp_W_G) +
        beta3*(simdata$Z - simdata$cexp_Z_G) + 
        rep(phi, dim(simdata)[1]) + rnorm(dim(simdata)[1], 0, sqrt(var_epsilon))
      
      
      #remove duplicated times
      #if subject has two of the same times, remove one of them.
      is.dup <- duplicated(simdata[, c(1,2)])#find which are duplicates (by ID and time)
      if(sum(is.dup > 0)){
        simdata <- simdata[!is.dup,] #remove them
      }
      
      #make sure observation times >0
      simdata <- simdata[simdata$time>0, ]
      
      # filter censored observations
      censored <- which(simdata$time > censortime)
      
      if(length(censored) > 0){
        # remove censored observations
        simdata <- simdata[-censored,]
      }
      
      simdatafull <- rbind(simdatafull, simdata)
      
    }
    
    
    
    id = id + 1
    
  }
  
  simdata <- simdatafull[-1,] #remove empty first row
  
  
  # order dataset by ID then time
  simdata <- simdata[with(simdata, order(id, time)), ]
  
  
  #include a variable counting observation number, to be used for lagging time for Surv function
  simdata$obsnumber <- with(simdata, ave(id, id, FUN = seq_along))
  
  # #create lagged time variable
  simdata$time.lag <- simdata$time[c(nrow(simdata),1:(nrow(simdata)-1))]
  simdata$time.lag[simdata$obsnumber == 1] <- 0
  
  simdata$event <- 1
  
  
  #check the generated data
  numevents <- summary(tapply(simdata$event, simdata$id, sum)) 
  
  
  #also get data for individuals at baseline (all same here)
  baselinedata <- simdata[simdata$obsnumber == 1, ]
  newn <- dim(baselinedata)[1] #number of people after censoring etc
  
  out <- list(simdata, numevents, baselinedata, newn)
  names(out) <- c("simdata", "numevents", "baselinedata", "newn")
  
  return(out) 
    
}


simulateOneMIPI <- function(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau){
  # Simulates one instance of the simulation, obtaining estimates for beta1 under various weighting
  #
  # n: number of subjects
  # tau: maximum follow up
  #
  # Data generation model: Y_i(t)  = (2-t) + beta1G_i + beta2[W_i - E(w_i|G_i)] + 
  #                                     beta3[Z_i - E(Z_i|G_i)] + phi_i + epsilon_i(t)
  #
  # Marginal model: E(Y_i(t) | G_i) = E{(2-t) + beta1G_i + beta2[W_i - E(w_i|G_i)] + beta3(Z_i - E(Z_i|G_i)) | G_i]
  #                                 = (2-t) + beta1G_i + betaE{2[W_i - W(W_i|G_i) | G_i] + beta3E{Z_i - E(Z_i|G_i)|G_i}
  #                                 = (2-t) + beta1G_i 
  #  
  # Observation intensity: lambda_i(t) = eta_i*exp(gamma1G_i + gamma2W_i + gamma3Z_i)
  # Treatment assignment: P(G_i = 1) = expit(alpha0 + alpha1W_i) 
  #
  # where
  #
  # G: Treatment assignment indicator
  # W: Treatment confounder
  # Z: Observation times confounder
  # Y: Time-varying outcome
  # t: observation time
  # phi: random effect
  # epsilon: random error, function of time
  
  
  
  
  singlerun <- gendataMIPI(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau)
  simdata <- singlerun$simdata
  numevents <- singlerun$numevents
  baselinedata <- singlerun$baselinedata
  newn <- singlerun$newn
  
  
  beta1_naive <- summary(geeglm(y ~  G + offset(2 - time) - 1, id = id, data = simdata))$coef[1,1] #ATE is biased
  
  
  # IIW Weights
  gamma.hat <- coxph(Surv(time.lag, time, event) ~ G + W + Z - 1, data = simdata)$coef
  delta.hat <- coxph(Surv(time.lag, time, event) ~ G - 1, data = simdata)$coef
    
  iiw <- 1/exp(cbind(simdata$G, simdata$W, simdata$Z)%*%gamma.hat)
  iiw_stab <- exp(cbind(simdata$G)%*%delta.hat)/exp(cbind(simdata$G, simdata$W, simdata$Z)%*%gamma.hat)
    
  beta1_iiw <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = iiw))$coef[1,1]
  beta1_iiw_stab <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = iiw_stab))$coef[1,1]
  
  # IPW weights
  psmod <- glm(G ~ W, family = binomial(link = "logit"), data = simdata)
  ps <- expit(predict(psmod))
  prGmod <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG <- expit(predict(prGmod))
  
  ipw <- 1/ps*simdata$G+1/(1-ps)*(1-simdata$G)
  ipw_stab <- prG/ps*simdata$G+(1-prG)/(1-ps)*(1-simdata$G)
  
  
  beta1_ipw <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = ipw))$coef[1,1]
  beta1_ipw_stab <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = ipw_stab))$coef[1,1]
  
  mipiw <- ipw*iiw
  mipiw_stab <- ipw_stab*iiw_stab
  beta1_mipiw <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = mipiw))$coef[1,1]
  beta1_mipiw_stab <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = mipiw_stab))$coef[1,1]
  
  
  
  out <- list(beta1,  beta1_naive,  beta1_iiw,  beta1_iiw_stab, 
              beta1_ipw, beta1_ipw_stab, 
              beta1_mipiw, beta1_mipiw_stab )
  names(out) <- c('beta1',  'beta1_naive',  'beta1_iiw',  'beta1_iiw_stab', 
                  'beta1_ipw', 'beta1_ipw_stab', 'beta1_mipiw', 'beta1_mipiw_stab' )
  return(out)
}



simulateResultsMIPI<-  function(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau){
  # Simulates N instances of the given scheme 
  #
  # N: number of simulation runs
  # n: number of subjects
  # tau: maximum follow up
  #
  # Data generation model: Y_i(t)  = (2-t) + beta1G_i + beta2[W_i - E(w_i|G_i)] + 
  #                                     beta3[Z_i - E(Z_i|G_i)] + phi_i + epsilon_i(t)
  #
  # Marginal model: E(Y_i(t) | G_i) = E{(2-t) + beta1G_i + beta2[W_i - E(w_i|G_i)] + beta3(Z_i - E(Z_i|G_i)) | G_i]
  #                                 = (2-t) + beta1G_i + betaE{2[W_i - W(W_i|G_i) | G_i] + beta3E{Z_i - E(Z_i|G_i)|G_i}
  #                                 = (2-t) + beta1G_i 
  #  
  # Observation intensity: lambda_i(t) = eta_i*exp(gamma1G_i + gamma2W_i + gamma3Z_i)
  # Treatment assignment: P(G_i = 1) = expit(alpha0 + alpha1W_i) 
  #
  # where
  #
  # G: Treatment assignment indicator
  # W: Treatment confounder
  # Z: Observation times confounder
  # Y: Time-varying outcome
  # t: observation time
  # phi: random effect
  # epsilon: random error, function of time
  
  
  #
  
  results_beta1 <- matrix(data = NA, nrow = N, ncol = 7)
  
  for(i in 1:N){
    if(i%%100 == 0){print(i)}
    simrun <- simulateOneMIPI(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau)
    results_beta1[i,] <- c(simrun$beta1_naive, simrun$beta1_iiw, simrun$beta1_iiw_stab,
                           simrun$beta1_ipw, simrun$beta1_ipw_stab,
                           simrun$beta1_mipiw, simrun$beta1_mipiw_stab) 
  }
  
  
  
  bias_beta1 <- round(apply(results_beta1, FUN =  mean, MARGIN = 2)  - beta1 ,3)
  names(bias_beta1) <- c("Naive", "IIW", "IIW Stabilized", "IPW", "IPW Stabilized", "MIPIW", "MIPIW Stabilized")
 
  var_beta1 <- round(apply(results_beta1, FUN = var, MARGIN = 2), 3)
  names(var_beta1) <- c("Naive", "IIW", "IIW Stabilized", "IPW", "IPW Stabilized", "MIPIW", "MIPIW Stabilized")
  
  mse_beta1 <- round(apply((results_beta1 - beta1)^2, FUN = mean, MARGIN = 2), 3)
  names(mse_beta1) <- c("Naive", "IIW", "IIW Stabilized", "IPW", "IPW Stabilized", "MIPIW", "MIPIW Stabilized")
  
  naive_beta1 <- c(bias_beta1[1], var_beta1[1], mse_beta1[1])
  names(naive_beta1) <- c("Bias", "Var", "MSE")

  IIW_beta1 <- c(bias_beta1[2], var_beta1[2], mse_beta1[2])
  names(IIW_beta1) <- c("Bias", "Var", "MSE")

  IIWstab_beta1 <- c(bias_beta1[3], var_beta1[3], mse_beta1[3])
  names(IIWstab_beta1) <- c("Bias", "Var", "MSE")
  
  IPW_beta1 <- c(bias_beta1[4], var_beta1[4], mse_beta1[4])
  names(IPW_beta1) <- c("Bias", "Var", "MSE")
  
  IPWstab_beta1 <- c(bias_beta1[5], var_beta1[5], mse_beta1[5])
  names(IPWstab_beta1) <- c("Bias", "Var", "MSE")
  
  MIPIW_beta1 <- c(bias_beta1[6], var_beta1[6], mse_beta1[6])
  names(MIPIW_beta1) <- c("Bias", "Var", "MSE")
  
  MIPIWstab_beta1 <- c(bias_beta1[7], var_beta1[7], mse_beta1[7])
  names(MIPIWstab_beta1) <- c("Bias", "Var", "MSE")

  

  out <- list(bias_beta1, var_beta1, mse_beta1, 
              naive_beta1,  IIW_beta1,  IIWstab_beta1, IPW_beta1, IPWstab_beta1, 
              MIPIW_beta1, MIPIWstab_beta1)
  
  names(out) <- c('bias_beta1', 'var_beta1', 'mse_beta1', 
                  'naive_beta1',  'IIW_beta1',  'IIWstab_beta1', 'IPW_beta1', 'IPWstab_beta1', 
                  'MIPIW_beta1', 'MIPIWstab_beta1')
  return(out)
  
}


simulateALLMIPI <- function(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, alpha0, alpha1vec, tau){
  #N: number of simulation runs
  #n: vector of sample sizes
  #beta1: coefficient for Xi(t) in logistic outcome model
  #beta2: vector of coefficients to consider for outcome generation model
  #gamma1, gamma2 parameters for intensity for Xi(t) and Zi, respectively
  #tau: maximum follow-up time
  
  #This function aggregates simulation results for varying n and beta2
  
  
  
  resultsmat <- matrix(NA, nrow = 1, ncol = 15)
  resultsmat_stabilized <- matrix(NA, nrow = 1, ncol = 15)
  
  i = 1
  
  
  for(gamma2 in gamma2vec){
    for(gamma3 in gamma3vec){
      for(alpha1 in alpha1vec){
      
        print(paste("Now on gamma2 = ", gamma2, ", gamma3 =", gamma3, ", alpha1 =  ", alpha1, sep = ""))
        result <- simulateResultsMIPI(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau)
        
        resultsmat <- rbind(resultsmat, c(gamma2,  gamma3, alpha1,
                                          result$naive_beta1, result$IIW_beta1, 
                                          result$IPW_beta1, result$MIPIW_beta1))
        
        resultsmat_stabilized <- rbind(resultsmat_stabilized, c(gamma2,  gamma3, alpha1,
                                          result$naive_beta1, result$IIWstab_beta1, 
                                          result$IPWstab_beta1, result$MIPIWstab_beta1))
        i <- i + 1
    
      }
    }
  }
  
  
  
  
  resultsmat<- resultsmat[-1,]
  resultsmat_stabilized <- resultsmat_stabilized[-1,]
  
  colnames(resultsmat) <- c( "gamma2",  "gamma3", "alpha1", "Bias", "Var" , "MSE" , "Bias", "Var" , "MSE" , "Bias", "Var" , "MSE", "Bias", "Var" , "MSE" )
  colnames(resultsmat_stabilized) <- c( "gamma2", "gamma3", "alpha1", "Bias", "Var" , "MSE" , "Bias", "Var" , "MSE" , "Bias", "Var" , "MSE", "Bias", "Var" , "MSE" )
  
  out <- list(resultsmat, resultsmat_stabilized)
  names(out) <- c("resultsmat", "resultsmat_stab")
  return(out)
}
  
  