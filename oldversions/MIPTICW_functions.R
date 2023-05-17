# The purpose of this sensitivity analysis is to look at the performance of MIPTCIW

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
# Censoring hazard: h(t) = lambda_c*exp(kappa1G + kappa2 W1 + kappa3 W2 + kappa4Z1 + kappa5Z2)
#
# tau: maximum follow up time, independent censoring.
# N: number of simulations


require(survival)
require(geepack)
require(nleqslv)
require(knitr)
require(kableExtra)

expit <- function(x){ return(exp(x)/(1+exp(x)))}


#~~~~~~~~~~~~~ Data Generation ~~~~~~~~~~~~~~~~#

gendataMIPI_timevar_5cov_censoring <- function(n, beta1, beta2, beta3, beta4, beta5,
                                               gamma1, gamma2, gamma3, gamma4, gamma5,
                                               kappa1, kappa2, kappa3, kappa4, kappa5,
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
  # Censoring hazard: h(t) = lambda_c*exp(kappa1G + kappa2 W1 + kappa3 W2 + kappa4Z1 + kappa5Z2)
  #
  # tau: maximum follow up time, independent censoring.
  # N: number of simulations
  
  mu_Z_G0 <- 2
  var_Z_G0 <- 1 
  mu_Z_G1 <- 0
  var_Z_G1 <- 0.5
  
  var_phi <- 0.25
  var_epsilon <- 1
  
  lambdac <- 0.1
  
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
    # Censoring hazard: h(t) = h_0*exp(kappa1G_i + kappa2W_i + kappa3Z_i)
    U <- runif(1, 0, 1)
    censortime <- -log(U)/(lambdac*exp(kappa1*G + kappa2*W1 + kappa3*W2 + kappa4*Z1 + kappa5*Z2))
    # censortime <- runif(1, tau/2, tau)
    
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

simulateOneMIPI_timevar_5cov_censoring<- function(n, beta1, beta2, beta3, beta4, beta5,
                                                  gamma1, gamma2, gamma3, gamma4, gamma5,
                                                  kappa1, kappa2, kappa3, kappa4, kappa5,
                                                  alpha0, alpha1, alpha2, tau){
  
  # Simulates one instance of the simulation, obtaining estimates for beta1 under various weighting
  # procedures
  
  
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
  # Censoring hazard
  #
  # tau: maximum follow up time, independent censoring.
  # N: number of simulations
  
  
  
  singlerun <- gendataMIPI_timevar_5cov_censoring(n, beta1, beta2, beta3, beta4, beta5,
                                                  gamma1, gamma2, gamma3, gamma4, gamma5,
                                                  kappa1, kappa2, kappa3, kappa4, kappa5,
                                                  alpha0, alpha1, alpha2, tau)
  
  simdata <- singlerun$simdata #only observed data
  numevents <- singlerun$numevents
  newn <- dim(simdata[simdata$obsnumber == 1,])[1]
  numcensored <- singlerun$numcensored
  
  
  
  #IPCW Weights

  ipcw <- tryCatch(
      expr = {
        ipcwmod <- coxph(Surv(time.lag, time, censored) ~ G + W1 + W2 + Z1 + Z2, data = simdata)
        ipcwfit <- summary(survfit(ipcwmod), times = simdata$time.lag)
        ipcw <- 1/ipcwfit$surv
    },
    error = function(e){
      ipcw = rep(1, dim(simdata)[1])
    },
    warning = function(w){
      ipcw = rep(1, dim(simdata)[1])
    }
  )

  #filter on only the observed data (pre censoring) for the rest of the weight calculations
  observedindices <- which(simdata$visit == 1)
  simdata <- simdata[observedindices,]
  ipcwweights <- ipcw[observedindices]
  
  beta1_ipcw <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = ipcwweights))$coef[1,1]
  
  
  
  
  
  
  beta1_naive <- summary(geeglm(y ~  G + offset(2 - time) - 1, id = id, data = simdata))$coef[1,1] #ATE is biased
  
  
  # IIW Weights
  gamma.hat <- coxph(Surv(time.lag, time, visit) ~ G + tt(W1) +W2 +Z1 + tt(Z2) - 1, tt = function(x, t, ...) x*log(t), 
                     data = simdata, timefix = FALSE)$coef
  delta.hat <- coxph(Surv(time.lag, time, visit) ~ G - 1, data = simdata, timefix = FALSE)$coef
  
  iiw <- 1/exp(cbind(simdata$G, simdata$W1*log(simdata$time), simdata$W2, simdata$Z1, simdata$Z2*log(simdata$time))%*%gamma.hat)
  # iiw_stab <- exp(cbind(simdata$G)%*%delta.hat)/exp(cbind(simdata$G, simdata$W*log(simdata$time), simdata$Z)%*%gamma.hat)
  
  beta1_iiw <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = iiw))$coef[1,1]
  # beta1_iiw_stab <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = iiw_stab))$coef[1,1]
  
  # IPW weights
  psmod <- glm(G ~ W1 + W2, family = binomial(link = "logit"), data = simdata)
  ps <- expit(predict(psmod))
  prGmod <- glm(G ~ 1, family = binomial(link = "logit"), data = simdata)
  prG <- expit(predict(prGmod))
  
  ipw <- 1/ps*simdata$G+1/(1-ps)*(1-simdata$G)
  # ipw_stab <- prG/ps*simdata$G+(1-prG)/(1-ps)*(1-simdata$G)
  
  beta1_ipw <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = ipw))$coef[1,1]
  # beta1_ipw_stab <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = ipw_stab))$coef[1,1]
  
  

  # MIPIW weights
  
  mipiw <- ipw*iiw
  # mipiw_stab <- ipw_stab*iiw_stab
  beta1_mipiw <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = mipiw))$coef[1,1]
  # beta1_mipiw_stab <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = mipiw_stab))$coef[1,1]
  # 
  #MIPICW Weights
  mipticw <- ipw*iiw*ipcwweights
  beta1_mipticw <- summary(glm(y ~ G + offset(2-time) - 1, data=simdata, weights = mipticw))$coef[1,1]

  
  out <- list(beta1,  beta1_naive,  beta1_iiw, 
              beta1_ipw, beta1_ipcw, 
              beta1_mipiw, beta1_mipticw)
  names(out) <- c('beta1',  'beta1_naive',  'beta1_iiw',  
                  'beta1_ipw', 'beta1_ipcw ',
                  'beta1_mipiw', 'beta1_mipticw')
  return(out)
}



simulateResultsMIPI_timevar_5cov_censoring <-  function(N, n, beta1, beta2, beta3, beta4, beta5,
                                                        gamma1, gamma2, gamma3, gamma4, gamma5,
                                                        kappa1, kappa2, kappa3, kappa4, kappa5,
                                                        alpha0, alpha1, alpha2, tau){
  # Simulates N instances of the given scheme 
  
  results_beta1 <- matrix(data = NA, nrow = N, ncol = 6)
  
  for(i in 1:N){
    if(i%%100 == 0){print(i)}
    simrun <- simulateOneMIPI_timevar_5cov_censoring(n, beta1, beta2, beta3,beta4, beta5,
                                                     gamma1, gamma2, gamma3, gamma4, gamma5,
                                                     kappa1, kappa2, kappa3, kappa4, kappa5,
                                                     alpha0, alpha1, alpha2, tau)
    
    results_beta1[i,] <- c(simrun$beta1_naive, simrun$beta1_iiw, 
                           simrun$beta1_ipw,
                           simrun$beta1_ipcw,
                           simrun$beta1_mipiw,
                           simrun$beta1_mipticw) 
  }
  
  
  
  bias_beta1 <- round(apply(results_beta1, FUN =  mean, MARGIN = 2)  - beta1 ,3)
  names(bias_beta1) <- c("Naive", "IIW", "IPTW", "IPCW", "MIPTIW",  "MIPTICW")
  
  var_beta1 <- round(apply(results_beta1, FUN = var, MARGIN = 2), 3)
  names(var_beta1) <- c("Naive", "IIW", "IPTW", "IPCW", "MIPTIW",  "MIPTICW")
  
  mse_beta1 <- round(apply((results_beta1 - beta1)^2, FUN = mean, MARGIN = 2), 3)
  names(mse_beta1) <- c("Naive", "IIW", "IPTW", "IPCW", "MIPTIW",  "MIPTICW")
  
  
  
  
  naive_beta1 <- c(bias_beta1[1], var_beta1[1], mse_beta1[1])
  names(naive_beta1) <- c("Bias", "Var", "MSE")
  
  
  IIW_beta1 <- c(bias_beta1[2], var_beta1[2], mse_beta1[2])
  names(IIW_beta1) <- c("Bias", "Var", "MSE")
  
  
  IPTW_beta1 <- c(bias_beta1[3], var_beta1[3], mse_beta1[3])
  names(IPTW_beta1) <- c("Bias", "Var", "MSE")
  
  IPCW_beta1 <- c(bias_beta1[4], var_beta1[4], mse_beta1[4])
  names(IPCW_beta1) <- c("Bias", "Var", "MSE")
  
  MIPTIW_beta1 <- c(bias_beta1[5], var_beta1[5], mse_beta1[5])
  names(MIPTIW_beta1) <- c("Bias", "Var", "MSE")
  
  
  MIPTICW_beta1 <- c(bias_beta1[6], var_beta1[6], mse_beta1[6])
  names(MIPTICW_beta1) <- c("Bias", "Var", "MSE")
  
  
  out <- list(bias_beta1, var_beta1, mse_beta1, 
              naive_beta1,  IIW_beta1,  IPTW_beta1, IPCW_beta1,
              MIPTIW_beta1, MIPTICW_beta1, results_beta1)
  
  names(out) <- c('bias_beta1', 'var_beta1', 'mse_beta1', 
                  'naive_beta1',  'IIW_beta1',  'IPTW_beta1', 
                  "IPCW_beta1", 
                  'MIPTIW_beta1', 'MIPTICW_beta1', "fullresults_beta1")
  return(out)
  
}



simulateALLMIPI_timevar_5cov_censoring <- function(N, n, beta1, beta2, beta3, beta4, beta5,
                                                   gamma1, gamma2, gamma3, gamma4, gamma5,
                                                   kappa1, kappa2, kappa3, kappa4, kappa5,
                                                   alpha0, alpha1, alpha2, tau){
  
  # Simulates various scenarios for all possible combinations of the ___vec  parameters
  #N: number of simulation runs for each scenario
  
  
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
  # tau: maximum follow up time, independent censoring.
  # N: number of simulations
  
  
  
  
  resultsmat <- matrix(NA, nrow = 1, ncol = 16)
  fulldatalist <- list()
  fulldatalistnames <- c()
  
  i = 1
  
  #lambda_i(t) = eta_i*sqrt(t)/2*exp(gamma1G_i + gamma2W_i1(t) + gamma3W_i2 + gamma4Z_i1 + gamma5Z_i2(t))
  for(kappa2 in kappa2vec){
    for(kappa3 in kappa3vec){
      for(kappa4 in kappa4vec){
        for(kappa5 in kappa5vec){
          
          
          
          print(paste("Now on kappa2 = ", kappa2, 
                      ", kappa3 =  ", kappa3, 
                      ", kappa4  = ", kappa4,
                      ", kappa5 = ", kappa5,
                      sep = ""))
          result <- simulateResultsMIPI_timevar_5cov_censoring(N, n, beta1, beta2, beta3, beta4, beta5, 
                                                               gamma1, gamma2, gamma3, gamma4, gamma5,
                                                               kappa1, kappa2, kappa3, kappa4, kappa5,
                                                               alpha0, alpha1, alpha2, tau)
          
          resultsmat <- rbind(resultsmat, c(kappa2, kappa3, kappa4, kappa5,
                                            result$naive_beta1, result$IIW_beta1, 
                                            result$IPW_beta1, result$MIPIW_beta1))
          
          listname <- paste("fulldata_kappa2_", kappa2, "_kappa3_", kappa3, "_kappa4_", 
                            kappa4, "_kappa5_", kappa5, sep = "")
          
          fulldatalist[[i]] <- result$fullresults_beta1
          fulldatalistnames <- c(fulldatalistnames, listname)
          
          i <- i + 1
          
        }
      }
    }
  }
  
  
  #names(fulldatalist) <- fulldatalistnames
  
  
  resultsmat<- resultsmat[-1,]
  
  colnames(resultsmat) <- c( "kappa2", "kappa3", "kappa4", "kappa5", 
                             "Bias", "Var" , "MSE", 
                             "Bias", "Var" , "MSE", 
                             "Bias", "Var" , "MSE", 
                             "Bias", "Var" , "MSE")
  
  out <- list(resultsmat, fulldatalist)
  names(out) <- c("resultsmat", "fulldatalist")
  return(out)
}








