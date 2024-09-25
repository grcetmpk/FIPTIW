
# D: Treatment assignment indicator
# W: Treatment confounder (also generate time-var version W(t) for obs time model)
# Z: Observation times confounder
# Y: Time-varying outcome
# t: observation time
# phi: random effect
# epsilon: random error, function of time

# Data generation model: Y_i(t)  = (2-t) + beta1D_i + beta2[W_i - E(w_i|D_i)] + 
#                                     beta3[Z_i - E(Z_i|D_i)] + phi_i + epsilon_i(t)
#
# Marginal model: E(Y_i(t) | D_i) = E{(2-t) + beta1D_i + beta2[W_i(t) - E(w_i(t)|D_i)] + beta3(Z_i - E(Z_i|D_i)) | D_i]
#                                 = (2-t) + beta1D_i + betaE{2[W_i(t) - E(W_i(t)|D_i) | D_i] + beta3E{Z_i - E(Z_i|D_i)|D_i}
#                                 = (2-t) + beta1D_i 
#
# Observation intensity: lambda_i(t) = eta_i*exp(gamma1D_i + gamma2W_i(t) + gamma3Z_i)
# Treatment assignment: P(D_i = 1) = expit(alpha0 + alpha1W_i) [W not W(t)]
#
#
# tau: maximum follow up time
# censoring time C can either be informative or noninformative
# N: number of simulations
#


require(survival)
require(geepack)
require(nleqslv)
require(knitr)
require(kableExtra)
require(dplyr)
require(doParallel)
require(splines)
require(geeM)
require(tidyr)
require(latex2exp)
require(sandwich)

expit <- function(x){ return(exp(x)/(1+exp(x)))}



#~~~~~~~~~~~~~ Data Generation ~~~~~~~~~~~~~~~~#
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



gendata_FIPTIW_continuous <- function(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, censinform = F, 
                           eta1 = NULL, eta2 = NULL, eta3 = NULL){

  # n: number of subjects
  # tau: maximum follow up
  # censinform: set to T to generate censoring times that are informative
  #
  # Data generation model: Y_i(t)  = (2-t) + beta1D_i + beta2[W_i - E(w_i|D_i)] + 
  #                                     beta3[Z_i - E(Z_i|D_i)] + phi_i + epsilon_i(t)
  #
  # Marginal model: E(Y_i(t) | D_i) = E{(2-t) + beta1D_i + beta2[W_i - E(w_i|D_i)] + beta3(Z_i - E(Z_i|D_i)) | D_i]
  #                                 = (2-t) + beta1D_i + betaE{2[W_i - W(W_i|D_i) | D_i] + beta3E{Z_i - E(Z_i|D_i)|D_i}
  #                                 = (2-t) + beta1D_i 
  #  
  # Observation intensity: lambda_i(t) = eta_i*exp(gamma1D_i + gamma2W_i(t) + gamma3Z_i)
  # Treatment assignment: P(D_i = 1) = expit(alpha0 + alpha1W_i) 
  # Censoring: C_i ~ Unif(tau/2, tau) unless censinform = T, then C_i gen from 0.1texp(eta1D_i + eta2W_i(t) + eta3Z_i)
  # where
  #
  # D: Treatment assignment indicator
  # W: Treatment confounder
  # Z: Observation times confounder
  # Y: Time-varying outcome
  # t: observation time
  # phi: random effect
  # epsilon: random error, function of time
  #
  #
  # simdata_full is used for the estimation of the observation intensity
  # simdata_observed is used for the estimation of the IPTW weights and outcome model. 
  
  
  
  mu_Z_D0 <- 2
  var_Z_D0 <- 1 
  mu_Z_D1 <- 0
  var_Z_D1 <- 0.5
  
  var_phi <- 1.25
  var_epsilon <- 1
  
  
  
  
  # generate subjects one by one and then combine
  id = 1
  simdatafull <- data.frame(matrix(NA, ncol = 13))
  colnames(simdatafull) <- c("id", "time", "D", "W", "Wt", "Z", "cexp_Wt_D", "cexp_Z_D",  "y",  "censortime", "observed", "nevents", "censored")
  
  potentialtimes <- seq(0.01, tau, by = 0.01)
  
  while(id <= n){
    
    ## Generate Covariates

    
    # generate treatment confounder w ~ N(0,1)
    W <- runif(1, 0, 1)
    
    # generate treatment assignment
    prD <- expit(alpha0 + alpha1*W)
    D <- rbinom(1, 1, prD)
    
    # generate observation times confounder Z ~ N( mu_Z_D0, var_Z_D0) if D = 0 and 
    # N( mu_Z_D1, var_Z_D1) if D = 1
    
    Z <- ifelse(D == 0, rnorm(1, mu_Z_D0, sqrt(var_Z_D0)), rnorm(1, mu_Z_D1, sqrt(var_Z_D1)))
    
    
    # generate random effect
    phi <- rnorm(1, 0, var_phi)
    
    
    ## Generate Observation Times for all possible times
    
    simdata <- data.frame(cbind(rep(id, length(potentialtimes)), potentialtimes, rep(D, length(potentialtimes)), 
                                rep(W, length(potentialtimes)), W*log(potentialtimes), rep(Z, length(potentialtimes))))
                          
    colnames(simdata) <- c("id", "time", "D", "W", "Wt", "Z")
    
    
    ## Generate outcome
    # need conditional expectation and variance of Z | X
    
    simdata$cexp_Wt_D <- 0.5*log(simdata$time) #generated from N(0,1) independently from D
    simdata$cexp_Z_D <- ifelse(simdata$D == 0, mu_Z_D0, mu_Z_D1)
    
    
    simdata$y <- (2 - simdata$time) + beta1*simdata$D + beta2*(simdata$Wt - simdata$cexp_Wt_D) +
      beta3*(simdata$Z - simdata$cexp_Z_D) + 
      rep(phi, dim(simdata)[1]) + rnorm(dim(simdata)[1], 0, sqrt(var_epsilon))
    
    
    # simulate censoring time
    if(censinform == T){
      U <- runif(1, 0, 1)
      censortime <- sqrt(2/0.1*(-log(U)*exp(-1*eta1*D+-1*eta2*W+-1*eta3*Z)))
    }else{
      censortime <- runif(1, tau/2, tau)
    }
    
    simdata$censortime <- round(censortime,2)
    
    # calculate eta
    eta <- rgamma(1, shape = 100, scale = 0.01) #eta is for random effect with mean 1 and sd 0.1 
    #(obs times within subject correlated if sd !0)
    
    
    lambdafun <- function(t) eta*sqrt(t)/2*exp(gamma1*D + gamma2*W*log(t) + gamma3*Z)
    
    # generate observation times via thinning
    obstimes <- getObsTimesViaThinning(lambdafun, tau)
    obstimes <- round(obstimes,2)
    
    simdata$observed <- ifelse(simdata$time %in% obstimes, 1, 0)
    
    simdata$y <- ifelse(simdata$observed == 0, NA, simdata$y) #only record Y if we have a visit/observation
    
  
    # calculate the number of events for individual i
    simdata$nevents <- length(obstimes)
 
    
    # filter censored observations
    simdata$censored <- ifelse(simdata$time > simdata$censortime, 1, 0)
    simdata <- simdata %>% filter(censored == 0)
    
      
    simdatafull <- rbind(simdatafull, simdata)
      

  
    
    id = id + 1
    
  }
  
  simdata <- simdatafull[-1,] #remove empty first row
  
  
  # order dataset by ID then time
  simdata <- simdata[with(simdata, order(id, time)), ]
  
 
  simdata$time.lag <- simdata$time - 0.01

  
  
  ### OBSERVED DATA
  
  simdata_observed <- simdata %>% filter(censored == 0) %>% filter(observed == 1)
  
  
  #include a variable counting observation number, to be used for lagging time for Surv function
  simdata_observed$obsnumber <- with(simdata_observed, ave(id, id, FUN = seq_along))
  

  #check the generated data
  numevents <- summary(tapply(simdata_observed$observed, simdata_observed$id, sum)) 
  
  
  out <- list(simdata, numevents)
  names(out) <- c("simdata", "numevents")
  
  return(out) 
  
}



gendata_FIPTIW_binary <- function(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, censinform = F, 
                                      eta1 = NULL, eta2 = NULL, eta3 = NULL){
  
  mu_Z_D0 <- 2
  var_Z_D0 <- 1 
  mu_Z_D1 <- 0
  var_Z_D1 <- 0.5
  
  var_phi <- 1.25
  var_epsilon <- 1
  
  
  
  
  # generate subjects one by one and then combine
  id = 1
  simdatafull <- data.frame(matrix(NA, ncol = 16))
  colnames(simdatafull) <- c("id", "time", "D", "W", "Wt", "Z", "cexp_Wt_D", "cexp_Z_D", "cvar_Wt_D", "cvar_Z_D", "M",
                             "y",  "censortime", "observed", "nevents", "censored")
  
  potentialtimes <- seq(0.01, tau, by = 0.01)
  
  while(id <= n){
    
    ## Generate Covariates
    
    
    # generate treatment confounder w ~ N(0,1)
    W <- runif(1, 0, 1)
    
    # generate treatment assignment
    prD <- expit(alpha0 + alpha1*W)
    D <- rbinom(1, 1, prD)
    
    # generate observation times confounder Z ~ N( mu_Z_D0, var_Z_D0) if D = 0 and 
    # N( mu_Z_D1, var_Z_D1) if D = 1
    
    Z <- ifelse(D == 0, rnorm(1, mu_Z_D0, sqrt(var_Z_D0)), rnorm(1, mu_Z_D1, sqrt(var_Z_D1)))
    
    
    # generate random effect
    phi <- rnorm(1, 0, var_phi)
    
    
    ## Generate Observation Times for all possible times
    
    simdata <- data.frame(cbind(rep(id, length(potentialtimes)), potentialtimes, rep(D, length(potentialtimes)), 
                                rep(W, length(potentialtimes)), W*log(potentialtimes), rep(Z, length(potentialtimes))))
    
    colnames(simdata) <- c("id", "time", "D", "W", "Wt", "Z")
    
    
    ## Generate outcome
    # need conditional expectation and variance of Z | X
    
    simdata$cexp_Wt_D <- 0.5*log(simdata$time) #generated from N(0,1) independently from D
    simdata$cexp_Z_D <- ifelse(simdata$D == 0, mu_Z_D0, mu_Z_D1)
    simdata$cvar_Wt_D <- 1/12*(log(simdata$time)^2)
    simdata$cvar_Z_D <- ifelse(simdata$D == 0, var_Z_D0, var_Z_D1)
    
    #claculate M
    simdata$M <- sqrt(beta2^2*simdata$cvar_Wt_D + beta3^2*simdata$cvar_Z_D + var_epsilon +var_phi)/1.7
    
    simdata$y <- ifelse((2-simdata$time)*simdata$M - beta2*simdata$cexp_Wt_D - beta3*simdata$cexp_Z_D + #f0*
      beta1*simdata$M*simdata$D + beta2*simdata$Wt + beta3*simdata$Z + phi + rnorm(dim(simdata)[1], 0, var_epsilon) > 0, 1, 0)
    

    
    # simulate censoring time
    if(censinform == T){
      U <- runif(1, 0, 1)
      censortime <- sqrt(2/0.1*(-log(U)*exp(-1*eta1*D+-1*eta2*W+-1*eta3*Z)))
    }else{
      censortime <- runif(1, tau/2, tau)
    }
    
    simdata$censortime <- round(censortime,2)
    
    # calculate eta
    eta <- rgamma(1, shape = 100, scale = 0.01) #eta is for random effect with mean 1 and sd 0.1 
    #(obs times within subject correlated if sd !0)
    
    
    lambdafun <- function(t) eta*sqrt(t)/2*exp(gamma1*D + gamma2*W*log(t) + gamma3*Z)
    
    # generate observation times via thinning
    obstimes <- getObsTimesViaThinning(lambdafun, tau)
    obstimes <- round(obstimes,2)
    
    simdata$observed <- ifelse(simdata$time %in% obstimes, 1, 0)
    
    simdata$y <- ifelse(simdata$observed == 0, NA, simdata$y) #only record Y if we have a visit/observation
    
    
    # calculate the number of events for individual i
    simdata$nevents <- length(obstimes)
    
    
    # filter censored observations
    simdata$censored <- ifelse(simdata$time > simdata$censortime, 1, 0)
    simdata <- simdata %>% filter(censored == 0)
    
    
    simdatafull <- rbind(simdatafull, simdata)
    
    
    
    
    id = id + 1
    
  }
  
  simdata <- simdatafull[-1,] #remove empty first row
  
  
  # order dataset by ID then time
  simdata <- simdata[with(simdata, order(id, time)), ]
  
  
  simdata$time.lag <- simdata$time - 0.01
  
  
  
  ### OBSERVED DATA
  
  simdata_observed <- simdata %>% filter(censored == 0) %>% filter(observed == 1)
  
  
  #include a variable counting observation number, to be used for lagging time for Surv function
  simdata_observed$obsnumber <- with(simdata_observed, ave(id, id, FUN = seq_along))
  
  
  #check the generated data
  numevents <- summary(tapply(simdata_observed$observed, simdata_observed$id, sum)) 
  
  
  out <- list(simdata, numevents)
  names(out) <- c("simdata", "numevents")
  
  return(out) 
  
}





#~~~~~~~~~~~~~~~~~~ Reg. Simulation Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

simulateOneFIPTIW <- function(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                              censinform = F, eta1 = NULL, eta2 = NULL, eta3 = NULL,  usesplines = T, ncutpts = NULL){
  # Simulates one instance of the simulation, obtaining estimates for beta1 under various weighting
  # IIW uses stabilized weights
  
  
  if(outcome == "continuous"){
    
    singlerun <- gendata_FIPTIW_continuous(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau,
                                           censinform, eta1, eta2, eta3)
    
    family = gaussian(link = "identity")
    
  }else if(outcome == "binary"){
    singlerun <- gendata_FIPTIW_binary(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau,
                                           censinform, eta1, eta2, eta3)
    family = binomial(link = "logit")
  }else{
    return(print("insufficient outcome type"))
  }

  simdata <- singlerun$simdata
  numevents <- singlerun$numevents
  
  # filter to have only data when y is observed 
  simdata_obs <- simdata%>%filter(!is.na(y))
  
  
  if(usesplines == T){
    
    if(ncutpts == 2){
      #find tertiles for cubic spline to estimate the time-var intercept
      terti<-quantile(0:tau , c(0.3333, 0.66666), type = 1) 
    }else if(ncutpts == 1){
      terti<-quantile(0:tau , c(0.5), type = 1) 
    }else{
      return(print(paste("unsuitable number of cutp points for spline")))
    }
    
  }
  
 

  
  #### NAIVE ESTIMATION (beta1 and se(beta1))
 
  if(usesplines == T){
    naivemod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id))
    beta1_naive <- naivemod$coef[2,1]
    se_beta1_naive <- naivemod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    naive_95CI_lower <- beta1_naive - 1.96*se_beta1_naive
    naive_95CI_upper <- beta1_naive + 1.96*se_beta1_naive
    
    beta1_naive_covered <- ifelse(naive_95CI_lower <= beta1 && beta1 <= naive_95CI_upper, TRUE, FALSE)
    
  }else{
    naivemod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id))
    beta1_naive <- naivemod$coef[[1]]
    se_beta1_naive <- naivemod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    naive_95CI_lower <- beta1_naive - 1.96*se_beta1_naive
    naive_95CI_upper <- beta1_naive + 1.96*se_beta1_naive
    
    beta1_naive_covered <- ifelse(naive_95CI_lower <= beta1 && beta1 <= naive_95CI_upper, TRUE, FALSE)
    
  }
  
  
  
  #### IIW ESTIMATION
  
  # IIW Weights
  gamma.hat <- coxph(Surv(time.lag, time, observed) ~ D + Wt + Z - 1, 
                     data = simdata)$coef
  delta.hat <- coxph(Surv(time.lag, time, observed) ~ D - 1, data = simdata)$coef
  

  iiw <- exp(cbind(simdata_obs$D)%*%delta.hat)/exp(cbind(simdata_obs$D, 
                                    simdata_obs$Wt, simdata_obs$Z)%*%gamma.hat)
  
  if(usesplines == T){
    iiwmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id, weights = iiw))
    beta1_iiw <- iiwmod$coef[2,1]
    se_beta1_iiw <- iiwmod$coef[2,1] #robust SE
    
    
    #see if CI contains true value or not
    iiw_95CI_lower <- beta1_iiw - 1.96*se_beta1_iiw
    iiw_95CI_upper <- beta1_iiw + 1.96*se_beta1_iiw
    
    beta1_iiw_covered <- ifelse(iiw_95CI_lower <= beta1 && beta1 <= iiw_95CI_upper, TRUE, FALSE)
    
  }else{
    iiwmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id, weights = iiw))
    beta1_iiw <- iiwmod$coef[[1]]
    se_beta1_iiw <- iiwmod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    iiw_95CI_lower <- beta1_iiw - 1.96*se_beta1_iiw
    iiw_95CI_upper <- beta1_iiw + 1.96*se_beta1_iiw
    
    beta1_iiw_covered <- ifelse(iiw_95CI_lower <= beta1 && beta1 <= iiw_95CI_upper, TRUE, FALSE)
    
  }
  
 
  
  
  
  
  #####  IPTW estimation
  
  psmod <- glm(D ~ W, family = binomial(link = "logit"), data = simdata)
  ps <- expit(predict(psmod, newdata = simdata_obs))
  prDmod <- glm(D ~ 1, family = binomial(link = "logit"), data = simdata)
  prD <- expit(summary(prDmod)$coef[1,1])
  
  ipw <- prD/ps*simdata_obs$D + (1-prD)/(1-ps)*(1-simdata_obs$D)
  #ipw <- 1/ps*simdata_obs$D+1/(1-ps)*(1-simdata_obs$D)

  if(usesplines == T){
    ipwmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id, weights = ipw))
    beta1_ipw <- ipwmod$coef[2,1]
    se_beta1_ipw <- ipwmod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    ipw_95CI_lower <- beta1_ipw - 1.96*se_beta1_ipw
    ipw_95CI_upper <- beta1_ipw + 1.96*se_beta1_ipw
    
    beta1_ipw_covered <- ifelse(ipw_95CI_lower <= beta1 && beta1 <= ipw_95CI_upper, TRUE, FALSE)
    
    
  }else{
    ipwmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id, weights = ipw))
    beta1_ipw <- ipwmod$coef[[1]]
    se_beta1_ipw <- ipwmod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    ipw_95CI_lower <- beta1_ipw - 1.96*se_beta1_ipw
    ipw_95CI_upper <- beta1_ipw + 1.96*se_beta1_ipw
    
    beta1_ipw_covered <- ifelse(ipw_95CI_lower <= beta1 && beta1 <= ipw_95CI_upper, TRUE, FALSE)
    
    
  }
  


  ####  FIPTIW estimation
  fiptiw <- ipw*iiw

  if(usesplines == T){
    fiptiwmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id, weights = fiptiw))
    beta1_fiptiw <- fiptiwmod$coef[2,1]
    se_beta1_fiptiw <- fiptiwmod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    fiptiw_95CI_lower <- beta1_fiptiw - 1.96*se_beta1_fiptiw
    fiptiw_95CI_upper <- beta1_fiptiw + 1.96*se_beta1_fiptiw
    
    beta1_fiptiw_covered <- ifelse(fiptiw_95CI_lower <= beta1 && beta1 <= fiptiw_95CI_upper, TRUE, FALSE)
    
  }else{
    fiptiwmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id, weights = fiptiw))
    beta1_fiptiw <- fiptiwmod$coef[[1]]
    se_beta1_fiptiw <- fiptiwmod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    fiptiw_95CI_lower <- beta1_fiptiw - 1.96*se_beta1_fiptiw
    fiptiw_95CI_upper <- beta1_fiptiw + 1.96*se_beta1_fiptiw
    
    beta1_fiptiw_covered <- ifelse(fiptiw_95CI_lower <= beta1 && beta1 <= fiptiw_95CI_upper, TRUE, FALSE)
    
  }
  

  
  out <- list(beta1,  beta1_naive, se_beta1_naive, beta1_naive_covered, 
              beta1_iiw, se_beta1_iiw, beta1_iiw_covered,
              beta1_ipw, se_beta1_ipw, beta1_ipw_covered, 
              beta1_fiptiw, se_beta1_fiptiw, beta1_fiptiw_covered)
  names(out) <- c('beta1',  'beta1_naive', 'se_beta1_naive', 'beta1_naive_covered', 
                  'beta1_iiw', 'se_beta1_iiw', 'beta1_iiw_covered',
                  'beta1_ipw', 'se_beta1_ipw', 'beta1_ipw_covered', 
                  'beta1_fiptiw', 'se_beta1_fiptiw', 'beta1_fiptiw_covered')
  return(out)
}





simulateResultsFIPTIW <-  function(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
                                alpha0, alpha1, tau, outcome,
                                censinform = F, eta1 = NULL, eta2 = NULL, eta3 = NULL,
                                usesplines = T, ncutpts =  NULL, outputfulldatalist = T,
                                inParallel = F, nclusters = NULL){


  if(inParallel == TRUE){
    registerDoParallel(nclusters)
    
    results<- foreach(i = 1:N, .combine = rbind, .export = c("expit", "simulateOneFIPTIW", 
                                                  "geeglm", "coxph", "Surv", "gendata_FIPTIW_continuous",
                                                  "gendata_FIPTIW_binary", "bs", 
                                                  "getObsTimesViaThinning"), 
                      .packages = c("dplyr", "geepack")) %dopar% {
      simrun <- simulateOneFIPTIW(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                                              censinform, eta1, eta2, eta3, usesplines, ncutpts)
      
      c(simrun$beta1_naive, simrun$se_beta1_naive, simrun$beta1_naive_covered,
        simrun$beta1_iiw, simrun$se_beta1_iiw, simrun$beta1_iiw_covered,
        simrun$beta1_ipw, simrun$se_beta1_ipw, simrun$beta1_ipw_covered,
        simrun$beta1_fiptiw, simrun$se_beta1_fiptiw, simrun$beta1_fiptiw_covered)
    }
    
    
  } else {
    
    results <- matrix(data = NA, nrow = N, ncol = 12)
    
    for(i in 1:N){
      if(i%%100 == 0){print(paste("Sim num:", i))}
      simrun <- simulateOneFIPTIW(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                                  censinform, eta1, eta2, eta3, usesplines, ncutpts)
      results[i,] <- c(simrun$beta1_naive, simrun$se_beta1_naive, as.numeric(simrun$beta1_naive_covered),
                             simrun$beta1_iiw, simrun$se_beta1_iiw, simrun$beta1_iiw_covered,
                             simrun$beta1_ipw, simrun$se_beta1_ipw, simrun$beta1_ipw_covered,
                             simrun$beta1_fiptiw, simrun$se_beta1_fiptiw, simrun$beta1_fiptiw_covered)
    }
  }

  results_beta1 <- results[, c(1, 4, 7, 10)]
  
  
  colnames(results_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW")
  
  bias_beta1 <- round(apply(results_beta1, FUN =  mean, MARGIN = 2)  - beta1 ,3)  
  names(bias_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW")
  
  var_beta1 <- round(apply(results_beta1, FUN = var, MARGIN = 2), 3)
  names(var_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW")
  
  mse_beta1 <- round(apply((results_beta1 - beta1)^2, FUN = mean, MARGIN = 2), 3)
  names(mse_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW")
  
  
  avgSE_beta1 <- results[, c(2, 5, 8, 11)]
  avgSE_beta1 <- round(apply(avgSE_beta1, FUN =  mean, MARGIN = 2), 3)
  names(avgSE_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW")
  
  
  
  coverage_beta1 <- results[, c(3, 6, 9, 12)]
  coverage_beta1 <- round(apply(coverage_beta1, FUN =  mean, MARGIN = 2), 3)
  names(coverage_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW")
  
  
  
  ## ## ##
  
  naive_beta1 <- c(bias_beta1[1], var_beta1[1], mse_beta1[1], avgSE_beta1[1], coverage_beta1[1])
  names(naive_beta1) <- c("Bias", "Var", "MSE", "Avg SE", "Coverage")
  
  IIW_beta1 <- c(bias_beta1[2], var_beta1[2], mse_beta1[2], avgSE_beta1[2], coverage_beta1[2])
  names(IIW_beta1) <- c("Bias", "Var", "MSE", "Avg SE", "Coverage")
  
  IPW_beta1 <- c(bias_beta1[3], var_beta1[3], mse_beta1[3], avgSE_beta1[3], coverage_beta1[3])
  names(IPW_beta1) <- c("Bias", "Var", "MSE", "Avg SE", "Coverage")
  
  FIPTIW_beta1 <- c(bias_beta1[4], var_beta1[4], mse_beta1[4], avgSE_beta1[4], coverage_beta1[4])
  names(FIPTIW_beta1) <- c("Bias", "Var", "MSE", "Avg SE", "Coverage")
  
  
  
  if(outputfulldatalist == TRUE){
    out <- list(bias_beta1, var_beta1, mse_beta1, avgSE_beta1, coverage_beta1,
                naive_beta1,  IIW_beta1, IPW_beta1, 
                FIPTIW_beta1, results_beta1)
    
    names(out) <- c('bias_beta1', 'var_beta1', 'mse_beta1', 'avgSE_beta1', 'coverage_beta1',
                   'naive_beta1',  'IIW_beta1', 'IPW_beta1', 
                   'FIPTIW_beta1', 'results_beta1')
  }else{
    out <- list(bias_beta1, var_beta1, mse_beta1, avgSE_beta1, coverage_beta1,
                naive_beta1,  IIW_beta1, IPW_beta1, 
                FIPTIW_beta1)
    
    names(out) <- c('bias_beta1', 'var_beta1', 'mse_beta1', 'avgSE_beta1', 'coverage_beta1',
                    'naive_beta1',  'IIW_beta1', 'IPW_beta1', 
                    'FIPTIW_beta1')
  }

  return(out)
  
}


simulateALLFIPTIW <- function(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
                                    alpha0, alpha1vec, tau, outcome, 
                                    censinform = F, eta1 = NULL, eta2 = NULL, eta3= NULL,
                                    usesplines = T, ncutpts = NULL, outputfulldatalist = F,
                                    inParallel = F, nclusters = NULL){

  
  
  resultsmat <- matrix(NA, nrow = 1, ncol = 8)
  resultsmat[1,] <- c("gamma1", "gamma2", "alpha1", "", "Naive", "IIW", "IPTW", "FIPTIW")
  fulldatalist <- list()
  fulldatalistnames <- c()
  
  i = 1
  
  
  for(gamma2 in gamma2vec){
    for(gamma3 in gamma3vec){
      for(alpha1 in alpha1vec){

        
        print(paste("Now on gamma2 = ", gamma2, ", gamma3 =", gamma3, ", alpha1 =  ", alpha1, sep = ""))
        result <- simulateResultsFIPTIW(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
                                        alpha0, alpha1, tau, outcome,
                                        censinform, eta1, eta2, eta3,
                                        usesplines, ncutpts, outputfulldatalist,
                                        inParallel, nclusters)
        
        thisobs <- rbind(c(gamma2, gamma3, alpha1, "Bias:", result$naive_beta1[[1]], result$IIW_beta1[[1]], result$IPW_beta1[[1]], result$FIPTIW_beta1[[1]]),
                         c("", "", "", "MSE:", result$naive_beta1[[3]], result$IIW_beta1[[3]], result$IPW_beta1[[3]], result$FIPTIW_beta1[[3]]),
                         c("", "", "", "Avg SE:", result$naive_beta1[[4]], result$IIW_beta1[[4]], result$IPW_beta1[[4]], result$FIPTIW_beta1[[4]]),
                         c("", "", "", "Coverage:", result$naive_beta1[[5]], result$IIW_beta1[[5]], result$IPW_beta1[[5]], result$FIPTIW_beta1[[5]]))
        
        resultsmat <- rbind(resultsmat, thisobs)
        
        
        
        listname <- paste("fulldata_gamma2_",gamma2, "_gamma3_", gamma3, "_alpha1_", alpha1, sep = "")
        
        if(outputfulldatalist == TRUE){
          fulldatalist[[i]] <- result$results_beta1
          fulldatalistnames <- c(fulldatalistnames, listname)
        }

        
        i <- i + 1
        
      }
    }
  }
  
  if(outputfulldatalist == TRUE){
    names(fulldatalist) <- fulldatalistnames
  }

  
  
  resultsmat<- resultsmat[-1,]
  # resultsmat_stabilized <- resultsmat_stabilized[-1,]
  
  colnames(resultsmat) <- c("gamma2",  "gamma3", "alpha1", "", "Naive", "IIW" , "IPTW", "FIPTIW")

  if(outputfulldatalist == TRUE){
    out <- list(resultsmat, fulldatalist)
    names(out) <- c("resultsmat", "fulldatalist")
  }else{
    out <- list(resultsmat)
    names(out) <- c("resultsmat")
  }

  return(out)
}






#########################################
############ IPCW simulation ############
#########################################


simulateOneFIPTICW <- function(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                               censinform = F, eta1 = NULL, eta2 = NULL, eta3 = NULL,  usesplines = F, ncutpts = NULL){
  # Simulates one instance of the simulation, obtaining estimates for beta1 under various weighting
  # IIW uses stabilized weights
  
  
  
  # Simulates one instance of the simulation, obtaining estimates for beta1 under various weighting
  # IIW uses stabilized weights
  
  
  if(outcome == "continuous"){
    
    singlerun <- gendata_FIPTIW_continuous(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau,
                                           censinform, eta1, eta2, eta3)
    
    family = gaussian(link = "identity")
    
  }else if(outcome == "binary"){
    singlerun <- gendata_FIPTIW_binary(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau,
                                       censinform, eta1, eta2, eta3)
    family = binomial(link = "logit")
  }else{
    return(print("insufficient outcome type"))
  }
  
  simdata <- singlerun$simdata
  numevents <- singlerun$numevents
  
  # filter to have only data when y is observed 
  simdata_obs <- simdata%>%filter(!is.na(y))
  simdata$censorednow <- ifelse(simdata$time >= simdata$censortime, 1, 0)
  
  
  if(usesplines == T){
    
    if(ncutpts == 2){
      #find tertiles for cubic spline to estimate the time-var intercept
      terti<-quantile(0:tau , c(0.3333, 0.66666), type = 1) 
    }else if(ncutpts == 1){
      terti<-quantile(0:tau , c(0.5), type = 1) 
    }else{
      return(print(paste("unsuitable number of cutp points for spline")))
    }
    
  }
  
  
  
  
  #### NAIVE ESTIMATION (beta1 and se(beta1))
  
  if(usesplines == T){
    naivemod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id))
    beta1_naive <- naivemod$coef[2,1]
    se_beta1_naive <- naivemod$coef[2,2] #robust SE
    
  }else{
    naivemod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id))
    beta1_naive <- naivemod$coef[[1]]
    se_beta1_naive <- naivemod$coef[[2]] #robust SE
    
  }

  
  #see if CI contains true value or not
  naive_95CI_lower <- beta1_naive - 1.96*se_beta1_naive
  naive_95CI_upper <- beta1_naive + 1.96*se_beta1_naive
  
  beta1_naive_covered <- ifelse(naive_95CI_lower <= beta1 && beta1 <= naive_95CI_upper, TRUE, FALSE)
  
  
  
  
  #### IIW ESTIMATION
  
  # IIW Weights
  gamma.hat <- coxph(Surv(time.lag, time, observed) ~ D + Wt + Z - 1, 
                     data = simdata)$coef
  delta.hat <- coxph(Surv(time.lag, time, observed) ~ D - 1, data = simdata)$coef
  
  
  iiw <- exp(cbind(simdata_obs$D)%*%delta.hat)/exp(cbind(simdata_obs$D, 
                                                         simdata_obs$Wt, simdata_obs$Z)%*%gamma.hat)
  
  if(usesplines == T){
    iiwmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id, weights = iiw))
    beta1_iiw <- iiwmod$coef[2,1]
    se_beta1_iiw <- iiwmod$coef[2,2] #robust SE
    
    
  }else{
    iiwmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id, weights = iiw))
    beta1_iiw <- iiwmod$coef[[1]]
    se_beta1_iiw <- iiwmod$coef[[2]] #robust SE
    
    
  }
  

  #see if CI contains true value or not
  iiw_95CI_lower <- beta1_iiw - 1.96*se_beta1_iiw
  iiw_95CI_upper <- beta1_iiw + 1.96*se_beta1_iiw
  
  beta1_iiw_covered <- ifelse(iiw_95CI_lower <= beta1 && beta1 <= iiw_95CI_upper, TRUE, FALSE)
  
  
  
  
  
  #####  IPTW estimation
  
  psmod <- glm(D ~ W, family = binomial(link = "logit"), data = simdata)
  ps <- expit(predict(psmod, newdata = simdata_obs))
  prDmod <- glm(D ~ 1, family = binomial(link = "logit"), data = simdata)
  prD <- expit(summary(prDmod)$coef[1,1])
  
  ipw <- prD/ps*simdata_obs$D + (1-prD)/(1-ps)*(1-simdata_obs$D)
  #ipw <- 1/ps*simdata_obs$D+1/(1-ps)*(1-simdata_obs$D)
  
  if(usesplines == T){
    ipwmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id, weights = ipw))
    beta1_ipw <- ipwmod$coef[2,1]
    se_beta1_ipw <- ipwmod$coef[2,2] #robust SE
    
  }else{
    ipwmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id, weights = ipw))
    beta1_ipw <- ipwmod$coef[[1]]
    se_beta1_ipw <- ipwmod$coef[[2]] #robust SE
    
  }
  

  
  #see if CI contains true value or not
  ipw_95CI_lower <- beta1_ipw - 1.96*se_beta1_ipw
  ipw_95CI_upper <- beta1_ipw + 1.96*se_beta1_ipw
  
  beta1_ipw_covered <- ifelse(ipw_95CI_lower <= beta1 && beta1 <= ipw_95CI_upper, TRUE, FALSE)
  
  
  
  ####  FIPTIW estimation
  fiptiw <- ipw*iiw
  
  if(usesplines == T){
    fiptiwmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id, weights = fiptiw))
    beta1_fiptiw <- fiptiwmod$coef[2,1]
    se_beta1_fiptiw <- fiptiwmod$coef[2,2] #robust SE
  }else{
    fiptiwmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id, weights = fiptiw))
    beta1_fiptiw <- fiptiwmod$coef[[1]]
    se_beta1_fiptiw <- fiptiwmod$coef[[2]] #robust SE
  }
  

  
  
  #see if CI contains true value or not
  fiptiw_95CI_lower <- beta1_fiptiw - 1.96*se_beta1_fiptiw
  fiptiw_95CI_upper <- beta1_fiptiw + 1.96*se_beta1_fiptiw
  
  beta1_fiptiw_covered <- ifelse(fiptiw_95CI_lower <= beta1 && beta1 <= fiptiw_95CI_upper, TRUE, FALSE)
  


  # IPCW
  
  ipcwmod <- coxph(Surv(time.lag, time, censorednow) ~ D + W + Z - 1, data = simdata)
  ipcwfit <- summary(survfit(ipcwmod), times = simdata_obs$time.lag)
  ipcw <- 1/ipcwfit$surv
  
  
  #FIPTICW
  fipticw <- ipw*iiw*ipcw
  
  
  if(usesplines == T){
    fipticwmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id, weights = fipticw))
    
    beta1_fipticw <- fipticwmod$coef[2,1]
    se_beta1_fipticw <- fipticwmod$coef[2,2] #robust SE
    
  }else{
    fipticwmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id, weights = fipticw))
    
    beta1_fipticw <- fipticwmod$coef[[1]]
    se_beta1_fipticw <- fipticwmod$coef[[2]] #robust SE
  }

  
  
  #see if CI contains true value or not
  fipticw_95CI_lower <- beta1_fipticw - 1.96*se_beta1_fipticw
  fipticw_95CI_upper <- beta1_fipticw + 1.96*se_beta1_fipticw
  
  beta1_fipticw_covered <- ifelse(fipticw_95CI_lower <= beta1 && beta1 <= fipticw_95CI_upper, TRUE, FALSE)
  
  
  
  
  
  out <- list(beta1,  beta1_naive, se_beta1_naive, beta1_naive_covered, 
                     beta1_iiw, se_beta1_iiw, beta1_iiw_covered,
                     beta1_ipw, se_beta1_ipw, beta1_ipw_covered, 
                     beta1_fiptiw, se_beta1_fiptiw, beta1_fiptiw_covered,
                     beta1_fipticw, se_beta1_fipticw, beta1_fipticw_covered)
  
  names(out) <- c('beta1',  'beta1_naive', 'se_beta1_naive', 'beta1_naive_covered', 
                  'beta1_iiw', 'se_beta1_iiw', 'beta1_iiw_covered',
                  'beta1_ipw', 'se_beta1_ipw', 'beta1_ipw_covered', 
                  'beta1_fiptiw', 'se_beta1_fiptiw', 'beta1_fiptiw_covered',
                  'beta1_fipticw', 'se_beta1_fipticw', 'beta1_fipticw_covered')

  return(out)
}


simulateResultsFIPTICW <-  function(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
                                    alpha0, alpha1, tau, outcome,
                                    censinform = F, eta1 = NULL, eta2 = NULL, eta3 = NULL,
                                    usesplines = F, ncutpts =  NULL, outputfulldatalist = T,
                                    inParallel = F, nclusters = NULL){
  # Simulates N instances of the given scheme 
  # computes in parallel if inParallel = F
  #
  # N: number of simulation runs
  # n: number of subjects
  # tau: maximum follow up
  #
  
  if(inParallel == TRUE){
    registerDoParallel(nclusters)
    
    results<- foreach(i = 1:N, .combine = rbind, 
                      .export = c("expit", "simulateOneFIPTIW_continuous", "simulateOneFIPTIW_binary",
                      "geeglm", "coxph", "Surv", "gendata_FIPTIW_continous","gendata_FIPTIW_binary", "bs", "getObsTimesViaThinning"), 
                      .packages = c("dplyr")) %dopar% {
                      simrun <- simulateOneFIPTICW(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                                                               censinform, eta1, eta2, eta3,  usesplines, ncutpts)
                      
                     c(simrun$beta1_naive, simrun$se_beta1_naive, simrun$beta1_naive_covered,
                     simrun$beta1_iiw, simrun$se_beta1_iiw, simrun$beta1_iiw_covered,
                     simrun$beta1_ipw, simrun$se_beta1_ipw, simrun$beta1_ipw_covered,
                     simrun$beta1_fiptiw, simrun$se_beta1_fiptiw, simrun$beta1_fiptiw_covered,
                    simrun$beta1_fipticw, simrun$se_beta1_fipticw, simrun$beta1_fipticw_covered)
                                                             }
    
    
  } else {
    
    results <- matrix(data = NA, nrow = N, ncol = 15)
    
    for(i in 1:N){
      if(i%%100 == 0){print(paste("Sim num:", i))}
      simrun <- simulateOneFIPTICW(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                                   censinform, eta1, eta2, eta3,  usesplines, ncutpts)
      results[i,] <- c(simrun$beta1_naive, simrun$se_beta1_naive, simrun$beta1_naive_covered,
                       simrun$beta1_iiw, simrun$se_beta1_iiw, simrun$beta1_iiw_covered,
                       simrun$beta1_ipw, simrun$se_beta1_ipw, simrun$beta1_ipw_covered,
                       simrun$beta1_fiptiw, simrun$se_beta1_fiptiw, simrun$beta1_fiptiw_covered,
                       simrun$beta1_fipticw, simrun$se_beta1_fipticw, simrun$beta1_fipticw_covered
                       )
    }
  }
  
  fullresults_beta1 <- results[, c(1, 4, 7, 10, 13)]
  
  colnames(fullresults_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW", "FIPTICW")
  
  bias_beta1 <- round(apply(fullresults_beta1, FUN =  mean, MARGIN = 2)  - beta1 ,3)  
  names(bias_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW", "FIPTICW")
  
  var_beta1 <- round(apply(fullresults_beta1, FUN = var, MARGIN = 2), 3)
  names(var_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW", "FIPTICW")
  
  mse_beta1 <- round(apply((fullresults_beta1 - beta1)^2, FUN = mean, MARGIN = 2), 3)
  names(mse_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW", "FIPTICW")
  
  
  avgSE_beta1 <- results[, c(2, 5, 8, 11, 14)]
  avgSE_beta1 <- round(apply(avgSE_beta1, FUN =  mean, MARGIN = 2), 3)
  names(avgSE_beta1) <-c("Naive", "IIW", "IPW", "FIPTIW", "FIPTICW")
  
  
  
  coverage_beta1 <- results[, c(3, 6, 9, 12, 15)]
  coverage_beta1 <- round(apply(coverage_beta1, FUN =  mean, MARGIN = 2), 3)
  names(coverage_beta1) <- c("Naive", "IIW", "IPW", "FIPTIW", "FIPTICW")
  
  
  
  ## ## ##
  
  naive_beta1 <- c(bias_beta1[1], var_beta1[1], mse_beta1[1], avgSE_beta1[1], coverage_beta1[1])
  names(naive_beta1) <- c("Bias", "Var", "MSE", "Avg SE", "Coverage")
  
  IIW_beta1 <- c(bias_beta1[2], var_beta1[2], mse_beta1[2], avgSE_beta1[2], coverage_beta1[2])
  names(IIW_beta1) <- c("Bias", "Var", "MSE", "Avg SE", "Coverage")
  
  IPW_beta1 <- c(bias_beta1[3], var_beta1[3], mse_beta1[3], avgSE_beta1[3], coverage_beta1[3])
  names(IPW_beta1) <- c("Bias", "Var", "MSE", "Avg SE", "Coverage")
  
  FIPTIW_beta1 <- c(bias_beta1[4], var_beta1[4], mse_beta1[4], avgSE_beta1[4], coverage_beta1[4])
  names(FIPTIW_beta1) <- c("Bias", "Var", "MSE", "Avg SE", "Coverage")
  
  FIPTICW_beta1 <- c(bias_beta1[5], var_beta1[5], mse_beta1[5], avgSE_beta1[5], coverage_beta1[5])
  names(FIPTICW_beta1) <- c("Bias", "Var", "MSE", "Avg SE", "Coverage")
  
  
  if(outputfulldatalist == TRUE){
    out <- list(bias_beta1, var_beta1, mse_beta1, avgSE_beta1, coverage_beta1,
                naive_beta1,  IIW_beta1, IPW_beta1, 
                FIPTIW_beta1, FIPTICW_beta1, fullresults_beta1)
    
    names(out) <- c('bias_beta1', 'var_beta1', 'mse_beta1', 'avgSE_beta1', 'coverage_beta1',
                    'naive_beta1',  'IIW_beta1', 'IPW_beta1', 
                    'FIPTIW_beta1', 'FIPTICW_beta1', 'fullresults_beta1')
  }else{
    out <- list(bias_beta1, var_beta1, mse_beta1, avgSE_beta1, coverage_beta1,
                naive_beta1,  IIW_beta1, IPW_beta1, 
                FIPTIW_beta1, FIPTICW_beta1)
    
    names(out) <- c('bias_beta1', 'var_beta1', 'mse_beta1', 'avgSE_beta1', 'coverage_beta1',
                    'naive_beta1',  'IIW_beta1', 'IPW_beta1', 
                    'FIPTIW_beta1', 'FIPTICW_beta1')
  }
  
  return(out)

}


simulateALLFIPTICW_CENS<- function(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
                                   alpha0, alpha1, tau, outcome,
                                   censinform = F, eta1vec, eta2vec, eta3vec,
                                   usesplines = F, ncutpts =  NULL, outputfulldatalist = T,
                                   inParallel = F, nclusters = NULL){
  #N: number of simulation runs
  #n: vector of sample sizes
  #beta1: coefficient for Xi(t) in logistic outcome model
  #beta2: vector of coefficients to consider for outcome generation model
  #gamma1, gamma2 parameters for intensity for Xi(t) and Zi, respectively
  #tau: maximum follow-up time
  
  #This function aggregates simulation results for varying n and beta2
  
  
  
  resultsmat <- matrix(NA, nrow = 1, ncol = 28)
  fulldatalist <- data.frame(matrix(NA, ncol = 8))
  colnames(fulldatalist) <- c("eta1", "eta2", "eta3", "Naive", "IIW", "IPTW", "FIPTIW", "FIPTICW")
  

  
  i = 1
  
  
  for(eta1 in eta1vec){
    for(eta2 in eta2vec){
      for(eta3 in eta3vec){
        
        
        print(paste("Now on eta1 = ", eta1, ", eta2 =", eta2, ", eta3 =  ", eta3, sep = ""))
        result <- simulateResultsFIPTICW(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3,
                                         alpha0, alpha1, tau, outputfulldatalist, censinform,
                                         eta1, eta2, eta3, ncutpts = ncutpts, outcome = outcome)
        
        resultsmat <- rbind(resultsmat, c(eta1, eta2, eta3,
                                          result$naive_beta1, result$IIW_beta1,
                                          result$IPW_beta1, result$FIPTIW_beta1, result$FIPTICW_beta1))
        
        
        listname <- paste("fulldata_eta1_",eta1, "_eta2_", eta2, "_eta3_", eta3, sep = "")
        
        if(outputfulldatalist == TRUE){
          fulldata <- cbind(rep(eta1, N), rep(eta2, N), rep(eta3, N), result$fullresults_beta1)
          colnames(fulldata) <- c("eta1", "eta2", "eta3", "Naive", "IIW", "IPTW", "FIPTIW", "FIPTICW")
          
          fulldatalist <- rbind(fulldatalist, fulldata)
        }
        
        
        i <- i + 1
        
      }
    }
  }
  

  
  resultsmat<- resultsmat[-1,]
  
  colnames(resultsmat) <- c( "eta1",  "eta2", "eta3", "Bias", "Var", "MSE", "AvgSE", "Coverage",
                             "Bias", "Var", "MSE", "AvgSE", "Coverage",
                             "Bias", "Var", "MSE", "AvgSE", "Coverage",
                             "Bias", "Var", "MSE", "AvgSE", "Coverage",
                             "Bias", "Var", "MSE", "AvgSE", "Coverage")
  
  if(outputfulldatalist == TRUE){
    fulldatalist <- fulldatalist[-1, ] #remove NA row
    out <- list(resultsmat, fulldatalist)
    names(out) <- c("resultsmat", "fulldatalist")
  }else{
    out <- list(resultsmat)
    names(out) <- c("resultsmat")
  }
  
  return(out)
}




makeBoxPlotsATE <- function(atemat, eta2vec){
  
  atemat <- gather(atemat, "Method", "ATE", -c(eta1, eta2, eta3))
  atemat$Method <- factor(atemat$Method, levels = c("Naive", "IIW", "IPTW", "FIPTIW", "FIPTICW"))
  
  ymin <- min(atemat$ATE)
  ymax <- max(atemat$ATE)
  # ymin <- -2
  # ymax <- 5
  
  
  #setup plot for legend
  plot <- ggplot(atemat, aes(x = factor(eta3), y = ATE, fill = factor(Method))) +
    geom_boxplot(width = 0.5)+
    theme(legend.position = "right", legend.text=element_text(size=16), legend.title=element_text(size=18, face = "bold")) +
    scale_fill_discrete(name = "Weighting Method")

  suppressWarnings(legend <- cowplot::get_legend(plot))
  
  
  
  # plot by eta2
  
  for(i in 1:length(eta2vec)){
    
    eta2 <- eta2vec[i]
    
    currentatemat <- atemat[atemat$eta2==eta2, ]
    currentatemat$Method <- factor(currentatemat$Method, levels = c("Naive", "IIW", "IPTW", "FIPTIW", "FIPTICW"))
    
    plot <- ggplot(currentatemat, aes(x = factor(eta3), y = ATE, fill = factor(Method))) +
      geom_boxplot(width = 0.5)+
      theme(legend.position = "none", axis.text=element_text(size=14),
            axis.title=element_text(size=16,face="bold")) +
      scale_fill_discrete(name = "Weighting Method") +
      xlab(TeX("$\\eta_3$")) + 
      ylab(TeX("$\\beta_1$ (ATE)")) + 
      ggtitle(TeX(paste0("$\\eta_2=$", eta2))) +
      geom_hline(yintercept = 0.5, col = "red", linetype = "dashed") + 
      ylim(ymin, ymax)
    
    assign(paste0("p", i), plot)
    

  }
  
  
  finalplot <- ggarrange(p1, p2, p3, legend, ncol = 2, nrow = 2)
  return(finalplot)
  
}







#~~~~~~~~~~~~~~~~~ Variable Selection Function (for IIW only) ~~~~~~~~~~~~~~~~~~#


simulateOneVarSel <- function(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                              usesplines = F, ncutpts = NULL){

  
  
  if(outcome == "continuous"){
    
    singlerun <- gendata_FIPTIW_continuous(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau,
                                           censinform = F, eta1 = NULL, eta2 = NULL, eta3= NULL)
    
    family = gaussian(link = "identity")
    
  }else if(outcome == "binary"){
    singlerun <- gendata_FIPTIW_binary(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau,
                                       censinform = F, eta1 = NULL, eta2 = NULL, eta3= NULL)
    family = binomial(link = "logit")
  }else{
    return(print("insufficient outcome type"))
  }
  
  simdata <- singlerun$simdata
  numevents <- singlerun$numevents
  
  # filter to have only data when y is observed 
  simdata_obs <- simdata%>%filter(!is.na(y))
  
  
  if(usesplines == T){
    
    if(ncutpts == 2){
      #find tertiles for cubic spline to estimate the time-var intercept
      terti<-quantile(0:tau , c(0.3333, 0.66666), type = 1) 
    }else if(ncutpts == 1){
      terti<-quantile(0:tau , c(0.5), type = 1) 
    }else{
      return(print(paste("unsuitable number of cutp points for spline")))
    }
    
  }
  
  
  
  
  #### NAIVE ESTIMATION (beta1 and se(beta1))
  
  if(usesplines == T){
    naivemod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id))
    
    beta1_naive <- naivemod$coef[2,1]
    se_beta1_naive <- naivemod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    naive_95CI_lower <- beta1_naive - 1.96*se_beta1_naive
    naive_95CI_upper <- beta1_naive + 1.96*se_beta1_naive
    
    beta1_naive_covered <- ifelse(naive_95CI_lower <= beta1 && beta1 <= naive_95CI_upper, TRUE, FALSE)
    
  }else{
    naivemod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id))
    
    beta1_naive <- naivemod$coef[[1]]
    se_beta1_naive <- naivemod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    naive_95CI_lower <- beta1_naive - 1.96*se_beta1_naive
    naive_95CI_upper <- beta1_naive + 1.96*se_beta1_naive
    
    beta1_naive_covered <- ifelse(naive_95CI_lower <= beta1 && beta1 <= naive_95CI_upper, TRUE, FALSE)
    
  }

  
  
  

  ## estimate delta hat (used for all IIW models)
  delta.hat <-coxph(Surv(time.lag, time, observed) ~ D - 1, data = simdata)$coef

  
  ## estimate gamma hat, intensity, beta1
  
  #####  intensity using D, W, Z
  gamma.hat_DWZ <-   gamma.hat <- coxph(Surv(time.lag, time, observed) ~ D + Wt + Z - 1,  data = simdata)$coef
  iiw_DWZ <- exp(cbind(simdata_obs$D)%*%delta.hat)/exp(cbind(simdata_obs$D, simdata_obs$Wt, simdata_obs$Z)%*%gamma.hat)

  if(usesplines == T){
    iiw_DWZmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, 
                                 family = family, id = id, weights = iiw_DWZ))
    beta1_iiw_DWZ <- iiw_DWZmod$coef[2,1]
    se_beta1_iiw_DWZ <- iiw_DWZmod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    iiw_DWZ_95CI_lower <- beta1_iiw_DWZ - 1.96*se_beta1_iiw_DWZ
    iiw_DWZ_95CI_upper <- beta1_iiw_DWZ + 1.96*se_beta1_iiw_DWZ
    
    beta1_iiw_DWZ_covered <- ifelse(iiw_DWZ_95CI_lower <= beta1 && beta1 <= iiw_DWZ_95CI_upper, TRUE, FALSE)
    
  }else{
    iiw_DWZmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, 
                                 id = id, weights = iiw_DWZ))
    beta1_iiw_DWZ <- iiw_DWZmod$coef[[1]]
    se_beta1_iiw_DWZ <- iiw_DWZmod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    iiw_DWZ_95CI_lower <- beta1_iiw_DWZ - 1.96*se_beta1_iiw_DWZ
    iiw_DWZ_95CI_upper <- beta1_iiw_DWZ + 1.96*se_beta1_iiw_DWZ
    
    beta1_iiw_DWZ_covered <- ifelse(iiw_DWZ_95CI_lower <= beta1 && beta1 <= iiw_DWZ_95CI_upper, TRUE, FALSE)
    
  }
  

  

  
  # intensity using D,W
  gamma.hat_DW <-   gamma.hat <- coxph(Surv(time.lag, time, observed) ~ D + Wt  - 1,  data = simdata)$coef
  iiw_DW <- exp(cbind(simdata_obs$D)%*%delta.hat)/exp(cbind(simdata_obs$D, simdata_obs$Wt)%*%gamma.hat)
  
  if(usesplines == T){
    iiw_DWmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, 
                                family = family, id = id, weights = iiw_DW))
    
    beta1_iiw_DW <- iiw_DWmod$coef[2,1]
    se_beta1_iiw_DW <- iiw_DWmod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    iiw_DW_95CI_lower <- beta1_iiw_DW - 1.96*se_beta1_iiw_DW
    iiw_DW_95CI_upper <- beta1_iiw_DW + 1.96*se_beta1_iiw_DW
    
    beta1_iiw_DW_covered <- ifelse(iiw_DW_95CI_lower <= beta1 && beta1 <= iiw_DW_95CI_upper, TRUE, FALSE)
    
  }else{
    iiw_DWmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, 
                                id = id, weights = iiw_DW))
    
    beta1_iiw_DW <- iiw_DWmod$coef[[1]]
    se_beta1_iiw_DW <- iiw_DWmod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    iiw_DW_95CI_lower <- beta1_iiw_DW - 1.96*se_beta1_iiw_DW
    iiw_DW_95CI_upper <- beta1_iiw_DW + 1.96*se_beta1_iiw_DW
    
    beta1_iiw_DW_covered <- ifelse(iiw_DW_95CI_lower <= beta1 && beta1 <= iiw_DW_95CI_upper, TRUE, FALSE)
    
  }

  
  # intensity using D, Z
  gamma.hat_DZ <-   gamma.hat <- coxph(Surv(time.lag, time, observed) ~ D + Z - 1, data = simdata)$coef
  iiw_DZ <- exp(cbind(simdata_obs$D)%*%delta.hat)/exp(cbind(simdata_obs$D,  simdata_obs$Z)%*%gamma.hat)
  
  if(usesplines == T){
    
    iiw_DZmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, 
                                family = family, id = id, weights = iiw_DZ))
    beta1_iiw_DZ <- iiw_DZmod$coef[2,1]
    se_beta1_iiw_DZ <- iiw_DZmod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    iiw_DZ_95CI_lower <- beta1_iiw_DZ - 1.96*se_beta1_iiw_DZ
    iiw_DZ_95CI_upper <- beta1_iiw_DZ + 1.96*se_beta1_iiw_DZ
    
    beta1_iiw_DZ_covered <- ifelse(iiw_DZ_95CI_lower <= beta1 && beta1 <= iiw_DZ_95CI_upper, TRUE, FALSE)
    
    
  }else{
    iiw_DZmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, 
                                id = id, weights = iiw_DZ))
    beta1_iiw_DZ <- iiw_DZmod$coef[[1]]
    se_beta1_iiw_DZ <- iiw_DZmod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    iiw_DZ_95CI_lower <- beta1_iiw_DZ - 1.96*se_beta1_iiw_DZ
    iiw_DZ_95CI_upper <- beta1_iiw_DZ + 1.96*se_beta1_iiw_DZ
    
    beta1_iiw_DZ_covered <- ifelse(iiw_DZ_95CI_lower <= beta1 && beta1 <= iiw_DZ_95CI_upper, TRUE, FALSE)
    
  }
  

  
  
  # intensity using  W, Z
  gamma.hat_WZ <-   gamma.hat <- coxph(Surv(time.lag, time, observed) ~ Wt + Z - 1,  data = simdata)$coef
  iiw_WZ <- exp(cbind(simdata_obs$D)%*%delta.hat)/exp(cbind( simdata_obs$Wt, simdata_obs$Z)%*%gamma.hat)
  
  if(usesplines == T){
    iiw_WZmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), 
                                data = simdata_obs, family = family, id = id, weights = iiw_WZ))
    
    beta1_iiw_WZ <- iiw_WZmod$coef[2,1]
    se_beta1_iiw_WZ <- iiw_WZmod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    iiw_WZ_95CI_lower <- beta1_iiw_WZ - 1.96*se_beta1_iiw_WZ
    iiw_WZ_95CI_upper <- beta1_iiw_WZ + 1.96*se_beta1_iiw_WZ
    
    beta1_iiw_WZ_covered <- ifelse(iiw_WZ_95CI_lower <= beta1 && beta1 <= iiw_WZ_95CI_upper, TRUE, FALSE)
    
  }else{
    iiw_WZmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, 
                                id = id, weights = iiw_WZ))
    
    beta1_iiw_WZ <- iiw_WZmod$coef[[1]]
    se_beta1_iiw_WZ <- iiw_WZmod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    iiw_WZ_95CI_lower <- beta1_iiw_WZ - 1.96*se_beta1_iiw_WZ
    iiw_WZ_95CI_upper <- beta1_iiw_WZ + 1.96*se_beta1_iiw_WZ
    
    beta1_iiw_WZ_covered <- ifelse(iiw_WZ_95CI_lower <= beta1 && beta1 <= iiw_WZ_95CI_upper, TRUE, FALSE)
    
  }

  
  
  # intensity using D
  gamma.hat_D <-   gamma.hat <- coxph(Surv(time.lag, time, observed) ~ D  - 1,  data = simdata)$coef
  iiw_D <- exp(cbind(simdata_obs$D)%*%delta.hat)/exp(cbind(simdata_obs$D)%*%gamma.hat)
  
  if(usesplines == T){
    
    iiw_Dmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, 
                               family = family, id = id, weights = iiw_D))
    beta1_iiw_D <- iiw_Dmod$coef[2,1]
    se_beta1_iiw_D <- iiw_Dmod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    iiw_D_95CI_lower <- beta1_iiw_D - 1.96*se_beta1_iiw_D
    iiw_D_95CI_upper <- beta1_iiw_D + 1.96*se_beta1_iiw_D
    
    beta1_iiw_D_covered <- ifelse(iiw_D_95CI_lower <= beta1 && beta1 <= iiw_D_95CI_upper, TRUE, FALSE)
  }else{
    
    iiw_Dmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, 
                               id = id, weights = iiw_D))
    beta1_iiw_D <- iiw_Dmod$coef[[1]]
    se_beta1_iiw_D <- iiw_Dmod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    iiw_D_95CI_lower <- beta1_iiw_D - 1.96*se_beta1_iiw_D
    iiw_D_95CI_upper <- beta1_iiw_D + 1.96*se_beta1_iiw_D
    
    beta1_iiw_D_covered <- ifelse(iiw_D_95CI_lower <= beta1 && beta1 <= iiw_D_95CI_upper, TRUE, FALSE)
  }
  

  
  # intensity using W
  gamma.hat_W <-   gamma.hat <- coxph(Surv(time.lag, time, observed) ~ Wt - 1,  data = simdata)$coef
  iiw_W <- exp(cbind(simdata_obs$D)%*%delta.hat)/exp(cbind(simdata_obs$Wt)%*%gamma.hat)
  
  if(usesplines == T){
    
    iiw_Wmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, 
                               family = family, id = id, weights = iiw_W))
    beta1_iiw_W <- iiw_Wmod$coef[2,1]
    se_beta1_iiw_W <- iiw_Wmod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    iiw_W_95CI_lower <- beta1_iiw_W - 1.96*se_beta1_iiw_W
    iiw_W_95CI_upper <- beta1_iiw_W + 1.96*se_beta1_iiw_W
    
    beta1_iiw_W_covered <- ifelse(iiw_W_95CI_lower <= beta1 && beta1 <= iiw_W_95CI_upper, TRUE, FALSE)
    
  }else{
    iiw_Wmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, 
                               id = id, weights = iiw_W))
    beta1_iiw_W <- iiw_Wmod$coef[[1]]
    se_beta1_iiw_W <- iiw_Wmod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    iiw_W_95CI_lower <- beta1_iiw_W - 1.96*se_beta1_iiw_W
    iiw_W_95CI_upper <- beta1_iiw_W + 1.96*se_beta1_iiw_W
    
    beta1_iiw_W_covered <- ifelse(iiw_W_95CI_lower <= beta1 && beta1 <= iiw_W_95CI_upper, TRUE, FALSE)
    
  }
  

  
  # intensity using Z
  gamma.hat_Z <-   gamma.hat <- coxph(Surv(time.lag, time, observed) ~  Z - 1,  data = simdata)$coef
  iiw_Z <- exp(cbind(simdata_obs$D)%*%delta.hat)/exp(cbind( simdata_obs$Z)%*%gamma.hat)
  if(usesplines == T){
    iiw_Zmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs,
                               family = family, id = id, weights = iiw_Z))
    beta1_iiw_Z <- iiw_Zmod$coef[2,1]
    se_beta1_iiw_Z <- iiw_Zmod$coef[2,2] #robust SE
    
    
    #see if CI contains true value or not
    iiw_Z_95CI_lower <- beta1_iiw_Z - 1.96*se_beta1_iiw_Z
    iiw_Z_95CI_upper <- beta1_iiw_Z + 1.96*se_beta1_iiw_Z
    
    beta1_iiw_Z_covered <- ifelse(iiw_Z_95CI_lower <= beta1 && beta1 <= iiw_Z_95CI_upper, TRUE, FALSE)
    
  }else{
    iiw_Zmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, 
                               id = id, weights = iiw_Z))
    beta1_iiw_Z <- iiw_Zmod$coef[[1]]
    se_beta1_iiw_Z <- iiw_Zmod$coef[[2]] #robust SE
    
    
    #see if CI contains true value or not
    iiw_Z_95CI_lower <- beta1_iiw_Z - 1.96*se_beta1_iiw_Z
    iiw_Z_95CI_upper <- beta1_iiw_Z + 1.96*se_beta1_iiw_Z
    
    beta1_iiw_Z_covered <- ifelse(iiw_Z_95CI_lower <= beta1 && beta1 <= iiw_Z_95CI_upper, TRUE, FALSE)
    
  }
  

  
  out <- list(beta1, beta1_naive, se_beta1_naive, beta1_naive_covered,
              beta1_iiw_DWZ, se_beta1_iiw_DWZ, beta1_iiw_DWZ_covered,
              beta1_iiw_DW, se_beta1_iiw_DW, beta1_iiw_DW_covered,
              beta1_iiw_DZ, se_beta1_iiw_DZ, beta1_iiw_DZ_covered,
              beta1_iiw_WZ, se_beta1_iiw_WZ, beta1_iiw_WZ_covered,
              beta1_iiw_D, se_beta1_iiw_D, beta1_iiw_D_covered,
              beta1_iiw_W, se_beta1_iiw_W, beta1_iiw_W_covered,
              beta1_iiw_Z, se_beta1_iiw_Z, beta1_iiw_Z_covered)
  
  
  names(out) <- c('beta1', 'beta1_naive', 'se_beta1_naive', 'beta1_naive_covered',
                  'beta1_iiw_DWZ', 'se_beta1_iiw_DWZ', 'beta1_iiw_DWZ_covered',
                  'beta1_iiw_DW', 'se_beta1_iiw_DW', 'beta1_iiw_DW_covered',
                  'beta1_iiw_DZ', 'se_beta1_iiw_DZ', 'beta1_iiw_DZ_covered',
                  'beta1_iiw_WZ', 'se_beta1_iiw_WZ', 'beta1_iiw_WZ_covered',
                  'beta1_iiw_D', 'se_beta1_iiw_D', 'beta1_iiw_D_covered',
                  'beta1_iiw_W', 'se_beta1_iiw_W', 'beta1_iiw_W_covered',
                  'beta1_iiw_Z', 'se_beta1_iiw_Z', 'beta1_iiw_Z_covered')
  return(out)
}




simulateResultsVarSel <-  function(N,n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                                   usesplines = F, ncutpts = NULL,
                                   inParallel = F, nclusters = NULL){
  # Gets required results from N simulation runs. We are interested in Beta1 and beta2

  if(inParallel == T){
    
    registerDoParallel(nclusters)
    
    results <- foreach(i = 1:N, .combine = rbind, .export = c("expit", "simulateOneFIPTIW", 
                                                                    "geeglm", "coxph", "Surv", "gendata_FIPTIW_continuous", 'gendata_FIPTIW_binary', 
                                                                    "simulateOneVarSel", "bs", "getObsTimesViaThinning"), 
                       .packages = c("dplyr", "geepack")) %dopar% {
                                                                      
              simrun <- simulateOneVarSel(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                                          usesplines, ncutpts)
              c(simrun$beta1_naive,  simrun$se_beta1_naive, simrun$beta1_naive_covered,
                simrun$beta1_iiw_D, simrun$se_beta1_iiw_D, simrun$beta1_iiw_D_covered,
                simrun$beta1_iiw_W, simrun$se_beta1_iiw_W, simrun$beta1_iiw_W_covered,
                simrun$beta1_iiw_Z, simrun$se_beta1_iiw_Z, simrun$beta1_iiw_Z_covered,
                simrun$beta1_iiw_DW, simrun$se_beta1_iiw_DW, simrun$beta1_iiw_DW_covered,
                simrun$beta1_iiw_DZ, simrun$se_beta1_iiw_DZ, simrun$beta1_iiw_DZ_covered,
                simrun$beta1_iiw_WZ, simrun$se_beta1_iiw_WZ, simrun$beta1_iiw_WZ_covered,
                simrun$beta1_iiw_DWZ,simrun$se_beta1_iiw_DWZ, simrun$beta1_iiw_DWZ_covered)
        }
    
  } else {
    
    results <- matrix(data = NA, nrow = N, ncol = 24)
    
    
    for(i in 1:N){
      if(i%%100 == 0){print(i)}
      simrun <- simulateOneVarSel(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                                  usesplines, ncutpts)
      results[i,] <- c(simrun$beta1_naive,  simrun$se_beta1_naive, simrun$beta1_naive_covered,
                             simrun$beta1_iiw_D, simrun$se_beta1_iiw_D, simrun$beta1_iiw_D_covered,
                             simrun$beta1_iiw_W, simrun$se_beta1_iiw_W, simrun$beta1_iiw_W_covered,
                             simrun$beta1_iiw_Z, simrun$se_beta1_iiw_Z, simrun$beta1_iiw_Z_covered,
                             simrun$beta1_iiw_DW, simrun$se_beta1_iiw_DW, simrun$beta1_iiw_DW_covered,
                             simrun$beta1_iiw_DZ, simrun$se_beta1_iiw_DZ, simrun$beta1_iiw_DZ_covered,
                             simrun$beta1_iiw_WZ, simrun$se_beta1_iiw_WZ, simrun$beta1_iiw_WZ_covered,
                             simrun$beta1_iiw_DWZ,simrun$se_beta1_iiw_DWZ, simrun$beta1_iiw_DWZ_covered)
      
    }
    
  }
 
  
  results_beta1 <- results[, c(1,4,7, 10, 13, 16, 19, 22)]
  
  #~~~~ Results for Beta1 ~~~~#
  
  bias_beta1 <- round(apply(results_beta1, FUN =  mean, MARGIN = 2)  - beta1 ,3)
  names(bias_beta1) <- c("Naive", "D", "W", "Z", "D,W", "D,Z", "W,Z", "D,W,Z")

  
  var_beta1 <- round(apply(results_beta1, FUN = var, MARGIN = 2), 3)
  names(var_beta1) <- c("Naive", "D", "W", "Z", "D,W", "D,Z", "W,Z", "D,W,Z")

  
  mse_beta1 <- round(apply((results_beta1 - beta1)^2, FUN = mean, MARGIN = 2),3)
  names(mse_beta1) <- c("Naive", "D", "W", "Z", "D,W", "D,Z", "W,Z", "D,W,Z")


  
  avgSE_beta1 <- results[, c(2, 5, 8, 11, 14, 17, 20, 23)]
  avgSE_beta1 <- round(apply(avgSE_beta1, FUN =  mean, MARGIN = 2), 3)
  names(avgSE_beta1) <- c("Naive", "D", "W", "Z", "D,W", "D,Z", "W,Z", "D,W,Z")
  
  
  coverage_beta1 <- results[, c(3, 6, 9, 12, 15, 18, 21, 24)]
  coverage_beta1 <- round(apply(coverage_beta1, FUN =  mean, MARGIN = 2), 3)
  names(coverage_beta1) <- c("Naive", "D", "W", "Z", "D,W", "D,Z", "W,Z", "D,W,Z")
  
  
  
  out <- list(bias_beta1, var_beta1, mse_beta1, avgSE_beta1, coverage_beta1)
  names(out) <- c('bias_beta1', 'var_beta1', 'mse_beta1', 'avgSE_beta1', 'coverage_beta1')
  return(out)
  
}



simulateALLVarSel <- function(N,n, beta1, beta2vec, beta3, gamma1, gamma2vec, gamma3, alpha0, alpha1, tau, outcome,
                              usesplines = F, ncutpts = NULL,
                              inParallel = F, nclusters = NULL){
  #N: number of simulation runs
  #This function aggregates simulation results for varying n and beta02
  
  
  
  resultsmat_beta01 <- matrix(NA, nrow = 1, ncol = 11)
  
  i = 1
  
  for(gamma2 in gamma2vec){
    for(beta2 in beta2vec){
      print(paste("Now on gamma02 = ", gamma2, " and beta2 = ", beta2, sep = ""))
      result <- simulateResultsVarSel(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau, outcome,
                                      usesplines, ncutpts,
                                      inParallel, nclusters)
      
      resultsmat_beta01 <- rbind(resultsmat_beta01, 
                                 c(gamma2, beta2, rep(" ", 9)),
                                 c("", "", "Bias:", result$bias_beta1), 
                                 c("", "", "Variance:", result$var_beta1), 
                                 c("", '', "MSE:", result$mse_beta1),
                                 c("", '', "Avg SE:", result$avgSE_beta1),
                                 c("", '', "Coverage:", result$coverage_beta1))
      
      
      i <- i + 1
    }
  }
  
  resultsmat_beta01 <- resultsmat_beta01[-1,]
  
  
  out <- list(resultsmat_beta01)
  names(out) <- c("resultsmat_beta01")
  return(out)
}








#############################################################
############ Weight Trimming Simulation #####################
#############################################################



simulateOneWeightTrimming <- function( n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
                                       alpha0, alpha1, tau, usesplines, ncutpts = 2, outcome){
  #weights are stabilized
  

  
  if(outcome == "continuous"){
    
    singlerun <- gendata_FIPTIW_continuous(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau,
                                           censinform = F, eta1 = NULL, eta2 = NULL, eta3= NULL)
    
    family = gaussian(link = "identity")
    
  }else if(outcome == "binary"){
    singlerun <- gendata_FIPTIW_binary(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, alpha0, alpha1, tau,
                                       censinform = F, eta1 = NULL, eta2 = NULL, eta3= NULL)
    family = binomial(link = "logit")
  }else{
    return(print("insufficient outcome type"))
  }
  
  simdata <- singlerun$simdata
  numevents <- singlerun$numevents
  
  # filter to have only data when y is observed 
  simdata_obs <- simdata%>%filter(!is.na(y))
  
  
  if(usesplines == T){
    
    if(ncutpts == 2){
      #find tertiles for cubic spline to estimate the time-var intercept
      terti<-quantile(0:tau , c(0.3333, 0.66666), type = 1) 
    }else if(ncutpts == 1){
      terti<-quantile(0:tau , c(0.5), type = 1) 
    }else{
      return(print(paste("unsuitable number of cutp points for spline")))
    }
    
  }
  
  
  
 
  
  # IIW Weights
  delta.hat <-coxph(Surv(time.lag, time, observed) ~ D - 1, data = simdata)$coef
  gamma.hat <- coxph(Surv(time.lag, time, observed) ~ D + Wt + Z - 1, data = simdata)$coef
  iiw <- exp(cbind(simdata_obs$D)%*%delta.hat)/exp(cbind(simdata_obs$D, simdata_obs$Wt, simdata_obs$Z)%*%gamma.hat)

  
  
  
  # IPTW weights
  
  psmod <- glm(D ~ W, family = binomial(link = "logit"), data = simdata)
  ps <- expit(predict(psmod, newdata = simdata_obs))
  prDmod <- glm(D ~ 1, family = binomial(link = "logit"), data = simdata)
  prD <- expit(summary(prDmod)$coef[1,1])
  iptw <- prD/ps*simdata_obs$D + (1-prD)/(1-ps)*(1-simdata_obs$D)
  
  
  fiptiw <- iiw*iptw
  
  
  
  #### NO TRIMMING
  
  if(usesplines == T){
    notrimmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, family = family, id = id, weights = fiptiw))
    beta1_notrim <- notrimmod$coef[2,1]
    se_beta1_notrim <- notrimmod$coef[2,2] #robust SE
    
  }else{
    notrimmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id, weights = fiptiw))
    beta1_notrim <- notrimmod$coef[[1]]
    se_beta1_notrim <- notrimmod$coef[[2]] #robust SE
    
  }
  

  
  #see if CI contains true value or not
  notrim_95CI_lower <- beta1_notrim - 1.96*se_beta1_notrim
  notrim_95CI_upper <- beta1_notrim + 1.96*se_beta1_notrim
  
  beta1_notrim_covered <- ifelse(notrim_95CI_lower <= beta1 && beta1 <= notrim_95CI_upper, TRUE, FALSE)
  
  
  

  
  beta1_trimmedfirst_vec <- c()
  beta1_trimmedafter_vec <- c()
  
  covered_trimmedfirst_vec <- c()
  covered_trimmedafter_vec <- c()
  
  
  thres_iiw <- quantile(iiw, seq(0.5, 1, by = 0.01))
  thres_iptw <- quantile(iptw, seq(0.5, 1, by = 0.01))
  thres_fiptiw <- quantile(fiptiw, seq(0.5, 1, by = 0.01))
  
  
  
  for(i in 1:51){
    # for each percentile from 50 to 100 (no trimming)
    iiw_trimmed <- iiw
    iptw_trimmed <- iptw
    
    # trim the IIW and IPTW weights individually, multiply them
    iiw_trimmed[iiw_trimmed > thres_iiw[i]] <- thres_iiw[i]
    iptw_trimmed[iptw_trimmed > thres_iptw[i]] <- thres_iptw[i]
    
    fiptiw_trimmedfirst <-  iiw_trimmed*iptw_trimmed
    
    
    
    # trim the FIPTIW weights after mult
    fiptiw_trimmedafter <- fiptiw
    fiptiw_trimmedafter[fiptiw_trimmedafter > thres_fiptiw[i]] <- thres_fiptiw[i]
    
    
    ###### Calculate beta and coverage for trimming individually first 
    
    if(usesplines == T){
      trimmedfirstmod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, 
                                        family = family, id = id, weights = fiptiw_trimmedfirst))
      beta1_trimmedfirst <- trimmedfirstmod$coef[2,1]
      se_beta1_trimmedfirst <- trimmedfirstmod$coef[2,2] #robust SE
      
    }else{
      trimmedfirstmod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id, weights = fiptiw_trimmedfirst))
      beta1_trimmedfirst <- trimmedfirstmod$coef[[1]]
      se_beta1_trimmedfirst <- trimmedfirstmod$coef[[2]] #robust SE
    }
    

    
    
    #see if CI contains true value or not
    trimmedfirst_95CI_lower <- beta1_trimmedfirst - 1.96*se_beta1_trimmedfirst
    trimmedfirst_95CI_upper <- beta1_trimmedfirst + 1.96*se_beta1_trimmedfirst
    
    beta1_trimmedfirst_covered <- ifelse(trimmedfirst_95CI_lower <= beta1 && beta1 <= trimmedfirst_95CI_upper, TRUE, FALSE)
    
    
    
    ###### Calculate beta and coverage for trimming after multiplying fiptiw weights
    
    
    if(usesplines == T){
      trimmedaftermod <- summary(geeglm(y ~ D + bs(simdata_obs$time,degree=3,knots=c(terti)), data = simdata_obs, 
                                        family = family, id = id, weights = fiptiw_trimmedafter))
      beta1_trimmedafter <- trimmedaftermod$coef[2,1]
      se_beta1_trimmedafter <- trimmedaftermod$coef[2,2] #robust SE
    }else{
      trimmedaftermod <- summary(geeglm(y ~ D + offset(2-time) -1 , data = simdata_obs, family = family, id = id, weights = fiptiw_trimmedafter))
      beta1_trimmedafter <- trimmedaftermod$coef[[1]]
      se_beta1_trimmedafter <- trimmedaftermod$coef[[2]] #robust SE
    }
    

    
    
    #see if CI contains true value or not
    trimmedafter_95CI_lower <- beta1_trimmedafter - 1.96*se_beta1_trimmedafter
    trimmedafter_95CI_upper <- beta1_trimmedafter + 1.96*se_beta1_trimmedafter
    
    beta1_trimmedafter_covered <- ifelse(trimmedafter_95CI_lower <= beta1 && beta1 <= trimmedafter_95CI_upper, TRUE, FALSE)
    
    
    
    
    beta1_trimmedfirst_vec <- c(beta1_trimmedfirst_vec, beta1_trimmedfirst)
    beta1_trimmedafter_vec <- c(beta1_trimmedafter_vec, beta1_trimmedafter)
    covered_trimmedfirst_vec <- c(covered_trimmedfirst_vec, beta1_trimmedfirst_covered)
    covered_trimmedafter_vec <- c(covered_trimmedafter_vec, beta1_trimmedafter_covered)
    
    
    
    
    }
    
  names(beta1_trimmedfirst_vec) <- c(50:100)
  names(beta1_trimmedafter_vec) <- c(50:100)
  names(covered_trimmedfirst_vec) <- c(50:100)
  names(covered_trimmedafter_vec) <- c(50:100)
  
  proportion_raw_greater_than_5_fiptiw <- sum(fiptiw > 5)/length(fiptiw)
  proportion_raw_greater_than_10_fiptiw <- sum(fiptiw > 10)/length(fiptiw)
  proportion_raw_greater_than_20_fiptiw <- sum(fiptiw > 20)/length(fiptiw)
  maxweight_fiptiw <- max(fiptiw)
  
  proportion_raw_greater_than_5_iiw <- sum(iiw > 5)/length(iiw)
  proportion_raw_greater_than_10_iiw <- sum(iiw > 10)/length(iiw)
  proportion_raw_greater_than_20_iiw <- sum(iiw > 20)/length(iiw)
  maxweight_iiw <- max(iiw)
  
  proportion_raw_greater_than_5_iptw <- sum(iptw > 5)/length(iptw)
  proportion_raw_greater_than_10_iptw <- sum(iptw > 10)/length(iptw)
  proportion_raw_greater_than_20_iptw <- sum(iptw > 20)/length(iptw)
  maxweight_iptw <- max(iptw)
  

  #find optimal cutpoint?
  
  out <- list(beta1_trimmedfirst_vec, beta1_trimmedafter_vec, covered_trimmedfirst_vec, covered_trimmedafter_vec,
              proportion_raw_greater_than_5_fiptiw, proportion_raw_greater_than_10_fiptiw, 
              proportion_raw_greater_than_20_fiptiw, proportion_raw_greater_than_5_iiw, 
              proportion_raw_greater_than_10_iiw, 
              proportion_raw_greater_than_20_iiw, proportion_raw_greater_than_5_iptw, 
              proportion_raw_greater_than_10_iptw, 
              proportion_raw_greater_than_20_iptw,
              maxweight_fiptiw, maxweight_iiw, maxweight_iptw)
  names(out) <- c("TrimmedFirst", "TrimmedAfterMult", 'covered_trimmedfirst_vec', 'covered_trimmedafter_vec',
                  'proportion_raw_greater_than_5_fiptiw', 'proportion_raw_greater_than_10_fiptiw', 
                  'proportion_raw_greater_than_20_fiptiw', 'proportion_raw_greater_than_5_iiw', 
                  'proportion_raw_greater_than_10_iiw', 
                  'proportion_raw_greater_than_20_iiw', 'proportion_raw_greater_than_5_iptw', 
                  'proportion_raw_greater_than_10_iptw', 'proportion_raw_greater_than_20_iptw',
                  'maxweight_fiptiw', 'maxweight_iiw', 'maxweight_iptw')
  
  return(out)
  
  }





simulateResultsWeightTrimming <- function(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
                                           alpha0, alpha1, tau, usesplines, ncutpts, outcome){
  #for a given scheme, simulates N realizations and aggregates results
  
  resultsmat_trimmedfirst <- matrix(data = NA, nrow = N, ncol = 51)
  resultsmat_trimmedaftermult <- matrix(data = NA, nrow = N, ncol = 51)
  
  coveredmat_trimmedfirst <- matrix(data = NA, nrow = N, ncol = 51)
  coveredmat_trimmedaftermult <- matrix(data = NA, nrow = N, ncol = 51)
  
  other_fiptiw <- matrix(data = NA, nrow = N, ncol = 4)
  other_iiw <- matrix(data = NA, nrow = N, ncol = 4)
  other_iptw <- matrix(data = NA, nrow = N, ncol = 4)
  
  
  
  for(i in 1:N){
    if(i %% 100 == 0){print(i)}
    
    simrun <- simulateOneWeightTrimming(n = n, beta1 = beta1, beta2 = beta2, beta3 = beta3, gamma1 = gamma1, 
                                        gamma2 = gamma2, gamma3 = gamma3, 
                                     alpha0 = alpha0, alpha1 = alpha1, tau = tau,usesplines = usesplines, ncutpts = ncutpts, outcome = outcome)
    resultsmat_trimmedfirst[i,] <- simrun$TrimmedFirst
    resultsmat_trimmedaftermult[i,] <- simrun$TrimmedAfterMult
    coveredmat_trimmedfirst[i,] <- simrun$covered_trimmedfirst_vec
    coveredmat_trimmedaftermult[i,] <- simrun$covered_trimmedafter_vec
    
    
    
    other_fiptiw[i,] <- c(simrun$proportion_raw_greater_than_5_fiptiw, simrun$proportion_raw_greater_than_10_fiptiw, 
                                simrun$proportion_raw_greater_than_20_fiptiw, simrun$maxweight_fiptiw)
    other_iiw[i,] <- c(simrun$proportion_raw_greater_than_5_iiw, simrun$proportion_raw_greater_than_10_iiw, 
                         simrun$proportion_raw_greater_than_20_iiw, simrun$maxweight_iiw)
    other_iptw[i,] <- c(simrun$proportion_raw_greater_than_5_iptw, simrun$proportion_raw_greater_than_10_iptw, 
                         simrun$proportion_raw_greater_than_20_iptw, simrun$maxweight_iptw)
  }
  
  
  
  colnames(resultsmat_trimmedfirst) <- c(50:100)
  colnames(resultsmat_trimmedaftermult) <- c(50:100)
  colnames(coveredmat_trimmedfirst) <- c(50:100)
  colnames(coveredmat_trimmedaftermult) <- c(50:100)
  
  
  
  bias_beta1_trimmedfirst <- round(apply(resultsmat_trimmedfirst, FUN =  mean, MARGIN = 2)  - beta1, 3)
  bias_beta1_trimmedaftermult <- round(apply(resultsmat_trimmedaftermult, FUN =  mean, MARGIN = 2)  - beta1 , 3)
  
  coverage_beta1_trimmedfirst <- round(apply(coveredmat_trimmedfirst, FUN =  mean, MARGIN = 2) , 3)
  coverage_beta1_trimmedaftermult <- round(apply(coveredmat_trimmedaftermult, FUN =  mean, MARGIN = 2) , 3)
  
  var_beta1_trimmedfirst <- round(apply(resultsmat_trimmedfirst, FUN = var, MARGIN = 2), 3)
  var_beta1_trimmedaftermult <- round(apply(resultsmat_trimmedaftermult, FUN = var, MARGIN = 2), 3)
  
  mse_beta1_trimmedfirst <- round(apply((resultsmat_trimmedfirst - beta1)^2, FUN = mean, MARGIN = 2), 3)
  mse_beta1_trimmedaftermult <- round(apply((resultsmat_trimmedaftermult - beta1)^2, FUN = mean, MARGIN = 2), 3)
  
  
  
  biasmat <- data.frame(cbind(c(rep("TrimmedFirst", 51), rep("TrimmedAfter", 51)), 
                       as.numeric(c(50:100, 50:100)),
                       c(bias_beta1_trimmedfirst, bias_beta1_trimmedaftermult)))
  
  colnames(biasmat) <- c("Method", "Percentile", "Bias")
  biasmat$Percentile <- as.numeric(biasmat$Percentile)
  biasmat$Bias <- as.numeric(biasmat$Bias)
  
  
  coveragemat <- data.frame(cbind(c(rep("TrimmedFirst", 51), rep("TrimmedAfter", 51)), 
                       c(50:100, 50:100),
                       c(coverage_beta1_trimmedfirst, coverage_beta1_trimmedaftermult)))
  colnames(coveragemat) <- c("Method", "Percentile", "Coverage")
  coveragemat$Percentile <- as.numeric(coveragemat$Percentile)
  coveragemat$Coverage <- as.numeric(coveragemat$Coverage)
  
  
  varmat <- data.frame(cbind(c(rep("TrimmedFirst", 51), rep("TrimmedAfter", 51)), 
                       c(50:100, 50:100),
                       c(var_beta1_trimmedfirst, var_beta1_trimmedaftermult)))
  colnames(varmat) <- c("Method", "Percentile", "Variance")
  varmat$Percentile <- as.numeric(varmat$Percentile)
  varmat$Variance <- as.numeric(varmat$Variance)
  
  
  
  msemat <- data.frame(cbind(c(rep("TrimmedFirst", 51), rep("TrimmedAfter", 51)), 
                       as.numeric(c(50:100, 50:100)),
                       c(mse_beta1_trimmedfirst, mse_beta1_trimmedaftermult)))
  colnames(msemat) <- c("Method", "Percentile", "MSE")
  msemat$Percentile <- as.numeric(msemat$Percentile)
  msemat$MSE <- as.numeric(msemat$MSE)
  
  
  
 
  
  
  othermat_fiptiw <- round(apply(other_fiptiw, FUN =  mean, MARGIN = 2), 4)
  names(othermat_fiptiw) <- c("mean_propweightgreaterthan5", "mean_propweightgreaterthan10",
                       "mean_propweightgreaterthan20", "mean_maxweight")
  othermat_iiw <- round(apply(other_iiw, FUN =  mean, MARGIN = 2), 4)
  names(othermat_iiw) <- c("mean_propweightgreaterthan5", "mean_propweightgreaterthan10",
                             "mean_propweightgreaterthan20", "mean_maxweight")
  othermat_iptw <- round(apply(other_iptw, FUN =  mean, MARGIN = 2), 4)
  names(othermat_iptw) <- c("mean_propweightgreaterthan5", "mean_propweightgreaterthan10",
                             "mean_propweightgreaterthan20", "mean_maxweight")
  
  
  out <- list(biasmat, varmat, msemat, coveragemat,
              othermat_fiptiw, othermat_iiw, othermat_iptw)
  
  names(out) <- c('biasmat', 'varmat', 'msemat', 'coveragemat',
                  'othermat_fiptiw', 'othermat_iiw', 'othermat_iptw')
  
  return(out)
}


plotWeightTrimmingResults <- function(results){
  
  biasmat <- results$biasmat
  varmat <- results$varmat
  msemat <- results$msemat
  coveragemat <- results$coveragemat
  
  
  # bias plot
    p1 <- ggplot(biasmat, aes(x = Percentile, y = Bias, color = Method, shape = Method)) +
    geom_line() +
    geom_point(size = 0.8) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed")+
    theme(legend.title=element_blank())+
    xlab("Threshold (Percentile)") +
    ylab("Bias")+
    ylim(0.0, 0.4) +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10))

  
  ### Variance plots
  
  p2 <- ggplot(varmat, aes(x = Percentile, y = Variance, color = Method, shape = Method)) + 
    geom_line() + 
    geom_point(size = 0.8) + 
    theme(legend.title=element_blank())+
    xlab("Threshold (Percentile)")+
    ylim(0.1, 0.2)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10))
  
  ### MSE plots
  
  
  p3 <- ggplot(msemat, aes(x = Percentile, y = MSE, color = Method, shape = Method)) + 
    geom_line() + 
    geom_point(size = 0.8)+ 
    theme(legend.title=element_blank()) +
    xlab("Threshold (Percentile)")+
    ylim(0.1, 0.8)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10))
  
  p4 <- ggplot(coveragemat, aes(x = Percentile, y = Coverage, color = Method, shape = Method)) + 
    geom_line() + 
    geom_point(size = 0.8)+ 
    theme(legend.title=element_blank()) +
    xlab("Threshold (Percentile)")+
    ylim(0.70, 0.97)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10)) +
    geom_hline(yintercept = 0.95, color = "red", linetype = "dashed")
  
  
  #combinedplots <- p1 + p2 + p3 + p4 & theme(legend.position = "bottom")
  
  out <- list(p1, p2, p3, p4)
  names(out) <- c("biasplot", "varplot", "mseplot", "coverageplot")
  
  return(out)
  
  
}


compareMinimums <- function(results){
  

  
  ### determine minimum bias, var, mse, etc
  
  #bias closest to zero
  biasmat <- results$biasmat
  msemat <- results$msemat
  coveragemat <- results$coveragemat
  
  minbias <- biasmat[which(abs(biasmat$Bias) == min(abs(biasmat$Bias))),] 
  bias95 <- biasmat[biasmat$Percentile==95, ]
  bias100 <- biasmat[biasmat$Percentile==100, ]
  
  
  minmse <- msemat[which(abs(msemat$MSE) == min(abs(msemat$MSE))),] 
  mse95 <- msemat[msemat$Percentile==95, ]
  mse100 <- msemat[msemat$Percentile==100, ]
  
  cov95 <- coveragemat[coveragemat$Percentile==95, ]
  cov100 <- coveragemat[coveragemat$Percentile==100, ]
  

  
  out <- list(results$othermat_iptw, results$othermat_iiw, results$othermat_fiptiw,
              minbias, bias95, bias100,
              minmse, mse95, mse100,
              cov95, cov100)
  names(out) <- c("dist_iptw", "dist_iiw","dist_fiptiw",
              "minbias", "bias95", "bias100",
              "minmse", "mse95", "mse100",
              "cov95", "cov100")
  return(out)
  
}

