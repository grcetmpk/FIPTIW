#old code


# 
# gendataMIPI_trimming_timevartrt <- function(n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                  alpha0, alpha1, tau){
#   
#   
#   # data generation for the thinning algorithm
#   #consider discretized times from 0 to tau by 0.01
#   # don't use thinning algorithm
#   
#   
#   # D: Treatment assignment indicator (may be time-varying if timevartrt = T)
#   # W(t): Treatment confounder = W*0.5*t
#   # Z: Observation times confounder
#   # Y: Time-varying outcome
#   # t: observation time
#   # phi: random effect
#   # epsilon: random error, function of time
#   
#   # Data generation model: Y_i(t)  = (2-t) + beta1D_i(t) + beta2[W_i(t) - E(W_i(t)|D_i(t))] + 
#   #                                     beta3[Z_i(t) - E(Z_i(t)|D_i(t))] + phi_i + epsilon_i(t)
#   #
#   #
#   # Observation intensity: lambda_i(t) = eta_i*exp(gamma1D_i(t) + gamma2W_i(t) + gamma3Z_i(t))
#   # Treatment assignment: P(D_i(t) = 1) = expit(alpha0 + alpha1W_i(t)) 
#   #
#   #
#   # tau: maximum follow up time, independent censoring.
#   # N: number of simulations
#   #
#   
#   
#   
#   mu_Z_D0 <- 2
#   var_Z_D0 <- 1
#   mu_Z_D1 <- 0
#   var_Z_D1 <- 0.5
#   
#   var_phi <- 0.25
#   var_epsilon <- 1
#   
#   
#   
#   #generate subjects one by one and then combine
#   id = 1
#   simdatafull <- data.frame(matrix(NA, ncol = 13))
#   colnames(simdatafull) <- c("id", "time", "D", "W", "Z", "phi", "censortime", "intensity",
#                              "observed", "nevents", "cexp_W_D", "cexp_Z_D",  "y")
#   
#   while(id <= n){
#     
#     time = seq(0, tau, by = 0.01)
#     
#     ## Denerate Covariates
#     
#     # generate treatment confounder w ~ N(0,1)
#     W <- runif(1, 0, 1)*0.5*time
#     
#     
#     # generate treatment assignment
#     prD <- expit(alpha0 + alpha1*W)
#     D <- rbinom(length(prD), 1, prD)
#     
#     # generate observation times confounder Z ~ N( mu_Z_D0, var_Z_D0) if D = 0 and
#     # N( mu_Z_D1, var_Z_D1) if D = 1
#     
#     Z <- ifelse(D == 0, rnorm(1, mu_Z_D0, sqrt(var_Z_D0)), rnorm(1, mu_Z_D1, sqrt(var_Z_D1)))
#     
#     
#     # generate random effect
#     phi <- rnorm(1, 0, var_phi)
#     
#     
#     ## Denerate Observation Times
#     
#     # simulate censoring time
#     censortime <- runif(1, tau/2, tau)
#     
#     # calculate eta
#     eta <- rgamma(1, shape = 100, scale = 0.01) #eta is for random effect with mean 1 and sd 0.1
#     #(obs times within subject correlated if sd !0)
#     
#     intensity <- eta*sqrt(time)/2*exp(gamma1*D + gamma2*W + gamma3*Z)
#     propprob <- 0.1*intensity
#     
#     observed <- rbinom(length(time), 1, propprob)
#     
#     #force everyone to have baseline observation
#     observed[1] <- 1
#     
#   
#     simdata <- data.frame(cbind(rep(id, length(time)), time, D, W,  Z, phi,
#                                 censortime, intensity, observed))
#     colnames(simdata)[1] <- "id"
#     
#     simdataobserved <- simdata[simdata$observed == 1, ]
#     simdataobserved <- simdataobserved[simdataobserved$time <= censortime, ]
#     
#     if(nrow(simdataobserved) > 0){
#       
#       # calculate the number of events for individual i
#       nevents <- nrow(simdataobserved)
#       
#       simdataobserved$nevents <- nevents
#       
#       ## Denerate outcome
# 
#       simdataobserved$cexp_W_D <- 0.5*0.5*simdataobserved$time #W generated from N(0,1) independently from D
#       simdataobserved$cexp_Z_D <- ifelse(simdataobserved$D == 0, mu_Z_D0, mu_Z_D1)
#       
#       
#       simdataobserved$y <- (2 - simdataobserved$time) + beta1*simdataobserved$D + 
#         beta2*(simdataobserved$W*0.5*simdataobserved$time - simdataobserved$cexp_W_D) +
#         beta3*(simdataobserved$Z - simdataobserved$cexp_Z_D) +
#         rep(phi, dim(simdataobserved)[1]) + rnorm(dim(simdataobserved)[1], 0, sqrt(var_epsilon))
#       
#       # 
#       # #remove duplicated time
#       # #if subject has two of the same time, remove one of them.
#       # is.dup <- duplicated(simdataobserved[, c(1,2)])#find which are duplicates (by ID and time)
#       # if(sum(is.dup > 0)){
#       #   simdataobserved <- simdataobserved[!is.dup,] #remove them
#       # }
#       # 
#       # #make sure observation times >0
#       # simdata <- simdata[simdata$time>=0, ]
#       
#       # filter censored observations
#       # censored <- which(simdata$time > censortime)
#       # 
#       # if(length(censored) > 0){
#       #   # remove censored observations
#       #   simdata <- simdata[-censored,]
#       # }
#       
#       simdatafull <- rbind(simdatafull, simdataobserved)
#       
#     }
#     
#     
#     
#     id = id + 1
#     
#   }
#   
#   #remove NA row
#   simdata <- simdatafull[-1,]
#   
#   # order dataset by ID then time
#   simdata <- simdata[with(simdata, order(id, time)), ]
#   
#   
#   #include a variable counting observation number, to be used for lagging time for Surv function
#   simdata$obsnumber <- with(simdata, ave(id, id, FUN = seq_along))
#   
#   # #create lagged time variable
#   simdata$time.lag <- simdata$time[c(nrow(simdata),1:(nrow(simdata)-1))]
#   simdata$time.lag[simdata$obsnumber == 1] <- -0.01
#   
#   simdata$event <- 1
#   
#   
#   #check the generated data
#   numevents <- summary(tapply(simdata$event, simdata$id, sum))
#   
#   
#   #also get data for individuals at baseline (all same here)
#   baselinedata <- simdata[simdata$obsnumber == 1, ]
#   newn <- dim(baselinedata)[1] #number of people after censoring etc
#   
#   out <- list(simdata, numevents, baselinedata, newn)
#   names(out) <- c("simdata", "numevents", "baselinedata", "newn")
#   
#   return(out)
#   
# }

