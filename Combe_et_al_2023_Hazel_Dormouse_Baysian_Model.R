#########################################################################
# IPM - Integrated Population Model Dormouse data            			#
#											#
# by											#
# Fraser J. Combe and W Edwin Harris				#
#########################################################################
nyears<- # number of years capture study
popcount<-#number of unique individuals
J <- # Number of offspring 
R <- # Number of surveyed broods
#A function to make a m-array from Capture History data for adults and juveniles.  Capture #history for adults (CHA) and juveniles (CHJ) are separately combined into this m-array.  #For more information on creating m-array see Kery and Shaub (2012)
#so-called m-array
marray<-function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  for (t in 1:n.occasions){ # Calc # of released individuals each time period
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}
# Run marray function on CH
marrayA <- marray(CHA)
marrayJ<-marray(CHJ)
#Join 2 arrays together both Adults and Juveniles
m<-rbind(marrayJ,marrayA


# Specify model in JAGS and BUGS language
sink("ipm.dormice.jags")
cat("
    model {
    #------------------------------------------------------------
    #  Integrated population model
    #  - Age structured model with 2 age classes: 
    #    Juveniles under 2 months old Adults >2months
    #  - Age at first breeding = 1 year
    #  - All vital rates are assumed to be time-dependent (random)
    #-------------------------------------------------------------
    
    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    # Initial population sizes
    n1 ~ dnorm(100, 0.0001)I(0,)           # Adult individuals
    nadSurv ~ dnorm(100, 0.0001)I(0,)      # Adults >= 2 years
    nadimm ~ dnorm(100, 0.0001)I(0,)       # Immigrants
    N1[1] <- round(n1)
    NadSurv[1] <- round(nadSurv)
    Nadimm[1] <- round(nadimm)
    
    # Mean demographic parameters (on appropriate scale)
    l.mphij ~ dnorm(0, 0.0001)
    l.mphia ~ dnorm(0, 0.0001)
    l.mfec ~ dnorm(0, 0.0001)
    l.mim ~ dnorm(0, 0.0001)
    l.p ~ dnorm(0, 0.0001)
    
    # Precision of standard deviations of temporal variability
    sig.phij ~ dunif(0, 10)
    tau.phij <- pow(sig.phij, -2)
    sig.phia ~ dunif(0, 10)
    tau.phia <- pow(sig.phia, -2)
    sig.fec ~ dunif(0, 10)
    tau.fec <- pow(sig.fec, -2)
    sig.im ~ dunif(0, 10)
    tau.im <- pow(sig.im, -2)
sig.obs~ dunif(0.5,50)
tau.obs<- pow(sig.obs,-2)
    
    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)	
    epsilon.phia[t] ~ dnorm(0, tau.phia)
    epsilon.fec[t] ~ dnorm(0, tau.fec)
    epsilon.im[t] ~ dnorm(0, tau.im)
    }
    
    #-------------------------
    # 2. Constrain parameters
    #-------------------------
    for (t in 1:(nyears-1)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]  # Juv. apparent survival
    logit(phia[t]) <- l.mphia + epsilon.phia[t]  # Adult apparent survival
    log(f[t]) <- l.mfec + epsilon.fec[t]         # Productivity
    log(omega[t]) <- l.mim + epsilon.im[t]       # Immigration
    logit(p[t]) <- l.p                           # Recapture probability
    }
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))   # Mean juvenile dormice survival probability
    mphia <- exp(l.mphia)/(1+exp(l.mphia))   # Mean adult dormice survival probability
    mfec <- exp(l.mfec)                      # Mean productivity/fecundity
    mim <- exp(l.mim)                        # Mean immigration rate
    
    # Population growth rate Lamda
    for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / Ntot[t]
    logla[t] <- log(lambda[t])
    }
    mlam <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))   # Geometric mean
    
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model see Fig.2)
    # 4.1.1 System process
    for (t in 2:nyears){
    mean1[t] <- 0.5 * f[t-1] * phij[t-1] * Ntot[t-1]
    N1[t] ~ dpois(mean1[t])
    NadSurv[t] ~ dbin(phia[t-1], Ntot[t-1])
    mpo[t] <- Ntot[t-1] * omega[t-1]
    Nadimm[t] ~ dpois(mpo[t])
    }
    
    # 4.1.2 Observation process
    for (t in 1:nyears){
    Ntot[t] <- NadSurv[t] + Nadimm[t] + N1[t]
    y[t] ~ dpois(Ntot[t])
    }
    
    # 4.2 Likelihood for capture-recapture data: CJS model (2 age classes Adult of Juvenile)
    # Multinomial likelihood
    for (t in 1:(nyears-1)){
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t])
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-1)){
    q[t] <- 1-p[t]
    # Main diagonal
    pr.j[t,t] <- phij[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.j[t,j] <- phij[t]*prod(phia[(t+1):j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults
    for (t in 1:(nyears-1)){
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    
    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t] * f[t]
    }
    }
    ",fill = TRUE)
sink()


# Bundle data to run Jags model
jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = popcount, 
                  J = J, R = R, r.j = rowSums(marray.j), r.a = rowSums(marray.a))
# Run the model in JAGS from R using jagsUI package speeds it up --------------------------
library(jagsUI)
# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), 
                    l.mfec = rnorm(1, 0.2, 0.5), l.mim = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), 
                    sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1,10), 
                    sig.im = runif(1, 0.1, 10), n1 = round(runif(1, 1, 50), 0), 
                    nadSurv = round(runif(1, 5, 50), 0), 
                    nadimm = round(runif(1, 1, 50), 0))}
# Parameters monitored
parameters <- c("phij", "phia", "f", "omega", "p", "lambda", "mphij", "mphia", "mfec",
                "mim", "mlam", "sig.phij", "sig.phia", "sig.fec", "sig.im", "N1", "NadSurv", "Nadimm", "Ntot")

# MCMC settings
ni <- 400000
nt <- 10
nb <- 50000
nc <- 3


# Call JAGS from R (long run 20 mins)
ipm <- jags(jags.data, inits, parameters, "ipm.dormice.jags", n.chains = nc, n.thin = nt,
                   n.iter = ni, n.burnin = nb, parallel=TRUE) #Add parallel=TRUE to increase speed

