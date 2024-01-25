model{
  ## PRIORS
  # hyperpriors
  mu.lpsi1 <- logit(mean.psi1)
  mean.psi1 ~ dunif(0,1)
  tau.lpsi1 <- pow(sd.lpsi1, -2)
  sd.lpsi1 ~ dunif(0,10)
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dunif(0,1)
  tau.lphi <- pow(sd.lphi, -2)
  sd.lphi ~ dunif(0,10)
  mu.lgamma <- logit(mean.gamma)
  mean.gamma ~ dunif(0,1)
  tau.lgamma <- pow(sd.lgamma, -2)
  sd.lgamma ~ dunif(0, 10)
  mu.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 10)
  # priors for random species-level effects
  for(k in 1:nspec) {   
    logit(psi1[k]) <- lpsi1[k]
    lpsi1[k] ~ dnorm(mu.lpsi1, tau.lpsi1)
    logit(phi[k]) <- lphi[k]
    lphi[k] ~ dnorm(mu.lphi, tau.lphi)
    logit(gamma[k]) <- lgamma[k]
    lgamma[k] ~ dnorm(mu.lgamma, tau.lgamma)
    logit(p[k]) <- lp[k]
    lp[k] ~ dnorm(mu.lp, tau.lp)
  }
  # ECOLOGICAL SUBMODEL
  for (i in 1:nsite) {
    for(k in 1:nspec) {
      z[i,1,k] ~ dbern(psi1[k])
      for (t in 2:nyear) {
        z[i,t,k] ~ dbern(z[i,t-1,k] * phi[k] + (1-z[i,t-1,k]) * gamma[k])
      }
    }
  }
  # OBSERVATION SUBMODEL
  for (i in 1:nsite) {
    for (k in 1:nspec) {
      for (j in 1:nsurvey) {
        for (t in 1:nyear) {
          y[i,j,t,k] ~ dbern(z[i,t,k] * p[k])
        }
      }
    }
  }
  # DERIVED PARAMETERS
  # Number of occupied sites & population occupancy
  for (k in 1:nspec) {
    n.occ[1,k] <- sum(z[,1,k])    # Number of occupied sites
    psi[1,k] <- psi1[k]           # Population occupancy
    for (t in 2:nyear) {
      n.occ[t,k] <-sum(z[,t,k])
      psi[t,k] <- psi[t-1,k] * phi[k] + (1-psi[t-1,k] * gamma[k])
    }
  }
}