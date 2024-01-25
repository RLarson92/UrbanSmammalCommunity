model{
  ### PRIORS
  ## Hyperpriors (overall average distributions to pull species-specific responses from)
  # for abundance
  mu.beta0 ~ dnorm(0, 0.1)
  tau.beta0 ~ dgamma(1,1)
  sig.beta0 <- 1 / sqrt(tau.beta0)
  mu.beta1 ~ dnorm(0, 0.1)
  tau.beta1 ~ dgamma(1,1)
  sig.beta1 <- 1 / sqrt(tau.beta1)
  mu.beta2 ~ dnorm(0, 0.1)
  tau.beta2 ~ dgamma(1,1)
  sig.beta2 <- 1 / sqrt(tau.beta2)
  mu.beta3 ~ dnorm(0, 0.1)
  tau.beta3 ~ dgamma(1,1)
  sig.beta3 <- 1 / sqrt(tau.beta3)
  # for captures
  mu.alpha0 ~ dnorm(0, 0.1)
  tau.alpha0 ~ dgamma(1,1)
  sig.alpha0 <- 1 / sqrt(tau.alpha0)
  alpha1 ~ dnorm(0, 0.1)                  # technically not hyperpriors? but here for convenience
  alpha2 ~ dnorm(0, 0.1)
  alpha3 ~ dnorm(0, 0.1) 
  # for recruitment
  gamma ~ dnorm(0, 0.1)
  # for persistence
  mu.phi ~ dnorm(0, 0.1)
  tau.phi ~ dgamma(1,1)
  sig.phi <- 1 / sqrt(tau.phi)
  ## Species-Specific Priors
  # species-specific coefficients
  for (i in 1:nspec) {
    # for abundance
    beta0[i] ~ dnorm(mu.beta0, tau.beta0)
    beta1[i] ~ dnorm(mu.beta1, tau.beta1)
    beta2[i] ~ dnorm(mu.beta2, tau.beta2)
    beta3[i] ~ dnorm(mu.beta3, tau.beta3)
    # for captures
    alpha0[i] ~ dnorm(mu.alpha0, tau.alpha0)
    # for persistence
    phi[i] ~ dnorm(mu.phi, tau.phi)
  }
  
  ### MODEL
  ## Season 1
  # Latent State (Abundance)
  for (i in 1:nspec) {
    for (j in 1:nsite) {
      N[i,j,1] ~ dpois(lambda[i,j,1])
      log(lambda[i,j,1]) <- beta0[i] + beta1[i]*PC1[j] + beta2[i]*PC2[j] + beta3[i]*contag[j]
  # Observation State (Captures)
      for (k in 1:nNight) {
        logit(p[i,j,k,1]) <- alpha0[i] + alpha1*moon[j,k,1] + alpha2*jDate[j,k,1] + alpha3*effort[j,k,1]
        y[i,j,k,1] ~ dbin(p[i,j,k,1], N[i,j,1])
      }
  ## Season > 1
  # Latent State
      for (t in 2:nseason) {
        N[i,j,t] ~ dpois(lambda[i,j,t])
        log(lambda[i,j,t]) <- beta0[i] + beta1[i]*PC1[j] + beta2[i]*PC2[j] + beta3[i]*contag[j] + gamma*season[t] + phi[i]*N[i,j,t-1]
  # Observation State
        for (k in 1:nNight) {
          logit(p[i,j,k,t]) <- alpha0[i] + alpha1*moon[j,k,t] + alpha2*jDate[j,k,t] + alpha3*effort[j,k,t]
          y[i,j,k,t] ~ dbin(p[i,j,k,t], N[i,j,t])
        }
      }
    }
  }
}