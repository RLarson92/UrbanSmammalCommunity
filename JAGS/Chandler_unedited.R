model{
  int.mu.beta0 ~ dnorm(0, 0.1)
  int.mu.beta1 ~ dnorm(0, 0.1)
  int.mu.beta2 ~ dnorm(0, 0.1)
  int.mu.beta3 ~ dnorm(0, 0.1)
  int.mu.beta4 ~ dnorm(0, 0.1)
  x.mu.beta0 ~ dnorm(0, 0.1)
  x.mu.beta1 ~ dnorm(0, 0.1)
  x.mu.beta2 ~ dnorm(0, 0.1)
  x.mu.beta3 ~ dnorm(0, 0.1)
  x.mu.beta4 ~ dnorm(0, 0.1)
  sig.beta0 ~ dunif(0, 5)
  tau.beta0 <- 1/(sig.beta0 * sig.beta0)
  sig.beta1 ~ dunif(0, 5)
  tau.beta1 <- 1/(sig.beta1 * sig.beta1)
  sig.beta2 ~ dunif(0, 5)
  tau.beta2 <- 1/(sig.beta2 * sig.beta2)
  sig.beta3 ~ dunif(0, 5)
  tau.beta3 <- 1/(sig.beta3 * sig.beta3)
  # mu.beta4 ~ dnorm(0, 0.1)
  sig.beta4 ~ dunif(0, 5)
  tau.beta4 <- 1/(sig.beta4 * sig.beta4)
  mu.logitp ~ dnorm(0, 0.1)
  sig.logitp ~ dunif(0, 5)
  tau.logitp <- 1/(sig.logitp * sig.logitp)
  alpha1 ~ dnorm(0, 0.1)
  alpha2 ~ dnorm(0, 0.1)
  alpha3 ~ dnorm(0, 0.1)
  alpha4 ~ dnorm(0, 0.1)
  alpha5 ~ dnorm(0, 0.1)
  for (g in 1:5) { # 5 forest-dependence classes
    # expected abundance as a function of forest-dependence class (ordinal)
    mu.beta0[g] <- int.mu.beta0 + x.mu.beta0 * (g-1)
    mu.beta1[g] <- int.mu.beta1 + x.mu.beta1 * (g-1)
    mu.beta2[g] <- int.mu.beta2 + x.mu.beta2 * (g-1)
    mu.beta3[g] <- int.mu.beta3 + x.mu.beta3 * (g-1)
    mu.beta4[g] <- int.mu.beta4 + x.mu.beta4 * (g-1)
  }
  for (i in 1:148) {
    # species-specific coefficients
    beta0[i] ~ dnorm(mu.beta0[depindex[i]], tau.beta0)
    beta1[i] ~ dnorm(mu.beta1[depindex[i]], tau.beta1)
    beta2[i] ~ dnorm(mu.beta2[depindex[i]], tau.beta2)
    beta3[i] ~ dnorm(mu.beta3[depindex[i]], tau.beta3)
    beta4[i] ~ dnorm(mu.beta4[depindex[i]], tau.beta4)
    alpha0[i] ~ dnorm(mu.logitp, tau.logitp)
    for (j in 1:32) {
      log(lambda[i,j]) <- beta0[i] + 
        beta1[i]*ioc[j] + 
        beta2[i]*shade[j] + 
        beta3[i]*secondary[j] +
        beta4[i]*distDivide[j]
      logit(p[i,j,1]) <- alpha0[i] +
        alpha1*nethours[j,1] +
        alpha2*wind[j,1] +
        alpha3*ioc[j] +
        alpha4*shade[j] +
        alpha5*secondary[j]
      logit(p[i,j,2]) <- alpha0[i] +
        alpha1*nethours[j,2] +
        alpha2*wind[j,2] +
        alpha3*ioc[j] +
        alpha4*shade[j] +
        alpha5*secondary[j]
      logit(p[i,j,3]) <- alpha0[i] +
        alpha1*nethours[j,3] +
        alpha2*wind[j,3] +
        alpha3*ioc[j] +
        alpha4*shade[j] +
        alpha5*secondary[j]
      N[i,j] ~ dpois(lambda[i,j])
      y[i,j,1] ~ dbin(p[i,j,1], N[i,j])
      N2[i,j] <- N[i,j] - y[i,j,1]
      y[i,j,2] ~ dbin(p[i,j,2], N2[i,j])
      N3[i,j] <- N2[i,j] - y[i,j,2]
      y[i,j,3] ~ dbin(p[i,j,3], N3[i,j])
    }
    # total abundance in each habitat type
    Nh[i,1] <- N[i,] %*% shade
    Nh[i,2] <- N[i,] %*% ioc
    Nh[i,3] <- N[i,] %*% secondary
    Nh[i,4] <- N[i,] %*% primary
  }
}