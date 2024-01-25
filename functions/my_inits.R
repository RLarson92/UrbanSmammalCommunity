# initial starting values for each chain
my_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      mu.beta0 = rnorm(1),
      tau.beta0 = rgamma(1,1,1),
      mu.beta1 = rnorm(1),
      tau.beta1 = rgamma(1,1,1),
      mu.beta2 = rnorm(1),
      tau.beta2 = rgamma(1,1,1),
      mu.beta3 = rnorm(1),
      tau.beta3 = rgamma(1,1,1),
      mu.alpha0 = rnorm(1),
      tau.alpha0 = rgamma(1,1,1),
      alpha1 = rnorm(1),
      alpha2 = rnorm(1),
      alpha3 = rnorm(1),
      gamma = rnorm(1),
      mu.phi = rnorm(1),
      tau.phi = rgamma(1,1,1),
      beta0 = rnorm(data_list$nspec),
      beta1 = rnorm(data_list$nspec),
      beta2 = rnorm(data_list$nspec),
      beta3 = rnorm(data_list$nspec),
      alpha0 = rnorm(data_list$nspec),
      phi = rnorm(data_list$nspec),
      N = array(10, 
                dim = c(data_list$nspec, data_list$nsite, data_list$nseason)
                 ),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(
    switch(chain,
           "1" = gen_list(chain),
           "2" = gen_list(chain),
           "3" = gen_list(chain),
           "4" = gen_list(chain),
           "5" = gen_list(chain),
           "6" = gen_list(chain),
           "7" = gen_list(chain),
           "8" = gen_list(chain)
    )
  )
}