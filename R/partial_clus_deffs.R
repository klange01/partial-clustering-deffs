### Design effects based on GEEs for partially clustered trials
### Author: Kylie Lange (kylie.lange@adelaide.edu.au)

##### HELPER FUNCTIONS #####

# Calculate the total number of obs from the number of clusters
# Input: a list of the number of clusters of sizes 1..K
# Output: total number of observations, M
numclus_to_N = function(M_k){
  N_k     = M_k %>% imap(~ .y*.x)    # N_k = k*M_k
  N       = sum(unlist(N_k))
  return(N)
}

# Calculate the proportions of observations in clusters of each size
# Input: a list of the number of clusters of size 1..K
# Output: a list of the proportions gamma_k
numclus_to_gamma = function(M_k){
  N     = numclus_to_N(M_k)
  gamma = M_k %>% imap(~ .y*.x/N)   # gamma_k = k*M_k/N
  return(gamma)
}

##### DESIGN EFFECTS #####

# Inputs: randomisation method ('cluster' or 'ind'), rho, list of proportions gamma_k
# Output: design effect

# Independence GEE for a continuous outcome
deff_cont_indgee = function(rand_method, rho, gamma) {
  if (rand_method == "cluster") {
    sum1 = gamma %>% imap(~ .y*.x)  # k*gamma_k
    deff = 1 + rho*(sum(unlist(sum1)) - 1)
  }
  else if (rand_method == "ind") {
    deff = 1
  }
  return(deff)
}

# Exchangeable GEE for a continuous outcome
deff_cont_exchgee = function(rand_method, rho, gamma) {
  if (rand_method == "cluster") {
    sum1 = gamma %>% imap(~ (1/(1 + (.y-1)*rho))*.x)  # 1/(1+(k-1)rho)*gamma_k
    deff = 1 / sum(unlist(sum1))
  }
  else if (rand_method == "ind") {
    sum1 = gamma %>% imap(~ (1 + (.y-2)*rho)/((1 - rho)*(1 + (.y-1)*rho))*.x)  # (1+(k-2)rho)/((1-rho)*(1+(k-1)rho))*gamma_k
    deff = 1 / sum(unlist(sum1))
  }
  return(deff)
}

# Independence GEE with a logit link for a binary outcome 
deff_logit_indgee = function(rand_method, rho, gamma, pi_C, pi_I) {
  if (rand_method == "cluster") {
    sum1 = gamma %>% imap(~ (.y-1)*.x)   # (k-1)*gamma_k
    deff = 1 + rho*(sum(unlist(sum1)))
  }
  else if (rand_method == "ind") {
    sum1 = gamma %>% imap(~ (.y-1)*.x)   # (k-1)*gamma_k
    pi_term = (sqrt(pi_I*pi_C*(1-pi_I)*(1-pi_C))) / (pi_I*(1-pi_I)+pi_C*(1-pi_C))
    deff = 1 + rho*((0.5 - pi_term)*sum(unlist(sum1)))
  }
  return(deff)
}

# Exchangeable GEE with a logit link for a binary outcome 
deff_logit_exchgee = function(rand_method, rho, gamma, pi_C, pi_I) {
  if (rand_method == "cluster") {
    sum1 = gamma %>% imap(~ 1/(1 + (.y-1)*rho)*.x)   # 1/(1+(k-1)rho)*gamma_k
    deff = 1 / sum(unlist(sum1))
  }
  else if (rand_method == "ind") {
    sum1 = gamma %>% imap(~ 1/(1 + (.y-1)*rho)*.x)       # 1/(1+(k-1)rho)*gamma_k
    sum2 = gamma %>% imap(~ (.y-1)/(1 + (.y-1)*rho)*.x)  # (k-1)/(1+(k-1)rho)*gamma_k
    pi_term = (sqrt(pi_I*pi_C*(1-pi_I)*(1-pi_C))) / (pi_I*(1-pi_I)+pi_C*(1-pi_C))
    deff = (sum(unlist(sum1)) + (0.5 - pi_term)*rho/(1-rho)*sum(unlist(sum2))) / 
      (sum(unlist(sum1))**2 + rho/(1-rho)*sum(unlist(sum1))*sum(unlist(sum2)))
  }
  return(deff)
}

# Independence GEE with a log link for a binary outcome 
deff_log_indgee = function(rand_method, rho, gamma, pi_C, pi_I) {
  if (rand_method == "cluster") {
    sum1 = gamma %>% imap(~ (.y-1)*.x)   # (k-1)*gamma_k
    deff = 1 + rho*(sum(unlist(sum1)))
  }
  else if (rand_method == "ind") {
    sum1 = gamma %>% imap(~ (.y-1)*.x)   # (k-1)*gamma_k
    pi_term = (sqrt(pi_I*pi_C*(1-pi_I)*(1-pi_C))) / (pi_I*(1-pi_C)+pi_C*(1-pi_I))
    deff = 1 + rho*((0.5 - pi_term)*sum(unlist(sum1)))
  }
  return(deff)
}

# Exchangeable GEE with a log link for a binary outcome 
deff_log_exchgee = function(rand_method, rho, gamma, pi_C, pi_I) {
  if (rand_method == "cluster") {
    sum1 = gamma %>% imap(~ 1/(1 + (.y-1)*rho)*.x)   # 1/(1+(k-1)rho)*gamma_k
    deff = 1 / sum(unlist(sum1))
  }
  else if (rand_method == "ind") {
    sum1 = gamma %>% imap(~ 1/(1 + (.y-1)*rho)*.x)       # 1/(1+(k-1))*gamma_k
    sum2 = gamma %>% imap(~ (.y-1)/(1 + (.y-1)*rho)*.x)  # (k-1)/(1+(k-1)rho)*gamma_k
    pi_term = (sqrt(pi_I*pi_C*(1-pi_I)*(1-pi_C))) / (pi_I*(1-pi_C)+pi_C*(1-pi_I))
    deff = (sum(unlist(sum1)) + (0.5 - pi_term)*rho/(1-rho)*sum(unlist(sum2))) / 
      (sum(unlist(sum1))**2 + rho/(1-rho)*sum(unlist(sum1))*sum(unlist(sum2)))
  }
  return(deff)
}


# Calculate the expected design effect, effective N and expected power for GEEs for a continuous or binary outcome
# Inputs are a list of trial design parameters, including:
# the number of clusters of each cluster size, the ICC, the randomisation method, and treatment effect
# Expected power = power of a Wald test for independent data for the specified effect size and the effective sample size
# Effective sample size = number of observations / expected deff

# Independence GEE for a continuous outcome
calculate_cont_indgee_deff = function(params) {
  # expected deff
  M_k   = c(params$num_clus1, params$num_clus2, params$num_clus3, params$num_clus4)
  gamma = numclus_to_gamma(M_k)
  deff  = deff_cont_indgee(rand_method = params$rand_method, rho = params$ICC, gamma = gamma)
  # expected power
  N         = numclus_to_N(M_k)
  N_eff     = N / deff
  power_exp = pwrss.t.2means(mu1=params$b0, mu2=params$b1, sd1=1, sd2=1, kappa=1, n2=N_eff/2, alpha=0.05, verbose=FALSE)$power
  res       = data.frame(deff, N_eff, power_exp)
  return(res)
}

# Exchangeable GEE for a continuous outcome
calculate_cont_exchgee_deff = function(params) {
  # expected deff
  M_k   = c(params$num_clus1, params$num_clus2, params$num_clus3, params$num_clus4)
  gamma = numclus_to_gamma(M_k)
  deff  = deff_cont_exchgee(rand_method = params$rand_method, rho = params$ICC, gamma = gamma)
  # expected power
  N         = numclus_to_N(M_k)
  N_eff     = N / deff
  power_exp = pwrss.t.2means(mu1=params$b0, mu2=params$b1, sd1=1, sd2=1, kappa=1, n2=N_eff/2, alpha=0.05, verbose=FALSE)$power
  res       = data.frame(deff, N_eff, power_exp)
  return(res)
}

# Independence GEE with a logit link for a binary outcome 
calculate_logit_indgee_deff = function(params) {
  # expected deff
  M_k   = c(params$num_clus1, params$num_clus2, params$num_clus3, params$num_clus4)
  gamma = numclus_to_gamma(M_k)
  deff  = deff_logit_indgee(rand_method = params$rand_method, rho = params$ICC, gamma = gamma, pi_C = params$pC, pi_I = params$pI)
  # expected power
  N         = numclus_to_N(M_k)
  N_eff     = N / deff
  power_exp = pwrss.z.logreg(p0=params$pC, p1=params$pI, alpha=0.05, distribution="bernoulli", n=N_eff, verbose=FALSE)$power
  res       = data.frame(deff, N_eff, power_exp)
  return(res)
}

# Exchangeable GEE with a logit link for a binary outcome 
calculate_logit_exchgee_deff = function(params) {
  # expected deff
  M_k   = c(params$num_clus1, params$num_clus2, params$num_clus3, params$num_clus4)
  gamma = numclus_to_gamma(M_k)
  deff  = deff_logit_exchgee(rand_method = params$rand_method, rho = params$ICC, gamma = gamma, pi_C = params$pC, pi_I = params$pI)
  # expected power
  N         = numclus_to_N(M_k)
  N_eff     = N / deff
  power_exp = pwrss.z.logreg(p0=params$pC, p1=params$pI, alpha=0.05, distribution="bernoulli", n=N_eff, verbose=FALSE)$power
  res       = data.frame(deff, N_eff, power_exp)
  return(res)
}

# Independence GEE with a log link for a binary outcome 
calculate_log_indgee_deff = function(params) {
  # expected deff
  M_k   = c(params$num_clus1, params$num_clus2, params$num_clus3, params$num_clus4)
  gamma = numclus_to_gamma(M_k)
  deff  = deff_log_indgee(rand_method = params$rand_method, rho = params$ICC, gamma = gamma, pi_C = params$pC, pi_I = params$pI)
  # expected power
  N         = numclus_to_N(M_k)
  N_eff     = N / deff
  power_exp = pwrss.z.2props(p1=params$pC, p2=params$pI, alpha=0.05, kappa=1, n2=N_eff/2, verbose=FALSE)$power
  res       = data.frame(deff, N_eff, power_exp)
  return(res)
}

# Exchangeable GEE with a log link for a binary outcome 
calculate_log_exchgee_deff = function(params) {
  # expected deff
  M_k   = c(params$num_clus1, params$num_clus2, params$num_clus3, params$num_clus4)
  gamma = numclus_to_gamma(M_k)
  deff  = deff_log_exchgee(rand_method = params$rand_method, rho = params$ICC, gamma = gamma, pi_C = params$pC, pi_I = params$pI)
  # expected power
  N         = numclus_to_N(M_k)
  N_eff     = N / deff
  power_exp = pwrss.z.2props(p1=params$pC, p2=params$pI, alpha=0.05, kappa=1, n2=N_eff/2, verbose=FALSE)$power
  res       = data.frame(deff, N_eff, power_exp)
  return(res)
}
