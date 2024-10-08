## Validation of GEE design effects for partially clustered trials
## Kylie Lange (kylie.lange@adelaide.edu.au)

##### GENERATING DATA #####

# Generate a random starting seed based on system time
generate_seed = function() {
  my_seed = as.integer(Sys.time()) %% 100000
  return(my_seed)
}

# Create a dataframe for a specified trial design
# Under cluster randomisation, treatment groups are balanced within cluster size.
# Under individual randomisation, treatment groups are balanced within cluster size,
# but the number of clusters of each within-cluster allocation type is not fixed.
# Inputs: randomisation method, number of clusters of each cluster size
# Outputs: a dataframe of observation and cluster ids and treatment allocation
create_trial = function(rand_method, num_clus1, num_clus2, num_clus3, num_clus4) {
  
  k        = c(1, 2, 3, 4)                                 # cluster sizes
  M_k      = c(num_clus1, num_clus2, num_clus3, num_clus4) # numbers of clusters of each size
  num_clus = sum(unlist(M_k))                              # total number of clusters
  N_k      = M_k %>% imap(~ .y*.x)                         # number of obs within each cluster size
  
  # create the clusters
  clusters = rep(k, times = M_k)
  dat = data.frame(cluster_id = seq(1:num_clus), clus_size = clusters)
  
  # create the treatment indicator variable
  if (rand_method == "cluster") {
    tmt = rep(c(0,1), length = num_clus)   # cluster level assignment is balanced by cluster size as the clusters are in size order
    dat$tmt = tmt
    dat = dat %>% uncount(clus_size, .remove = FALSE)   # expand the clusters into observation level data
    
  } else if (rand_method == "ind") {
    dat = dat %>% uncount(clus_size, .remove = FALSE)   # expand the clusters into observation level data
    tmt = unlist(lapply(N_k, function(x) (sample(rep(c(0, 1), length = x)))))  # generate balanced treatment assignments for each cluster size
    dat$tmt = tmt
  }
  
  dat = dat %>%
    group_by(cluster_id) %>%
    mutate(withinclus_id = row_number(cluster_id),
           t = sum(tmt)) %>%
    ungroup() %>% 
    mutate(obs_id = seq(1:nrow(.))) %>% 
    select(obs_id, cluster_id, withinclus_id, clus_size, t, tmt)
  return(dat)
}

# Create a dataframe for a specified trial design
# A relaxed version of create_trial() where the assumption
# of treatment balance within cluster size is not enforced.
# The assumption of overall treatment balance is still enforced.
# Inputs: randomisation method, number of clusters of each cluster size
# Outputs: a dataframe of observation and cluster ids and treatment allocation
create_trial_relaxed = function(rand_method, num_clus1, num_clus2, num_clus3, num_clus4) {
  
  k        = c(1, 2, 3, 4)                                 # cluster sizes
  M_k      = c(num_clus1, num_clus2, num_clus3, num_clus4) # numbers of clusters of each size
  num_clus = sum(unlist(M_k))                              # total number of clusters
  N_k      = M_k %>% imap(~ .y*.x)                         # number of obs within each cluster size
  
  # create the clusters
  clusters = rep(k, times = M_k)
  dat = data.frame(cluster_id = seq(1:num_clus), clus_size = clusters)
  
  # create the treatment indicator variable
  if (rand_method == "cluster") {
    tmt = sample(rep(c(0,1), length = num_clus))        # cluster level assignment with overall treatment balance
    dat$tmt = tmt
    dat = dat %>% uncount(clus_size, .remove = FALSE)   # expand the clusters into observation level data
    
  } else if (rand_method == "ind") {
    dat = dat %>% uncount(clus_size, .remove = FALSE)   # expand the clusters into observation level data
    tmt = sample(rep(c(0,1), length = nrow(dat)))       # individual level assignment with overall treatment balance
    dat$tmt = tmt
  }
  
  dat = dat %>%
    group_by(cluster_id) %>%
    mutate(withinclus_id = row_number(cluster_id),
           t = sum(tmt)) %>%
    ungroup() %>% 
    mutate(obs_id = seq(1:nrow(.))) %>% 
    select(obs_id, cluster_id, withinclus_id, clus_size, t, tmt)
  return(dat)
}

# Generate a single simulated dataset with a continuous outcome
# Inputs: simulation parameters
# Output: a dataframe
generate_cont_data = function(params) {

  # create dataframe of observations and the treatment variable
  dat = create_trial(params$rand_method, params$num_clus1, params$num_clus2, params$num_clus3, params$num_clus4)
  num_clus   = n_distinct(dat$cluster_id)
  clus_sizes = dat %>% summarise(clus_sizes=first(clus_size), .by=cluster_id) %>% select(clus_sizes)

  # random variables
  sd_alpha = sqrt(params$ICC)
  sd_eps   = sqrt(1 - params$ICC)
  alpha_j  = rnorm(num_clus, 0, sd_alpha)
  alpha_j  = rep(alpha_j, times = unlist(clus_sizes))
  eps_ij   = rnorm(nrow(dat), 0, sd_eps)
  dat      = data.frame(dat, alpha_j, eps_ij)

  # response variable
  dat = dat %>% mutate(Y = params$b0 + params$b1*tmt + alpha_j + eps_ij)
  
  # attach sim params
  dat = dat %>% 
    mutate(rand_method = params$rand_method,
           M1          = params$num_clus1,
           M2          = params$num_clus2,
           M3          = params$num_clus3,
           M4          = params$num_clus4,
           ICC         = params$ICC,
           trueb0      = params$b0,
           trueb1      = params$b1)
  return(dat)
}

# Gather the prevalence values for a cluster into one vector, as required for cBern()
gather_p = function(clusterdat){
  p = unlist(as.vector(clusterdat[2:5]))
  p = p[!is.na(p)]
}

# Generate a single simulated dataset with a binary outcome
# Data generated using method of Jiang et al, 2021 (doi: 10.1080/00031305.2020.1816213)
# as implemented in the CorBin package (https://cran.r-project.org/web/packages/CorBin/index.html)
# Inputs: simulation parameters
# Output: a dataframe
generate_binary_data = function(params) {
  
  # create dataframe of observations and treatment variable
  dat = create_trial(params$rand_method, params$num_clus1, params$num_clus2, params$num_clus3, params$num_clus4)
  num_clus   = n_distinct(dat$cluster_id)
  clus_sizes = dat %>% summarise(clus_sizes=first(clus_size), .by=cluster_id) %>% select(clus_sizes)
  
  # p = expected prevalence of outcome
  dat = dat %>% mutate(p = case_when(tmt == 0 ~ params$pC,
                                     tmt == 1 ~ params$pI))

  # create list of prevalence values for each cluster
  datw = dat %>% 
    select(c(cluster_id, withinclus_id, p)) %>% 
    pivot_wider(names_from = withinclus_id, values_from = p)
  p_vec = apply(datw, 1, gather_p)
  
  # generate outcome for each cluster
  res   = lapply(p_vec, cBern, n = 1, rho = params$ICC, type = "exchange")
  resdf = bind_rows(lapply(res, as.data.frame)) # one row per cluster
  resdf = resdf %>% mutate(cluster_id = row_number())
  
  # pivot to one row per observation and merge into dat
  ydat = resdf %>% 
    rename("1" = "V1", "2" = "V2", "3" = "V3", "4" = "V4") %>% 
    pivot_longer(cols           = c(1, 2, 3, 4),
                 names_to       = 'withinclus_id',
                 values_to      = 'Y',
                 values_drop_na = TRUE) %>% 
    mutate_if(is.character, as.numeric)
  dat = dat %>% 
    left_join(ydat, c("cluster_id", "withinclus_id"))
  
  # attach sim params
  dat = dat %>% 
    mutate(rand_method = params$rand_method,
           M1          = params$num_clus1,
           M2          = params$num_clus2,
           M3          = params$num_clus3,
           M4          = params$num_clus4,
           ICC         = params$ICC,
           truepC      = params$pC,
           truepI      = params$pI)
  return(dat)
}

# Generate multiple simulated datasets for a single scenario
# Inputs: simulation parameters for a single scenario
# Outputs: a list of datasets created by generate_data() and a matrix of random seeds for each dataset
generate_multiple_datasets = function(sims, params, outcome) {
  # initialise
  nsims = 0
  simulated_data = list()
  states = matrix(ncol=626, nrow=sims)
  # set the generate data function depending on outcome type
  if (outcome == "cont") {
    generate_data = function() { generate_cont_data(params) }
  }
  else if (outcome == "binary"){
    generate_data = function() { generate_binary_data(params) }
  }
  # simulation loop
  while (nsims < sims) {
    states[nsims + 1, ] = .Random.seed
    simulated_data[[nsims + 1]] = generate_data()
    nsims = nsims + 1
    }
  simulated_data = Map(cbind, simulated_data, index=seq_along(simulated_data))
  list(data=simulated_data, states=states)
}

# Run simulation for one scenario
# Inputs: a dataframe containing the simulation parameters; number of simulated datasets to create per scenario; type of outcome (cont/binary)
# Outputs: a list of lists of datasets created by generate_data()
# Data and random seeds are saved to file
simulate = function(params, nsims, outcome) {
  # prep params
  params = params %>% mutate_at(c("num_clus1", "num_clus2", "num_clus3", "num_clus4"), as.integer)
  print('Simulating data with params:')
  print(params)
  # simulate the datasets
  result = generate_multiple_datasets(nsims, params, outcome)
  data = result$data
  data = Map(cbind, scenario_id=params$scenario_id, data)
  # save the results
  saveRDS(data, file=paste0("simdata/data_", params$scenario_id, ".rds"))
  saveRDS(result$states, file=paste0("simdata/states_", params$scenario_id, ".rds"))
  return(data)
#  }
}

###### MODEL FITTING ######

# Fit models for partially clustered data with a continuous outcome to a dataset
# Inputs: a dataset as produced by generate_data()
# Output: a dataframe of simulation parameters and model estimates, one row per analysis
# Analyses: geeglm, lm
fit_cont_models = function(dat) {
  # logging
  print(dat[1,c("scenario_id", "index")]) 
  print(Sys.time())
  # sim params
  d = data.frame(scenario_id = dat$scenario_id[1],
                 dataset     = dat$index[1],
                 rand_method = dat$rand_method[1],
                 nobs        = nrow(dat),
                 M1          = dat$M1[1],
                 M2          = dat$M2[1],
                 M3          = dat$M3[1],
                 M4          = dat$M4[1],
                 ICC         = dat$ICC[1],
                 trueb0      = dat$trueb0[1],
                 trueb1      = dat$trueb1[1],
                 stringsAsFactors = FALSE, row.names = NULL)

  # independence gee via geepack::geeglm
  m1 = geeglm(Y ~ tmt, id=cluster_id, family=gaussian, corstr="independence", data=dat)
  s1 = summary(m1)
  d1 = data.frame(model  = "gee-ind",
                  b0     = s1$coefficients["(Intercept)", "Estimate"],
                  se_b0  = s1$coefficients["(Intercept)", "Std.err"],
                  b1     = s1$coefficients["tmt", "Estimate"],
                  se_b1  = s1$coefficients["tmt", "Std.err"],
                  var_b1 = m1$geese$vbeta[2,2],
                  p_b1   = s1$coefficients["tmt", "Pr(>|W|)"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d1 = bind_cols(d, d1)

  # exchangeable gee via geepack::geeglm
  m2 = geeglm(Y ~ tmt, id=cluster_id, family=gaussian, corstr="exchangeable", data=dat)
  s2 = summary(m2)
  d2 = data.frame(model    = "gee-exch",
                  b0       = s2$coefficients["(Intercept)", "Estimate"],
                  se_b0    = s2$coefficients["(Intercept)", "Std.err"],
                  b1       = s2$coefficients["tmt", "Estimate"],
                  se_b1    = s2$coefficients["tmt", "Std.err"],
                  var_b1   = m2$geese$vbeta[2,2],
                  p_b1     = s2$coefficients["tmt", "Pr(>|W|)"],
                  alpha    = s2$corr["alpha", "Estimate"],
                  alpha_se = s2$corr["alpha", "Std.err"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d2 = bind_cols(d, d2)

  # linear regression assuming independence
  m3 = lm(Y ~ tmt, data=dat)
  s3 = summary(m3)
  d3 = data.frame(model    = "linreg",
                  b0       = s3$coefficients["(Intercept)", "Estimate"],
                  se_b0    = s3$coefficients["(Intercept)", "Std. Error"],
                  b1       = s3$coefficients["tmt", "Estimate"],
                  se_b1    = s3$coefficients["tmt", "Std. Error"],
                  var_b1   = vcov(m3)["tmt", "tmt"],
                  sd_resid = sd(s3$residuals),
                  p_b1     = s3$coefficients["tmt", "Pr(>|t|)"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d3 = bind_cols(d, d3)
  result = bind_rows(d1, d2, d3)
  return(result)
}

# Fit models for partially clustered data with a binary outcome to a dataset
# Inputs: a dataset as produced by generate_data(); link function to use (logit or log)
# Output: a dataframe of simulation parameters and model estimates, one row per analysis
# Analyses: geeglm, glm
fit_binary_link = function(dat, my_link = "logit") {
  # logging
  print(dat[1,c("scenario_id", "index")]) 
  print(Sys.time())
  # sim params
  d = data.frame(scenario_id = dat$scenario_id[1],
                 dataset     = dat$index[1],
                 rand_method = dat$rand_method[1],
                 nobs        = nrow(dat),
                 M1          = dat$M1[1],
                 M2          = dat$M2[1],
                 M3          = dat$M3[1],
                 M4          = dat$M4[1],
                 ICC         = dat$ICC[1],
                 truepC      = dat$truepC[1],
                 truepI      = dat$truepI[1],
                 stringsAsFactors = FALSE, row.names = NULL)
  
  # independence gee via geepack::geeglm
  m1 = geeglm(Y ~ tmt, id=cluster_id, family=binomial(link = my_link), corstr="independence", data=dat)
  s1 = summary(m1)
  d1 = data.frame(model  = "gee-ind",
                  link   = my_link,
                  b0     = s1$coefficients["(Intercept)", "Estimate"],
                  se_b0  = s1$coefficients["(Intercept)", "Std.err"],
                  b1     = s1$coefficients["tmt", "Estimate"],
                  exp_b1 = exp(s1$coefficients["tmt", "Estimate"]),
                  se_b1  = s1$coefficients["tmt", "Std.err"],
                  var_b1 = m1$geese$vbeta[2,2],
                  p_b1   = s1$coefficients["tmt", "Pr(>|W|)"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d1 = bind_cols(d, d1)
  
  # exchangeable gee via geepack::geeglm
  m2 = geeglm(Y ~ tmt, id=cluster_id, family=binomial(link = my_link), corstr="exchangeable", data=dat)
  s2 = summary(m2)
  d2 = data.frame(model    = "gee-exch",
                  link     = my_link,
                  b0       = s2$coefficients["(Intercept)", "Estimate"],
                  se_b0    = s2$coefficients["(Intercept)", "Std.err"],
                  b1       = s2$coefficients["tmt", "Estimate"],
                  exp_b1   = exp(s2$coefficients["tmt", "Estimate"]),
                  se_b1    = s2$coefficients["tmt", "Std.err"],
                  var_b1   = m2$geese$vbeta[2,2],
                  p_b1     = s2$coefficients["tmt", "Pr(>|W|)"],
                  alpha    = s2$corr["alpha", "Estimate"],
                  alpha_se = s2$corr["alpha", "Std.err"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d2 = bind_cols(d, d2)
  
  # standard regression assuming independence via glm
  m3 = glm(Y ~ tmt, family=binomial(link = my_link), data=dat)
  s3 = summary(m3)
  d3 = data.frame(model    = "glm",
                  link     = my_link,
                  b0       = s3$coefficients["(Intercept)", "Estimate"],
                  se_b0    = s3$coefficients["(Intercept)", "Std. Error"],
                  b1       = s3$coefficients["tmt", "Estimate"],
                  exp_b1   = exp(s3$coefficients["tmt", "Estimate"]),
                  se_b1    = s3$coefficients["tmt", "Std. Error"],
                  var_b1   = vcov(m3)["tmt", "tmt"],
                  p_b1     = s3$coefficients["tmt", "Pr(>|z|)"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d3 = bind_cols(d, d3)
  result = bind_rows(d1, d2, d3)
  return(result)
}

# Fit models for partially clustered data with a binary outcome to a dataset
# Calls fit_binary_link() with each of the logit and log links
# Inputs: a dataset as produced by generate_data()
# Output: a dataframe of simulation parameters and model estimates, one row per analysis (rsimsum compatible)
fit_binary_models = function(dat) {
  res_logit = fit_binary_link(dat, my_link = "logit")
  res_log   = fit_binary_link(dat, my_link = "log")
  res = bind_rows(res_logit, res_log)
  return(res)
}
  

###### REPORTING RESULTS ######

# Pivot a dataframe of summary results from long (one row per model) to wide
# The resulting wide datasets are passed to the table functions below

convert_wide = function(res){
  res_wide = res %>%
    pivot_wider(names_from = model, values_from = c(deff_exp, deff_obs_med, deff_reldiff_med, pwr_exp, pwr_obs, pwr_diff)) %>% 
    rename(deff_obs_ind      = "deff_obs_med_gee-ind",
           deff_exp_ind      = "deff_exp_gee-ind",
           deff_reldiff_ind  = "deff_reldiff_med_gee-ind",
           deff_obs_exch     = "deff_obs_med_gee-exch",
           deff_exp_exch     = "deff_exp_gee-exch",
           deff_reldiff_exch = "deff_reldiff_med_gee-exch",
           pwr_obs_ind       = "pwr_obs_gee-ind",
           pwr_exp_ind       = "pwr_exp_gee-ind",
           pwr_diff_ind      = "pwr_diff_gee-ind",
           pwr_obs_exch      = "pwr_obs_gee-exch",
           pwr_exp_exch      = "pwr_exp_gee-exch",
           pwr_diff_exch     = "pwr_diff_gee-exch") %>% 
    mutate(numclus = M1 + M2 + M3 + M4,
           gamma1 = M1 / numclus,
           gamma2 = M2 / numclus,
           gamma3 = M3 / numclus,
           gamma4 = M4 / numclus,
           clusdist = paste0("(", gamma1, ", ", gamma2, ", ", gamma3, ", ", gamma4, ")")) %>%
    arrange(rand_method, nobs, scenario_id) %>% 
    relocate(rand_method, nobs, clusdist, ICC,
             deff_obs_ind, deff_exp_ind, deff_reldiff_ind, deff_obs_exch, deff_exp_exch, deff_reldiff_exch,
             pwr_obs_ind, pwr_exp_ind, pwr_diff_ind, pwr_obs_exch, pwr_exp_exch, pwr_diff_exch)
  return(res_wide)
}

# Creates a table of observed and expected design effects

deff_table = function(dat, numsim = numsim, caption = caption) {
  tab = flextable(dat[,1:10]) %>% 
    set_caption(caption) %>%
    set_header_labels(rand_method       = "Randomisation method",
                      nobs              = "Sample size",
                      clusdist          = "Proportion of clusters of size 1-4",
                      deff_obs_ind      = "Observed DEFF*",
                      deff_exp_ind      = "Expected DEFF",
                      deff_reldiff_ind  = "Relative difference (%)*",
                      deff_obs_exch     = "Observed DEFF*",
                      deff_exp_exch     = "Expected DEFF",
                      deff_reldiff_exch = "Relative difference (%)*") %>% 
    add_header_row(colwidths = c(4, 3, 3), values = c("", "GEE independence", "GEE exchangeable")) %>% 
    align(align = "center", part = "header") %>%
    vline(j = c(4, 7)) %>%
    hline(i = 8) %>% 
    line_spacing(space=0.5) %>% 
    add_footer_lines(paste0("*Median value across ", numsim, " simulated datasets")) %>% 
    add_footer_lines("GEE = generalised estimating equation, DEFF = design effect, ICC = intracluster correlation coefficient") %>% 
    fontsize(size=9) %>%
    fontsize(part = "header", size = 10) %>% 
    fontsize(part="footer", size=9) %>% 
    colformat_double(j = c("deff_obs_ind", "deff_exp_ind", "deff_reldiff_ind",
                           "deff_obs_exch", "deff_exp_exch", "deff_reldiff_exch"), digits = 2)
  return(tab)
}

# Creates a table of observed and expected power results

power_table = function(dat, caption = caption){
  tab = flextable(dat[,c(1:4, 11:16)]) %>% 
    set_caption(caption) %>%
    set_header_labels(rand_method   = "Randomisation method",
                      nobs          = "Sample size",
                      clusdist      = "Proportion of clusters of size 1-4",
                      pwr_obs_ind   = "Observed power",
                      pwr_exp_ind   = "Expected power",
                      pwr_diff_ind  = "Difference",
                      pwr_obs_exch  = "Observed power",
                      pwr_exp_exch  = "Expected power",
                      pwr_diff_exch = "Difference") %>% 
    add_header_row(colwidths = c(4, 3, 3), values = c("", "GEE independence", "GEE exchangeable")) %>% 
    align(align = "center", part = "header") %>%
    vline(j = c(4, 7)) %>%
    hline(i = 8) %>% 
    line_spacing(space=0.5) %>% 
    add_footer_lines("GEE = generalised estimating equation, ICC = intracluster correlation coefficient") %>% 
    fontsize(size=9) %>%
    fontsize(part = "header", size = 10) %>% 
    fontsize(part="footer", size=9) %>% 
    colformat_double(j = c("pwr_obs_ind", "pwr_exp_ind", "pwr_diff_ind",
                           "pwr_obs_exch", "pwr_exp_exch", "pwr_diff_exch"), digits = 2)
  return(tab)
}
