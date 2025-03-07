## GEE design effects for partially clustered trials
## Kylie Lange (kylie.lange@adelaide.edu.au)

# Validation of theoretically derived design effects

## CONTINOUS OUTCOME

# Simulation of partially clustered data
# Outcome: continuous (linear link)
# Treatment group allocation (2 groups): cluster randomisation, individual randomisation

# Sensitivity analysis where the total number of cluster is fixed,
# but the number of clusters of each size are randomly generated from multinomial distribution.

library(tidyverse)
library(readxl)
library(pwrss)
library(lme4)
library(geepack)
library(rsimsum)
library(rmarkdown)
library(flextable)
library(officer)

source("R/partial_clus_functions.R")
source("R/partial_clus_deffs.R")

# simulation parameters
parameters = readxl::read_xlsx("R/scenarios_cont_multinom.xlsx")
# calculate the expected number of clusters of each size
parameters = parameters %>% 
  mutate(num_clus1 = ceiling(num_clus * delta_1),
         num_clus2 = ceiling(num_clus * delta_2),
         num_clus3 = ceiling(num_clus * delta_3),
         num_clus4 = ceiling(num_clus * delta_4))
param_list = parameters %>% pmap(data.frame)  #convert to a list so can pass the parameters to map() later on


# calculate the independence gee expected deff for each scenario
deff_exp_indgee  = map(param_list, calculate_cont_indgee_deff) %>% 
  reduce(., rbind.data.frame) %>% 
  setNames(c("deff_exp", "N_eff", "power_exp")) %>% 
  mutate(scenario_id = seq(1:nrow(.)),
         model = "gee-ind")
# calculate the exchangeable gee expected deff for each scenario
deff_exp_exchgee = map(param_list, calculate_cont_exchgee_deff) %>% 
  reduce(., rbind.data.frame) %>% 
  setNames(c("deff_exp", "N_eff", "power_exp")) %>% 
  mutate(scenario_id = seq(1:nrow(.)),
         model = "gee-exch")
deff_exp = bind_rows(deff_exp_indgee, deff_exp_exchgee)
deff_exp = deff_exp %>%
  left_join(parameters, "scenario_id") %>% 
  select(c(1,2,3,4,5,6))
deff_plot = ggplot(deff_exp, aes(x=model, y=deff_exp, fill=rand_method)) +
  geom_dotplot(binaxis='y', stackdir='center')

# random seed
s = format(generate_seed(), scientific=FALSE)
set.seed(s)
saveRDS(s, file="simdata/cont-multinom/seed.rds")

# generate the datasets
sims=10000
datalist = map(param_list, simulate_randomclusters, nsims=sims, outcome="cont")
saveRDS(datalist, file="simdata/cont-multinom/datalist.rds")

# analyse each dataset
model_results = modify_depth(datalist, .depth=2, ~fit_cont_models(.x))
model_results = bind_rows(map(model_results, bind_rows))

# use rsimsum::dropbig() to identify non-converged models
model_results = dropbig(data = model_results, estvarname = "b1", se = "se_b1", by = "scenario_id")

# add the expected deffs to the results
model_results = model_results %>%
  left_join(deff_exp, c("scenario_id", "model"), suffix=c("", ".y")) %>% 
  select(-rand_method.y)
saveRDS(model_results, file="simdata/cont-multinom/model_results.rds")

# only report results from models that converged
model_results_conv = model_results %>% filter(.dropbig == FALSE)

# calculate the observed design effects
summ = model_results_conv  %>% 
  group_by(scenario_id, dataset) %>% 
  mutate(indp_var = var_b1[model == "linreg"]) %>% 
  ungroup %>%
  filter(model != "linreg") %>% 
  mutate(deff_obs      = var_b1 / indp_var,
         deff_reldiff  = (deff_obs-deff_exp)/deff_exp*100)

# summaries for each scenario

model_res = model_results %>% 
  filter(model == "gee-ind" | model == "gee-exch") %>% 
  group_by(scenario_id, model) %>%
  summarise(conv    = sum(.dropbig)/sims*100,
            npd     = sum(abs(alpha) > 1)/sims*100,
            estICC  = median(alpha),
            .groups = "drop")

res = summ %>% 
  group_by(scenario_id, model) %>% 
  summarise(rand_method       = first(rand_method),
            nclus             = first(M1) + first(M2) + first(M3) + first(M4),
            ICC               = first(ICC),
            deff_exp          = first(deff_exp),
            deff_obs_med      = median(deff_obs),
            deff_reldiff_med  = median(deff_reldiff),
            pwr_exp           = first(power_exp)*100,
            pwr_obs           = sum(p_b1 < 0.05)/n()*100,
            pwr_diff          = pwr_obs - pwr_exp,
            .groups           = "drop")

# generate tables

model_res_wide = model_res %>%
  pivot_wider(names_from = model, values_from = c(conv, npd, estICC)) %>% 
  rename(conv_ind  = "conv_gee-ind",
         conv_exch = "conv_gee-exch",
         npd_exch  = "npd_gee-exch",
         estICC    = "estICC_gee-exch") %>% 
  select(-c("npd_gee-ind", "estICC_gee-ind")) %>% 
  left_join(parameters, by="scenario_id")

model_table = flextable(model_res_wide[,1:6]) %>% 
  colformat_double(j = c("estICC"), digits = 2) %>% 
  fontsize(size=9)

res_wide = convert_wide_sens(res)
tableS10 = deff_table(res_wide, numsim = sims,
                    caption = "Supplementary Table 10: Observed and expected design effects for a continuous outcome,
                    with a varying number of clusters of each size")
tableS11 = power_table(res_wide, caption = "Supplementary Table 2: Observed and expected power for a continuous outcome,
                       with a varying number of clusters of each size")

# write output file of results

read_docx() %>%
  body_add_gg(deff_plot) %>%
  body_add_break() %>% 
  body_add_flextable(tableS10) %>%
  body_add_break() %>% 
  body_add_flextable(tableS11) %>%
  body_add_break() %>% 
  body_add_flextable(model_table) %>%
  print(target="output/deff_results_cont_multinom.docx")
