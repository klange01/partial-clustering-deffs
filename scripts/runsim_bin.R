## Validation of GEE design effects for partially clustered trials
## Kylie Lange (kylie.lange@adelaide.edu.au)

# Validation of theoretically derived design effects

## BINARY OUTCOME

# Simulation of partially clustered data
# Outcome: binary (logit link), binary (log link)
# Treatment group allocation (2 groups): cluster randomisation, individual randomisation

rm(list = ls())

library(tidyverse)
library(readxl)
library(pwrss)
library(lme4)
library(geepack)
library(rsimsum)
library(CorBin)
library(rmarkdown)
library(flextable)
library(officer)

source("R/partial_clus_functions.R")
source("R/partial_clus_deffs.R")

# simulation parameters
parameters = readxl::read_xlsx("R/scenarios_bin.xlsx")
param_list = parameters %>% pmap(data.frame)  #convert to a list so can pass the parameters to map() later on

# calculate the independence gee expected deff for each scenario
# logit link
deff_exp_indgee_logit  = map(param_list, calculate_logit_indgee_deff) %>% 
  reduce(., rbind.data.frame) %>% 
  setNames(c("deff_exp", "N_eff", "power_exp")) %>% 
  mutate(scenario_id = seq(1:nrow(.)),
         model = "gee-ind",
         link  = "logit")
# log link
deff_exp_indgee_log  = map(param_list, calculate_log_indgee_deff) %>% 
  reduce(., rbind.data.frame) %>% 
  setNames(c("deff_exp", "N_eff", "power_exp")) %>% 
  mutate(scenario_id = seq(1:nrow(.)),
         model = "gee-ind",
         link  = "log")
# calculate the exchangeable gee expected deff for each scenario
# logit link
deff_exp_exchgee_logit = map(param_list, calculate_logit_exchgee_deff) %>% 
  reduce(., rbind.data.frame) %>%  
  setNames(c("deff_exp", "N_eff", "power_exp")) %>% 
  mutate(scenario_id = seq(1:nrow(.)),
         model = "gee-exch",
         link  = "logit")
# log link
deff_exp_exchgee_log = map(param_list, calculate_log_exchgee_deff) %>% 
  reduce(., rbind.data.frame) %>% 
  setNames(c("deff_exp", "N_eff", "power_exp")) %>% 
  mutate(scenario_id = seq(1:nrow(.)),
         model = "gee-exch",
         link  = "log")

deff_exp_logit = bind_rows(deff_exp_indgee_logit, deff_exp_exchgee_logit)
deff_exp_logit = deff_exp_logit %>%
  left_join(parameters, "scenario_id") %>% 
  select(c(1,2,3,4,5,6,7))
deff_plot_logit = ggplot(deff_exp_logit, aes(x=model, y=deff_exp, fill=rand_method)) +
  geom_dotplot(binaxis='y', stackdir='center') +
  labs(title = 'Design effects for a binary outcome with a logit link')
deff_exp_log = bind_rows(deff_exp_indgee_log, deff_exp_exchgee_log)
deff_exp_log = deff_exp_log %>%
  left_join(parameters, "scenario_id") %>% 
  select(c(1,2,3,4,5,6,7))
deff_plot_log = ggplot(deff_exp_log, aes(x=model, y=deff_exp, fill=rand_method)) +
  geom_dotplot(binaxis='y', stackdir='center') +
  labs(title = 'Design effects for a binary outcome with a log link')
deff_exp = bind_rows(deff_exp_logit, deff_exp_log)

# random seed
s = format(generate_seed(), scientific=FALSE)
set.seed(s)
saveRDS(s, file="simdata/binary/seed.rds")

# generate the datasets
sims=10000
datalist = map(param_list, simulate, nsims=sims, outcome="binary")
saveRDS(datalist, file="simdata/binary/datalist.rds")

# analyse each dataset
model_results = modify_depth(datalist, .depth=2, ~fit_binary_models(.x))
model_results = bind_rows(map(model_results, bind_rows))

# use rsimsum::dropbig() to identify non-converged models
model_results = dropbig(data = model_results, estvarname = "b1", se = "se_b1", by = "scenario_id")

# add the expected deffs to the results
model_results = model_results %>%
  left_join(deff_exp, c("scenario_id", "model", "link")) %>% 
  select(-"rand_method.y") %>% 
  rename(rand_method = rand_method.x)
saveRDS(model_results, file="simdata/binary/model_results.rds")

# only report results from models that converged
model_results_conv = model_results %>% filter(.dropbig == FALSE)

# calculate the observed design effects
summ = model_results_conv  %>% 
  group_by(scenario_id, dataset, link) %>% 
  mutate(indp_var = var_b1[model == "glm"]) %>% 
  ungroup %>%
  filter(model != "glm") %>% 
  mutate(deff_obs = var_b1 / indp_var,
         deff_reldiff  = (deff_obs-deff_exp)/deff_exp*100)

# summaries for each scenario

model_res = model_results %>% 
  filter(model == "gee-ind" | model == "gee-exch") %>% 
  group_by(scenario_id, model, link) %>%
  summarise(conv    = sum(.dropbig)/sims*100,
            npd     = sum(abs(alpha) > 1)/sims*100,
            estICC  = median(alpha),
            .groups = "drop")

res = summ %>% 
  group_by(scenario_id, model, link) %>% 
  summarise(rand_method       = first(rand_method),
            nobs              = first(nobs),
            ICC               = first(ICC),
            M1                = first(M1),
            M2                = first(M2),
            M3                = first(M3),
            M4                = first(M4),
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
  left_join(parameters, by="scenario_id") %>% 
  arrange(desc(link))

model_table = flextable(model_res_wide[,1:6]) %>% 
  colformat_double(j = c("estICC"), digits = 2) %>% 
  fontsize(size=9)

res_logit = res %>% filter(link == "logit")
res_wide_logit = convert_wide(res_logit)
table5 = deff_table(res_wide_logit, numsim = sims,
                    caption = "Table 5: Observed and expected design effects for a binary outcome with a logit link")
table6 = power_table(res_wide_logit, caption = "Table 6: Observed and expected power for a binary outcome with a logit link")

res_log = res %>% filter(link == "log")
res_wide_log = convert_wide(res_log)
tableS2 = deff_table(res_wide_log, numsim = sims,
                    caption = "Supp Table 2: Observed and expected design effects for a binary outcome with a log link")
tableS3 = power_table(res_wide_log, caption = "Supp Table 3: Observed and expected power for a binary outcome with a log link")

# write output file of results
read_docx() %>%
  body_add_gg(deff_plot_logit) %>%
  body_add_gg(deff_plot_log) %>% 
  body_add_break() %>% 
  body_add_flextable(table5) %>%
  body_add_break() %>% 
  body_add_flextable(table6) %>%
  body_add_break() %>% 
  body_add_flextable(tableS2) %>%
  body_add_break() %>% 
  body_add_flextable(tableS3) %>%
  body_add_break() %>% 
  body_add_flextable(model_table) %>%
  print(target="output/deff_results_bin.docx")
