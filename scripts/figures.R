## GEE design effects for partially clustered trials
## Kylie Lange (kylie.lange@adelaide.edu.au)

# Figures of the theoretical design effect functions

library(tidyverse)
library(patchwork)

## Functions

cont_deffs = function(params) {
  rho = params$rho
  gamma_k = list(params$gamma1, params$gamma2, params$gamma3, params$gamma4)
  
  ind_clus_sum1 = gamma_k %>% imap(~ .y*.x)  # k*gamma_k
  deff_ind_clus = 1 + rho*(sum(unlist(ind_clus_sum1)) - 1)
  
  exch_clus_sum1 = gamma_k %>% imap(~ (1/(1 + (.y-1)*rho))*.x)  # 1/(1+(k-1)rho)*gamma_k
  deff_exch_clus = 1 / sum(unlist(exch_clus_sum1))
  
  deff_ind_ind = 1
  
  exch_ind_sum1 = gamma_k %>% imap(~ (1 + (.y-2)*rho)/((1 - rho)*(1 + (.y-1)*rho))*.x)  # (1+(k-2)rho)/((1-rho)*(1+(k-1)rho))*gamma_k
  deff_exch_ind = 1 / sum(unlist(exch_ind_sum1))
  
  res = data.frame(params, deff_ind_clus, deff_exch_clus, deff_ind_ind, deff_exch_ind)
  return(res)
}

logit_deffs = function(params) {
  rho = params$rho
  gamma_k = list(params$gamma1, params$gamma2, params$gamma3, params$gamma4)
  pi_I = params$pi_I
  pi_C = params$pi_C
  
  ind_clus_sum1 = gamma_k %>% imap(~ (.y-1)*.x)   # (k-1)*gamma_k
  deff_ind_clus = 1 + rho*(sum(unlist(ind_clus_sum1)))
  
  exch_clus_sum1 = gamma_k %>% imap(~ 1/(1 + (.y-1)*rho)*.x)   # 1/(1+(k-1)rho)*gamma_k
  deff_exch_clus = 1 / sum(unlist(exch_clus_sum1))
  
  ind_ind_sum1 = gamma_k %>% imap(~ (.y-1)*.x)   # (k-1)*gamma_k
  deff_ind_ind = 1 + rho*((0.5 - sqrt(pi_I*pi_C*(1-pi_I)*(1-pi_C))/(pi_I*(1-pi_I)+pi_C*(1-pi_C)))*sum(unlist(ind_ind_sum1)))
  
  exch_ind_sum1 = gamma_k %>% imap(~ 1/(1 + (.y-1)*rho)*.x)       # 1/(1+(k-1)rho)*gamma_k
  exch_ind_sum2 = gamma_k %>% imap(~ (.y-1)/(1 + (.y-1)*rho)*.x)  # (k-1)/(1+(k-1)rho)*gamma_k
  deff_exch_ind = (sum(unlist(exch_ind_sum1)) + (0.5 - sqrt(pi_I*pi_C*(1-pi_I)*(1-pi_C))/(pi_I*(1-pi_I)+pi_C*(1-pi_C)))*rho/(1-rho)*sum(unlist(exch_ind_sum2))) / 
    (sum(unlist(exch_ind_sum1))**2 + rho/(1-rho)*sum(unlist(exch_ind_sum1))*sum(unlist(exch_ind_sum2)))
  
  res = data.frame(params, deff_ind_clus, deff_exch_clus, deff_ind_ind, deff_exch_ind)
  return(res)
}

log_deffs = function(params) {
  rho = params$rho
  gamma_k = list(params$gamma1, params$gamma2, params$gamma3, params$gamma4)
  pi_I = params$pi_I
  pi_C = params$pi_C
  
  ind_clus_sum1 = gamma_k %>% imap(~ (.y-1)*.x)   # (k-1)*gamma_k
  deff_ind_clus =1 + rho*(sum(unlist(ind_clus_sum1)))
  
  exch_clus_sum1 = gamma_k %>% imap(~ 1/(1 + (.y-1)*rho)*.x)   # 1/(1+(k-1)rho)*gamma_k
  deff_exch_clus = 1 / sum(unlist(exch_clus_sum1))
  
  ind_ind_sum1 = gamma_k %>% imap(~ (.y-1)*.x)   # (k-1)*gamma_k
  deff_ind_ind = 1 + rho*((0.5 - sqrt(pi_I*pi_C*(1-pi_I)*(1-pi_C))/(pi_I*(1-pi_C)+pi_C*(1-pi_I)))*sum(unlist(ind_ind_sum1)))
  
  exch_ind_sum1 = gamma_k %>% imap(~ 1/(1 + (.y-1)*rho)*.x)       # 1/(1+(k-1))*gamma_k
  exch_ind_sum2 = gamma_k %>% imap(~ (.y-1)/(1 + (.y-1)*rho)*.x)  # (k-1)/(1+(k-1)rho)*gamma_k
  deff_exch_ind = (sum(unlist(exch_ind_sum1)) + (0.5 - sqrt(pi_I*pi_C*(1-pi_I)*(1-pi_C))/(pi_I*(1-pi_C)+pi_C*(1-pi_I)))*rho/(1-rho)*sum(unlist(exch_ind_sum2))) / 
    (sum(unlist(exch_ind_sum1))**2 + rho/(1-rho)*sum(unlist(exch_ind_sum1))*sum(unlist(exch_ind_sum2)))
  
  res = data.frame(params, deff_ind_clus, deff_exch_clus, deff_ind_ind, deff_exch_ind)
  return(res)
}

prep_data = function(dat, link="linear"){
  dat = dat %>% pmap(data.frame) #convert to a list so can pass the parameters to map()
  
  if (link == "linear") {
    deffs = map(dat, cont_deffs)
  }
  else if (link == "logit") {
    deffs = map(dat, logit_deffs)
  }
  else if (link == "log") {
    deffs = map(dat, log_deffs)
  }
  
  deffs = deffs %>%
    reduce(., rbind.data.frame) %>% 
    pivot_longer(cols = starts_with("deff"), names_to = "model", names_prefix = "deff_", values_to = "deff")
  deffs$model = factor(deffs$model, levels = c("ind_clus", "ind_ind", "exch_clus", "exch_ind"))
  return(deffs)
}

plot_de = function(data, x, y, group, xlab=NULL, ylab=NULL, title=NULL, subtitle=NULL, caption=NULL, ylim=c(NA,NA)) {
  gg = ggplot(data = data, aes(x = x, y = y, group = group, colour = group)) +
    geom_line(aes(linetype = group), linewidth = 0.8) +
    scale_linetype_manual(name = "",
                          labels = c("GEE-ind, Cluster rand", "GEE-ind, Individual rand",
                                     "GEE-exch, Cluster rand", "GEE-exch, Individual rand"),
                          values=c("dashed", "solid", "dashed", "solid")) +
    scale_color_manual(name = "",
                       labels = c("GEE-ind, Cluster rand", "GEE-ind, Individual rand",
                                  "GEE-exch, Cluster rand", "GEE-exch, Individual rand"),
                       values=c("#E69F00", "#E69F00", "#56B4E9", "#56B4E9")) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle, caption = caption) +
    scale_y_continuous(limits = c(ylim[1], ylim[2])) +
    theme_bw() +
    theme(legend.key.size = unit(1.5, "cm"), legend.text = element_text(size = 9),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))
  return(gg)
}

## Figure 1: Continuous outcome - varying ICC

# Panel A
rho = seq(0, 0.995, by=0.0002)
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
data1A = data.frame(rho, gamma1, gamma2, gamma3, gamma4)
data1A_deffs = prep_data(data1A, link = "linear")
plot1A = plot_de(data1A_deffs, x=data1A_deffs$rho, y=data1A_deffs$deff, group=data1A_deffs$model, ylim=c(0, 2.5),
                 xlab="ICC", ylab="Design effect", subtitle="A", caption="")

# Panel B
rho = seq(0, 0.995, by=0.0002)
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
data1B = data.frame(rho, gamma1, gamma2, gamma3, gamma4)
data1B_deffs = prep_data(data1B, link = "linear")
plot1B = plot_de(data1B_deffs, x=data1B_deffs$rho, y=data1B_deffs$deff, group=data1B_deffs$model, ylim=c(0, 2.5),
                 xlab="ICC", ylab="Design effect", subtitle="B", caption="")

# title = "Figure 1: Design effects for continuous outcomes by randomisation method and GEE working correlation structure"
fig1 = plot1A + plot1B + plot_layout(guides = "collect") +
  plot_annotation(title = "") & theme(legend.position="bottom")
ggsave("output/continuous.tiff",height=8.5,width=9.75,unit="in",dpi=800)

## Supp Figure 1: Logit link - varying ICC

# Panel A
rho = seq(0, 0.995, by=0.0002)
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
pi_C = 0.4
pi_I = 0.3
dataS1A = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I)
dataS1A_deffs = prep_data(dataS1A, link = "logit")
plotS1A = plot_de(dataS1A_deffs, x=dataS1A_deffs$rho, y=dataS1A_deffs$deff, group=dataS1A_deffs$model, ylim=c(0, 2.5),
                 xlab="ICC", ylab="Design effect", subtitle="A", caption="")

# Panel B
rho = seq(0, 0.995, by=0.0002)
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
pi_C = 0.4
pi_I = 0.3
dataS1B = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I)
dataS1B_deffs = prep_data(dataS1B, link = "logit")
plotS1B = plot_de(dataS1B_deffs, x=dataS1B_deffs$rho, y=dataS1B_deffs$deff, group=dataS1B_deffs$model, ylim=c(0, 2.5),
                 xlab="ICC", ylab="Design effect", subtitle="B", caption="")

# title = "Supp Figure 1: Design effects for binary outcomes with a logit link by randomisation method and GEE working correlation structure"
figS1 = plotS1A + plotS1B + plot_layout(guides = "collect") +
  plot_annotation(title = "") & theme(legend.position="bottom")
ggsave("output/logit-ICC.tiff",height=8.5,width=9.75,unit="in",dpi=800)

## Figure 2: Logit link - varying OR when pC=0.4

# Panel A
rho = 0.2
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
pi_C = 0.4
pi_I = seq(0.005, 0.995, by=0.0002)
or = (pi_I / (1 - pi_I)) / (pi_C / (1 - pi_C))
data2A = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, or)
data2A_deffs = prep_data(data2A, link = "logit")
data2A_plot = data2A_deffs %>% filter(or <= 5)
plot2A = plot_de(data2A_plot, x=data2A_plot$or, y=data2A_plot$deff, group=data2A_plot$model, ylim=c(0.3, 2.25),
                 xlab="Odds ratio", ylab="Design effect", subtitle="A", caption="")

# Panel B
rho = 0.2
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
pi_C = 0.4
pi_I = seq(0.005, 0.995, by=0.0002)
or = (pi_I / (1 - pi_I)) / (pi_C / (1 - pi_C))
data2B = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, or)
data2B_deffs = prep_data(data2B, link = "logit")
data2B_plot = data2B_deffs %>% filter(or <= 5)
plot2B = plot_de(data2B_plot, x=data2B_plot$or, y=data2B_plot$deff, group=data2B_plot$model, ylim=c(0.3, 2.25),
                 xlab="Odds ratio", ylab="Design effect", subtitle="B", caption="")

# Panel C
rho = 0.8
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
pi_C = 0.4
pi_I = seq(0.005, 0.995, by=0.0002)
or = (pi_I / (1 - pi_I)) / (pi_C / (1 - pi_C))
data2C = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, or)
data2C_deffs = prep_data(data2C, link = "logit")
data2C_plot = data2C_deffs %>% filter(or <= 5)
plot2C = plot_de(data2C_plot, x=data2C_plot$or, y=data2C_plot$deff, group=data2C_plot$model, ylim=c(0.3, 2.25),
                 xlab="Odds ratio", ylab="Design effect", subtitle="C", caption="")

# Panel D
rho = 0.8
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
pi_C = 0.4
pi_I = seq(0.005, 0.995, by=0.0002)
or = (pi_I / (1 - pi_I)) / (pi_C / (1 - pi_C))
data2D = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, or)
data2D_deffs = prep_data(data2D, link = "logit")
data2D_plot = data2D_deffs %>% filter(or <= 5)
plot2D = plot_de(data2D_plot, x=data2D_plot$or, y=data2D_plot$deff, group=data2D_plot$model, ylim=c(0.3, 2.25),
                 xlab="Odds ratio", ylab="Design effect", subtitle="D", caption="")

# title = "Figure 2: Design effects for binary outcomes with a logit link by randomisation method and GEE working correlation structure"
fig2 = plot2A + plot2B + plot2C + plot2D + plot_layout(guides = "collect") +
  plot_annotation(title = "") & theme(legend.position="bottom")
ggsave("output/logit-OR04.tiff",height=8.5,width=9.75,unit="in",dpi=800)

## Supp Figure 2: Logit link - varying OR when pC=0.1

# Panel A
rho = 0.2
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
pi_C = 0.1
pi_I = seq(0.005, 0.995, by=0.0002)
or = (pi_I / (1 - pi_I)) / (pi_C / (1 - pi_C))
dataS2A = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, or)
dataS2A_deffs = prep_data(dataS2A, link = "logit")
dataS2A_plot = dataS2A_deffs %>% filter(or <= 5)
plotS2A = plot_de(dataS2A_plot, x=dataS2A_plot$or, y=dataS2A_plot$deff, group=dataS2A_plot$model, ylim=c(0.3, 2.25),
                  xlab="Odds ratio", ylab="Design effect", subtitle="A", caption="")

# Panel B
rho = 0.2
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
pi_C = 0.1
pi_I = seq(0.005, 0.995, by=0.0002)
or = (pi_I / (1 - pi_I)) / (pi_C / (1 - pi_C))
dataS2B = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, or)
dataS2B_deffs = prep_data(dataS2B, link = "logit")
dataS2B_plot = dataS2B_deffs %>% filter(or <= 5)
plotS2B = plot_de(dataS2B_plot, x=dataS2B_plot$or, y=dataS2B_plot$deff, group=dataS2B_plot$model, ylim=c(0.3, 2.25),
                  xlab="Odds ratio", ylab="Design effect", subtitle="B", caption="")

# Panel C
rho = 0.8
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
pi_C = 0.1
pi_I = seq(0.005, 0.995, by=0.0002)
or = (pi_I / (1 - pi_I)) / (pi_C / (1 - pi_C))
dataS2C = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, or)
dataS2C_deffs = prep_data(dataS2C, link = "logit")
dataS2C_plot = dataS2C_deffs %>% filter(or <= 5)
plotS2C = plot_de(dataS2C_plot, x=dataS2C_plot$or, y=dataS2C_plot$deff, group=dataS2C_plot$model, ylim=c(0.3, 2.25),
                  xlab="Odds ratio", ylab="Design effect", subtitle="C", caption="")

# Panel D
rho = 0.8
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
pi_C = 0.1
pi_I = seq(0.005, 0.995, by=0.0002)
or = (pi_I / (1 - pi_I)) / (pi_C / (1 - pi_C))
dataS2D = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, or)
dataS2D_deffs = prep_data(dataS2D, link = "logit")
dataS2D_plot = dataS2D_deffs %>% filter(or <= 5)
plotS2D = plot_de(dataS2D_plot, x=dataS2D_plot$or, y=dataS2D_plot$deff, group=dataS2D_plot$model, ylim=c(0.3, 2.25),
                  xlab="Odds ratio", ylab="Design effect", subtitle="D", caption="")

# title = "Supp Figure 2: Design effects for binary outcomes with a logit link by randomisation method and GEE working correlation structure"
figS2 = plotS2A + plotS2B + plotS2C + plotS2D + plot_layout(guides = "collect") +
  plot_annotation(title = "") & theme(legend.position="bottom")
ggsave("output/logit-OR01.tiff",height=8.5,width=9.75,unit="in",dpi=800)

## Supp Figure 3: Log link - varying ICC

# Panel A
rho = seq(0, 0.995, by=0.0002)
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
pi_C = 0.4
pi_I = 0.3
dataS3A = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I)
dataS3A_deffs = prep_data(dataS3A, link = "log")
plotS3A = plot_de(dataS3A_deffs, x=dataS3A_deffs$rho, y=dataS3A_deffs$deff, group=dataS3A_deffs$model, ylim=c(0, 2.5),
                 xlab="ICC", ylab="Design effect", subtitle="A", caption="")

# Panel B
rho = seq(0, 0.995, by=0.0002)
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
pi_C = 0.4
pi_I = 0.3
dataS3B = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I)
dataS3B_deffs = prep_data(dataS3B, link = "log")
plotS3B = plot_de(dataS3B_deffs, x=dataS3B_deffs$rho, y=dataS3B_deffs$deff, group=dataS3B_deffs$model, ylim=c(0, 2.5),
                 xlab="ICC", ylab="Design effect", subtitle="B", caption="")

# title = "Supp Figure 3: Design effects for binary outcomes with a log link by randomisation method and GEE working correlation structure"
figS3 = plotS3A + plotS3B + plot_layout(guides = "collect") +
  plot_annotation(title = "") & theme(legend.position="bottom")
ggsave("output/log-ICC.tiff",height=8.5,width=9.75,unit="in",dpi=800)

## Supp Figure 4: Log link - varying RR when pC=0.4

# Panel A
rho = 0.2
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
pi_C = 0.4
pi_I = seq(0.005, 0.995, by=0.0002)
rr = pi_I / pi_C 
dataS4A = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, rr)
dataS4A_deffs = prep_data(dataS4A, link = "log")
dataS4A_plot = dataS4A_deffs %>% filter(rr <= 2)
plotS4A = plot_de(dataS4A_deffs, x=dataS4A_deffs$rr, y=dataS4A_deffs$deff, group=dataS4A_deffs$model, ylim=c(0.3, 2.25),
                 xlab="Relative risk", ylab="Design effect", subtitle="A", caption="")

# Panel B
rho = 0.2
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
pi_C = 0.4
pi_I = seq(0.005, 0.995, by=0.0002)
rr = pi_I / pi_C 
dataS4B = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, rr)
dataS4B_deffs = prep_data(dataS4B, link = "log")
dataS4B_plot = dataS4B_deffs %>% filter(rr <= 2)
plotS4B = plot_de(dataS4B_plot, x=dataS4B_plot$rr, y=dataS4B_plot$deff, group=dataS4B_plot$model, ylim=c(0.3, 2.25),
                 xlab="Relative risk", ylab="Design effect", subtitle="B", caption="")

# Panel C
rho = 0.8
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
pi_C = 0.4
pi_I = seq(0.005, 0.995, by=0.0002)
rr = pi_I / pi_C 
dataS4C = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, rr)
dataS4C_deffs = prep_data(dataS4C, link = "log")
dataS4C_plot = dataS4C_deffs %>% filter(rr <= 2)
plotS4C = plot_de(dataS4C_plot, x=dataS4C_plot$rr, y=dataS4C_plot$deff, group=dataS4C_plot$model, ylim=c(0.3, 2.25),
                 xlab="Relative risk", ylab="Design effect", subtitle="C", caption="")

# Panel D
rho = 0.8
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
pi_C = 0.4
pi_I = seq(0.005, 0.995, by=0.0002)
rr = pi_I / pi_C 
dataS4D = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, rr)
dataS4D_deffs = prep_data(dataS4D, link = "log")
dataS4D_plot = dataS4D_deffs %>% filter(rr <= 2)
plotS4D = plot_de(dataS4D_plot, x=dataS4D_plot$rr, y=dataS4D_plot$deff, group=dataS4D_plot$model, ylim=c(0.3, 2.25),
                 xlab="Relative risk", ylab="Design effect", subtitle="D", caption="")

# title = "Supp Figure 4: Design effects for binary outcomes with a log link by randomisation method and GEE working correlation structure"
figS4 = plotS4A + plotS4B + plotS4C + plotS4D + plot_layout(guides = "collect") +
  plot_annotation(title = "") &
  theme(legend.position="bottom")
ggsave("output/log-RR04.tiff",height=8.5,width=9.75,unit="in",dpi=800)

## Supp Figure 5: Log link - varying RR when pC=0.1

# Panel A
rho = 0.2
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
pi_C = 0.1
pi_I = seq(0.005, 0.995, by=0.0002)
rr = pi_I / pi_C 
dataS5A = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, rr)
dataS5A_deffs = prep_data(dataS5A, link = "log")
dataS5A_plot = dataS5A_deffs %>% filter(rr <= 2)
plotS5A = plot_de(dataS5A_plot, x=dataS5A_plot$rr, y=dataS5A_plot$deff, group=dataS5A_plot$model, ylim=c(0.3, 2.25),
                 xlab="Relative risk", ylab="Design effect", subtitle="A", caption="")

# Panel B
rho = 0.2
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
pi_C = 0.1
pi_I = seq(0.005, 0.995, by=0.0002)
rr = pi_I / pi_C 
dataS5B = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, rr)
dataS5B_deffs = prep_data(dataS5B, link = "log")
dataS5B_plot = dataS5B_deffs %>% filter(rr <= 2)
plotS5B = plot_de(dataS5B_plot, x=dataS5B_plot$rr, y=dataS5B_plot$deff, group=dataS5B_plot$model, ylim=c(0.3, 2.25),
                 xlab="Relative risk", ylab="Design effect", subtitle="B", caption="")

# Panel C
rho = 0.8
gamma1 = 0.70
gamma2 = 0.15
gamma3 = 0.10
gamma4 = 0.05
pi_C = 0.1
pi_I = seq(0.005, 0.995, by=0.0002)
rr = pi_I / pi_C 
dataS5C = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, rr)
dataS5C_deffs = prep_data(dataS5C, link = "log")
dataS5C_plot = dataS5C_deffs %>% filter(rr <= 2)
plotS5C = plot_de(dataS5C_plot, x=dataS5C_plot$rr, y=dataS5C_plot$deff, group=dataS5C_plot$model, ylim=c(0.3, 2.25),
                 xlab="Relative risk", ylab="Design effect", subtitle="C", caption="")

# Panel D
rho = 0.8
gamma1 = 0.25
gamma2 = 0.25
gamma3 = 0.25
gamma4 = 0.25
pi_C = 0.1
pi_I = seq(0.005, 0.995, by=0.0002)
rr = pi_I / pi_C 
dataS5D = data.frame(rho, gamma1, gamma2, gamma3, gamma4, pi_C, pi_I, rr)
dataS5D_deffs = prep_data(dataS5D, link = "log")
dataS5D_plot = dataS5D_deffs %>% filter(rr <= 2)
plotS5D = plot_de(dataS5D_plot, x=dataS5D_plot$rr, y=dataS5D_plot$deff, group=dataS5D_plot$model, ylim=c(0.3, 2.25),
                 xlab="Relative risk", ylab="Design effect", subtitle="D", caption="")

# title = "Supp Figure 5: Design effects for binary outcomes with a log link by randomisation method and GEE working correlation structure"
figS5 = plotS5A + plotS5B + plotS5C + plotS5D + plot_layout(guides = "collect") +
  plot_annotation(title = "") & theme(legend.position="bottom")
ggsave("output/log-RR01.tiff",height=8.5,width=9.75,unit="in",dpi=800)
