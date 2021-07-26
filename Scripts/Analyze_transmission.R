#### Preamble ####
# Analyze_transmission.R
#
# ---
### title: Script B for Objective 2: effects of seasonal network differences on pathogen transmission
# author: Marie Gilbertson
# date: "08/30/2020"
#---
# 
# What this code does:
# 1. Analyzes differences between wet and dry season epidemics

#### Clear Environment ####
remove(list=ls())

##### Set seed #####
set.seed(3268)

# Load libraries
library(ggplot2)
library(dplyr)
library(stringr) # for "str_split"
library(viridis)
library(ggpubr)



#### load data ####
# load parameter sets
param.sets <- get(load("Data/simulation_paramsets.Rdata"))

# set which parameter sets to analyze
param.ids <- c(1:420)

# empty dataframe to store results
all.results <- NULL
summ.stats <- data.frame(matrix(nrow = length(param.ids), ncol = 7))
colnames(summ.stats) <- c("param.set", "m.dur.time", "m.total.i", "m.final.i", 
                          "m.max.i", "m.max.it_s", "prop.fail")

# loop through and load all results
for(i in 1:length(param.ids)){
  print(i)
  param.set.num <- param.ids[i]
  
  full.name <- paste("Output/Simulation_Results/full set results_paramset ", param.set.num, ".Rdata", sep = "")
  set_results <- get(load(full.name))
  
  all.results <- rbind(all.results, set_results)
  
  #### store summary stats ####
  summ.stats$param.set[i] <- set_results$param.set[1]
  
  # proportion failed epidemics; if only 1 individual is ever infected, consider it a "failed" epidemic
  bad.draws <- sum(set_results$num.failed)
  only.1 <- nrow(set_results[set_results$total.i<=(1/33),])
  # number of bad draws + number only infecting 1 / simulations + number of bad draws
  summ.stats$prop.fail[i] <- (bad.draws + only.1)/(nrow(set_results) + bad.draws)
    
  ## now, of only the "successful" epidemics:
  set_results_succ <- set_results[set_results$total.i>(1/33),]
  # mean duration
  summ.stats$m.dur.time[i] <- mean(set_results_succ$dur.time)
  # mean total ever infected
  summ.stats$m.total.i[i] <- mean(set_results_succ$total.i)
  # mean final infected
  summ.stats$m.final.i[i] <- mean(set_results_succ$final.i)
  # mean max infected
  summ.stats$m.max.i[i] <- mean(set_results_succ$max.i)
  # mean start point of max infected
  summ.stats$m.max.it_s[i] <- mean(set_results_succ$max.it_s)
}

param.sets$model.type <- NA
param.sets$gamma <- NA
for(i in 1:nrow(param.sets)){
  temp.params <- param.sets[i,]
  model.type_gamma <- unlist(str_split(temp.params$model.type_gamma, "_"))
  param.sets$model.type[i] <- model.type_gamma[1]  # type of disease process; options are "SI", "SIR", or "SIS"
  
  if(model.type_gamma[1]=="SI"){
    param.sets$gamma[i] <- 0
  }else{
    param.sets$gamma[i] <- as.numeric(model.type_gamma[2])
  }
}

param.sets$param.set <- paste("set_", param.sets$set.num, sep = "")
summ.stats <- left_join(summ.stats, param.sets, by = "param.set")
all.results <- left_join(all.results, param.sets, by = "param.set")


weight.data <- subset(summ.stats, summ.stats$weights==T)
binary.data <- subset(summ.stats, summ.stats$weights==F)

w.si <- subset(weight.data, weight.data$model.type=="SI")
b.si <- subset(binary.data, binary.data$model.type=="SI")

w.sis <- subset(weight.data, weight.data$model.type=="SIS")
b.sis <- subset(binary.data, binary.data$model.type=="SIS")
  
w.sir <- subset(weight.data, weight.data$model.type=="SIR")
b.sir <- subset(binary.data, binary.data$model.type=="SIR")

### save these datasets for making manuscript figures ###
# save(w.si, file = "Output/Manuscript_Figures/data_for_plotting/w.si.Rdata")
# save(b.si, file = "Output/Manuscript_Figures/data_for_plotting/b.si.Rdata")
# save(w.sis, file = "Output/Manuscript_Figures/data_for_plotting/w.sis.Rdata")
# save(b.sis, file = "Output/Manuscript_Figures/data_for_plotting/b.sis.Rdata")
# save(w.sir, file = "Output/Manuscript_Figures/data_for_plotting/w.sir.Rdata")
# save(b.sir, file = "Output/Manuscript_Figures/data_for_plotting/b.sir.Rdata")


#### AUTOMATED GENERATION OF RESULTS HEATMAPS ####
model.types <- c("w.si", "b.si" , "w.sir", "b.sir", "w.sis", "b.sis")

for(i in 1:length(model.types)){
  
  temp.model <- model.types[i]
  print(temp.model)
  
  
  # read in data
  name <- paste("Output/Manuscript_Figures/data_for_plotting/", temp.model, ".Rdata", sep = "")
  temp.data <- get(load(name))
  
  temp.data.d <- subset(temp.data, temp.data$season=="Dry")
  temp.data.w <- subset(temp.data, temp.data$season=="Wet")
  temp.data_diff <- temp.data.d
  temp.data_diff[,c("m.dur.time", "m.total.i","prop.fail")] <- temp.data.d[,c("m.dur.time", "m.total.i","prop.fail")] - temp.data.w[,c("m.dur.time", "m.total.i","prop.fail")]
  
  
  
  #### MEAN DURATION ####
  if(temp.model %in% c("w.sir", "b.sir", "w.sis", "b.sis")){

    g1_st <- ggplot(temp.data, aes(x = c.rate, y = prob)) + geom_tile(aes(fill = m.dur.time), colour = "white") +
      scale_fill_viridis(discrete=F, option = "B", direction = 1, limits = c(1, 26)) + 
      ylab("Probability of transmission") +
      xlab("Edge weight scaling") 
    
    
    gg1_st <- g1_st + facet_grid(gamma ~ season) +
      theme(text = element_text(size=15)) + 
      theme(axis.text.x = element_text(size = 10),
            legend.position = "bottom") +
      labs(fill = "Mean duration (weeks)")
    gg1_st
    
    # differences
    g2_st <- ggplot(temp.data_diff, aes(x = c.rate, y = prob)) + geom_tile(aes(fill = m.dur.time), colour = "white") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, breaks = c(-2, 0, 2, 4), limits = c(-2.6, 4.4)) +
      ylab("Probability of transmission") +
      xlab("Edge weight scaling") 
    
    
    gg2_st <- g2_st + facet_grid(gamma~.) +
      theme(text = element_text(size=15)) + 
      theme(axis.text.x = element_text(size = 10),
            legend.position = "bottom") +
      labs(fill = "Dry mean difference (weeks)")
    gg2_st
    
    
    g <- ggarrange(gg1_st, gg2_st, labels = c("A", "B"), nrow = 1, widths = c(1.5, 1))
    g
    
    dur.name <- paste("Output/Manuscript_Figures/", temp.model,"_duration_full figure.jpeg", sep = "")
    ggsave(g, file = dur.name, width = 11, height = 11, units = "in", dpi = 300)
  }
  
  #### MEAN TOTAL INFECTED ####
  g1_st <- ggplot(temp.data, aes(x = c.rate, y = prob)) + geom_tile(aes(fill = m.total.i), colour = "white") +
    scale_fill_viridis(discrete=F, option = "B", direction = 1, limits = c(0.06, 1)) + 
    ylab("Probability of transmission") +
    xlab("Edge weight scaling") 
  
  gg1_st <- g1_st + facet_grid(gamma ~ season) +
    theme(text = element_text(size=15)) + 
    theme(axis.text.x = element_text(size = 10),
          legend.position = "bottom", legend.key.width = unit(1, "cm")) +
    labs(fill = "Mean total proportion infected")
  gg1_st
  
  
  g2_st <- ggplot(temp.data_diff, aes(x = c.rate, y = prob)) + geom_tile(aes(fill = m.total.i), colour = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, breaks = c(-0.08, 0, 0.08, 0.16), limits = c(-0.08, 0.2)) +
    ylab("Probability of transmission") +
    xlab("Edge weight scaling") 
  
  gg2_st <- g2_st + facet_grid(gamma~.) +
    theme(text = element_text(size=15)) + 
    theme(axis.text.x = element_text(size = 10),
          legend.position = "bottom", legend.key.width = unit(1, "cm")) +
    labs(fill = "Dry mean difference")
  gg2_st
  
  g <- ggarrange(gg1_st, gg2_st, labels = c("A", "B"), nrow = 1, widths = c(1.5, 1))
  g
  
  
  if(temp.model %in% c("w.si", "b.si")){
    total.name <- paste("Output/Manuscript_Figures/", temp.model,"_total infected_full figure.jpeg", sep = "")
    ggsave(g, file = total.name, width = 11, height = 5, units = "in", dpi = 300)
  }else{
    total.name <- paste("Output/Manuscript_Figures/", temp.model,"_total infected_full figure.jpeg", sep = "")
    ggsave(g, file = total.name, width = 11, height = 11, units = "in", dpi = 300)
  }
  
  
  #### MEAN PROPORTION FAILED ####
  
  g1_st <- ggplot(temp.data, aes(x = c.rate, y = prob)) + geom_tile(aes(fill = prop.fail), colour = "white") +
    scale_fill_viridis(discrete=F, option = "B", direction = -1, limits = c(0, 1)) + 
    ylab("Probability of transmission") +
    xlab("Edge weight scaling") 
  
  gg1_st <- g1_st + facet_grid(gamma ~ season) +
    theme(text = element_text(size=15)) + 
    theme(axis.text.x = element_text(size = 10),
          legend.position = "bottom", legend.key.width = unit(1, "cm")) +
    labs(fill = "Proportion failed")
  gg1_st
  
  
  
  
  g2_st <- ggplot(temp.data_diff, aes(x = c.rate, y = prob)) + geom_tile(aes(fill = prop.fail), colour = "white") +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, breaks = c(-0.2, -0.1, 0, 0.1), limits = c(-0.22, 0.12)) +
    ylab("Probability of transmission") +
    xlab("Edge weight scaling") 
  
  gg2_st <- g2_st + facet_grid(gamma~.) +
    theme(text = element_text(size=15)) + 
    theme(axis.text.x = element_text(size = 10),
          legend.position = "bottom") + 
    labs(fill = "Dry difference")
  gg2_st
  
  
  g <- ggarrange(gg1_st, gg2_st, labels = c("A", "B"), nrow = 1, widths = c(1.5, 1))
  g
  
  if(temp.model %in% c("w.si", "b.si")){
    fail.name <- paste("Output/Manuscript_Figures/", temp.model,"_prop failed_full figure.jpeg", sep = "")
    ggsave(g, file = fail.name, width = 11, height = 5, units = "in", dpi = 300)
  }else{
    fail.name <- paste("Output/Manuscript_Figures/", temp.model,"_prop failed_full figure.jpeg", sep = "")
    ggsave(g, file = fail.name, width = 11, height = 11, units = "in", dpi = 300)
  }
  
}

