#### Preamble ####
# Transmission_on_networks.R
#
# ---
### title: Script A for Objective 2: effects of seasonal network differences on pathogen transmission
# author: Marie Gilbertson
# date: "08/30/2020"
#---
# 
# What this code does:
# 1. Simulates transmission on networks

#### Clear Environment ####
remove(list=ls())

#### Load libraries ####
library(fitdistrplus)
library(ergm)
library(intergraph)
library(stringr) # for "str_split"
library(reshape2)
library(ggplot2)


#### Load external function ####
# network and disease simulation functions
source("Scripts/sim.SO.network.R")
source("Scripts/sim.dz.R")
source("Scripts/props.affected.R")

#### Set seed ####
set.seed(6252)



#### DETERMINE NETWORK STRUCTURE FOR SIMULATIONS ####

#### read in previous network analysis results
# load seasons vector
seasons <- c("Wet_1996", "Dry_1996", "Wet_1997", "Dry_1997", "Wet_1998", "Dry_1998", "Wet_1999", "Dry_1999", "Wet_2000",
             "Dry_2000", "Wet_2001", "Dry_2001", "Wet_2002", "Dry_2004", "Wet_2005", "Dry_2005", "Wet_2006", "Dry_2006",
             "Dry_2002", "Wet_2003", "Dry_2003", "Wet_2004")


# choose UDOI cutoff
co.UDOI <- 0


#### loop to load SNA data and fit network structure distributions ####
summ.stats <- data.frame(matrix(nrow = length(seasons), ncol = 8))
colnames(summ.stats) <- c("season", "net.size", "net.den", "nd.s1", "nd.s2", "str.shape", "str.rate", "isos")

for(j in 1:length(seasons)){
  print(j)

  ndl.results.filename <- paste("Output/SNA_results/", seasons[j], "_", co.UDOI, "_UDOI_nodelevel_results.Rdata", sep="")
  node.data <- get(load(ndl.results.filename))

  # fit beta distribution to normalized degree distribution
  norm.deg <- node.data$norm.deg[node.data$norm.deg>0]
  y1 <- fitdist(norm.deg, "beta")


  # fit gamma distribution to strength distribution (for those nodes with strength>0)
  str <- node.data$std.str[node.data$std.str>0]
  y2 <- fitdist(str, "gamma")

  
  # load network-level data
  nw.results.filename <- paste("Output/SNA_results/", seasons[j], "_", co.UDOI, "_UDOI_networklevel_results.Rdata", sep="")
  nw.data <- get(load(nw.results.filename))

  # save results
  summ.stats$season[j] <- seasons[j]
  summ.stats$net.size[j] <- nw.data$net.size
  summ.stats$net.den[j] <- nw.data$dens
  summ.stats$nd.s1[j] <- y1$estimate[1]
  summ.stats$nd.s2[j] <- y1$estimate[2]
  summ.stats$str.shape[j] <- y2$estimate[1]
  summ.stats$str.rate[j] <- y2$estimate[2]
  summ.stats$isos[j] <- nw.data$n.isolates
}


summ.stats$season.only <- as.factor(substr(summ.stats$season, start = 1, stop = 3))
wet.stats <- subset(summ.stats, summ.stats$season.only=="Wet")
dry.stats <- subset(summ.stats, summ.stats$season.only=="Dry")


# fit poisson to isolates for each season-type
# a random draw from these distributions will give the number of isolates per simulated network
y3 <- fitdist(wet.stats$isos, "pois")

y4 <- fitdist(dry.stats$isos, "pois")

### Store mean distribution results
mean.stats <- data.frame(matrix(nrow = 2, ncol = ncol(summ.stats)-1))
colnames(mean.stats) <- c(colnames(summ.stats[c(2:7,9)]), "iso.lambda")

mean.stats[1,] <- c(colMeans(wet.stats[,c(2:7)]), paste(wet.stats$season.only[1]), y3$estimate)
mean.stats[2,] <- c(colMeans(dry.stats[,c(2:7)]), paste(dry.stats$season.only[1]), y4$estimate)



#### loop to load edgelist data and fit UDOI edge weight distributions ####
### determine which distribution consistently fits UDOI data the best
AIC.results <- data.frame(matrix(nrow = length(seasons), ncol = 4))
colnames(AIC.results) <- c("season", "gamma", "exp", "lnorm")

### empty dataframe to store results
udoi.dists <- data.frame(matrix(nrow = length(seasons), ncol = 3))
colnames(udoi.dists) <- c("season", "udoi.shape", "udoi.rate")


for(j in 1:length(seasons)){
  print(j)

  ### read in current edgelist
  edgelist.filename <- paste("Data/UDOI_edgelists/", seasons[j], "_edgelist.Rdata", sep="")
  el <- get(load(edgelist.filename))

  ### subset edgelist by UDOI cutoff
  el2 <- subset(el, el$UDOI>co.UDOI)

  ### fit test distributions to UDOI
  z1 <- fitdist(el2$UDOI, "gamma") # gamma distribution consistently fits the best
  z2 <- fitdist(el2$UDOI, "exp")
  z3 <- fitdist(el2$UDOI, "lnorm")

  ### store AICs for each test distribution
  AIC.results$season[j] <- seasons[j]
  AIC.results$gamma[j] <- z1$aic
  AIC.results$exp[j] <- z2$aic
  AIC.results$lnorm[j] <- z3$aic


  ### store parameters for gamma distribution
  udoi.dists$season[j] <- seasons[j]
  udoi.dists$udoi.shape[j] <- z1$estimate[1]
  udoi.dists$udoi.rate[j] <- z1$estimate[2]

}


## check which distribution provides lowest AICs across seasons
long <- melt(AIC.results, id.var = "season")
p <- ggplot(long, aes(x = as.factor(season), y = value, colour = variable)) + geom_point()
p
# gamma consistently has lowest AIC, so use gamma distribution for UDOI

## extract means for distributions
udoi.dists$season.only <- as.factor(substr(udoi.dists$season, start = 1, stop = 3))
wet.udoi <- subset(udoi.dists, udoi.dists$season.only=="Wet")
dry.udoi <- subset(udoi.dists, udoi.dists$season.only=="Dry")

mean.stats$udoi.shape[mean.stats$season.only=="Wet"] <- mean(wet.udoi$udoi.shape)
mean.stats$udoi.rate[mean.stats$season.only=="Wet"] <- mean(wet.udoi$udoi.rate)
mean.stats$udoi.shape[mean.stats$season.only=="Dry"] <- mean(dry.udoi$udoi.shape)
mean.stats$udoi.rate[mean.stats$season.only=="Dry"] <- mean(dry.udoi$udoi.rate)


## save "mean.stats" object for use in simulations
save(mean.stats, file = "Output/distribution stats for network sims_UDOI 0.Rdata")






#### TRANSMISSION SIMULATIONS ####

# load parameter sets
param.sets <- get(load("Data/simulation_paramsets.Rdata"))

# set parameters to apply to all sets/simulations
nsims <- 100
duration.yrs <- 0.5   # number of years disease simulations should last

# set which parameter sets to run
params.to.run <- seq(1, 420)


#### start of "q" loop for running multiple parameter sets ####
for(q in 1:length(params.to.run)){
  
  iter <- params.to.run[q]
  
  #### set parameter set number and pull parameters ####
  # ensures select correct parameter set, even when running discontinuous parameter sets, etc
  real.setnos <- param.sets$set.num
  param.set.num <- real.setnos[iter]
  print(paste("iter ", iter, "_set ", param.set.num, sep = ""))
  
  
  # pull parameters
  params <- param.sets[param.sets$set.num==param.set.num,]
  
  
  #### Set seed ####
  sim.seed <- 57245+param.set.num
  set.seed(sim.seed) # ensure unique but reproducible seeding
  
  
  #### assign network simulation parameters ####

  # set network and disease simulation parameters
  sim.season <- params$season # options are "Wet" or "Dry"
  co.UDOI <- params$co.UDOI # cutoff for UDOI for loading "mean stats" object
  
  # load distribution data
  mean.stats.filename <- paste("Output/distribution stats for network sims_UDOI ", co.UDOI, ".Rdata", sep = "")
  mean.stats <- get(load(mean.stats.filename))
  
  
  #### assign disease transmission parameters ####
  model.type_gamma <- unlist(str_split(params$model.type_gamma, "_"))
  model.type <- model.type_gamma[1]  # type of disease process; options are "SI", "SIR", or "SIS"
  gamma <- as.numeric(model.type_gamma[2])        # recovery rate (1/duration of infection) in weeks; does not apply for SI model.types
  
  prob <- params$prob         # probability of transmission, given infectious contact
  c.rate <- params$c.rate       # weekly probability of a contact
  weights <- params$weights        # should UDOI-based edge weights be used to scale contact probability? T or F
  
  
  
  # prep dataframe for storing simulation results for this parameter set
  full.sims.results <- data.frame(sim.num = numeric(),
                                 dur.time = numeric(),
                                 total.i = numeric(),
                                 total.r = numeric(),
                                 final.s = numeric(),
                                 final.i = numeric(),
                                 max.i = numeric(),
                                 max.it_s = numeric(),
                                 max.it_e = numeric(),
                                 num.failed = numeric()
                                )
  
  
  for(z in 1:nsims){  
    print(paste("sim number ", z, sep = ""))
    sim.num <- z
    
    #### SIMULATE NETWORK ####
    g.sim <- sim.SO.network(sim.season = sim.season,
                            net.size = 33,
                            mean.stats = mean.stats,
                            method = method
                            )
    
    # convert to igraph object
    g <- asIgraph(g.sim)

    

    #### SIMULATE DISEASE ON SIMULATED NETWORKS ####
    suppressMessages(library(igraph)) # must detach igraph in network simulation function, so reload here
    
    m.list <- sim.dz(g = g,
                     model.type = model.type,
                     duration.yrs = duration.yrs,
                     prob = prob,                      
                     c.rate = c.rate,                  
                     gamma = gamma,
                     weights = weights
                     )
    
    
    
    #### process disease simulation results ####
    
    ## convert results to population proportions
    m.new <- m.list[[1]]
    
    p.results <- props_affected(m.new=m.new)
    
    
    # save results
    g.name <- paste("Output/Simulation_Results/g network_paramset ", param.set.num, "_sim ",sim.num, ".Rdata", sep = "")
    save(g, file = g.name)
    
    m.name <- paste("Output/Simulation_Results/mlist data_paramset ", param.set.num, "_sim ", sim.num, ".Rdata", sep = "")
    save(m.list, file = m.name)
    
    mnew.name <- paste("Output/Simulation_Results/mnew data_paramset ", param.set.num, "_sim ", sim.num, ".Rdata", sep = "")
    save(m.new, file = mnew.name)
    
    presults.name <- paste("Output/Simulation_Results/presults data_paramset ", param.set.num, "_sim ", sim.num, ".Rdata", sep = "")
    save(p.results, file = presults.name)
    
    
    
    
    #### Store primary outcomes of interest #####
    
    # make sure pop.size is loaded
    pop.size <- ncol(m.new)
    
    # duration of outbreak = time point at which there are 0 individuals with status = 1
    dur.data <- p.results
    dur.time <- (which(dur.data$prop.i==0)[1])-1 # first instance where proportion infectious = 0 (-1 because first entry is time=0)
    # (will give "NA" if there is no time at which prop.i = 0)
    dur.time <- ifelse(is.na(dur.time), duration.yrs*52, dur.time)
    
    # total proportion of individuals EVER infected (individuals that were ever status 1); 
    # note: this means s+i+r does not necessarily = pop.size
    total.i <- sum(apply(m.new==1, 2, any))/pop.size
    
    
    # total proportion recovered (individuals that were ever status 2)
    total.r <- sum(apply(m.new==2, 2, any))/pop.size
    
    
    # total proportion still susceptible at end of outbreak
    last.time <- tail(m.new, 1)
    final.s <- length(which(last.time==0))/pop.size
    
    
    # total still infectious at end of simulation
    # most relevant to outbreaks that are ongoing at end of simulation
    final.i <- length(which(last.time==1))/pop.size
    
    # max proportion infected at one time and the associated time point
    max.i <- max(p.results$prop.i)
    max.it <- p.results$time[which(p.results$prop.i==max(p.results$prop.i))]
    max.it <- c(head(max.it,1), tail(max.it,1))
    
    
    # number of "failed" epidemics (number that initiated with an isolate)
    # (-1 is because initiate.dur gave the number of iterations before a non-isolate was selected
    # meaning that the number of times an isolate was selected = initiate.dur-1)
    num.failed <- m.list[[2]]-1
    
    
    
    
    temp.results <- data.frame(sim.num = sim.num,
                               dur.time = dur.time,
                               total.i = total.i,
                               total.r = total.r,
                               final.s = final.s,
                               final.i = final.i,
                               max.i = max.i,
                               max.it_s = max.it[1],
                               max.it_e = max.it[2],
                               num.failed = num.failed
                              )
    
    
    full.sims.results[z,] <- temp.results

  }
  

  full.sims.results$param.set <- paste("set_", param.set.num, sep = "")
  # Save final results
  full.name <- paste("Simulation_Results/full set results_paramset ", param.set.num, ".Rdata", sep = "")
  save(full.sims.results, file = full.name)

}

