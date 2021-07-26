# props affected
#
#========================================================	
# ---
### title: Proportions affected in outbreak
# author: Marie Gilbertson
# date: "08/31/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function that processes outbreak results; adapted from FIV-->FeLV project

props_affected <- function(m.new = m.new){
  
  #### prep data ####
  # NOTE: read in "m.new" from "post_process_outbreak_data_031320.R" script
  # new dataframe of m.new to work with (and conserve original m.new data)
  p.data <- m.new
  
  pop.size <- ncol(p.data)
  
  
  # make data frame for population results
  p.results <- data.frame(
    time=numeric(nrow(m.new)),
    pop.size=numeric(nrow(m.new)),
    prop.s=numeric(nrow(m.new)),
    prop.i=numeric(nrow(m.new)),
    prop.r=numeric(nrow(m.new))
  )
  
  p.results$time <- seq(0,nrow(m.new)-1)
  
  
  #### calculate population size per time step ####
  # population size per time, subtracting all dead and not-born-yet individuals
  p.results$pop.size <- pop.size
  
  
  
  
  #### calculate propotion of population in each disease category per time step ####
  
  p.results$prop.s <- apply(p.data==0, 1, sum)/p.results$pop.size
  p.results$prop.i <- apply(p.data==1, 1, sum)/p.results$pop.size
  p.results$prop.r <- apply(p.data==2, 1, sum)/p.results$pop.size
  
  
  #### output data results ####
  
  return(p.results)
  
}