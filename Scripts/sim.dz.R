# sim.dz.R
#
#========================================================	
# ---
### title: Function for simulating transmission on networks
# author: Marie Gilbertson
# date: "08/30/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Simulates disease transmission on networks
### from Jenny's network SIR code





sim.dz <- function(g = g,                          # network on which to simulate (should be an igraph object)
                 model.type = model.type,          # SI, SIR, or SIS
                 duration.yrs = duration.yrs,      # number of years to run simulation
                 prob = prob,                      # probability of transmission, given infectious contact
                 c.rate = c.rate,                  # weekly probability of contact; to be scaled by edge weight
                 gamma = gamma,                    # recovery rate (1/duration of infection) for binomial distribution draws
                 weights = weights                 # should UDOI-based edge weights be used to scale contact probability?
                 ) {
  
  
  # assign pop.size value 
  pop.size <- length(V(g))
  

  
  ############# Network model ########################
  
  if(weights==T){
    ER<-get.adjacency(g, attr = "weight") #make an adjacency network; convert graph to an adjacency matrix or an edge list
  }else if(weights==F){
    ER<-get.adjacency(g, attr = NULL)
  }
  
  netwb <- ER
  
  
  #### assign dz model parameters ####
  prob <- prob # probability of transmission given infectious contact
  duration <- 52*duration.yrs # 52 * number of years to run simulation
  

  
  #### run the disease simulation ####
  time<-1
  gp <- ncol(netwb)    # number of columns in network    
  
  # randomly select individual initiating infection in the population
  # for efficiency, only keep instances in which infection initiates in a non-isolate
  # but keep record of how many iterations before a non-isolate is selected
  
  initiate <- c()
  repeat {
    
    stat<-sample(c(numeric(gp-1),1));   #choose first node to get infected. This shuffles the order of one 1 and the rest zeros (to make bigger network, need to add more 0's)
    #sample function can be used to return a random permutation of a vector
    #numeric creates a real vector of the specified length. The elements of the vector are all 0.
    #c combines values into a vector or list
    
    initiate <- c(initiate, which(stat==1))
    
    # exit loop when infection initiates in a non-isolate
    if (any(netwb[which(stat==1),]>0) & any(netwb[,which(stat==1)]>0)) break
  }
  
  # "initiate.dur" gives the number of iterations before a non-isolate was selected
  # in other words, the number of "failed" epidemics + 1
  initiate.dur <- length(initiate)
  
  
  m <- matrix(0,duration+1,gp)  # make empty matrix for number of time steps long by number of individuals wide (this will be the progression of infection per time step per indiv.) 52 plus 1 more for time=0
  m[time,] <- stat
  
  # create matrix to accommodate duration of infection assignments
  n <- matrix(0,duration*10,gp) # need this to be extra long to accommodate duration of infection past end of monitoring
  
  # assign duration of infection for first infected individual
  if(model.type=="SIR" | model.type=="SIS"){
    recovery.all <- c()
    repeat {
      # repeatedly draw until recovery = 1, recording number of draws until that point
      recovery <- rbinom(1,1,gamma)
      recovery.all <- c(recovery.all, recovery)
      
      # exit loop when death = 1
      if (recovery > 0) break
    }
    dur.1 <- length(recovery.all)
    n[c(1:dur.1),which(stat==1)] <- 1
  
  }else if(model.type=="SI"){
    n[c(1:nrow(n)),which(stat==1)] <- 1
  }
  
  # optional
  write(c(0,stat),file="test_sim.txt",ncolumns=(gp+1),append=F) #writes initial conditions to file
  
  
  for (l in 1:duration) {  
    # print(l)
    #0 means susceptible; 1 means infectious; 2 recovered
    
    # for each (weekly) time step
    statc <- stat  #assign stat to copy of stat so can change 'statc' while still working off of 'stat'
    

    for (i in (which(stat==0))){     #for every susceptible individual in original stat...
      
      if(rbinom(1,1,c.rate)==1){    #does that individual have a contact this time step? If so...
        for (j in (1:gp)[-i]){          #for each j in 1 to number groups, with exclusion of itself (because that is always 0)
          
          if ((rbinom(1,1,netwb[i,j])==1) & ((stat[j]==1)&(rbinom(1,1,prob)==1)) ) #is there an edge to an infected individual and does a transmission event take place?  (By looking down column of that indiv in netw for a 1)  
          {
            
            statc[i] <- 1
            
            # assign duration of infection based on model type
              if(model.type=="SIR" | model.type=="SIS"){ # SIR and SIS become non-infectious at some point
                recovery.all <- c()
                repeat {
                  # repeatedly draw until recovery = 1, recording number of draws until that point
                  recovery <- rbinom(1,1,gamma)
                  recovery.all <- c(recovery.all, recovery)
                  
                  # exit loop when death = 1
                  if (recovery > 0) break
                }
                dur.temp <- length(recovery.all)
                n[c((time+1):(time+dur.temp)),i] <- 1
                
                }else if(model.type=="SI"){ # in SI model, all infectious are infectious forever
                  n[c((time+1):nrow(n)),i] <- 1
              }
            break # assign outcome of infection to the copy and get out of loop
            }
            
          }
        }
      }

    

    # if model.type is SIR or SIS, determine which individuals recover from infection or become re-susceptible
    if(model.type=="SIR" | model.type=="SIS"){
      for (i in (which(stat==1)))     #for every infectious individual in original stat
      {
        if (n[time+1,i]==0) # has that individual reached the duration of their infectiousness?
        {
          
          if(model.type=="SIR"){ # if SIR, they transition to state 2 (recovered)
          
            statc[i] <- 2

            # then assign that death duration in the n matrix
            n[c((time+1):nrow(n)),i] <- 2 
          
          }else if(model.type=="SIS"){ # if SIS, they transition back to state 0 (susceptible)
            
            statc[i] <- 0
            # don't need to update n matrix because already 0
            
          }
          
        }
      }
    }
    
    
    
    stat <- statc
    m[time+1,] <- stat
    write(c(time,m[time+1,]),file="test_sim.txt",ncolumns=(gp+1),append=TRUE)
    time<-(time+1)
  }
  
  m.list <- list(m, initiate.dur)
  names(m.list) <- c("m", "initiate.dur")
  return(m.list)
  
}