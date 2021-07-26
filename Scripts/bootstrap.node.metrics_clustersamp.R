####  Preamble	####
# bootstrap.node.metrics_clustersamp.R
#
# ---
### title: Cluster sampling bootstrapping
# author: Marie Gilbertson
# date: "08/30/2020"
#---
# 
# What this code does:
# 1. Cluster-level bootstrap of panther node-level network metrics

bootstrap.node.metrics_clustersamp <- function(nsims = nsims, 
                                               metric = metric, 
                                               coefs = coefs,
                                               dataset = all.node.level,
                                               co.UDOI = co.UDOI,
                                               progress = T
                                               ){
  
  # set up progress bar
  if(progress){
    pb <- txtProgressBar(min = 0,
                         max = nsims/10,
                         style = 3,
                         width = nsims/10, # Needed to avoid multiple printings
                         char = "=") 
  }
  
  # set up output
  betas<- data.frame(matrix(NA,	nrow	=	nsims,	ncol	=	length(coefs))) # matrix to hold estimates of betas
  colnames(betas) <- coefs # label columns
  
  metric <- metric
  
  # Loop to repeat bootstrapping and extract betas from each linear regression
  cluster.sizes <- data.frame(table(all.node.level$node.id))
  colnames(cluster.sizes) <- c("node.id", "cluster.size")
  
  cluster.nums <- data.frame(table(table(all.node.level$node.id)))
  colnames(cluster.nums) <- c("cluster.size", "freq")
  
  for(j in 1:nsims){
    if(progress){
      setTxtProgressBar(pb, j/10)
    }
    
    bootDat <- NULL
    # create bootstrap sample data, sampling based on cluster size
    for (i in 1:nrow(cluster.nums)){
      temp.cluster.size <- cluster.nums$cluster.size[i]
      
      temp.clusters <- subset(cluster.sizes, cluster.sizes$cluster.size==temp.cluster.size)
      
      bootIDs <- data.frame(node.id = sample(x = temp.clusters$node.id, size = cluster.nums$freq[i], replace = TRUE)) # Take a sample with replacement
      temp.bootDat <- merge(bootIDs, dataset, by = "node.id", type = "left") # Use this to sample from original data by beach
      
      bootDat <- rbind(bootDat, temp.bootDat)
      
    }
    if(nrow(dataset) != nrow(bootDat)){
      print("Warning! Bootstrap sample is not the right size!")
    }
    
    # run the linear regression on the latest bootstrapped dataset
    lm.temp <- lm(formula(paste(metric, "~ log(terr.km) + season.only", sep = "")), data=bootDat) # Linear regression from this new sample
    betas[j,] <- coef(lm.temp) # extract intercept and slope estimates
    
  }
  
  lmfit <- lm(formula(paste(metric, "~ log(terr.km) + season.only", sep = "")), data=dataset)
  
  
  #### Optional plotting ####
  ### NOTE: would need to create appropriate directories for this plotting code to run appropriately 
  # # Plot distribution of bootstrap estimates for intercept (beta0)
  # plot.name <- paste("bootstrap_plots/cluster_bootstrap/", metric, "_", coefs[1], "_UDOI_", co.UDOI,".jpg", sep = "")
  # jpeg(plot.name)
  # hist(betas[,1], col="gray",xlab="", main=expression(paste("Sampling Distribution of ", hat(beta)[0])))
  # abline(v=coef(lmfit)[1], col = "red") # add original intercept estimate
  # dev.off()
  # 
  # # Plot distribution of bootstrap estimates for terr.km (beta1)
  # plot.name <- paste("bootstrap_plots/cluster_bootstrap/", metric, "_", coefs[2], "_UDOI_", co.UDOI,".jpg", sep = "")
  # jpeg(plot.name)
  # hist(betas[,2], col="gray",xlab="", main=expression(paste("Sampling Distribution of ", hat(beta)[1])))
  # abline(v=coef(lmfit)[2], col = "red") # add original beta1 estimate
  # dev.off()
  # 
  # # Plot distribution of bootstrap estimates for season (beta2)
  # plot.name <- paste("bootstrap_plots/cluster_bootstrap/", metric, "_", coefs[3], "_UDOI_", co.UDOI,".jpg", sep = "")
  # jpeg(plot.name)
  # hist(betas[,3], col="gray",xlab="", main=expression(paste("Sampling Distribution of ", hat(beta)[2])))
  # abline(v=coef(lmfit)[3], col = "red") # add original beta2 estimate
  # dev.off()
  
  # Bootstrap confidence intervals (e.g., using quantile function)
  quantiles <- data.frame(matrix(NA,	nrow	=	length(coefs),	ncol	=	3)) # matrix to hold estimates of betas
  colnames(quantiles) <- c("coef", "lower.quant", "upper.quant") # label columns
  quantiles$coef <- coefs
  
  quantiles[1,c("lower.quant", "upper.quant")] <- quantile(betas[,1], prob=c(0.025, 0.975)) # Bootstrap CI for beta0 
  quantiles[2,c("lower.quant", "upper.quant")] <- quantile(betas[,2], prob=c(0.025, 0.975)) # Bootstrap CI for beta1
  quantiles[3,c("lower.quant", "upper.quant")] <- quantile(betas[,3], prob=c(0.025, 0.975)) # Bootstrap CI for beta1
  
  
  results <- list(metric, betas, quantiles)
  names(results) <- c("metric", "betas", "quantiles")
  return(results)
}