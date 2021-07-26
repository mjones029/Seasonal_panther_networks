# Function for simulation networks based on observed wet and dry season spatial overlap networks

sim.SO.network <- function(sim.season = sim.season,
                           net.size = 33,
                           mean.stats = mean.stats,
                           method = method
                           ){

  # must unload igraph to assign edge attributes
  if( "tnet" %in% (.packages())) {
    detach("package:tnet", unload=TRUE) 
  }
  
  if( "igraph" %in% (.packages())) {
    detach("package:igraph", unload=TRUE) 
  }
  
  
  if(sim.season=="Wet"){
    temp.dists <- mean.stats[mean.stats$season.only=="Wet",]  
  }else{
    temp.dists <- mean.stats[mean.stats$season.only=="Dry",]
  }
  # make sure correct numeric formatting
  temp.dists[,c(1:6,8:10)] <- as.numeric(temp.dists[,c(1:6,8:10)])
  
  # if constraining by average density, set that here
  temp.dens <- mean(as.numeric(mean.stats$net.den))
  
  
  # set number of isolates
  isos <- rpois(1, lambda = temp.dists$iso.lambda)
  

  # draw normalized degree distribution for non-isolates
  sim_norm.deg.dist <- rbeta(n = net.size-isos, shape1 = temp.dists$nd.s1, shape2 = temp.dists$nd.s2)
  # add isolates and "un-normalize" degree by multiplying by (n-1) where n is the number of vertices (nodes) in the graph
  sim_deg.dist <- c(rep(0, isos), round(sim_norm.deg.dist*(net.size-isos-1)))
  # note: due to rounding, can end up with more isolates than those listed by "isos" object
  

  # simulate network based on target degree distribution
  targets <- as.vector(table(sim_deg.dist))
  rand.degs <- sort(unique(sim_deg.dist))
  
  # make sure simulated network is of reasonable density
  repeat{
    x <- network(net.size, directed = FALSE, density = temp.dists$net.den)
    g <- san(x~degree(rand.degs), target.stats = targets, constraints = ~edges) # tries to maintain mean degree
    
    # determine if network density is constrained by seasonal values or overall mean
    # if(network.density(g)>(0.75*temp.dists$net.den) & network.density(g)<(1.25*temp.dists$net.den)) break
    if(network.density(g)>(0.75*temp.dens) & network.density(g)<(1.25*temp.dens)) break
  }
  
  
  # assign edge weight based on UDOI distribution
  sim_udois <- rgamma(n = network.edgecount(g), shape = temp.dists$udoi.shape, rate = temp.dists$udoi.rate)
  scaled_sim_udois <- sim_udois/max(sim_udois)
  set.edge.attribute(g, "udoi", sim_udois, e=seq_along(g$mel))
  set.edge.attribute(g, "weight", scaled_sim_udois, e=seq_along(g$mel))
  
  return(g)

}