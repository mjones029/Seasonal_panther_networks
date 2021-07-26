# Seasonal_panther_networks
Code for analyzing differences in panther seasonal home range overlap networks, and consequences for pathogen transmission

This repository features code and data used in the manuscript, "Seasonal changes in spatial connectivity alter predicted epidemic dynamics in a solitary carnivore."  

The analysis pipeline follows the following steps:

1. Analysis of puma home range overlap networks in wet and dry seasons.
2. Simulation of wet and dry season networks, followed by transmission of a range of pathogens on said networks.
3. Analysis of transmission outcomes in wet versus dry seasons.  

**Note:** All R scripts can be found in the "Scripts" folder.

## 1. Analysis of home range overlap networks.  
This step uses social network analysis to quantify differences in puma home range overlap network connectivity in wet versus dry seasons.  

Analysis is performed using:  
1. **Analyze_networks.R**, which references external function:
    1. ***bootstrap.node.metrics_clustersamp.R***  

## 2. Simulation of networks and transmission
This step simulates new wet and dry season networks based on characteristics of the observed networks. After networks are simulated, transmission of a range of theoretical pathogens is simulated through said networks.

Simulations are performed using:  
1. **Transmission_on_networks.R**, which references external functions:
    1. ***sim.SO.network.R***
    2. ***sim.dz.R***
    3. ***props.affected.R***

## 3. Analysis of transmission outcomes
This step analyzes the results from step 2 (transmission simulations) and generates heat map figures to visualize differences between wet and dry season outbreaks.

Analysis is performed using:
1. **Analyze_transmission.R**


## Markdown
Note that analysis steps (steps 1 and 3) can also be viewed as RMarkdown files with output. These can be found in the "Markdown" folder.
