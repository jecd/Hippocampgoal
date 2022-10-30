# Hippocampgoal
[![DOI](https://zenodo.org/badge/406137918.svg)](https://zenodo.org/badge/latestdoi/406137918)

This repo contains the code and data to reproduce main and supplemental figures in [Crivelli-Decker et al., 2021](https://www.biorxiv.org/content/10.1101/2021.08.18.456881v2)
This is a living repository and will change in response to peer review 
There are still things that need to be adjusted, if you run into an please issue submit an issue on github

# Dependencies
- Matlab 2020a
- R > 3.6.0
- R Packages: ggplot2, plyr, dplyr, tidyr, afex, emmeans
- See r_env.txt for full list of version numbers and information about how to recreate environment

# OS 
- MacOS Catalina 10.15.7

# Contents 
- Colormaps: Code to repro colormaps in publication figures
- Data: contains data to reproduce main and supplemental figures
- Figures: Source figures from manuscript
- Scripts: Scripts to reproduce figures
- Stats: Statistics files created from cluster based permutation tests

# Instructions
- To reproduce figures 4, S3, S4 in the manuscript run **TR_TR_analyses.m**
- To reproduce SR simulations in figure S1 run **SR_simulation.m**
- To reproduce figures 2, S2 in the manuscript run **MixedModelPlots.R**
- Note paths will need to be changed to correct path on your local machine 