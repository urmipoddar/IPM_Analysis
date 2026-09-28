# IPM_Analysis: Post-fire pitch pine demography

## Non-technical summary

Pitch pine (*Pinus rigida*), the defining tree of the pine barrens of Long Island and New Jersey, is well known for being adapted to fire. But what happens to these tree populations after severe wildfires that kill almost all adult trees? This project follows nearly 5,000 pitch pine seedlings that came up after a severe wildfire in 1995 on Long Island, NY, tracking which ones survived and how they grew over 13 years, to understand how pitch pine populations respond to and recover from fire. Seedlings that emerged after the wildfire were tagged and observed until 2009, with survival, height, stem and crown diameter, cone production etc. being recorded at various censuses. Using this data, I am building an integral projection model (IPM), a type of model for predicting population change over time, to assess how these tree populations recovered after the 1995 fire.


## Technical summary

This repository contains the data-processing and modeling code for an integral projection model (IPM) of pitch pine (*Pinus rigida*) populations following stand-replacing wildfire in the Pine Barrens of Long Island, NY. *Pinus rigida* is a fire-adapted species, however, decades of fire suppression in the Long Island pine barrens resulted in severe wildfires in 1995. This resulted in high tree and seed mortality, and may have had long-term effects on population trajectories. Following these fires, emerging seedlings were tagged and censused 15 times between June and June 2009. Survival, height, stem and crown diameter, cone production, and cone serotiny were recorded at various censuses. Using this data, I am the parameterizing size-dependent survival, growth and reproduction submodels of an IPM to capture how this species' demography changes with plant height and time since fire. These will then be used to build an IPM to assess population recovery and project future population trajectories. 

