# blockdss

An R package to run block direct sequential simulation with R to model the spatial distribution of a
disease, as in Azevedo et al. 2020 (https://doi.org/10.1186/s12942-020-00221-5). For that purpose we will follow an example using COVID-19 data where
simulated risk maps, median risk map and risk uncertainty map are obtained in the end.

# data

A set of functions written in R will generate the required files and call dss.c.64.exe for block sequential
simulation. To run the code, you will need a disease dataset (data.frame) and a georeferenced grid with id region values at all simulation nodes (locations).
