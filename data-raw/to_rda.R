
# author: m ribeiro;
# date : 18-03-22; last revision : 18-03-2022
# email: manuel.ribeiro@tecnico.ulisboa.pt

# description:
# convert files .txt to .rda to include in folder 'data' of blockdss package.

ptdata <- read.table("./data-raw/ptdata.txt", header = T, dec = ".", sep = "\t")
save(ptdata, file ="data/ptdata.rda")


ptgrid <- read.table("./data-raw/ptgrid.txt", header = T, dec = ".", sep = "\t")
save(ptgrid, file ="data/ptgrid.rda")
