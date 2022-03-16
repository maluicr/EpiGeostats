ptdata <- read.table("data-raw/covid19pt.txt", header = T, sep = "\t", dec = ".")
save(ptdata, file="data/covid19pt.RData")
