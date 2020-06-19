#' Mizer extension aiming to introduce evolutionary processes to Mizer
#' Long story short, the algo add species from time to time to the existing ecosystem
#' First step is to set a mutation vector which tells mizer to stop at certain time so I can add species
#' TO REMEMBER: the saved step is not per dt but per time_step (so 10 times less) makes things easier but vary less
rm(list = ls())
require(mizerEvo)


#evolution param
folder <- file.path(tempdir(), "simTemp")
dir.create(folder)

params <- evoParams(no_sp = 5, RDD = "extinctionRDD")

sim <- evoProject(params = params, t_max = 300, mutation = 5, folder = folder)

plot(sim)


