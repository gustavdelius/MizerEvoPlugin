#' Set up parameters for an evolutionary model
#'
#' this will either set up new params as usual or convert Mizer params to MizerEvo params if given a mizer.params object (TODO)
#'
#' @param no_sp The number of species in the model.
#' @param min_w_inf The asymptotic size of the smallest species in the
#'   community. This will be rounded to lie on a grid point.
#' @param max_w_inf The asymptotic size of the largest species in the community.
#'   This will be rounded to lie on a grid point.
#' @param min_w The size of the the egg of the smallest species. This also
#'   defines the start of the community size spectrum.
#' @param max_w The largest size in the model. By default this is set to the
#'   largest asymptotic size `max_w_inf`. Setting it to something larger
#'   only makes sense if you plan to add larger species to the model later.
#' @param eta Ratio between maturity size and asymptotic size of a species.
#'   Ignored if `min_w_mat` is supplied. Default is 10^(-0.6),
#'   approximately 1/4.
#' @param min_w_mat The maturity size of the smallest species. Default value is
#'   \code{eta * min_w_inf}. This will be rounded to lie on a grid point.
#' @param no_w The number of size bins in the community spectrum. These bins
#'   will be equally spaced on a logarithmic scale. Default value is such that
#'   there are 50 bins for each factor of 10 in weight.
#' @param min_w_pp The smallest size of the resource spectrum. By default this
#'   is set to the smallest value at which any of the consumers can feed.
#' @param w_pp_cutoff The largest size of the resource spectrum. Default value
#'   is max_w_inf unless \code{perfect_scaling = TRUE} when it is Inf.
#' @param n Scaling exponent of the maximum intake rate.
#' @param p Scaling exponent of the standard metabolic rate. By default this is
#'   equal to the exponent `n`.
#' @param lambda Exponent of the abundance power law.
#' @param r_pp Growth rate parameter for the resource spectrum.
#' @param kappa Coefficient in abundance power law.
#' @param alpha The assimilation efficiency of the community.
#' @param ks Standard metabolism coefficient. If not provided, default will be
#'   calculated from critical feeding level argument `fc`.
#' @param fc Critical feeding level. Used to determine `ks` if it is not given
#'   explicitly.
#' @param h Maximum food intake rate.
#' @param beta Preferred predator prey mass ratio.
#' @param sigma Width of prey size preference.
#' @param f0 Expected average feeding level. Used to set `gamma`, the
#'   coefficient in the search rate. Ignored if `gamma` is given
#'   explicitly.
#' @param gamma Volumetric search rate. If not provided, default is determined
#'   by [get_gamma_default()] using the value of `f0`.
#' @param zeta ...
#' @param ext_mort_prop The proportion of the total mortality that comes from
#'   external mortality, i.e., from sources not explicitly modelled. A number in
#'   the interval [0, 1).
#' @param R_factor The factor such that \code{R_max = R_factor * R}, where `R_max`
#'   is the maximum reproduction rate allowed and `R` is the steady-state
#'   reproduction rate. Thus the larger `R_factor` the less the impact of the
#'   density-dependence.
#' @param gear_names The names of the fishing gears for each species. A
#'   character vector, the same length as the number of species.
#' @param knife_edge_size The minimum size at which the gear or gears select
#'   fish. A single value for each gear or a vector with one value for each
#'   gear.
#' @param RDD ...
#' @param egg_size_scaling If TRUE, the egg size is a constant fraction of the
#'   maximum size of each species. This fraction is \code{min_w / min_w_inf}. If
#'   FALSE, all species have the egg size `w_min`.
#' @param resource_scaling If TRUE, the carrying capacity for larger resource
#'   is reduced to compensate for the fact that fish eggs and larvae are
#'   present in the same size range.
#' @param perfect_scaling If TRUE then parameters are set so that the community
#'   abundance, growth before reproduction and death are perfect power laws. In
#'   particular all other scaling corrections are turned on.
#' @export
#' @return An object of type `MizerParams`
#' @family functions for setting up models
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' sim <- project(params, t_max = 5, effort = 0)
#' plotSpectra(sim)
#' }
evoParams <- function(no_sp = 11,
                      min_w_inf = 10,
                      max_w_inf = 10 ^ 4,
                      min_w = 10 ^ (-3),
                      max_w = max_w_inf,
                      eta = 10^(-0.6),
                      min_w_mat = min_w_inf * eta,
                      no_w = log10(max_w_inf / min_w) * 20 + 1,
                      min_w_pp = 1e-10,
                      w_pp_cutoff = min_w_mat,
                      n = 2 / 3,
                      p = n,
                      lambda = 2.05,
                      r_pp = 0.1,
                      kappa = 0.005,
                      alpha = 0.4,
                      h = 40,
                      beta = 100,
                      sigma = 1.3,
                      f0 = 0.6,
                      fc = 0.25,
                      ks = NA,
                      gamma = NA,
                      zeta = .2,
                      ext_mort_prop = 0,
                      R_factor = 4,
                      gear_names = "knife_edge_gear",
                      knife_edge_size = 1000,
                      RDD = "extinctionRDD",
                      egg_size_scaling = FALSE,
                      resource_scaling = FALSE,
                      perfect_scaling = FALSE) {

  params <- newTraitParams(
    no_sp = no_sp,
    min_w_inf = min_w_inf,
    max_w_inf = max_w_inf,
    min_w = min_w,
    max_w = max_w,
    eta = eta,
    min_w_mat = min_w_mat,
    no_w = no_w,
    min_w_pp = min_w_pp,
    w_pp_cutoff = w_pp_cutoff,
    n = n,
    p = p,
    lambda = lambda,
    r_pp = r_pp,
    kappa = kappa,
    alpha = alpha,
    h = h,
    beta = beta,
    sigma = sigma,
    f0 = f0,
    fc = fc,
    ks = ks,
    gamma = gamma,
    ext_mort_prop = ext_mort_prop,
    R_factor = R_factor,
    gear_names = gear_names,
    knife_edge_size = knife_edge_size,
    egg_size_scaling = egg_size_scaling,
    resource_scaling = resource_scaling,
    perfect_scaling = perfect_scaling)

  params@species_params$zeta <- zeta
  params@species_params$lineage <- as.factor(1:no_sp)
  params <- setReproduction(params, RDD = RDD)

  return(params)
}

# this function calculate the distance bewteen neighbouring values of a vector
neighbourDistance <- function(x)
{
  y <- vector("numeric", length = length(x))
  x <- c(0,x)
  for(i in 1:(length(x)-1))  y[i] <- x[i+1] - x[i]
  return(y)
}

#function that paste the different MizerSim objects
finalTouch <- function(saveFolder,params,t_max)
{
  # for now I am not going to be removing extinct species so the last sim will have the right species dimension
  # the time dimension will be the sum of the time dim of the runs
  # for more efficiency I am not going to load all the run at once but one by one
  # the biomass of all species but the resident stays constant for one time step when introducing a new species
  no_w <- length(params@w) # will have to this properly somewhere else
  no_w_pp <- length(params@w_full)
  no_phen <- dim(params@species_params)[1]

  # biomass | n_other is not accounted for at the moment cause I don't know what that is
  # not sure about the effort one yet as I need to do different fisheries scenarios to see how it behaves
  biomass <- array(NA, dim = c(t_max+1,no_phen,no_w), dimnames = list("time" = 1:(t_max+1), "sp" = 1:no_phen, "w" = params@w))
  biomassPP <- array(NA, dim = c(t_max+1,no_w_pp), dimnames = list("time" = 1:(t_max+1), "w" = params@w_full))
  effort <- array(0, dim = c(t_max+1,length(unique(params@gear_params$gear))), dimnames = list("time" = 1:(t_max+1), "gear" = unique(params@gear_params$gear)))

  sim_start = 1
  for(iRun in 1:length(dir(saveFolder)))
  {
    tempRun <- readRDS(paste(saveFolder,"/run",iRun,".rds",sep=""))
    biomass[(sim_start):(sim_start-1+dim(tempRun@n)[1]),1:dim(tempRun@n)[2],] <- tempRun@n
    biomassPP[(sim_start):(sim_start-1+dim(tempRun@n)[1]),] <- tempRun@n_pp
    effort[(sim_start):(sim_start-1+dim(tempRun@n)[1]),] <- tempRun@effort
    sim_start = sim_start-1+dim(tempRun@n)[1]
  }

  # reconstruct the mizer object; the last tempRun loaded contains the right @params

  tempRun@n <- biomass
  tempRun@effort <- effort
  tempRun@n_pp <- biomassPP
  tempRun@n_other <- matrix(0, nrow = (t_max+1),ncol = 1, dimnames = list("time" = 1:(t_max+1), "component"))



  return(tempRun)
}

#' Project
#'
#' @param params ...
#' @param t_max ...
#' @param mutation ...
#' @param saveFolder ...
#' @return ...
#' @export
#' @examples
#' \donttest{
#' saveFolder <- file.path(tempdir(), "simTemp")
#' dir.create(saveFolder)
#' params <- evoParams(no_sp = 5, RDD = "extinctionRDD")
#' sim <- evoProject(params = params, t_max = 300, mutation = 5,
#'                   saveFolder = saveFolder)
#' plot(sim)
#' }
evoProject <- function(params,t_max = 100, dt = 0.1, 
                       mutation = 2, trait = "w_mat",
                       saveFolder = file.path(tempdir(), "simTemp"), effort = 0)
{
  # Initialisation
  SpIdx <- unique(params@species_params$species)
  # check if saveFolder is ok
  if(!dir.exists(saveFolder)) dir.create(saveFolder)
  
  if(is.numeric(mutation))
  {
# users gives a mutation rate per species per year in %, it is converted into a matrix giving which species gets mutant at which time_step
  mutationPerSteps <- mutation#*dt
  t_mutation <- matrix(0,nrow = length(SpIdx), ncol = (t_max/1), dimnames = list("species" = SpIdx, "time" = 1:(t_max/1)))  
  for(iSpecies in SpIdx) # for each species
  {
    for(iTime in 1:(t_max/1)) # for each time step
    {
      if(mutationPerSteps >= sample(seq(0,100,.1), 1))
        t_mutation[iSpecies,iTime] <- 1
    }
  }
  t_event <- apply(t_mutation,2,sum) # t_mutation knows which species mutates and t_event knows which time
  t_event <- which(t_event >=1)
# t_event <- floor(t_event*dt) # need to go back to this dimension as project(that I cannot change) is going to divide by dt
# but I start to get issues
# if(t_event[1] == 0) t_event[1] = 1
#print(t_event)
#print(t_mutation)
  } else if (is.data.frame(mutation))
  {
    # species invastion case
    t_event <- mutation$time
    
  } else (stop("what do you want from me!? I said, numeric only!"))
  
  # use this if want specific number of mutation
  # t_event <- sort(sample(1:(t_max-1),mutation))
  # # corrected_t_event <- t_event+1# need to add 1 per run as they contain a 0 time step (should I get rid of it?)
  # # we know when the mutations will appear, now we need to calculate the different time intervals between these mutations
  t_max_vec <- neighbourDistance(x = c(t_event,t_max))
  cat(sprintf('%i mutants will be added\n', (length(t_max_vec)-1)))
  # print(t_mutation)
  # print(t_max_vec)

  mySim <- project(params, t_max = t_max_vec[1],progress_bar = F, effort = effort)
  if(length(t_max_vec) >1) # if there is at least one mutation planned
  {
  saveRDS(mySim,file= paste(saveFolder,"/run1.rds", sep = ""))
  for(iSim in 2:length(t_max_vec))
  {
    #print(iSim)
    if(is.numeric(mutation))
    {
    ## new mutant param
    # randomly take a resident
    resident <- as.character(which(t_mutation[,t_event[iSim-1]]>0)[1]) # not handling 2 new mutatns at the same time for now
    #print(resident)
    # resident <- as.character(sample(mySim@params@species_params$species, 1))

    # create a new species param

    newSp <- mySim@params@species_params[resident,] # get a copy of the resident
# switch to determine what happens to the new species
    #TODO  need to rewrite the specific cases properly
    switch(trait, # cases with specific names and default if users gives just a parameter df name as trait
           size = {
             # Trait = asymptotic size
             sd = as.numeric(mAmplitude *  resident_params["w_inf"]) # standard deviation
             mutant["w_inf"] <- resident_params["w_inf"] + rnorm(1, 0, sd) # change a bit the asymptotic size
             mutant["w_mat"] <- mutant["w_inf"] * eta # calculate from the new w_inf value
             mutant["z0"] <- z0pre * as.numeric(mutant["w_inf"]) ^ (n - 1) # if I don't put as.numeric I lose the name z0
             #cat(sprintf("Its size mutes slightly.\n"))
           },
           Beta = {
             # Trait = PPMR
             sd = as.numeric(mAmplitude *  resident_params["beta"]) # standard deviation
             mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd) # change a bit the PPMR
             while(mutant$beta < 10) 
             {print("need to reroll beta")
               mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd)
             }
             # calculate the new gamma
             alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)*   
               (pnorm(3 - (lambda - 2) * mutant$sigma) + pnorm(log(mutant$beta)/mutant$sigma + (lambda - 2) * mutant$sigma) - 1)
             
             # mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
             
             cat(sprintf("parent beta:%f, gamma:%f\n",resident_params$beta,resident_params$gamma))
             cat(sprintf("mutant beta:%f, gamma:%f\n",mutant$beta,mutant$gamma))
             # cat(sprintf("Its PPMR is:%f\n",mutant$beta))
           },
           Sigma = {
             # Trait = fedding kernel
             sd = as.numeric(mAmplitude *  resident_params["sigma"]) # standard deviation
             mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd) # change a bit the diet breadth
             while(mutant$sigma < .5 | mutant$sigma > 5) 
             {print("need to reroll sigma")
               mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd)}
             # calculate the new gamma
             alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)*   
               (pnorm(3 - (lambda - 2) * mutant$sigma) + pnorm(log(mutant$beta)/mutant$sigma + (lambda - 2) * mutant$sigma) - 1)
             mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
             #cat(sprintf("Its diet breadth mutes slightly.\n"))
           },
           predation = {
             # PPMR
             sd = as.numeric(mAmplitude *  resident_params["beta"]) # standard deviation
             mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd) # change a bit the PPMR
             while(mutant$beta < 10) 
             {print("need to reroll beta")
               mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd)}
             # feeding kernel
             sd = as.numeric(mAmplitude *  resident_params["sigma"]) # standard deviation
             mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd) # change a bit the diet breadth
             while(mutant$sigma < .5 | mutant$sigma >5) 
             {print("need to reroll sigma")
               mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd)}
             # recalculate gamma if necessary
             alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)*   
               (pnorm(3 - (lambda - 2) * mutant$sigma) + pnorm(log(mutant$beta)/mutant$sigma + (lambda - 2) * mutant$sigma) - 1)
             mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
             cat(sprintf("parent beta:%f, sigma:%f, gamma:%f\n",resident_params$beta,resident_params$sigma,resident_params$gamma))
             cat(sprintf("mutant beta:%f, sigma:%f, gamma:%f\n",mutant$beta,mutant$sigma,mutant$gamma))
           },
           eta = {
             # Trait = eta
             sd = as.numeric(mAmplitude *  resident_params["eta"]) # standard deviation
             mutant["eta"] <- resident_params["eta"] + rnorm(1, 0, sd) # change a bit eta
             if (mutant["eta"] >= 1) mutant["eta"] <- 0.95 # because yes it does happen
             mutant["w_mat"] <- mutant["w_inf"] * mutant["eta"] # update
             #cat(sprintf("Its w_mat is: %g\n",mutant["w_mat"]))
           },
           ed_int = {
             # Trait = ed_int
             sd = as.numeric(mAmplitude *  resident_params["ed_int"])
             mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd)
             while(mutant$ed_int < 2.5) mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd) # ed_int cannot go lower than 2.5
           },
           t_d = {
             # Trait = ed_int
             sd = as.numeric(mAmplitude *  resident_params["t_d"])
             mutant$t_d <- mutant$t_d + rnorm(1, 0, sd)
             cat(sprintf("Its name is %i and its trait value is %g\n", mutant$ecotype,mutant["t_d"]))
           },
           temperature = {
             sd = as.numeric(mAmplitude * resident_params["ed_int"])
             mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd)
             while(mutant$ed_int < 2.5) mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd) # ed_int cannot go lower than 2.5
             sd = as.numeric(mAmplitude * resident_params["t_d"])
             mutant$t_d <- mutant$t_d + rnorm(1, 0, sd)
           },
           all = {
             # Trait = asymptotic size
             sd = as.numeric(mAmplitude *  resident_params["w_inf"]) # standard deviation
             mutant["w_inf"] <- resident_params["w_inf"] + rnorm(1, 0, sd) # change a bit the asymptotic size
             mutant["w_mat"] <- mutant["w_inf"] * eta # calculate from the new w_inf value
             mutant["z0"] <- z0pre * as.numeric(mutant["w_inf"]) ^ (n - 1) # if I don't put as.numeric I lose the name z0
             # Trait = predation
             sd = as.numeric(mAmplitude *  resident_params["beta"]) # standard deviation
             mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd) # change a bit the PPMR
             sd = as.numeric(mAmplitude *  resident_params["sigma"]) # standard deviation
             mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd) # change a bit the diet breadth
             # calculate the new gamma
             alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)*   
               (pnorm(3 - (lambda - 2) * mutant$sigma) + pnorm(log(mutant$beta)/mutant$sigma + (lambda - 2) * mutant$sigma) - 1)
             mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
             #cat(sprintf("Its traits mute slightly.\n"))
           },
           {
             sd = as.numeric(mySim@params@species_params$zeta[as.numeric(resident)] * mySim@params@species_params[trait][as.numeric(resident),1])
             newSp[trait] <- abs(newSp[trait] + rnorm(1, 0, sd))

           })
    
    # set the abundance for all species to start a new project
    lastBiom <- mySim@n[dim(mySim@n)[1],,]
    #n_newSp <- rep(0,dim(mySim@n)[3])
    n_newSp = 0.05 * lastBiom[dimnames(mySim@n)$sp == resident,] # the initial abundance is 5% of the resident pop
    lastBiom[dimnames(mySim@n)$sp ==resident,]= lastBiom[dimnames(mySim@n)$sp == resident,] - 0.05*lastBiom[dimnames(mySim@n)$sp ==resident,] # Witdraw the abundance of the mutant from its parent (we're not talking about eggs here but different ecotype already present)
    
    
} else if (is.data.frame(mutation))
{
  newSp <- mutation[iSim-1,]
  
  lastBiom <- mySim@n[dim(mySim@n)[1],,]
  n_newSp <- rep(0,dim(mySim@n)[3])
  #need to create size spectrum of abundance from one value
  n0_mult = mutation$init_n_multiplier
  a = 0.35
 
    no_w <- length(params@w)
    initial_n <- array(NA, dim = c(1, no_w))
    # N = N0 * Winf^(2*n-q-2+a) * w^(-n-a)
    # Reverse calc n and q from intake_max and search_vol slots (could add get_n function)
    n <- (log(params@intake_max[,1] / params@species_params$h) / log(params@w[1]))[1]
    q <- (log(params@search_vol[,1] / params@species_params$gamma) / log(params@w[1]))[1]
    # Guessing at a suitable n0 value based on kappa - this was figured out using trial and error and should be updated
    if (is.null(n0_mult)) {
      lambda <- 2 + q - n
      kappa <- params@cc_pp[1] / (params@w_full[1]^(-lambda))
      n0_mult <- kappa / 1000
    }
    initial_n <- unlist(tapply(params@w, 1:no_w, function(wx,n0_mult,w_inf,a,n,q)
      n0_mult * w_inf^(2 * n - q - 2 + a) * wx^(-n - a),
      n0_mult = n0_mult, w_inf = mutation$w_inf[iSim-1], a=a, n=n, q=q))
    #set densities at w > w_inf to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_inf) w_inf<wx, w_inf=mutation$w_inf[iSim-1]))] <- 0
    # Also any densities at w < w_min set to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_min)w_min>wx, w_min=mutation$w_min[iSim-1]))] <- 0    
n_newSp <- t(initial_n)
}
    newSp$species <-factor(as.character(max(as.numeric(mySim@params@species_params$species))+1), levels = max(as.numeric(mySim@params@species_params$species))+1) # new species name but lineage stays the same

    init_n <- rbind(lastBiom,n_newSp) # this include the new mutant as last column
    names(dimnames(init_n)) <- c("sp","w")
    # print(newSp)
    rownames(init_n)[length(rownames((init_n)))] <- as.character(newSp$species) # update the name of the mutant accordingly

    params <- addSpecies(params = params, species_params = newSp, init_n= init_n)
    if(t_max_vec[iSim]>0) # happens if mutant appears at the last time step, makes the code crash
    mySim <- project(params, t_max = t_max_vec[iSim],progress_bar = F, effort = effort)
    
    saveRDS(mySim,file= paste(saveFolder,"/run",iSim,".rds", sep = ""))
  }
# return(list(saveFolder,params,t_max))
  sim <- finalTouch(saveFolder = saveFolder, params = params, t_max = t_max)
  unlink(saveFolder,recursive = T)

  return(sim)
  }
  return(mySim) # no mutation, just normal run
}

addSpecies <- function(params, species_params, interaction, defaultInteraction = 1, init_n) {
  # check validity of parameters ----
  # assert_that(is(params, "MizerParams"), # does not work in my R version
  #             is.data.frame(species_params))
  if (any(species_params$species %in% params@species_params$species)) {
    stop("You can not add species that are already there.")
  }
  no_old_sp <- nrow(params@species_params)
  old_sp <- 1:no_old_sp
  no_new_sp <- nrow(species_params)
  new_sp <- 1:no_new_sp + no_old_sp
  no_sp <- no_old_sp + no_new_sp
  if (missing(interaction)) {
    # keep existing interactions between old species and
    # set interactions involving new species to 1
    inter <- matrix(defaultInteraction, nrow = no_sp, ncol = no_sp)
    inter[old_sp, old_sp] <- params@interaction
  } else if (all(dim(interaction) == c(no_new_sp, no_new_sp))) {
    # keep existing interactions between old species,
    # set interactions involving an old and a new species to 1
    # and use supplied matrix for interaction among new species
    inter <- matrix(defaultInteraction, nrow = no_sp, ncol = no_sp)
    inter[old_sp, old_sp] <- params@interaction
    inter[new_sp, new_sp] <- interaction
  } else if (all(dim(interaction) != c(no_sp, no_sp))) {
    stop("interaction matrix has invalid dimensions.")
  } else {
    inter <- interaction
  }

  # combine species params ----

  # Move linecolour and linetype into species_params
  # params@species_params$linetype <-
  #   params@linetype[as.character(params@species_params$species)]
  # params@species_params$linecolour <-
  #   params@linecolour[as.character(params@species_params$species)]

  # Make sure that all columns exist in both data frames
  missing <- setdiff(names(params@species_params), names(species_params))
  species_params[missing] <- NA
  missing <- setdiff(names(species_params), names(params@species_params))
  params@species_params[missing] <- NA

  # add the new species (with parameters described by species_params),
  # to make a larger species_params dataframe.
  combi_species_params <- rbind(params@species_params, species_params,
                                stringsAsFactors = FALSE)
  # new params object ----
  # use dataframe and global settings from params to make a new MizerParams
  # object.
  # TODO going to need to check if params need to be updated with the new sp or not
  p <- newMultispeciesParams(
    combi_species_params,
    interaction = inter,
    min_w = min(params@w),
    max_w = max(params@w),
    min_w_pp = min(params@w_full),
    no_w = length(params@w),
    initial_effort = params@initial_effort,
    RDD = params@rates_funcs$RDD
  )
  # Use the same resource spectrum as params
  p@initial_n_pp <- params@initial_n_pp
  p@cc_pp <- params@cc_pp
  p@rr_pp <- params@rr_pp
  p@resource_dynamics <- params@resource_dynamics
  p@resource_params <- params@resource_params
  # Preserve comment
  comment(p) <- comment(params)

  # initial solution ----
  #TODO set initial_N of new sp as a parameters (like 5% of parent and such)
  # p@initial_n[old_sp, ] <- params@initial_n
  # p@initial_n[new_sp, ] <- init_n
  p@initial_n <- init_n
  p@A[old_sp] <- params@A
  # Use the same psi and mu_b as before for old species
  p@psi[old_sp, ] <- params@psi
  p@sc <- params@sc
  p@mu_b[old_sp, ] <- params@mu_b
  # we assume same background death for all species
  p@mu_b[new_sp, ] <- rep(params@mu_b[1, ], each = no_new_sp)

  # Turn off self-interaction among the new species, so we can determine the
  # growth rates, and death rates induced upon them by the pre-existing species
  # p@interaction[new_sp, new_sp] <- 0
  # mumu <- getMort(p)
  # gg <- getEGrowth(p)
  #
  # # Compute solution for new species
  # for (i in new_sp) {
  #   g <- gg[i, ]
  #   mu <- mumu[i, ]
  #   w_inf_idx <- sum(p@w < p@species_params$w_inf[i])
  #   idx <- p@w_min_idx[i]:(w_inf_idx - 1)
  #   if (any(g[idx] == 0)) {
  #     stop("Can not compute steady state due to zero growth rates for ",
  #          p@species_params$species[i])
  #   }
  #   p@initial_n[i, ] <- 0
  #   p@initial_n[i, p@w_min_idx[i]:w_inf_idx] <-
  #     c(1, cumprod(g[idx] / ((g + mu * p@dw)[idx + 1])))
  #
  #   # set low abundance ----
  #   # Normalise solution so that at its maximum it lies at 1/100 of the
  #   # Sheldon spectrum.
  #   # We look at the maximum of abundance times w^lambda
  #   # because that is always an increasing function at small size.
  #   idx <- which.max(p@initial_n[i, ] * p@w^p@resource_params$lambda)
  #   p@initial_n[i, ] <- p@initial_n[i, ] *
  #     p@resource_params$kappa * p@w[idx]^(-p@resource_params$lambda) / p@initial_n[i, idx] / 100
  #   p@A[i] <- sum(p@initial_n[i, ] * p@w * p@dw * p@maturity[i, ])
  # }
  #
  # if (any(is.infinite(p@initial_n))) {
  #   stop("Candidate steady state holds infinities.")
  # }
  # if (any(is.na(p@initial_n) | is.nan(p@initial_n))) {
  #   stop("Candidate steady state holds non-numeric values.")
  # }
  #
  # # Turn self interaction back on
  # p@interaction[new_sp, new_sp] <- inter[new_sp, new_sp]
  #
  # # Retune reproductive efficiencies of new species
  # p <- retune_erepro(p, p@species_params$species[new_sp])
  #
  return(p)
}

#' Stock recruitment relationship enabling extinction of species
#'
#' when a species reaches an abundance threshold its recruitment gets disabled, species is later removed
#' For now each species acts as its own
#' TODO at the moment I am using rdi to apply the threhold but it should be done on n directly
#'
#' @param rdi ...
#' @param species_params ...
#' @param ... ...
#' @export
extinctionRDD <- function(rdi, species_params, ...) {
  if (!("R_max" %in% names(species_params))) {
    stop("The R_max column is missing in species_params.")
  }
  rdiNormal = vector(mode = "numeric", length = length(rdi))
  names(rdi) <- species_params$lineage
  for (iSpecies in sort(unique(species_params$lineage))) # makes a vector of value from 0 to 1 showing the abundance proportion of each phenotypes within each species
  {
    rdiSp = rdi # save to manip
    rdiSp[which(names(rdi) != iSpecies)] = 0 # make everything but the targeted species to go 0 to have correct normalisation
    
    for (i in 1:length(rdiSp))
      # in case of NA
      if (is.na(rdiSp[i]) == TRUE)
        rdiSp[i] = 1e-30
    
    if (sum(rdiSp) != 0)
      rdiNormal = rdiNormal + rdiSp / sum(rdiSp)
  }
  r_maxN = species_params$R_max * rdiNormal # apply the scaling to rmax
  
  for (i in 1:length(r_maxN)) # do not want to divide by 0 so replacing the 0 value by the original rmax (does not matter as if there was a 0 value, it means that the rmax is going to be multiplied by 0)
    if (r_maxN[i] == 0)
      r_maxN[i] = 1
  
  rdd <- rdi / (1 + rdi/r_maxN)
  if(sum(which(rdd <= 1e-30))) rdd[which(rdd <= 1e-30)] <- 0 # if any of the rdd is under threshold, set it to 0
  return(rdd)
}




plotDynamics <- function(object, time_range = c(min(as.numeric(dimnames(object@n)$time)),max(as.numeric(dimnames(object@n)$time))), 
                         phenotype = TRUE, species = NULL, trait = NULL, SpIdx = NULL, print_it = T, returnData = F, save_it = F, 
                         nameSave = "Biomass.png", ylimit = c(NA,NA)){
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) # colorful gradient
  min_value <- 1e-30
  # get the phenotype biomass through time (need at least 2 time steps for now)
  biomass <- getBiomass(object)
  time_elements <- get_time_elements(object,time_range)
  biomass <- biomass[time_elements,]
  
  # find out SpIdx if not given / get rid of fully extinct species
  if (is.null(SpIdx))
    for (i in unique(object@params@species_params$lineage))
      if (sum(biomass[, i]) != 0) 
        SpIdx = c(SpIdx, i)
  
  # sum the biomass across species (no more phenotype identity)
  biomassSp = NULL
  biomassTemp = biomass
  colnames(biomassTemp) = object@params@species_params$lineage
  for (i in SpIdx)
  {
    biomassPhen = biomassTemp[,which(colnames(biomassTemp) == i)]
    if(!is.null(dim(biomassPhen))) biomassPhen = apply(biomassPhen,1,sum,na.rm =T)
    biomassSp = cbind(biomassSp,biomassPhen)
  }
  colnames(biomassSp) = SpIdx
  
  # Check if phenotypes are extinct in biomass as well
  spSub <- object@params@species_params$species[object@params@species_params$lineage %in% SpIdx]
  biomass <- biomass[,as.numeric(dimnames(biomass)$sp) %in% spSub]
  
  plotBiom <- function(x)
  {
    Biom <- melt(x) # melt for ggplot
    colnames(Biom) = c("time","phen","value")
    # create a species column
    Biom$sp = sapply(Biom$phen, function(x) object@params@species_params$lineage[x])
    return(Biom)
  }
  BiomSp <- plotBiom(biomassSp)
  BiomSp <- BiomSp[BiomSp$value >= min_value,]
  
  
  if (phenotype) # are we displaying the phenotypes or just species total?
  {
    
    BiomPhen <- plotBiom(biomass)
    BiomPhen <- plotBiom(biomass)
    if(!is.null(trait))
      BiomPhen$trait <- rep(trait, each = length(unique(BiomPhen$time)))
    BiomPhen <- BiomPhen[BiomPhen$value >= min_value,]
    
    p <- ggplot(BiomSp) +
      geom_line(aes(x = time, y = value, colour = as.factor(sp), group = sp), size = 1.2) +
      geom_line(data = BiomPhen, aes(x = time, y = value, colour = as.factor(sp), group = phen), alpha = 0.2) +
      scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(seq(-30,4,2)))) +
      scale_x_continuous(name = "Time in years") +
      labs(color='Species') +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.key = element_rect(fill = "white"))+
      scale_colour_manual(values=cbPalette)+ # colorblind
      ggtitle("Community biomass") 
    
    
    if (!is.null(species)) # Are we looking at one species in particular or all of them? Note: need to select only one species if we want to look at traits
    {
      # BiomPhen <- BiomPhen[BiomPhen$sp == species, ]
      # BiomSp <- BiomSp[BiomSp$sp == species, ]
      plotTitle <- paste("Species",species)
      
      p <- ggplot(filter(BiomPhen, sp == species)) +
        geom_line(aes(x = time, y = value, group = phen), alpha = .75) +
        scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(seq(-30,4,2)))) +
        scale_x_continuous(name = "Time in years") +
        # scale_colour_gradientn(colours=jet.colors(9), limits = c(NA,NA))+
        theme(panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),
              legend.key = element_rect(fill = "white"))+
        ggtitle(plotTitle) 
      
      if(!is.null(trait))
      {
        plotTitle <- paste("Trait of species",species)
        
        p <- ggplot(filter(BiomPhen, sp == species)) +
          geom_line(aes(x = time, y = value, colour = trait, group = trait), alpha = 1) +
          scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(seq(-30,4,2)))) +
          scale_x_continuous(name = "Time in years") +
          scale_colour_gradientn(colours=jet.colors(9), limits = c(NA,NA))+
          theme(panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),
                legend.key = element_rect(fill = "white"))+
          ggtitle(plotTitle) 
      }
    }
    
  } else {
    # just total biomass per species here
    p <- ggplot(BiomSp) +
      geom_line(aes(x = time, y = value, colour = as.factor(sp), group = sp), size = 1.2) +
      scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(-30:4))) +
      scale_x_continuous(name = "Time in years") +
      labs(color='Species') +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.key = element_rect(fill = "white"))+
      scale_colour_manual(values=cbPalette)+ # colorblind
      ggtitle("Community biomass") 
  }
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(list(BiomSp,BiomPhen)) else if(print_it) return(p)
}

plotSS <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), min_w =min(object@params@w)/100, ylim = c(NA,NA),
                   biomass = TRUE, print_it = TRUE, species = TRUE, community = FALSE, save_it = FALSE, nameSave = "SizeSpectrum.png", returnData = FALSE, ...){
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  # min_w = 0.001
  time_elements <- get_time_elements(object,time_range)
  spec_n <- apply(object@n[time_elements,,,drop=FALSE],c(2,3), mean)
  pkt_n <- apply(object@n_pp[time_elements,,drop=FALSE],2,mean)
  alg_n <- apply(object@n_aa[time_elements,,drop=FALSE],2,mean)
  ben_n <- apply(object@n_bb[time_elements,,drop=FALSE],2,mean)
  
  y_axis_name = "Abundance"
  if (biomass){
    spec_n <- sweep(spec_n,2,object@params@w,"*")
    pkt_n <- pkt_n * object@params@w_full
    alg_n <- alg_n * object@params@w_full
    ben_n <- ben_n * object@params@w_full
    y_axis_name = "Biomass"
  }
  # Make data.frame for plot
  plot_datSP <- data.frame(value = c(spec_n), Species = dimnames(spec_n)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)), bloodline = object@params@species_params$species)
  plot_datPkt <- data.frame(value = c(pkt_n), Species = "Phytoplankton", w = object@params@w_full)
  plot_datAlg <- data.frame(value = c(alg_n), Species = "Algae", w = object@params@w_full)
  plot_datBen <- data.frame(value = c(ben_n), Species = "Benthos", w = object@params@w_full)
  
  if(community) plot_datSP <- data.frame(value = apply(spec_n, 2, sum), w = object@params@w)
  
  else if (species)
  {
    dimnames(spec_n)$species = object@params@species_params$species
    SpIdx = unique(object@params@species_params$species)
    spec_sp = matrix(data = NA, ncol = dim(spec_n)[2], nrow = length(SpIdx), dimnames = list(as.character(SpIdx),dimnames(spec_n)$size))
    names(dimnames(spec_sp))=list("species","size")
    
    for (i in 1:dim(spec_sp)[1])
    {
      temp = spec_n # save to manip
      temp[which(rownames(spec_n) != i), ] = 0 # make everything but the targeted species to go 0 to have correct normalisation
      temp = apply(temp, 2, sum)
      spec_sp[i, ] = temp
    }
    plot_datSP <- data.frame(value = c(spec_sp), Species = dimnames(spec_sp)[[1]], w = rep(object@params@w, each=length(SpIdx)))
  }
  # lop off 0s in background and apply min_w
  plot_datSP <- plot_datSP[(plot_datSP$value > 0) & (plot_datSP$w >= min_w),]
  plot_datPkt <- plot_datPkt[(plot_datPkt$value > 0) & (plot_datPkt$w >= min_w),]
  plot_datAlg <- plot_datAlg[(plot_datAlg$value > 0) & (plot_datAlg$w >= min_w),]
  plot_datBen <- plot_datBen[(plot_datBen$value > 0) & (plot_datBen$w >= min_w),]
  #getPalette = colorRampPalette(brewer.pal(9, "Set1"))# increase the number of colors used
  
  if(community)
  {
    p <- ggplot(plot_datSP) + 
      geom_line(aes(x=w, y = value)) + 
      geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species),alpha = 0.5, color = "blue", size = 1.5) +
      geom_line(data = plot_datAlg, aes(x = w, y = value, group = Species),alpha = 0.5, color = "green", size = 1.5) +
      geom_line(data = plot_datBen, aes(x = w, y = value, group = Species),alpha = 0.5, color = "yellow", size = 1.5) +
      scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5)))+
      scale_y_continuous(name = "Abundance density in individuals.m^-3", limits = ylim, trans = "log10") +
      # labs(color='Species') +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.key = element_rect(fill = "white"))+
      scale_colour_manual(values=cbPalette)+ # colorblind
      ggtitle("Community spectrum") 
  }
  else if (species)
  {
    p <- ggplot(plot_datSP) + 
      geom_line(aes(x=w, y = value, colour = as.factor(Species), group = Species)) + 
      geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species),alpha = 0.5, color = "blue", size = 1.5) +
      geom_line(data = plot_datAlg, aes(x = w, y = value, group = Species),alpha = 0.5, color = "green", size = 1.5) +
      geom_line(data = plot_datBen, aes(x = w, y = value, group = Species),alpha = 0.5, color = "yellow", size = 1.5) +
      scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5)))+
      scale_y_continuous(name = "Abundance density in individuals.m^-3", limits = ylim, trans = "log10") +
      labs(color='Species') +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.key = element_rect(fill = "white"))+
      scale_colour_manual(values=cbPalette)+ # colorblind
      ggtitle("Size spectrum")
  }
  
  else
  {
    p <- ggplot(plot_datSP) + 
      geom_line(aes(x=w, y = value, colour = as.factor(bloodline), group = Species)) + 
      geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species),alpha = 0.5, color = "blue", size = 1.5) +
      geom_line(data = plot_datAlg, aes(x = w, y = value, group = Species),alpha = 0.5, color = "green", size = 1.5) +
      geom_line(data = plot_datBen, aes(x = w, y = value, group = Species),alpha = 0.5, color = "yellow", size = 1.5) +
      scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5)))+
      scale_y_continuous(name = "Abundance density in individuals.m^-3", limits = ylim, trans = "log10") +
      labs(color='Species') +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.key = element_rect(fill = "white"))+
      scale_colour_manual(values=cbPalette)+ # colorblind
      ggtitle("Size spectrum")
    
  }
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_datSP) else if(print_it) return(p)
}


#need to update below

plotFood <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, throughTime = F, start = 1000, every = 1000, 
                     print_it = T, returnData = F, save_it =F, nameSave = "Feeding.png"){
  
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  if (throughTime)
  {
    time_range = seq(start,max(as.numeric(dimnames(object@n)$time)),every)
    time_range = c(time_range,max(as.numeric(dimnames(object@n)$time))) # so it counts the last time step which is probably not even
    time_range = unique(time_range)
    feeding = array(data = NA, dim = c(length(unique(object@params@species_params$species)),100,length(time_range)),  
                    dimnames = list(as.character(unique(object@params@species_params$species)),object@params@w,time_range)) 
    Critfeeding = matrix(data=NA, nrow = length(time_range), ncol= 100, dimnames = list(time_range,object@params@w))
    for (i in time_range)
    {
      
      feed_time <- getFeedingLevel(object=object, time_range=i, drop=FALSE)#, ...) # get the feeding time
      feed <- apply(feed_time, c(2,3), mean) # average on the time frame
      
      Cfeed_time <- getCriticalFeedingLevel(object=object, time_range=i, drop=FALSE)#, ...) # get the critical feeding level
      Critfeed <- apply(Cfeed_time, c(2,3), mean) # average on the time frame
      Critfeed <- Critfeed[1,] # all rows the same
      
      dimnames(feed)$sp = object@params@species_params$species
      SpIdx = unique(object@params@species_params$species) # get the species names
      feed_sp = matrix(data = NA, ncol = dim(feed)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(feed)$w)) # prepare the new object
      names(dimnames(feed_sp))=list("species","size")
      
      for (j in SpIdx)
      {
        temp = feed # save to manip
        temp[which(rownames(feed) != j), ] = 0 # keep the ecotypes from the species only
        temp = apply(temp, 2, sum)
        temp = temp / length(which(rownames(feed)==j)) # do the mean (in 2 steps)
        feed_sp[which(rownames(feed_sp)==j), ] = temp
      }
      feeding[,,which(dimnames(feeding)[[3]] == i)] = feed_sp
      Critfeeding[which(dimnames(Critfeeding)[[1]] == i),] = Critfeed
    }
    a <- c(object@params@species_params$w_inf[SpIdx]) # to get vline of different col, need to create a data frame
    vlines <- data.frame(xint = a,grp = SpIdx)
    
    plot_dat = melt(feeding)
    colnames(plot_dat) = c("species","size","time","value")
    plot_crit = melt(Critfeeding)
    colnames(plot_crit) = c("time","size","value")
    p <- ggplot(plot_dat) + 
      geom_line(aes(x=size, y = value, colour = as.factor(species))) + 
      geom_line(data = plot_crit, aes(x = size, y = value), linetype = "dashed") +
      scale_x_log10(name = "Size") + 
      scale_y_continuous(name = "Feeding Level", lim=c(0,1))+
      geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
      facet_grid(time ~ .)+
      scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"))+
      
      ggtitle("Feeding level through time")
    
    if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
    
    if (returnData) return(list(plot_dat,plot_crit)) else if(print_it) return(p)
    
  }
  
  feed_time <- getFeedingLevel(object=object, time_range=time_range, drop=FALSE) #, ...) # get the feeding time
  feed <- apply(feed_time, c(2,3), mean) # average on the time frame
  
  Cfeed_time <- getCriticalFeedingLevel(object=object, time_range=time_range, drop=FALSE)#, ...) # get the critical feeding level
  Critfeed <- apply(Cfeed_time, c(2,3), mean) # average on the time frame
  Critfeed <- Critfeed[1,] # all rows the same
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(feed)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    feed_sp = matrix(data = NA, ncol = dim(feed)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(feed)$w)) # prepare the new object
    names(dimnames(feed_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = feed # save to manip
      temp[which(rownames(feed) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(feed)==i)) # do the mean (in 2 steps)
      feed_sp[which(rownames(feed_sp)==i), ] = temp
    }
    feed = feed_sp
  }
  
  a <- c(object@params@species_params$w_inf[SpIdx]) # to get vline of different col, need to create a data frame
  vlines <- data.frame(xint = a,grp = SpIdx)
  
  plot_dat <- data.frame(value = c(feed), species = dimnames(feed)[[1]], size = rep(object@params@w, each=length(dimnames(feed)[[1]])))
  
  name = paste("Feeding level at time",time_range,sep=" ")
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=size, y = value, colour = as.factor(species))) + 
    geom_hline(yintercept = Critfeed[1], linetype = "dashed", color = "red") +
    geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
    scale_x_log10(name = "Size", breaks = c(1 %o% 10^(-3:5)))  + 
    scale_y_continuous(name = "Feeding Level", lim=c(0,1))+
    scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
    theme(panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),
          legend.key = element_rect(fill = "white"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave,width = 18, height = 18,units = "cm" )
  
  if (returnData) return(list(plot_dat,Critfeed)) else if(print_it) return(p)
}

plotGrowth <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, print_it = T, returnData = F, save_it = F, ylim = c(NA,NA),
                       nameSave = "Growth.png",...){
  
  time_elements <- get_time_elements(object,time_range)
  growth_time <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    growth <- getEGrowth(object@params, n=n, n_pp = object@n_pp[x,])#,
                         #n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                         #intakeScalar = object@intTempScalar[,,x], metScalar = object@metTempScalar[,,x])
    return(growth)})
  
  #growth <- apply(growth_time, c(2,3), mean) # use this when I will have time_range on more than one time
  growth = growth_time
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(growth)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    growth_sp = matrix(data = NA, ncol = dim(growth)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(growth)$w)) # prepare the new object
    names(dimnames(growth_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = growth # save to manip
      temp[which(rownames(growth) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(growth)==i)) # do the mean (in 2 steps)
      growth_sp[which(rownames(growth_sp)==i), ] = temp
    }
    growth = growth_sp
  }
  
  name = paste("Growth level at time",time_range,sep=" ")
  plot_dat <- data.frame(value = c(growth), Species = dimnames(growth)[[1]], w = rep(object@params@w, each=length(dimnames(growth)[[1]])))
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-3:5))) + 
    scale_y_continuous(name = "instantaneous growth", trans ="log10", limits = ylim)+
    theme(legend.title=element_blank(),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}

plotStarvation <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, print_it = T, returnData = F, save_it =F, 
                           nameSave = "Starvation.png"){
  
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  # death_time <- getSmort(object=object, time_range=what_time, drop=FALSE)
  
  
  time_elements <- get_time_elements(object,time_range)
  death_time <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    starv <- getSMort(object@params, n=n, n_pp = object@n_pp[x,])#,n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                      #intakeScalar = object@intTempScalar[,,x], metScalar = object@metTempScalar[,,x])
    return(starv)})
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(death_time)[[1]] = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    death_sp = matrix(data = NA, ncol = dim(death_time)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(death_time)$w)) # prepare the new object
    names(dimnames(death_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = death_time # save to manip
      temp[which(rownames(death_time) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(death_time)==i)) # do the mean (in 2 steps)
      death_sp[which(rownames(death_sp)==i), ] = temp
    }
    death = death_sp
  }
  
  # a <- c(object@params@species_params$w_inf[1:9]) # to get vline of different col, need to create a data frame
  # vlines <- data.frame(xint = a,grp = c(1:9))
  
  plot_dat <- data.frame(value = c(death), species = dimnames(death)[[1]], size = rep(object@params@w, each=length(dimnames(death)[[1]])))
  
  name = paste("Starvation at time",time_range,sep=" ")
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=size, y = value, colour = as.factor(species))) + 
    #geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
    scale_x_log10(name = "Size", breaks = c(1 %o% 10^(-3:5)))  + 
    scale_y_continuous(name = "Instantaneous starvation mortality")+
    scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
    theme(panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),
          legend.key = element_rect(fill = "white"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}

plotScythe <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),print_it = TRUE, returnData = F, comments = T){
  
  # effort can be in 2 forms
  
  if(is.matrix(object@effort)) effort = object@effort[time_range,]
  else effort = object@effort[time_range]
  
  z <- getZ(object = object@params, n = object@n[time_range,,], n_pp = object@n_pp[time_range,],n_aa = object@n_aa[time_range,],n_bb = object@n_bb[time_range,],
            effort = effort, 
            intakeScalar = object@intTempScalar[,,time_range], metScalar = object@metTempScalar[,,time_range], morScalar = object@morTempScalar[,,time_range])
  dimnames(z)$prey = object@params@species_params$species
  #SpIdx = sort(unique(object@params@species_params$species)) # get the species names
  
  # need to get rid of the extinct species at that time in SpIdx
  a <- apply(object@n[time_range,,],1,sum)
  
  names(a) <- sapply(names(a), function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
  
  d <- rowsum(a, group = names(a))
  
  if (sum(d[,1] == 0)) 
  {
    
    d <- d[-which(d[,1] == 0),]
    SpIdx <- as.numeric(names(d))
    
    
  } else {SpIdx <- as.numeric(rownames(d))}
  
  z_sp = matrix(data = NA, ncol = dim(z)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(z)$w_prey)) # prepare the new object
  names(dimnames(z_sp))=list("prey","w_prey")
  
  for (i in SpIdx)
  {
    temp = z # save to manip
    temp[which(rownames(z) != i), ] = 0 # keep the ecotypes from the species only
    temp = apply(temp, 2, sum)
    temp = temp / length(which(rownames(z)==i)) # do the mean (in 2 steps)
    z_sp[which(rownames(z_sp)==i), ] = temp
  }
  z = z_sp
  
  name = paste("Total Mortality at time",time_range,sep=" ")
  
  plot_dat <- data.frame(value = c(z), Species = dimnames(z)[[1]], w = rep(object@params@w, each=length(dimnames(z)[[1]])))
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-3:5))) + 
    scale_y_continuous(name = "Mortality", lim=c(0,max(plot_dat$value))) +
    theme(legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(name)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}

plotSpawn <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, print_it = T, returnData = F, save_it = F, 
                      nameSave = "Spawn.png",...){
  
  time_elements <- get_time_elements(object,time_range)
  spawn_time <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    spawn <- getESpawning(object@params, n=n, n_pp = object@n_pp[x,],n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                          intakeScalar = object@intTempScalar[,,x], metScalar = object@metTempScalar[,,x])
    return(spawn)})
  
  #spawn <- apply(spawn_time, c(2,3), mean) # use this when I will have time_range on more than one time
  spawn = spawn_time
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(spawn)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    spawn_sp = matrix(data = NA, ncol = dim(spawn)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(spawn)$w)) # prepare the new object
    names(dimnames(spawn_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = spawn # save to manip
      temp[which(rownames(spawn) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(spawn)==i)) # do the mean (in 2 steps)
      spawn_sp[which(rownames(spawn_sp)==i), ] = temp
    }
    spawn = spawn_sp
  }
  
  a <- c(object@params@species_params$w_inf[1:9]) # to get vline of different col, need to create a data frame
  vlines <- data.frame(xint = a,grp = c(1:9))
  
  name = paste("Spawn level at time",time_range,sep=" ")
  plot_dat <- data.frame(value = c(spawn), Species = dimnames(spawn)[[1]], w = rep(object@params@w, each=length(dimnames(spawn)[[1]])))
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-3:5))) + 
    scale_y_continuous(name = "Energy allocated to spawning", trans ="log10")+
    geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
    theme(legend.title=element_blank(),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}

plotPredRate <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, ylim = c(NA,NA),
                         print_it = T, returnData = F, save_it = F, nameSave = "PredRate.png",...){
  
  time_range = max(as.numeric(dimnames(object@n)$time))
  time_elements <- get_time_elements(object,time_range)
  spawn_time <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    spawn <- getPredRate(object@params, n=n, n_pp = object@n_pp[x,],n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                         intakeScalar = object@intTempScalar[,,x])
    return(spawn)})
  
  #spawn <- apply(spawn_time, c(2,3), mean) # use this when I will have time_range on more than one time
  spawn = spawn_time
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(spawn)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    spawn_sp = matrix(data = NA, ncol = dim(spawn)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(spawn)$w)) # prepare the new object
    names(dimnames(spawn_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = spawn # save to manip
      temp[which(rownames(spawn) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(spawn)==i)) # do the mean (in 2 steps)
      spawn_sp[which(rownames(spawn_sp)==i), ] = temp
    }
    spawn = spawn_sp
  }
  
  a <- c(object@params@species_params$w_inf[1:9]) # to get vline of different col, need to create a data frame
  vlines <- data.frame(xint = a,grp = c(1:9))
  
  name = paste("Predation rate at time",time_range,sep=" ")
  plot_dat <- data.frame(value = c(spawn), Species = dimnames(spawn)[[1]], w = rep(object@params@w_full, each=length(dimnames(spawn)[[1]])))
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-7:5))) + 
    scale_y_continuous(name = "Potential death rate from predator", trans ="log10", limits = ylim)+
    geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
    theme(legend.title=element_blank(),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}
