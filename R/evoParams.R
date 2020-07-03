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
evoProject <- function(params,t_max = 100, mutation = 2, saveFolder = file.path(tempdir(), "simTemp"), effort = 0)
{
  # check if saveFolder is ok
  if(!dir.exists(saveFolder)) dir.create(saveFolder)
  
    t_event <- sort(sample(1:(t_max-1),mutation))
  # corrected_t_event <- t_event+1# need to add 1 per run as they contain a 0 time step (should I get rid of it?)
  # we know when the mutations will appear, now we need to calculate the different time intervals bewteen these mutations
  t_max_vec <- neighbourDistance(x = c(t_event,t_max))


  mySim <- project(params, t_max = t_max_vec[1],progress_bar = F, effort = effort)
  if(length(t_max_vec) >1) # if there is at least one mutation planned
  {
  saveRDS(mySim,file= paste(saveFolder,"/run1.rds", sep = ""))
  for(iSim in 2:length(t_max_vec))
  {
    ## new mutant param
    # randomly take a resident
    resident <- as.character(sample(mySim@params@species_params$species, 1))

    # create a new species param
    newSp <- mySim@params@species_params[resident,] # get a copy of the resident
    # Trait = maturation size
    sd = as.numeric(mySim@params@species_params$zeta[as.numeric(resident)] * mySim@params@species_params$w_mat[as.numeric(resident)])
    newSp$w_mat <- abs(newSp$w_mat + rnorm(1, 0, sd))

    newSp$species <-factor(as.character(max(as.numeric(mySim@params@species_params$species))+1), levels = max(as.numeric(mySim@params@species_params$species))+1) # new species name but lineage stays the same

    # Decide who the parent is

    # set the abundance for all species to start a new project
    lastBiom <- mySim@n[dim(mySim@n)[1],,]
    n_newSp <- rep(0,dim(mySim@n)[3])
    n_newSp = 0.05 * lastBiom[dimnames(mySim@n)$sp ==resident,] # the initial abundance is 5% of the resident pop
    lastBiom[dimnames(mySim@n)$sp ==resident,]= lastBiom[dimnames(mySim@n)$sp == resident,] - 0.05*lastBiom[dimnames(mySim@n)$sp ==resident,] # Witdraw the abundance of the mutant from its parent (we're not talking about eggs here but different ecotype already present)


    init_n <- rbind(lastBiom,n_newSp) # this include the new mutant as last column
    names(dimnames(init_n)) <- c("sp","w")
    rownames(init_n)[which(rownames(init_n) == "n_newSp")] <- as.character(newSp$species) # update the name of the mutant accordingly

    params <- addSpecies(params = params, species_params = newSp, init_n= init_n)
    mySim <- project(params, t_max = t_max_vec[iSim],progress_bar = F, effort = effort)
    saveRDS(mySim,file= paste(saveFolder,"/run",iSim,".rds", sep = ""))
  }

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
#' TODO make rmax apply at the species level and not inidividual phenotypes
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
  rdd <- rdi / (1 + rdi/species_params$R_max)
  if(sum(which(rdd <= 1e-30))) rdd[which(rdd <= 1e-30)] <- 0 # if any of the rdd is under threshold, set it to 0
  return(rdd)
}
