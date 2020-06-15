#' Set up parameters for an evolutionary model
#' this will either set up new params as usual or convert Mizer params to MizerEvo params if given a mizer.params object (TODO)
#'
#' This functions creates a `MizerParams` object describing a trait-based
#' model. This is a simplification of the general size-based model used in
#' `mizer` in which the species-specific parameters are the same for all
#' species, except for the asymptotic size, which is considered the most
#' important trait characterizing a species. Other parameters are related to the
#' asymptotic size. For example, the size at maturity is given by \code{w_inf *
#' eta}, where `eta` is the same for all species. For the trait-based model
#' the number of species is not important. For applications of the trait-based
#' model see Andersen & Pedersen (2010). See the `mizer` website for more
#' details and examples of the trait-based model.
#'
#' The function has many arguments, all of which have default values. Of
#' particular interest to the user are the number of species in the model and
#' the minimum and maximum asymptotic sizes.
#'
#' The characteristic weights of the smallest species are defined by
#' `min_w` (egg size), `min_w_mat` (maturity size) and
#' `min_w_inf` (asymptotic size). The asymptotic sizes of
#' the `no_sp` species
#' are logarithmically evenly spaced, ranging from `min_w_inf` to
#' `max_w_inf`.
#' Similarly the maturity sizes of the species are logarithmically evenly
#' spaced, so that the ratio `eta` between maturity size and asymptotic
#' size is the same for all species. If \code{egg_size_scaling = TRUE} then also
#' the ratio between asymptotic size and egg size is the same for all species.
#' Otherwise all species have the same egg size.
#'
#' In addition to setting up the parameters, this function also sets up an
#' initial condition that is close to steady state.
#'
#' Although the trait based model's steady state is often stable without
#' imposing additional density-dependence, the function can set a Beverton-Holt
#' type density-dependence that imposes a maximum for the reproduction rate that
#' is a multiple of the reproduction rate at steady state. That multiple is set
#' by the argument `R_factor`.
#'
#' The search rate coefficient `gamma` is calculated using the expected
#' feeding level, `f0`.
#'
#' The option of including fishing is given, but the steady state may lose its
#' natural stability if too much fishing is included. In such a case the user
#' may wish to include stabilising effects (like `R_factor`) to ensure the
#' steady state is stable. Fishing selectivity is modelled as a knife-edge
#' function with one parameter, `knife_edge_size`, which is the size at
#' which species are selected. Each species can either be fished by the same
#' gear (`knife_edge_size` has a length of 1) or by a different gear (the
#' length of `knife_edge_size` has the same length as the number of species
#' and the order of selectivity size is that of the asymptotic size).
#'
#' The resulting `MizerParams` object can be projected forward using
#' \code{project()} like any other `MizerParams` object. When projecting
#' the model it may be necessary to reduce `dt` below 0.1 to avoid any
#' instabilities with the solver. You can check this by plotting the biomass or
#' abundance through time after the projection.
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
                           RDD = "BevertonHoltRDD",
                           egg_size_scaling = FALSE,
                           resource_scaling = FALSE,
                           perfect_scaling = FALSE) {

  ## Check validity of parameters ----
  # assert_that(is.logical(egg_size_scaling),
  #             is.logical(resource_scaling),
  #             is.logical(perfect_scaling))
  if (ext_mort_prop >= 1 || ext_mort_prop < 0) {
    stop("ext_mort_prop must be a number between 0 and 1",
         " because it should be the proportion of the total mortality",
         " coming from sources other than predation.")
  }
  if (R_factor <= 1) {
    message("R_factor needs to be larger than 1. Setting R_factor = 1.01")
    R_factor <- 1.01
  }
  no_w <- round(no_w)
  if (no_w < log10(max_w_inf/min_w)*5) {
    no_w <- round(log10(max_w_inf / min_w) * 5 + 1)
    message(paste("Increased no_w to", no_w, "so that there are 5 bins ",
                  "for an interval from w and 10w."))
  }
  if (no_w > 10000) {
    message("Running a simulation with ", no_w,
            " size bins is going to be very slow.")
  }
  if (min_w <= 0) {
    stop("The smallest egg size min_w must be greater than zero.")
  }
  if (min_w_inf >= max_w_inf) {
    stop("The asymptotic size of the smallest species min_w_inf must be ",
         "smaller than the asymptotic size of the largest species max_w_inf")
  }
  if (min_w >= min_w_mat) {
    stop("The egg size of the smallest species min_w must be smaller than ",
         "its maturity size min_w_mat")
  }
  if (min_w_mat >= min_w_inf) {
    stop("The maturity size of the smallest species min_w_mat must be ",
         "smaller than its maximum size min_w_inf")
  }
  no_sp <- as.integer(no_sp)
  if (no_sp < 2) {
    stop("The number of species must be at least 2.")
  }
  if (!all(c(n, r_pp, lambda, kappa, alpha, h, beta, sigma, f0) > 0)) {
    stop("The parameters n, lambda, r_pp, kappa, alpha, h, beta, sigma ",
         "and f0, if supplied, need to be positive.")
  }
  if (!is.na(fc) && (fc < 0 || fc > f0)) {
    stop("The critical feeding level must lie between 0 and f0")
  }
  # Check gears
  if (length(knife_edge_size) > no_sp) {
    stop("knife_edge_size needs to be no longer than the number of species in the model")
  }
  if ((length(knife_edge_size) > 1) & (length(knife_edge_size) != no_sp)) {
    warning("Length of knife_edge_size is less than number of species so gear information is being recycled. Is this what you want?")
  }
  if ((length(gear_names) != 1) & (length(gear_names) != no_sp)) {
    stop("Length of gear_names argument must equal the number of species.")
  }

  if (perfect_scaling) {
    egg_size_scaling <- TRUE
    resource_scaling <- TRUE
    w_pp_cutoff <- Inf
    p <- n
  }

  ## Set grid points and characteristic sizes ----
  # in such a way that the sizes all line up with the grid and the species are
  # all equally spaced.

  # Divide the range from min_w to max_w into (no_w - 1) logarithmic bins of
  # log size dx so that the last bin starts at max_w
  min_x <- log10(min_w)
  max_x <- log10(max_w)
  dx <- (max_x - min_x) / (no_w - 1)
  x <- seq(min_x, by = dx, length.out = no_w)
  w <- 10 ^ x

  # Find index of nearest grid point to min_w_inf that is an integer multiple
  # of the species spacing away from max_w
  min_x_inf <- log10(min_w_inf)
  max_x_inf <- log10(max_w_inf)
  # bins_per_sp is the number of bins separating species
  bins_per_sp <- round((max_x_inf - min_x_inf) / (dx * (no_sp - 1)))
  min_i_inf <- no_w - (no_sp - 1) * bins_per_sp
  # Maximum sizes for all species
  w_inf <- w[seq(min_i_inf, by = bins_per_sp, length.out = no_sp)]

  # Find index of nearest grid point to min_w_mat
  min_x_mat <- log10(min_w_mat)
  min_i_mat <- round((min_x_mat - min_x) / dx) + 1
  # Maturity sizes for all species
  w_mat <- w[seq(min_i_mat, by = bins_per_sp, length.out = no_sp)]

  if (egg_size_scaling) {
    # Determine egg weights w_min for all species
    w_min_idx <- seq(1, by = bins_per_sp, length.out = no_sp)
    w_min <- w[w_min_idx]
  } else {
    w_min <- rep(min_w, no_sp)
    w_min_idx <- rep(1, no_sp)
  }

  ## Build Params Object ----
  erepro <- 0.1  # Will be changed later to achieve coexistence
  species_params <- data.frame(
    species = as.factor(1:no_sp),
    w_min = w_min,
    w_inf = w_inf,
    w_mat = w_mat,
    w_min_idx = w_min_idx,
    h = h,
    ks = ks,
    f0 = f0,
    fc = fc,
    beta = beta,
    sigma = sigma,
    z0 = 0,
    alpha = alpha,
    erepro = erepro,
    lineage = as.factor(1:no_sp), # species is going to be the name, lineage will be the common ancestor
    zeta = zeta,
    stringsAsFactors = FALSE
  )
  gear_params <- data.frame(
    gear = gear_names,
    species = species_params$species,
    sel_func = "knife_edge",
    knife_edge_size = knife_edge_size,
    catchability = 1,
    stringsAsFactors = FALSE
  )
  params <-
    newMultispeciesParams(
      species_params,
      gear_params = gear_params,
      min_w = min_w,
      no_w = no_w,
      max_w = max_w,
      w_pp_cutoff = max_w,
      lambda = lambda,
      kappa = kappa,
      n = n,
      p = p,
      min_w_pp = min_w_pp,
      r_pp = r_pp,
      RDD = RDD
    )

  w <- params@w
  dw <- params@dw
  w_full <- params@w_full
  ks <- params@species_params$ks[[1]]

  ## Construct steady state solution ----

  # Get constants for steady-state solution
  # Predation mortality rate coefficient
  mu0 <- get_power_law_mort(params)
  # Add backgound mortality rate
  mu0 <- mu0 / (1 - ext_mort_prop)
  hbar <- alpha * h * f0 - ks
  if (hbar < 0) {
    stop("The feeding level is not sufficient to maintain the fish.")
  }
  pow <- mu0 / hbar / (1 - n)
  if (pow < 1) {
    message("The ratio of death rate to growth rate is too small, leading to
                an accumulation of fish at their largest size.")
  }

  initial_n <- params@psi  # get array with correct dimensions and names
  initial_n[, ] <- 0
  mumu <- mu0 * w^(n - 1)  # Death rate
  i_inf <- min_i_inf  # index of asymptotic size
  i_min <- 1  # index of natural egg size
  for (i in 1:no_sp) {
    gg <- hbar * w^n * (1 - params@psi[i, ])  # Growth rate
    idx <- w_min_idx[i]:(i_inf - 2)
    # Steady state solution of the upwind-difference scheme used in project
    n_exact <- c(1, cumprod(gg[idx] / ((gg + mumu * dw)[idx + 1])))
    # Use the first species for normalisation
    if (i == 1) {
      dist_sp <- bins_per_sp * dx
      mult <- kappa /
        sum(n_exact * (w^(lambda - 1) * dw)[1:(min_i_inf - 1)]) *
        (10^(dist_sp*(1-lambda)/2) - 10^(-dist_sp*(1-lambda)/2)) /
        (1-lambda)
    }
    if (!egg_size_scaling) {
      n_exact <- n_exact / n_exact[i_min]
    }
    idxs <- w_min_idx[i]:(i_inf - 1)
    initial_n[i, idxs] <- n_exact * mult *
      (w_inf[1] / w_inf[i]) ^ lambda
    i_inf <- i_inf + bins_per_sp
    i_min <- i_min + bins_per_sp
  }

  # Calculate the community spectrum
  sc <- colSums(initial_n)
  params@sc <- sc

  ##  Setup resource ----
  if (resource_scaling) {
    resource_vec <- (kappa * w ^ (-lambda)) - sc
    # Cut off resource at w_pp_cutoff
    resource_vec[w >= w_pp_cutoff] <- 0
    if (any(resource_vec < 0)) {
      if (!perfect_scaling) {
        # Do not allow negative resource abundance
        message("Note: Negative resource abundance values overwritten with zeros")
        resource_vec[resource_vec < 0] <- 0
      } else {
        message("Note: Negative resource abundances")
      }
    }
    params@cc_pp[sum(params@w_full <= w[1]):length(params@cc_pp)] <-
      resource_vec
  }
  if (!perfect_scaling) {
    params@cc_pp[w_full >= w_pp_cutoff] <- 0
  }

  initial_n_pp <- params@cc_pp
  # The cc_pp factor needs to be higher than the desired steady state in
  # order to compensate for predation mortality
  m2_background <- getResourceMort(params, initial_n, initial_n_pp)
  params@cc_pp <- (params@rr_pp + m2_background ) * initial_n_pp / params@rr_pp

  ## Setup external death ----
  m2 <- getPredMort(params, initial_n, initial_n_pp)
  flag <- FALSE
  for (i in 1:no_sp) {
    # The steplike psi was only needed when we wanted to use the analytic
    # expression for the steady-state solution
    # params@psi[i,] <- (w / w_inf[i]) ^ (1 - n)
    # params@psi[i, w < (w_mat[i] - 1e-10)] <- 0
    # params@psi[i, w > (w_inf[i] - 1e-10)] <- 1
    params@mu_b[i,] <- mu0 * w ^ (n - 1) - m2[i, ]
    if (!perfect_scaling && any(params@mu_b[i,] < 0)) {
      params@mu_b[i, params@mu_b[i,] < 0] <- 0
    }
  }


  ## Set erepro to meet boundary condition ----
  rdi <- getRDI(params, initial_n, initial_n_pp)
  gg <- getEGrowth(params, initial_n, initial_n_pp)
  mumu <- getMort(params, initial_n, initial_n_pp)
  erepro_final <- 1:no_sp  # set up vector of right dimension
  for (i in (1:no_sp)) {
    gg0 <- gg[i, params@w_min_idx[i]]
    mumu0 <- mumu[i, params@w_min_idx[i]]
    DW <- params@dw[params@w_min_idx[i]]
    erepro_final[i] <- erepro *
      (initial_n[i, params@w_min_idx[i]] *
         (gg0 + DW * mumu0)) / rdi[i]
  }
  if (is.finite(R_factor)) {
    # erepro has been multiplied by a factor of (R_factor/(R_factor-1)) to
    # compensate for using Beverton Holt function.
    erepro_final <- (R_factor / (R_factor - 1)) * erepro_final
  }
  params@species_params$erepro <- erepro_final

  # Record abundance of fish and resource at steady state, as slots.
  params@initial_n <- initial_n
  params@initial_n_pp <- initial_n_pp
  # set rmax=fac*RDD
  # note that erepro has been multiplied by a factor of (R_factor/(R_factor-1)) to
  # compensate for using a Beverton Holt function
  params@species_params$R_max <-
    (R_factor - 1) * getRDI(params, initial_n, initial_n_pp)

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
finalTouch <- function(folder,params,t_max)
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
  for(iRun in 1:length(dir(folder)))
  {
    tempRun <- readRDS(paste(folder,"/run",iRun,".rds",sep=""))
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

evoProject <- function(params,t_max = 100, mutation = 2,saveFolder)
{
  t_event <- sort(sample(1:(t_max-1),mutation))
  # corrected_t_event <- t_event+1# need to add 1 per run as they contain a 0 time step (should I get rid of it?)
  # we know when the mutations will appear, now we need to calculate the different time intervals bewteen these mutations
  t_max_vec <- neighbourDistance(x = c(t_event,t_max))


  mySim <- project(params, t_max = t_max_vec[1])
  saveRDS(mySim,file= paste(folder,"/run1.rds", sep = ""))
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
    mySim <- project(params, t_max = t_max_vec[iSim])
    saveRDS(mySim,file= paste(folder,"/run",iSim,".rds", sep = ""))
  }

  sim <- finalTouch(folder = folder, params = params, t_max = t_max)

  return(sim)
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
    initial_effort = params@initial_effort
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

#' stock recruitment relationship enabling extinction of species
#' when a species reaches an abundance threshold its recruitment gets disabled, species is later removed
#' For now each species acts as its own
#' TODO make rmax apply at the species level and not inidividual phenotypes
#' TODO at the moment I am using rdi to apply the threhold but it should be done on n directly
extinctionRDD <- function(rdi, species_params, ...) {
  if (!("R_max" %in% names(species_params))) {
    stop("The R_max column is missing in species_params.")
  }
  rdd <- rdi / (1 + rdi/species_params$R_max)
  if(sum(which(rdd <= 1e-30))) rdd[which(rdd <= 1e-30)] <- 0 # if any of the rdd is under threshold, set it to 0
  return(rdd)
}
