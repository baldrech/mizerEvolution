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
#' @param updateParams ...
#' @param ea_int ...
#' @param ca_int ...
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
                      perfect_scaling = FALSE,
                      updateParams = NULL,
                      # interactionMatrix = NULL,
                      ea_int = 0, # temperature parameters
                      ca_int = 0
) {


  if(is.null(updateParams))
  {
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

    params@species_params$lineage <- as.factor(1:no_sp) # parameter to remember the origin species of the phenotypes
    params@species_params$name <- nameGenerator(no_sp)
    # temperature parameters
    params@species_params$ea_int <- ea_int
    params@species_params$ca_int <- ca_int

  } else {
    # if the species column has characters instead of numeric it becomes annoying
    new_sp_params <- updateParams@species_params
    speciesName <- new_sp_params$species
    new_sp_params$species <- as.factor(1:dim(new_sp_params)[1])

    params <- newMultispeciesParams(new_sp_params, interaction = updateParams@interaction, min_w_pp = updateParams@w_full[1])
    params <- setResource(params, r_pp = updateParams@resource_params$r_pp,
                          kappa = updateParams@resource_params$kappa,
                          lambda = updateParams@resource_params$lambda,
                          n = updateParams@resource_params$n,
                          w_pp_cutoff = updateParams@resource_params$w_pp_cutoff)
    # setPredKernel(params, pred_kernel = getPredKernel(updateParams))

    params@gear_params$species <- as.factor(params@gear_params$species) # for some reason the function above gives species as char
    params@species_params$lineage <- as.factor(1:dim(params@species_params)[1]) # parameter to remember the origin species of the phenotypes
    params@species_params$name <- speciesName


    #   if(!is.null(interactionMatrix))
    #   {
    #     dimnames(interactionMatrix) <- list(as.factor(1:dim(inter)[1]),as.factor(1:dim(inter)[2]))
    #     params <- setInteraction(params, interaction = interactionMatrix)
    #   }
  }

  params@species_params$zeta <- zeta # amplitude of lognorm distribution around trait when generating new phenotype
  params@species_params$pop <- 1 # when the phenotype entered the simulation
  params@species_params$ext <- F # when the phenotype left the simulation
  params <- setReproduction(params, RDD = RDD)

  # temperature parameters
  params@species_params$ea_int <- ea_int
  params@species_params$ca_int <- ca_int

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
  biomass[is.na(biomass)] <- 0 # don;t want any NA from mid-sim popping phenotypes
  tempRun@n <- biomass
  tempRun@effort <- effort
  tempRun@n_pp <- biomassPP
  # tempRun@n_other <- matrix(0, nrow = (t_max+1),ncol = 1, dimnames = list("time" = 1:(t_max+1), "component"))

  component_names <- names(params@other_dynamics)
  no_components <- length(component_names)
  no_t <- t_max +1

  list_n_other <- rep(list(NA), no_t * no_components)
  dim(list_n_other) <- c(no_t, no_components)
  dimnames(list_n_other) <- list(time = 1:(t_max+1),
                                 component = component_names)
  tempRun@n_other <- list_n_other

  # update the params for any extinct species
  for(iSpecies in tempRun@params@species_params$species)
  {
    itime <- dim(tempRun@n)[1]
    count = 0
    while(sum(tempRun@n[itime,iSpecies,])<=1e-30 && count < dim(tempRun@n)[1])
    {
      count = count +1
      # print("species")
      # print(iSpecies)
      # print(itime)
      itime = itime - 1
      tempRun@params@species_params$ext[as.numeric(iSpecies)] <- itime
      # print("extinction")
      # print(tempRun@params@species_params$ext)
    }
  }


  return(tempRun)
}

#' Project
#'
#' @param initCondition ...
#' @param params ...
#' @param t_max ...
#' @param dt ...
#' @param mutation ...
#' @param trait ...
#' @param initPool ...
#' @param initSpread ...
#' @param alien ...
#' @param alien_init_n ...
#' @param saveFolder ...
#' @param effort ...
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
evoProject <- function(params = NULL, initCondition = NULL, t_max = 100, dt = 0.1, # params and initCondition cannot be both null at same time, need to specify that
                       mutation = 2, trait = "w_mat", initPool = 0, initSpread = 5, alien = 0,
                       alien_init_n = NULL, trait_range = NULL,
                       saveFolder = file.path(tempdir(), "simTemp"), effort = 0)
{
  ## Check first if starting a new simulation or using a previous sim as parameter object
  if(is.null(initCondition))
  {
    ## If initPool is positive, we need to create a batch of new phenotypes at the start of the simulation.
    if (initPool > 0)
    {
      for (iSpecies in sort(unique(params@species_params$lineage)))
      {
        # Generate phenotypes pool
        for (iPhenotype in seq(1, initPool))
        {
          newSp <- params@species_params[params@species_params$species == iSpecies,] # perfect copy
          newSp$species <-factor(as.character(max(as.numeric(params@species_params$species))+1), levels = max(as.numeric(params@species_params$species))+1) # new species name but lineage stays the same

          switch(trait,
                 size = {
                   # Trait = asymptotic size
                   sd = as.numeric(mAmplitude * initSpread *  params@species_params[which(params@species_params$ecotype == iSpecies),]$w_inf)
                   newSp$w_inf <- abs(newSp$w_inf + rnorm(1, 0, sd)) # change a bit the asymptotic size
                   newSp$w_mat <- newSp$w_inf * eta # calculate from the new w_inf value
                   newSp$z0 <- z0pre * as.numeric(newSp$w_inf) ^ (n - 1) # if I don't put as.numeric I lose the name z0
                 },
                 Beta = {
                   # Trait = PPMR
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$beta)
                   newSp$beta <- abs(newSp$beta + rnorm(1, 0, sd)) # change a bit the PPMR
                   while(newSp$beta < 10) newSp$beta <-  abs(newSp$beta + rnorm(1, 0, sd))
                   alpha_e <- sqrt(2 * pi) * newSp$sigma * newSp$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * newSp$sigma ^ 2 / 2) *
                     (pnorm(3 - (lambda - 2) * newSp$sigma) + pnorm(log(newSp$beta)/newSp$sigma + (lambda - 2) * newSp$sigma) - 1)
                   # newSp$gamma <- h * f0 / (alpha_e * kappa * (1 - f0))
                 },
                 sigma = {
                   # Trait = fedding kernel
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$sigma)
                   newSp$sigma <- abs(newSp$sigma + rnorm(1, 0, sd)) # change a bit the diet breadth
                   while(newSp$sigma < .5 | newSp$sigma > 5)  newSp$sigma <-  abs(newSp$sigma + rnorm(1, 0, sd))
                   alpha_e <- sqrt(2 * pi) * newSp$sigma * newSp$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * newSp$sigma ^ 2 / 2)*
                     (pnorm(3 - (lambda - 2) * newSp$sigma) + pnorm(log(newSp$beta)/newSp$sigma + (lambda - 2) * newSp$sigma) - 1)
                   newSp$gamma <- h * f0 / (alpha_e * kappa * (1 - f0))
                 },
                 predation = {
                   # PPMR
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$beta)
                   newSp$beta <- abs(newSp$beta + rnorm(1, 0, sd)) # change a bit the PPMR
                   while(newSp$beta < 10) newSp$beta <-  abs(newSp$beta + rnorm(1, 0, sd))
                   # feeding kernel
                   sd = as.numeric(mAmplitude *  params@species_params[which(params@species_params$ecotype == iSpecies),]$sigma)
                   newSp$sigma <- abs(newSp$sigma + rnorm(1, 0, sd)) # change a bit the diet breadth
                   while(newSp$sigma < .5 | newSp$sigma >5 )  newSp$sigma <-  abs(newSp$sigma + rnorm(1, 0, sd))
                   # recalculate gamma if necessary
                   alpha_e <- sqrt(2 * pi) * newSp$sigma * newSp$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * newSp$sigma ^ 2 / 2)*
                     (pnorm(3 - (lambda - 2) * newSp$sigma) + pnorm(log(newSp$beta)/newSp$sigma + (lambda - 2) * newSp$sigma) - 1)
                   newSp$gamma <- h * f0 / (alpha_e * kappa * (1 - f0))
                 },
                 eta = {
                   # Trait = eta
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$eta)
                   newSp$eta <- abs(newSp$eta + rnorm(1, 0, sd)) # change a bit eta
                   newSp$w_mat <- newSp$w_inf * newSp$eta # update
                 },
                 ed_int = {
                   # Trait = ed_int
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$ed_int)
                   newSp$ed_int <- abs(newSp$ed_int + rnorm(1, 0, sd))
                   while(newSp$ed_int < 2.5) newSp$ed_int <- abs(newSp$ed_int + rnorm(1, 0, sd))
                 },
                 t_d = {
                   # Trait = ed_int
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$t_d)
                   newSp$t_d <- abs(newSp$t_d + rnorm(1, 0, sd))
                 },
                 temperature = {
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$ed_int)
                   newSp$ed_int <-abs( newSp$ed_int + rnorm(1, 0, sd))
                   while(newSp$ed_int < 2.5) newSp$ed_int <- abs(newSp$ed_int + rnorm(1, 0, sd)) # ed_int cannot go lower than 2.5
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$t_d)
                   newSp$t_d <- abs(newSp$t_d + rnorm(1, 0, sd))
                 },

                 all = {
                   # Trait = asymptotic size
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$w_inf)
                   newSp$w_inf <- abs(newSp$w_inf + rnorm(1, 0, sd)) # change a bit the asymptotic size
                   newSp$w_mat <- abs(newSp$w_inf * eta) # calculate from the new w_inf value
                   newSp$z0 <- z0pre * as.numeric(newSp$w_inf) ^ (n - 1) # if I don't put as.numeric I lose the name z0
                   # Trait = predation
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$beta)
                   newSp$beta <- abs(newSp$beta + rnorm(1, 0, sd)) # change a bit the PPMR
                   sd = as.numeric(mAmplitude * initSpread * params@species_params[which(params@species_params$ecotype == iSpecies),]$sigma)
                   newSp$sigma <- abs(newSp$sigma + rnorm(1, 0, sd)) # change a bit the diet breadth
                   # calculate the new gamma
                   alpha_e <- sqrt(2 * pi) * newSp$sigma * newSp$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * newSp$sigma ^ 2 / 2)*
                     (pnorm(3 - (lambda - 2) * newSp$sigma) + pnorm(log(newSp$beta)/newSp$sigma + (lambda - 2) * newSp$sigma) - 1)

                   newSp$gamma <- h * f0 / (alpha_e * kappa * (1 - f0))
                 },
                 {

                   sd = newSp$zeta * initSpread * as.numeric(newSp[trait])
                   newSp[trait] <- abs(newSp[trait] + rnorm(1, 0, sd))

                 })


          # update initial abundance matrix
          n_init <- params@initial_n
          n_newSp <- array(0,dim = c(1,dim(n_init)[2]), dimnames = list("sp" = newSp$species, "w" = dimnames(n_init)$w))
          n_init <- rbind(n_init,n_newSp) # this include the new mutant as last column
          names(dimnames(n_init)) <- c("sp","W") # getting the dimnames back

          params <- addSpecies(params = params, species_params = newSp, init_n= n_init)







        }
      }

      # redistribute the abundance of the phenotypes more randomly
      original_n <- params@initial_n[unique(params@species_params$lineage),] # getting the abundance of the original species, all the others are set to 0
      init_n <- params@initial_n # the one to fill with new values
      for (iSpecies in unique(params@species_params$lineage)) # for every species
      {
        #the total abundance is randomly distributed among all the phenotypes
        biomRandom<- runif(initPool+1,0,1)
        biomFrac <- biomRandom/sum(biomRandom) # that's the fraction to apply to the initial abundance to spread the biomass among phenotypes
        position = 0
        for(iPhen in which(params@species_params$lineage == iSpecies)) # which phen are within ispecies
        {
          position = position + 1
          init_n[iPhen,] <- original_n[iSpecies,] * biomFrac[position]
        }
      }


      params@initial_n <- init_n

    }
  }
  else
  { ## if initCondition is provided, using it's params object and last step of biomasses to start the new simulation
    ## TODO needs improvements, probably not importing everything at this stage
    # not using the built in stuff in the project function for now
    params <- initCondition@params
    params@initial_n <- initCondition@n[dim(initCondition@n)[1],,]
    params@initial_n_pp <- initCondition@n_pp[dim(initCondition@n_pp)[1],]

  }
  ## At this stage we have a param object, now need to complete the setup
  SpIdx <- unique(params@species_params$lineage) # Handy species index

  # The model creates short simulations and temporarily save them to free RAM space, combining them at the end
  if(!dir.exists(saveFolder)) dir.create(saveFolder)

  ## Evolution section, using mutation variable
  if(is.numeric(mutation))
  { # If mutation is numeric, it's a rate of mutation per species per time step per thousand.
    # New phenotypes issue from mutation are planned in advance in a matrix.
    # TODO this framework assumes stable initial ecosystem, maybe do case for extinctions?
    # mutationPerSteps <- mutation
    t_mutation <- matrix(0,nrow = length(SpIdx), ncol = (t_max/1), dimnames = list("species" = SpIdx, "time" = 1:(t_max/1)))
    for(iSpecies in SpIdx) # for each species
    {
      for(iTime in 1:(t_max/1)) # for each time step
      {
        if(mutation > sample(seq(0,100,.1), 1)) # random draw equivalent to a chance of event per thousand
          t_mutation[iSpecies,iTime] <- 1
      }
    }
    t_phen <- apply(t_mutation,2,sum) # t_mutation knows which species mutates and t_phen knows which time (with no species info)
    t_phen <- which(t_phen >=1)
  } else if (is.data.frame(mutation))
  { # If mutation is a dataframe, then it's a list of all the mutants that are going to be introduced during the simulation
    # TODO probably a better way to mimic t_mutation creation than the one below
    t_phen <- mutation$time
    t_mutation <- matrix(0,nrow = length(SpIdx), ncol = (t_max/1), dimnames = list("species" = SpIdx, "time" = 1:(t_max/1))) # feels like a remnant from the randome mutant generation, maybe can get rid of it down the line?
    for(iRow in 1:dim(mutation)[1]) t_mutation[mutation$species[iRow],mutation$time[iRow]] <- 1
  } else (stop("The mutation argument has to be numeric or dataframe. Try again."))

  ## Invasion section, using alien variable
  if(is.numeric(alien))
  { # If alien is numeric, it's a rate of invasion per time step per thousand
    # TODO t_invasion is built the same as t_mutation for code compatiblility below but it doesn't have to, shouldn't actually.
    t_invasion <- matrix(0,nrow = length(SpIdx), ncol = (t_max/1), dimnames = list("alien" = SpIdx , "time" = 1:(t_max/1)))
    for(iTime in 1:(t_max/1)) # for each time step
    {
      if(alien > sample(seq(0,100,.1), 1)) # random draw equivalent to a chance of event per thousand
        t_invasion[,iTime] <- 1
    }
    t_alien <- apply(t_invasion,2,sum) # t_invasion is just a t_mutation equivalent to build t_alien which contains invasion occurence
    t_alien <- which(t_alien >=1)
  } else if (is.data.frame(alien))
  { # If alien is a dataframe, then it's a list of all the invaders that are going to be introduced during the simulation
    t_alien <- alien$time
    # TODO probably a better way to mimic t_invasion creation than the one below
    t_invasion <- matrix(0,nrow = dim(alien)[1], ncol = (t_max/1), dimnames = list("alien" = alien$species , "time" = 1:(t_max/1)))
    for(iRow in 1:dim(alien)[1]) t_invasion[alien$species[iRow],alien$time[iRow]] <- 1
  } else (stop("The alien argument has to be numeric or dataframe. Try again."))

  ## Creating time frame
  # Need to combine all the events together which tells us when to stop and add species.
  t_event <- sort(unique(c(t_phen,t_alien)))
  # Time intervals between these events, used to know simulation time (t_max argument)
  t_max_vec <- neighbourDistance(x = c(t_event,t_max))
  # An event at the last time step cannot be handeled so removing it
  if(t_max_vec[length(t_max_vec)] == 0) t_max_vec <- t_max_vec[-length(t_max_vec)]

  ## Starting the first simulation, using core Mizer
  mySim <- project(params, t_max = t_max_vec[1],progress_bar = F, effort = effort)
  if(length(t_max_vec) >1) # if there is at least one event planned
  {
    saveRDS(mySim,file= paste(saveFolder,"/run1.rds", sep = ""))
    for(iSim in 2:length(t_max_vec))
    {
      lastBiom_updated <- FALSE #if lastBiom gets updated and put in mySim@params@init_n, this becomes true and further species addition use mySim@params@init_n instead of mySim@n[dim(mySim@n)[1],,]
      # print("iSim")
      # print(iSim)
      # print("time to add something")
      # print(t_event[iSim-1])
      # when Mizer stops it can come from new phen, species invasion, or both

      # new alien
      #       print("t_invasion")
      #       print(t_invasion)
      #       print(t_event[iSim-1])
      #       print(t_invasion[,t_event[iSim-1]])
      # print(sum(t_invasion[,t_event[iSim-1]]))
      if(sum(t_invasion[,t_event[iSim-1]])) # if t_invasion is positive it means there is an invasion
      {
        if(is.numeric(alien))
        {
          #produce alien
          print("alien trying to invade ecosystem")
          newSp <- alien_synthesis(trait_range = trait_range)
          # newSp <- mySim@params@species_params[sample(1:dim(mySim@params@species_params)[1],1),] # for now randomly select a species present and heavily change some parameters
          newSp$species <- factor(as.character(max(as.numeric(mySim@params@species_params$species))+1), levels = max(as.numeric(mySim@params@species_params$species))+1) # new species name
          # newSp$w_inf <- sample(10^(seq(log10(min(mySim@params@species_params$w_inf)), log10(max(mySim@params@species_params$w_inf)), length.out = 100 )),1) # randomly select a size between min and max, log balanced
          # newSp$w_mat <- .25*newSp$w_inf
          # newSp$ks <- abs(newSp$ks + rnorm(1, 0, newSp$zeta * newSp$ks))
          # newSp$beta <- abs(newSp$beta + rnorm(1, 0, newSp$zeta * newSp$beta * 5))
          # newSp$sigma <- abs(newSp$sigma + rnorm(1, 0, newSp$zeta * newSp$sigma * 3))
          # newSp$erepro <- abs(newSp$erepro + rnorm(1, 0, newSp$zeta * newSp$erepro))
          newSp$name <- nameGenerator()
          newSp$R_max <- resource_params(mySim@params)$kappa * newSp$w_inf^-1
          # TODO automatically fill below slots
          newSp$ea_int <- 0
          newSp$ca_int <- 0
          newSp$zeta <- .2
          newSp$fc <- .25
          newSp$f0 <- .6
          newSp$w_min_idx <- 1
          newSp$lineage <- factor(as.character(max(as.numeric(mySim@params@species_params$lineage))+1), levels = max(as.numeric(mySim@params@species_params$lineage))+1)
          newSp$pop <- t_event[iSim-1]
          newSp$ext <- FALSE
          # when creating alien randomly, might need to test their survivability first, in case of crazy values
        } else if (is.data.frame(alien))
        {
          alien
          newSp <- alien[(iSim-1),]
          # print("invader extracted from df")
        } else {print("something went wrong")}
        # adding initial biomass and stuff
        lastBiom <- mySim@n[dim(mySim@n)[1],,]
        lastBiom_updated <- TRUE
        # not really sure how much to put for now

        n_newSp <- rep(0,dim(mySim@n)[3])
        #need to create size spectrum of abundance from one value
        n0_mult = alien_init_n # (boost the initial abundance) | apparently the initial biomass of the invading species are really low
        a = 0.35
        no_w <- length(mySim@params@w)
        initial_n <- array(NA, dim = c(1, no_w))
        # N = N0 * Winf^(2*n-q-2+a) * w^(-n-a)
        # Reverse calc n and q from intake_max and search_vol slots (could add get_n function)
        n <- (log(mySim@params@intake_max[,1] / mySim@params@species_params$h) / log(mySim@params@w[1]))[1]
        q <- (log(mySim@params@search_vol[,1] / mySim@params@species_params$gamma) / log(mySim@params@w[1]))[1]
        # Guessing at a suitable n0 value based on kappa - this was figured out using trial and error and should be updated
        if (is.null(n0_mult)) {
          lambda <- 2 + q - n
          kappa <- mySim@params@cc_pp[1] / (mySim@params@w_full[1]^(-lambda))
          n0_mult <- kappa / 1000
        }

        initial_n <- unlist(tapply(mySim@params@w, 1:no_w, function(wx,n0_mult,w_inf,a,n,q)
          n0_mult * w_inf^(2 * n - q - 2 + a) * wx^(-n - a),
          n0_mult = n0_mult, w_inf = newSp$w_inf, a=a, n=n, q=q))

        # print("init_n")
        # print(initial_n)

        #set densities at w > w_inf to 0
        initial_n[unlist(tapply(mySim@params@w,1:no_w,function(wx,w_inf) w_inf<wx, w_inf=newSp$w_inf))] <- 0
        # Also any densities at w < w_min set to 0
        initial_n[unlist(tapply(mySim@params@w,1:no_w,function(wx,w_min)w_min>wx, w_min=newSp$w_min))] <- 0
        n_newSp <- t(initial_n)
        # print("n_newSp")
        # print(n_newSp)

        init_n <- rbind(lastBiom,n_newSp) # this include the new mutant as last column
        names(dimnames(init_n)) <- c("sp","w")
        # print("new invader")
        # print(newSp)
        rownames(init_n)[length(rownames((init_n)))] <- as.character(newSp$species) # update the name of the mutant accordingly


        mySim@params <- addSpecies(params = mySim@params, species_params = newSp, init_n= init_n)

        # print(mySim@params@species_params)
        # print("alien invaded the ecosystem")
      }


      # new phen
      if(sum(t_mutation[,t_event[iSim-1]])) # if it's positive it means at least one phen appears
      {
        resident_list <- as.character(which(t_mutation[,t_event[iSim-1]]>0))

        for(resident_lineage in resident_list) # loop to handle multiple species additions at once
        {
          # if(length(resident_list)>1)
          # {
          #   print("more than one phen going to be added")
          #   count <- which(resident_lineage == resident_list)
          #   cat(sprintf("phen number %i \n", count))
          #
          #
          # }


          if(is.numeric(mutation))
          {
            ## new mutant param
            # randomly take a resident
            # print("t_mutation")
            # print(t_mutation)
            # print("t_event")
            # print(t_event)
            # print("isim")
            # print(iSim)
            # print("t_max_vec")
            # print(t_max_vec)
            # print("resident")
            # print(which(t_mutation[,t_event[iSim-1]]>0))
            # print("new phenotype born")

            resident <- as.character(sample(mySim@params@species_params$species[which(mySim@params@species_params$lineage == resident_lineage)],1))
            # print(mySim@n[,resident,1:2])
            phenCount <- 0
            while(sum(mySim@n[t_max_vec[iSim-1],resident,]) < 2e-28) # if resident is too low abundance, cannot produce a new phenotype that would be under the extinction threshol (need to put it as variables)
            {
              resident <- as.character(sample(mySim@params@species_params$species[which(mySim@params@species_params$lineage == resident_lineage)],1))
              phenCount <- phenCount +1
              if(phenCount> length(which(mySim@params@species_params$lineage == resident_lineage))*2) # if selecting twice the number of phen randomly in the species did not yield a whorthy parent, take anohter sp
              {
                lineageList <- unique(mySim@params@species_params$lineage)
                lineageList <- lineageList[-as.numeric(resident_lineage)] # remove the extinct species from the list
                resident_lineage <- sample(lineageList,1)
                phenCount <- 0 # reset count
              }
            }

            # The more stuff go extinct, the slower it will get and the more chance it might fail.
            # TODO update extinction time whenever possible so can select subset of live phenotypes directly

            #print(resident)
            # resident <- as.character(sample(mySim@params@species_params$species, 1))

            # create a new species param

            newSp <- mySim@params@species_params[resident,] # get a copy of the resident
            newSp$pop <- t_event[iSim-1]
            # print(newSp)
            # switch to determine what happens to the new species
            #TODO  need to rewrite the specific cases properly
            switch(trait, # cases with specific names and default if users gives just a parameter df name as trait
                   size = {
                     # Trait = asymptotic size

                     sd = as.numeric(mySim@params@species_params$zeta[as.numeric(resident)] * mySim@params@species_params$w_inf[as.numeric(resident),1])
                     newSp$w_inf <- newSp$w_inf + rnorm(1, 0, sd)
                     # need to get eta and Z0pre from somewhere before this works
                     # newSp$w_mat <- newSp$w_inf * newSp$eta
                     # newSp$z0 <- z0pre * as.numeric(newSp$w_inf) ^ (n - 1)

                   },
                   # Beta = {
                   #   # Trait = PPMR
                   #   sd = as.numeric(mAmplitude *  resident_params["beta"]) # standard deviation
                   #   mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd) # change a bit the PPMR
                   #   while(mutant$beta < 10)
                   #   {print("need to reroll beta")
                   #     mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd)
                   #   }
                   #   # calculate the new gamma
                   #   alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)*
                   #     (pnorm(3 - (lambda - 2) * mutant$sigma) + pnorm(log(mutant$beta)/mutant$sigma + (lambda - 2) * mutant$sigma) - 1)
                   #
                   #   # mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
                   #
                   #   cat(sprintf("parent beta:%f, gamma:%f\n",resident_params$beta,resident_params$gamma))
                   #   cat(sprintf("mutant beta:%f, gamma:%f\n",mutant$beta,mutant$gamma))
                   #   # cat(sprintf("Its PPMR is:%f\n",mutant$beta))
                   # },
                   # Sigma = {
                   #   # Trait = fedding kernel
                   #   sd = as.numeric(mAmplitude *  resident_params["sigma"]) # standard deviation
                   #   mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd) # change a bit the diet breadth
                   #   while(mutant$sigma < .5 | mutant$sigma > 5)
                   #   {print("need to reroll sigma")
                   #     mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd)}
                   #   # calculate the new gamma
                   #   alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)*
                   #     (pnorm(3 - (lambda - 2) * mutant$sigma) + pnorm(log(mutant$beta)/mutant$sigma + (lambda - 2) * mutant$sigma) - 1)
                   #   mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
                   #   #cat(sprintf("Its diet breadth mutes slightly.\n"))
                   # },
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
                   # eta = {
                   #   # Trait = eta
                   #   sd = as.numeric(mAmplitude *  resident_params["eta"]) # standard deviation
                   #   mutant["eta"] <- resident_params["eta"] + rnorm(1, 0, sd) # change a bit eta
                   #   if (mutant["eta"] >= 1) mutant["eta"] <- 0.95 # because yes it does happen
                   #   mutant["w_mat"] <- mutant["w_inf"] * mutant["eta"] # update
                   #   #cat(sprintf("Its w_mat is: %g\n",mutant["w_mat"]))
                   # },
                   # ed_int = {
                   #   # Trait = ed_int
                   #   sd = as.numeric(mAmplitude *  resident_params["ed_int"])
                   #   mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd)
                   #   while(mutant$ed_int < 2.5) mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd) # ed_int cannot go lower than 2.5
                   # },
                   # t_d = {
                   #   # Trait = ed_int
                   #   sd = as.numeric(mAmplitude *  resident_params["t_d"])
                   #   mutant$t_d <- mutant$t_d + rnorm(1, 0, sd)
                   #   cat(sprintf("Its name is %i and its trait value is %g\n", mutant$ecotype,mutant["t_d"]))
                   # },
                   # temperature = {
                   #   sd = as.numeric(mAmplitude * resident_params["ed_int"])
                   #   mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd)
                   #   while(mutant$ed_int < 2.5) mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd) # ed_int cannot go lower than 2.5
                   #   sd = as.numeric(mAmplitude * resident_params["t_d"])
                   #   mutant$t_d <- mutant$t_d + rnorm(1, 0, sd)
                   # },
                   # all = {
                   #   # Trait = asymptotic size
                   #   sd = as.numeric(mAmplitude *  resident_params["w_inf"]) # standard deviation
                   #   mutant["w_inf"] <- resident_params["w_inf"] + rnorm(1, 0, sd) # change a bit the asymptotic size
                   #   mutant["w_mat"] <- mutant["w_inf"] * eta # calculate from the new w_inf value
                   #   mutant["z0"] <- z0pre * as.numeric(mutant["w_inf"]) ^ (n - 1) # if I don't put as.numeric I lose the name z0
                   #   # Trait = predation
                   #   sd = as.numeric(mAmplitude *  resident_params["beta"]) # standard deviation
                   #   mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd) # change a bit the PPMR
                   #   sd = as.numeric(mAmplitude *  resident_params["sigma"]) # standard deviation
                   #   mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd) # change a bit the diet breadth
                   #   # calculate the new gamma
                   #   alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)*
                   #     (pnorm(3 - (lambda - 2) * mutant$sigma) + pnorm(log(mutant$beta)/mutant$sigma + (lambda - 2) * mutant$sigma) - 1)
                   #   mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
                   #   #cat(sprintf("Its traits mute slightly.\n"))
                   # },
                   {
                     sd = as.numeric(mySim@params@species_params$zeta[as.numeric(resident)] * mySim@params@species_params[trait][as.numeric(resident),1])
                     newSp[trait] <- abs(newSp[trait] + rnorm(1, 0, sd))

                   })

            # set the abundance for all species to start a new project
            if(lastBiom_updated) # for all mutants after the first one being added, the initial value is already stored in the params instead of previous sim last time step
              lastBiom <- mySim@params@initial_n
            else lastBiom <- mySim@n[dim(mySim@n)[1],,]
            lastBiom_updated <- TRUE
            #n_newSp <- rep(0,dim(mySim@n)[3])
            n_newSp = 0.05 * lastBiom[dimnames(mySim@n)$sp == resident,] # the initial abundance is 5% of the resident pop
            lastBiom[dimnames(mySim@n)$sp ==resident,]= lastBiom[dimnames(mySim@n)$sp == resident,] - 0.05*lastBiom[dimnames(mySim@n)$sp ==resident,] # Withdraw the abundance of the mutant from its parent (we're not talking about eggs here but different ecotype already present)


          } else if (is.data.frame(mutation))
          {
            newSp <- mutation[iSim-1,]

            lastBiom <- mySim@n[dim(mySim@n)[1],,]
            lastBiom_updated <- TRUE
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
          # print("mutant generated")

          newSp$species <-factor(as.character(max(as.numeric(mySim@params@species_params$species))+1), levels = max(as.numeric(mySim@params@species_params$species))+1) # new species name but lineage stays the same

          init_n <- rbind(lastBiom,n_newSp) # this include the new mutant as last column
          names(dimnames(init_n)) <- c("sp","w")
          # print(newSp)
          rownames(init_n)[length(rownames((init_n)))] <- as.character(newSp$species) # update the name of the mutant accordingly

          # print("new phenotype implemented")

          mySim@params <- addSpecies(params = mySim@params, species_params = newSp, init_n= init_n)
        }
      }
      # print("mutant integrated in ecosystem, ready to project")
      if(t_max_vec[iSim]>0)
      {# happens if mutant appears at the last time step, makes the code crash | probably obsolete but does not hurt anyone
        mySim <- project(mySim@params, t_max = t_max_vec[iSim],progress_bar = F, effort = effort)
        # print("projected until next mutant, ready to save truncated run")
        saveRDS(mySim,file= paste(saveFolder,"/run",iSim,".rds", sep = ""))
      }
    }
    # return(list(saveFolder,mySim@params,t_max))
    print("Data handling")
    sim <- finalTouch(saveFolder = saveFolder, params = mySim@params, t_max = t_max)
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

  # fishing params | need to update them here as the default contruction functions will use wrong values / names
  params@gear_params <- rbind(params@gear_params,params@gear_params[species_params$lineage,]) #copy catchability of parent
  levels(params@gear_params$species) <- c(levels(params@gear_params$species),as.character(species_params$species)) # add new factor level (new mutant)
  params@gear_params$species[dim(params@gear_params)[1]] <- species_params$species # correct the mutant name

  # new params object ----
  # TODO check updated GUstav version for any mismatch with gear, effort, resource, etc
  # use dataframe and global settings from params to make a new MizerParams
  # object.
  # TODO going to need to check if params need to be updated with the new sp or not
  # print(3.5)
  # print(params@initial_effort)
  # print(params@gear_params)
  # print("catch addspecies")
  # print(params@catchability)
  p <- newMultispeciesParams(
    combi_species_params,
    interaction = inter,
    min_w = min(params@w),
    max_w = max(params@w),
    min_w_pp = min(params@w_full),
    no_w = length(params@w),
    initial_effort = params@initial_effort,
    gear_params = params@gear_params,
    # selectivity = params@selectivity,
    # catchability = params@catchability,
    RDD = params@rates_funcs$RDD,
    n = params@resource_params$n,
    r_pp = params@resource_params$r_pp,
    kappa = params@resource_params$kappa,
    lambda = params@resource_params$lambda,
    w_pp_cutoff = params@resource_params$w_pp_cutoff
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


#' function generating new alien species
#'
#' @description
#' Function uses a range of trait to produce alien
#' Default range of trait is based on North Sea
#' parameters but can be supplied by user
#'
#' @param trait_range a dataframe containing the
#' parameters to change and their accepted values.
#' Default is NULL, leading to the use of the NS_params
#' object
#' @param n An integer to know the number of aliens to
#' generate
#'
#' @export
#TODO add more default mizer params

alien_synthesis <- function(trait_range, n = 1){

  if(is.null(trait_range)) # using NS_params
  {
    # trait_range <- data.frame("trait" = c("w_inf", "betaS","betaL", "sigma", "k_vb", "ks","eta"),
    #                           "distribution" = c("lnorm","norm", "norm","norm","norm","norm","norm"),
    #                           "mean" = c(6.7648459,186.22222,243488.33,1.7166667,0.44408333,6.2893097,0.12071530),
    #                           "sd" = c(2.2524833,176.59733,144374.82,0.5742144,0.27293481,2.5737627,0.11480856))

    trait_range <- data.frame("trait" = c("w_inf", "beta", "sigma", "k_vb", "ks","eta"),
                              "distribution" = c("lnorm","uform", "uform","uform","uform","norm"),
                              "var1" = c(6.7648459,10,0.8,0.1,2.8,0.12071530),
                              "var2" = c(2.2524833,400000,3.2,1,13,0.11480856),row.names = 1)
  }

  for(iAlien in 1:n) # for n number of alien
  {
    # for(iRow in dim(trait_range)[1]) # go through each row of the trait df
    # {
    # hard to make it user fool proof so going to focus on default df for now

    # for beta, range is wide so most values are high and low values are rarely taken. Using logscale to circumvent

    w_inf <- sigma <- k_vb <- ks <- eta <- beta <- -1 # initialization for while loop

    while(w_inf <0 || w_inf >10000) w_inf <- rlnorm(1, trait_range["w_inf",]$var1, trait_range["w_inf",]$var2) #TODO remove size limit or make it a var
    while(beta <0) beta <- runif(1, log10(trait_range["beta",]$var1), log10(trait_range["beta",]$var2))
    while(sigma <0)    sigma <- runif(1, trait_range["sigma",]$var1, trait_range["sigma",]$var2)
    while(k_vb <0)    k_vb <- runif(1, trait_range["k_vb",]$var1, trait_range["k_vb",]$var2)
    while(ks <0)    ks <- runif(1, trait_range["ks",]$var1, trait_range["ks",]$var2)
    while(eta <0.027)    eta <- rnorm(1, trait_range["eta",]$var1, trait_range["eta",]$var2) # eta can get really small in NS_params, just making a threshold at the min value of the ecosystem

    w_mat <- w_inf * eta
    beta <- 10^beta

    species_df <- data.frame(
      "w_inf" = w_inf,
      "w_mat" = w_mat,
      "beta" = beta,
      "sigma" = sigma,
      "k_vb" = k_vb,
      "ks" = ks)
  }
  # print("alien embryo")
  # print(species_df)

  return(species_df)
}

#' name generator
#'
#' @description
#' This function randomly create a name composed of
#' 3 consonants and 2 vowels. It is used when randomly
#' generating invasive species
#'
#'
#' @export

nameGenerator <- function(iteration = 1){
  name_vec = NULL
  for(i in 1:iteration)
  {
    vowel <- letters[c(1,5,9,15,21,25)]
    consonant <- letters[-c(1,5,9,15,21,25)]
    CONSONANT <- LETTERS[-c(1,5,9,15,21,25)]

    firstLetter <- sample(CONSONANT,1)
    vowels <- sample(vowel,2,TRUE)
    consonants <- sample(consonant,2,TRUE)

    name <- paste0(firstLetter,vowels[1],consonants[1],vowels[2],consonants[2])
    name_vec <- c(name_vec,name)
  }

  return(name_vec)
}

#' tempFun is a function that takes temperature parameters (temperature, t_ref, t_d) and physiological parameters (w, Ea, c_a)
#' and returns a scalar based on the Padfield (2016) temperature performance equation. This will be upgraded in the future with
#' a deactivation portion with Anna Gardmark's equation once it's published
tempFun <- function(w, temperature, t_ref, Ea, c_a)
{
  k = 8.617332e-5 # Boltzmann constant
  # equation
  # (w^(c_a*(temperature-t_ref)))  *exp((-Ea/8.617332e-5)*((1/temperature) - (1/t_ref)))

  # converting to Kelvin from Celcius
  temperature <- temperature + 273
  t_ref <- t_ref + 273


  temperatureScalar <- t(sapply(w,FUN = function(x){x^(c_a*(temperature-(t_ref)))}) *exp((-Ea/k)*((1/temperature) - (1/(t_ref)))))
  return(temperatureScalar)
}


