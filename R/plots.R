#' Plot biomass variations through time of species and phenotypes
#'
#' @description
#'
#' @param object An object of class \linkS4class{MizerSim}
#' @param phenotype A boolean value to display the phenotype's
#' abundance on the plot. Default is TRUE
#' @param species Plot only a specific species and its phenotypes
#' Default is NULL
#' @param trait Display the trait value as color gradient. Works
#' only if `species` is not NULL
#' @param SpIdx Display a subset of species
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it A boolean value that displays comments. Default is TRUE
#' @param returnData A boolean value that controls the output. Either a dataframe
#' ready for ggplot (TRUE) or the plot itself (FALSE). Default is FALSE
#' @param save_it
#' @param nameSave
#' @param ylimit To edit the y axis scale on the plot
#'
#' @return A plot
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],
#' @examples
#' \donttest{
#' }

plotDynamics <- function(object,  phenotype = TRUE, species = NULL, trait = NULL, SpIdx = NULL,
                         time_range = c(min(as.numeric(dimnames(object@n)$time)),max(as.numeric(dimnames(object@n)$time))),
                         print_it = T, returnData = F, save_it = F,
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

            p <- ggplot(dplyr::filter(BiomPhen, sp == species)) +
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

                p <- ggplot(dplyr::filter(BiomPhen, sp == species)) +
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


#' Plot size spectrumms of species and phenotypes
#'
#' @description
#'
#' @inheritParams plotDynamics
#' @return A plot
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],
#' @examples
#' \donttest{
#' }

plotSS <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), min_w =min(object@params@w)/100, ylim = c(NA,NA),
                   biomass = TRUE, print_it = TRUE, species = TRUE, community = FALSE, save_it = FALSE, nameSave = "SizeSpectrum.png", returnData = FALSE, ...){

    cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
    # min_w = 0.001
    time_elements <- get_time_elements(object,time_range)
    spec_n <- apply(object@n[time_elements,,,drop=FALSE],c(2,3), mean)
    pkt_n <- apply(object@n_pp[time_elements,,drop=FALSE],2,mean)
    # alg_n <- apply(object@n_aa[time_elements,,drop=FALSE],2,mean)
    # ben_n <- apply(object@n_bb[time_elements,,drop=FALSE],2,mean)

    y_axis_name = "Abundance density in individuals.m^-3"
    if (biomass){
        spec_n <- sweep(spec_n,2,object@params@w,"*")
        pkt_n <- pkt_n * object@params@w_full
        # alg_n <- alg_n * object@params@w_full
        # ben_n <- ben_n * object@params@w_full
        y_axis_name = "Biomass in g.m^-3"
    }
    # Make data.frame for plot
    plot_datSP <- data.frame(value = c(spec_n), Species = dimnames(spec_n)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)), lineage = object@params@species_params$lineage)
    plot_datPkt <- data.frame(value = c(pkt_n), Species = "Phytoplankton", w = object@params@w_full)
    # plot_datAlg <- data.frame(value = c(alg_n), Species = "Algae", w = object@params@w_full)
    # plot_datBen <- data.frame(value = c(ben_n), Species = "Benthos", w = object@params@w_full)

    if(community) plot_datSP <- data.frame(value = apply(spec_n, 2, sum), w = object@params@w) else if (species)
    {
        dimnames(spec_n)$sp = object@params@species_params$lineage
        SpIdx = unique(object@params@species_params$lineage)
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
    # plot_datAlg <- plot_datAlg[(plot_datAlg$value > 0) & (plot_datAlg$w >= min_w),]
    # plot_datBen <- plot_datBen[(plot_datBen$value > 0) & (plot_datBen$w >= min_w),]
    #getPalette = colorRampPalette(brewer.pal(9, "Set1"))# increase the number of colors used

    if(community)
    {
        p <- ggplot(plot_datSP) +
            geom_line(aes(x=w, y = value)) +
            geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species),alpha = 0.5, color = "blue", size = 1.5) +
            # geom_line(data = plot_datAlg, aes(x = w, y = value, group = Species),alpha = 0.5, color = "green", size = 1.5) +
            # geom_line(data = plot_datBen, aes(x = w, y = value, group = Species),alpha = 0.5, color = "yellow", size = 1.5) +
            scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5)))+
            scale_y_continuous(name = y_axis_name, limits = ylim, trans = "log10") +
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
            # geom_line(data = plot_datAlg, aes(x = w, y = value, group = Species),alpha = 0.5, color = "green", size = 1.5) +
            # geom_line(data = plot_datBen, aes(x = w, y = value, group = Species),alpha = 0.5, color = "yellow", size = 1.5) +
            scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5)))+
            scale_y_continuous(name = y_axis_name, limits = ylim, trans = "log10") +
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
            geom_line(aes(x=w, y = value, colour = as.factor(lineage), group = Species)) +
            geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species),alpha = 0.5, color = "blue", size = 1.5) +
            # geom_line(data = plot_datAlg, aes(x = w, y = value, group = Species),alpha = 0.5, color = "green", size = 1.5) +
            # geom_line(data = plot_datBen, aes(x = w, y = value, group = Species),alpha = 0.5, color = "yellow", size = 1.5) +
            scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5)))+
            scale_y_continuous(name = y_axis_name, limits = ylim, trans = "log10") +
            labs(color='Species') +
            theme(panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"),
                  legend.key = element_rect(fill = "white"))+
            scale_colour_manual(values=cbPalette)+ # colorblind
            ggtitle("Size spectrum per phenotypes")

    }
    if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)

    if (returnData) return(plot_datSP) else if(print_it) return(p)
}


#' Plot feeding level
#'
#' @description
#'
#' @inheritParams plotDynamics
#' @return A plot
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],
#' @examples
#' \donttest{
#' }

plotevoFeeding <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, throughTime = F, start = 1000, every = 1000,
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

    Cfeed <- getCriticalFeedingLevel(object@params)#, ...) # get the critical feeding level

    if (species) # if I want to display species instead of ecotypes
    {
        dimnames(feed)$sp = object@params@species_params$lineage
        SpIdx = sort(unique(object@params@species_params$lineage)) # get the species names
        feed_sp = matrix(data = NA, ncol = dim(feed)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(feed)$w)) # prepare the new object
        Cfeed_sp = matrix(data = NA, ncol = dim(feed)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(feed)$w)) # prepare the new object

        names(dimnames(feed_sp))=list("species","size")
        names(dimnames(Cfeed_sp))=list("species","size")


        for (i in SpIdx)
        {
            #feeding level
            temp = feed # save to manip
            temp[which(rownames(feed) != i), ] = 0 # keep the ecotypes from the species only
            temp = apply(temp, 2, sum)
            temp = temp / length(which(rownames(feed)==i)) # do the mean (in 2 steps)
            feed_sp[which(rownames(feed_sp)==i), ] = temp

            #critical feeding
            temp = Cfeed # save to manip
            temp[which(rownames(Cfeed) != i), ] = 0 # keep the ecotypes from the species only
            temp = apply(temp, 2, sum)
            temp = temp / length(which(rownames(Cfeed)==i)) # do the mean (in 2 steps)
            Cfeed_sp[which(rownames(Cfeed_sp)==i), ] = temp


        }
        feed = feed_sp
        Cfeed = Cfeed_sp
    }

    a <- c(object@params@species_params$w_inf[SpIdx]) # to get vline of different col, need to create a data frame
    vlines <- data.frame(xint = a,grp = SpIdx)

    plot_dat <- data.frame(value = c(feed),critical = c(Cfeed), species = as.factor(dimnames(feed)[[1]]), size = rep(object@params@w, each=length(dimnames(feed)[[1]])))

    p <- ggplot(plot_dat) +
        geom_line(aes(x=size, y = value, colour = species)) +
        geom_line(aes(x=size, y = critical), linetype = "dashed", color = "red") +
        geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") +
        scale_x_log10(name = "Size", breaks = c(1 %o% 10^(-3:5)))  +
        scale_y_continuous(name = "Feeding Level", lim=c(0,1))+
        scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
        theme(panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),
              legend.key = element_rect(fill = "white"))+
        ggtitle(NULL)

    if(save_it) ggsave(plot = p, filename = nameSave,width = 18, height = 18, units = "cm" )

    if (returnData) return(list(plot_dat,Critfeed)) else if(print_it) return(p)
}


#' Plot growth
#'
#' @description
#'
#' @inheritParams plotDynamics
#' @return A plot
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],
#' @examples
#' \donttest{
#' }

plotevoGrowth <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, print_it = T, returnData = F, save_it = F, ylim = c(NA,NA),
                          nameSave = "Growth.png",...){

    time_elements <- get_time_elements(object,time_range)
    growth_time <- plyr::aaply(which(time_elements), 1, function(x){
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
        dimnames(growth)$sp = object@params@species_params$lineage
        SpIdx = sort(unique(object@params@species_params$lineage)) # get the species names
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

    # name = paste("Growth level at time",time_range,sep=" ")
    name = NULL
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


#' Plot starvation level
#'
#' @description
#'
#' @inheritParams plotDynamics
#' @return A plot
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],
#' @examples
#' \donttest{
#' }

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


#' Plot total mortality
#'
#' @description
#'
#' @inheritParams plotDynamics
#' @return A plot
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],
#' @examples
#' \donttest{
#' }

plotevoMortality <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),print_it = TRUE, returnData = F, comments = T){

    # effort can be in 2 forms

    if(is.matrix(object@effort)) effort = object@effort[time_range,]  else effort = object@effort[time_range]

    z <- getZ(object@params, n = object@n[time_range,,], n_pp = object@n_pp[time_range,],#n_aa = object@n_aa[time_range,],n_bb = object@n_bb[time_range,],
              effort = effort)#,
    #intakeScalar = object@intTempScalar[,,time_range], metScalar = object@metTempScalar[,,time_range], morScalar = object@morTempScalar[,,time_range])
    dimnames(z)$prey = object@params@species_params$lineage
    SpIdx = as.numeric(as.character(sort(unique(object@params@species_params$lineage)))) # get the species names

    # need to get rid of the extinct species at that time in SpIdx
    # a <- apply(object@n[time_range,,],1,sum)
    # names(a) <- sapply(names(a), function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
    # d <- rowsum(a, group = names(a))
    #
    # if (sum(d[,1] == 0))
    # {
    #   d <- d[-which(d[,1] == 0),]
    #   SpIdx <- as.numeric(names(d))
    # } else {SpIdx <- as.numeric(rownames(d))}

    z_sp = matrix(data = NA, ncol = dim(z)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(z)$w_prey)) # prepare the new object
    names(dimnames(z_sp))=list("prey","w_prey")

    for (iSpecies in SpIdx)
    {
        temp = z # save to manip
        temp[which(rownames(z) != iSpecies), ] = 0 # keep the ecotypes from the species only
        temp = apply(temp, 2, sum)
        temp = temp / length(which(rownames(z)==iSpecies)) # do the mean (in 2 steps)
        z_sp[which(rownames(z_sp)==iSpecies), ] = temp
    }
    z = z_sp

    # name = paste("Total Mortality at time",time_range,sep=" ")
    name = NULL

    plot_dat <- data.frame(value = c(z), Species = dimnames(z)[[1]], w = rep(object@params@w, each=length(dimnames(z)[[1]])))

    p <- ggplot(plot_dat) +
        geom_line(aes(x=w, y = value, colour = Species)) +
        scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-3:5))) +
        scale_y_continuous(name = "Mortality", lim=c(0,max(plot_dat$value))) +
        theme(legend.key = element_rect(fill = "white"),
              panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"))+
        ggtitle(NULL)

    if (returnData) return(plot_dat) else if(print_it) return(p)
}


#' Plot spawn
#'
#' @description
#'
#' @inheritParams plotDynamics
#' @return A plot
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],
#' @examples
#' \donttest{
#' }


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


#' Plot predation rate
#'
#' @description
#'
#' @inheritParams plotDynamics
#' @return A plot
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],
#' @examples
#' \donttest{
#' }

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


#' Plot predation mortality rate of each species against size
#'
#' After running a projection, plot the predation mortality rate of each species
#' by size. The mortality rate is averaged over the specified time range (a
#' single value for the time range can be used to plot a single time step).
#'
#'
#' @inheritParams plotDynamics
#' @return A plot
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],  [getPredMort()]
#' @examples
#' \donttest{
#' params <- suppressMessages(newMultispeciesParams(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotPredMort(sim)
#' plotPredMort(sim, time_range = 10:20)
#' }
plotPredMort <- function(object, species = NULL,
                         time_range,
                         highlight = NULL, ...) {
    if (is(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range  <- max(as.numeric(dimnames(object@n)$time))
        }
        params <- object@params
    } else {
        params <- validParams(object)
    }
    pred_mort <- getPredMort(object, time_range = time_range, drop = FALSE)
    # If a time range was returned, average over it
    if (length(dim(pred_mort)) == 3) {
        pred_mort <- apply(pred_mort, c(2, 3), mean)
    }

    # selector for desired species
    # Set species if missing to list of all non-background species
    if (is.null(species)) {
        species <- dimnames(params@initial_n)$sp[!is.na(params@A)]
    }
    # Need to keep species in order for legend
    species_levels <- c(as.character(params@species_params$species),
                        "Background", "Resource", "Total")
    pred_mort <- pred_mort[as.character(dimnames(pred_mort)[[1]]) %in% species, , drop = FALSE]
    plot_dat <- data.frame(value = c(pred_mort),
                           Species = factor(dimnames(pred_mort)[[1]],
                                            levels = species_levels),
                           w = rep(params@w, each = length(species)))
    p <- ggplot(plot_dat) +
        geom_line(aes(x = w, y = value, colour = Species,
                      linetype = Species, size = Species))

    linesize <- rep(0.8, length(params@linetype))
    names(linesize) <- names(params@linetype)
    linesize[highlight] <- 1.6
    p <- p +
        scale_x_continuous(name = "Size [g]", trans = "log10") +
        scale_y_continuous(name = "Predation mortality [1/year]",
                           limits = c(0, max(plot_dat$value))) +
        scale_colour_manual(values = params@linecolour) +
        scale_linetype_manual(values = params@linetype) +
        scale_size_manual(values = linesize)
    return(p)
}


#' Plot trait evolution through time
#'
#' @description
#'
#' @inheritParams plotDynamics
#' @return A plot
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],
#' @examples
#' \donttest{
#' }

plotevoTrait <- function(object, SpIdx = NULL, Normalisation = F, returnData = T, traitID = NULL)
{
    cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind

    if(is.null(SpIdx)) SpIdx <- unique(object@params@species_params$lineage)
    TimeMax <- dim(object@n)[1]
    if(is.null(traitID)) traitID <- "w_mat"

    Wmat <- object@params@species_params$w_mat[SpIdx] # maturation size to only select mature phenotypes

    # if (save_it & is.null(dir))  # if want to save but not specified where
    #   dir = paste(getwd(),"/temporary",sep = "")
    #
    # ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir)), FALSE) #create the file if it does not exists


    # Set-up the data right, we need the right apparition, extinction, traits, ... it is done dynamically mouahahah
    SumPar = object@params@species_params #shortcut
    # TT = cbind(SumPar$species,as.numeric(SumPar$ecotype),dt*SumPar$pop,dt*SumPar$extinct,SumPar[traitID]) # weird things happen without the as.numeric / *0.1 because dt / need to update this to have the traits as arguments
    TT = data.frame("Species" = SumPar$lineage, "Phenotype" =  as.numeric(SumPar$species),"Apparition" = SumPar$pop,
                    "Extinction" = SumPar$ext,SumPar[traitID]) # weird things happen without the as.numeric / *0.1 because dt / need to update this to have the traits as arguments

    # colnames(TT) = c("Species","Phenotype","Apparition","Extinction","Td","EdInt")
    # rownames(TT) = rownames(SumPar)
    # TT = TT[order(TT[,1],decreasing=FALSE),]
    # TT = as.data.frame(TT) # I use TT later in the graphs (its like an artefact)

    for (i in 1:dim(TT)[1]) if (TT$Extinction[i] == 0) TT$Extinction[i] = TimeMax # Fill the extinction value for the non-extinct species

    # Weighted mean of trait value
    # 1) matrix of summed abundance of mature individuals at each time step
    biomass = object@n
    # put 0 in object@n when w < w_mat
    for (iTime in 1:dim(biomass)[1]) # for each time step
    {
        for (iPhen in 1:dim(biomass)[2]) # for each ecotypes
        {
            w_lim = SumPar$w_mat[iPhen] # get the maturation size of the ecotype
            S <- numeric(length(object@params@w))
            S[sapply(w_lim, function(x) which.min(abs(x - object@params@w)))] <- 1 # find what w bin is the closest of the maturation size
            NoW_mat = which(S == 1) # what is the col number of that size bin
            biomass[iTime,iPhen,1:NoW_mat-1] <-0 # everything under this value become 0
        }
    }

    abundanceM = apply(biomass, c(1,2),sum) # sum the abundance left

    # 2) normalisation per species
    colnames(abundanceM) = TT$Species # phenotypes from same species have the same name
    abundanceM[is.na(abundanceM)] <- 0 # don't want any NA
    abundanceNormal = matrix(0,nrow = dim(abundanceM)[1], ncol = dim(abundanceM)[2])

    # I am getting rid of the species which went instinct at the begining and that went extinct without having mutants (no trait variation)
    # SpIdx = NULL
    # for (i in unique(SumPar$species))
    #   if (sum(abundanceM[,i]) != 0 & dim(SumPar[SumPar$species == i,])[1] != 1) # if not extinct at the beginning and more than one ecotype (for the apply)
    #     SpIdx = c(SpIdx,i)

    # SpIdx is annoying:
    # if (length(SpIdx) > length(unique(SumPar$species))) SpIdx = unique(SumPar$species) # ok so I might have species not even reaching this point so I'm short cutting spidx automatcaly

    for (iSpecies in SpIdx)
    {
        abundanceSp = abundanceM # save to manipulate
        abundanceSp[,which(colnames(abundanceM) != iSpecies)] = 0 # make everything but the targeted species to go 0 to have correct normalisation
        abundanceSp = sweep(abundanceSp,1,apply(abundanceSp,1,sum),"/") # normalise
        abundanceSp[is.nan(abundanceSp)] <-0 # when I divide by 0 I get nan
        abundanceNormal = abundanceNormal + abundanceSp # I just need to add them up to get the final matrix
    }
    # if (comments) cat(sprintf("Abundance normalised\n"))

    # Now I have normalised abundance, I need to apply the trait I want to plot on them


    plotList <- vector(mode = "list", length = length(traitID))
    dataList <- vector(mode = "list", length = length(traitID))

    # Trait

    for(iTrait in traitID)
    {
        plotStore <- list()
        dataStore <- NULL
        abundanceT = sweep(abundanceNormal,2,as.matrix(SumPar[iTrait]),"*") # I use the normalised abundance and multiply by the trait value

        # Calculate mean at each time step
        TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = max(as.numeric((unique(TT$Species)))),
                         dimnames = list(rownames(abundanceT),as.character(seq(1,max(as.numeric(unique(TT$Species))))))) #sort(unique(SumPar$species))[length(unique(SumPar$species))]
        names(dimnames(TotMean)) = list("time","species")

        for (i in SpIdx)
        {
            AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
            if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
            TotMean[,i] = AMean
        }

        # Calculate variance and standard deviation
        # it is the sum of the difference between value and mean squared and multiplied by the weight

        statData = list() # list with the stats of all species as I need to do this for each species separatly

        for (i in SpIdx)
        {
            meanSp = TotMean[,i] # take the mean of the species
            traitSp = dplyr::filter(TT, Species == i)[iTrait] # take the traits of the ecotypes in the species
            weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
            statMat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0,0), ncol = 5,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","sd","percentMean","percentSd"))) # initialise the matrix
            for (itime in 1:length(meanSp)) # for each time step
            {
                if (is.null(dim(weightSp))) {variance = sum(((traitSp-meanSp[itime])^2)*weightSp[itime])} else {variance = sum(((traitSp-meanSp[itime])^2)*weightSp[itime,])} # calculate the variance, condition if only one phen
                statMat[itime,3] <- variance
                statMat[itime,4] <- (meanSp[itime] - traitSp[1,])/traitSp[1,] # normalisation of mean
                statMat[itime,5] <- statMat[itime,3]/traitSp[1,] # normalised sd
            }
            statData[[i]] = as.data.frame(statMat) # put in the list
        }

        # I have the stats for every species, just need to plot now


        for (i in SpIdx)
        {
            stat = statData[[i]] # take the stats of the right species

            # phenotype = TT[TT$Species == i,c("Apparition","Extinction",iTrait)] # recup the traits time values
            # Phen = melt(phenotype,iTrait) # make the dataframe
            # #name = paste("Maturation size of species ",i, sep = "")
            #
            # #short cut the data frame when species does not reach end of simulation
            # # if (sum(which(Phen$value == TimeMax)) == 0) # if there is at least one value equal to the end of the sim it means that the species survived until then
            # #   stat = stat[-which(stat$mean == 0),] # first occurence of mean = 0, meaning dead
            #
            # name = paste("Species",i, sep = " ")
            #
            # # prepare the data for the ribbon
            # g1 <- ggplot(stat)+
            #   geom_smooth(aes(x = time, y = percentMean-percentSd)) +
            #   geom_smooth(aes(x= time, y = percentMean+percentSd))
            #
            # gg1 <- ggplot_build(g1)
            # dfRibbon <- data.frame(x = gg1$data[[1]]$x, ymin = gg1$data[[1]]$y, ymax = gg1$data[[2]]$y) #and extract the smooth data for the ribbon
            #
            # if (!Normalisation) stat$percentMean = stat$mean
            #
            # p = ggplot() +
            #   geom_smooth(data = stat, aes(x = time, y = percentMean)) +
            #   # geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.4)+
            #   # geom_hline(yintercept = 0, linetype = "dashed") +
            #   theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
            #         panel.grid.minor = element_line(colour = "grey92"), legend.position="none",
            #         legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
            #   scale_x_continuous(name = "Time", limits  = c(0,TimeMax)) +
            #   scale_y_continuous(name = name, limits = c(NA,NA)) + #, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
            #   ggtitle(paste("Trait =",iTrait))

            # data saving
            stat$species <- i
            dataStore <- rbind(dataStore,stat)

            plotStore[[i]] = NULL #p


        }

        plotList[[which(iTrait == traitID)]] <- plotStore
        dataList[[which(iTrait == traitID)]] <- dataStore
    }

    #summary plot
    plot_dat <- dataList[[1]]
    p<-  ggplot(plot_dat)+
        geom_smooth(aes(x=time,y=mean, group = as.factor(species), color = as.factor(species) )) +
        # geom_line(aes(x=time,y=PPMR + sd, group = as.factor(species), color = as.factor(species) ), linetype = "dashed", alpha = .5) +
        # geom_line(aes(x=time,y=PPMR - sd, group = as.factor(species), color = as.factor(species) ),linetype = "dashed", alpha = .5) +
        scale_x_continuous(name = "Time in years", limits = c(NA,NA))+
        scale_y_continuous(trans = "log10")+
        scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"), legend.position="bottom",
              strip.background = element_blank(),
              #legend.justification=c(1,1),
              legend.key = element_rect(fill = "white"))+
        guides(color = guide_legend(nrow=1)) +
        # geom_text(size = 5, data = dat_text, mapping = aes(x = - Inf, y = Inf, label = label), hjust = -0.75, vjust = 1.5) +
        ggtitle(NULL)



    if(returnData)  return(dataList) else {
        print(p)
        return(plotList)}

}
