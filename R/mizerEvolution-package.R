#' mizerEvolution: Extends the mizer package with evolution and species invasion dynamics
#'
#' @description
#' This is an extension package for the
#' [mizer package](https://sizespectrum.org/mizer/)
#' that contains additional features
#' contributed by the mizer community. These features can then be further
#' improved while being used by the community. Once matured, frequently-used
#' features can be moved to the core mizer package.
#'
#' @import mizer dplyr ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats pnorm runif time
#' @importFrom methods is
#' @importFrom grDevices colorRampPalette
#' @importFrom plyr aaply
#' @export
reshape2::melt

"_PACKAGE"

globalVariables(c("inter","value","sp","phen","Species","species","xint","grp","lineage","size","critical"))
