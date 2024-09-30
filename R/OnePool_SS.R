#' Steady State Function for the One-Pool First Order Soil Carbon Model with Linear Decay
#' 
#' 
#' 
#' @param parms 
#' Consists of a list containing numeric values 'const_inputs' and 'turnoverTime', where turnoverTime is in months and inputs can be either a number for static carbon inputs or a vector representing inputs over time.
#' 'ave_inputs' should generally be an average per month.
#' 
#' @return 
#' Returns unique values of the litter and soil C pools when the system is in equilibrium.
#' 
#' @examples 
#' (none yet)
#' 
#' 
OnePool_SS.fn <- function(parms){
  ans <- parms$ave_inputs * parms$turnoverTime
  return(c(Cumulative_Respiration = 0, 
           Soil_Carbon_Concentration = ans))
}