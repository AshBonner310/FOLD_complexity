#' Steady State Function for the Five-Pool First Order Soil Carbon Model with Linear Decay
#' 
#' 
#' 
#' 
#' @param parms 
#' Consists of a list containing 'inputs', 'input_to_struc', 'input_to_meta', 'input_to_fast', 'input_to_slow', 'input_to_passive', 'turnoverTime_meta', 'turnoverTime_fast', 'turnoverTime_struc', 'turnoverTime_slow', 'turnoverTime_passive','struc_to_fast', 'struc_to_slow', 'meta_to_fast', 'fast_to_slow', 'meta_to_struc', 'meta_to_slow', 'meta_to_passive', 'struc_to_meta', 'struc_to_passive', 'fast_to_meta', 'slow_to_meta', 'passive_to_meta', 'fast_to_struc', 'slow_to_struc', 'passive_to_struc', 'fast_to_slow', 'fast_to_passive', 'slow_to_fast', 'slow_to_passive', 'passive_to_fast', and'passive_to_slow'. 
#' 'inputs' can be either a number for static carbon inputs or a vector representing inputs (mass concentration in kg per m$^2$) over time (years).
#' All 'turnoverTime_pool' are expressed in years.
#' All 'inputs_to_pool' and all 'pool_to_pool' are unitless fractions representing either allocation or a portion of decomposition flow, respectively.
#' 
#' @return 
#' Returns unique values of the litter and soil C pools when the system is in equilibrium.
#' 
#' @examples 
#' (none yet)
#' 
#' 
FivePool_SS.fn <- function(parms){
  allocation_vector <- matrix(c(parms$input_to_meta, 
                                parms$input_to_struc, 
                                parms$input_to_fast, 
                                parms$input_to_slow, 
                                parms$input_to_passive),
                              nrow = 5)
  
  decay_matrix <- diag(x= 1/c(parms$turnoverTime_meta,
                              parms$turnoverTime_struc,
                              parms$turnoverTime_fast,
                              parms$turnoverTime_slow,
                              parms$turnoverTime_passive))
  
  
  transfer_matrix <- matrix(c(
    #metabolic to pool      #structural to pool     #fast to pool            #slow to pool           #passive to pool
    1,                      -parms$struc_to_meta,    -parms$fast_to_meta,     -parms$slow_to_meta,    -parms$passive_to_meta,
    -parms$meta_to_struc,   1,                       -parms$fast_to_struc,    -parms$slow_to_struc,   -parms$passive_to_struc,
    -parms$meta_to_fast,    -parms$struc_to_fast,    1,                       -parms$slow_to_fast,    -parms$passive_to_fast,
    -parms$meta_to_slow,    -parms$struc_to_slow,    -parms$fast_to_slow,     1,                      -parms$passive_to_slow,
    -parms$meta_to_passive, -parms$struc_to_passive, -parms$fast_to_passive,  -parms$slow_to_passive,  1
  ),
  nrow = 5,
  byrow=TRUE)
  
  ans <- solve(transfer_matrix %*% decay_matrix, parms$inputs * allocation_vector)
  
  return(c(Metabolic_Plant_Residue = ans[1],
           Structural_Plant_Residue = ans[2],
           Fast_Soil = ans[3],
           Slow_Soil = ans[4],
           Passive_Soil = ans[5]))
}