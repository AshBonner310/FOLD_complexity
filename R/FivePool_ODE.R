#' Five-Pool First Order Soil Carbon Model with Linear Decay
#'
#' @param t 
#' A sequence of time values, at whatever resolution desired. This model is designed to run in a yearly timestep.
#' 
#' @param y 
#' A vector containing six named elements, CO2, 2 litter pools (Lit_metabolic, Lit_structural) and 3 soil pools (Soil_fast, Soil_slow, Soil_passive), in that order, representing the initial value for an lsoda forward run/simulation. CO2 should be zero to start as this will track the Cumulative CO2 released.
#' 
#' @param parms 
#' Consists of a list containing 'ave_inputs', 'input_to_struc', 'input_to_meta', 'input_to_fast', 'input_to_slow', 'input_to_passive', 'turnoverTime_meta', 'turnoverTime_fast', 'turnoverTime_struc', 'turnoverTime_slow', 'turnoverTime_passive','struc_to_fast', 'struc_to_slow', 'meta_to_fast', 'fast_to_slow', 'meta_to_struc', 'meta_to_slow', 'meta_to_passive', 'struc_to_meta', 'struc_to_passive', 'fast_to_meta', 'slow_to_meta', 'passive_to_meta', 'fast_to_struc', 'slow_to_struc', 'passive_to_struc', 'fast_to_slow', 'fast_to_passive', 'slow_to_fast', 'slow_to_passive', 'passive_to_fast', and'passive_to_slow', 'inputs.fn'. 
#' 'inputs' can be either a number for static carbon inputs or a vector representing inputs (mass concentration in kg per m$^2$) over time (years).
#' All 'turnoverTime_pool' are expressed in years.
#' All 'inputs_to_pool' and all 'pool_to_pool' are unitless fractions representing either allocation or a portion of decomposition flow, respectively.
#' Finally, any parms required for inputs.fn should also be included.
#' 
#' @param rel_tol 
#' default relative tolerance is set to 1e-8
#'
#' @return 
#' Returns cumulative CO2 release and the current value of the litter and soil C pools.
#' 
#' @examples 
#' (none yet)
#' 
#' 
FivePool_ODE.fn <- function(t, #time of integration
                             y, #6-valued vector in the order of Cummulative Respiration, then pools in ascending order based on turnover time (descending based on decay rate - faster-decaying pools first, so Litter-Metabolic, Litter-Structural, Soil-Fast, Soil-Slow, Soil-Passive)
                             parms, #Parameter list conatining input_type, inputs (monthly average), all allocation of inputs into pools, and all transfer rates between the pools.
                             rel_tol = 1e-10 #our default relative tolerance
){ 
  
  pools <- matrix(y[2:6], nrow=5)
  allocation_vector <- matrix(c(parms$input_to_meta, 
                                parms$input_to_struc, 
                                parms$input_to_fast, 
                                parms$input_to_slow, 
                                parms$input_to_passive),
                              nrow = 5
  )
  
  decay_matrix <- diag(x= 1/c(parms$turnoverTime_meta,
                              parms$turnoverTime_struc,
                              parms$turnoverTime_fast,
                              parms$turnoverTime_slow,
                              parms$turnoverTime_passive
  )
  )
  
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
  
  
  if(! all(c('ave_inputs', 'inputs.fn',  'input_to_struc', 'input_to_meta', 'input_to_fast', 'input_to_slow', 'input_to_passive', 'turnoverTime_meta', 'turnoverTime_fast', 'turnoverTime_struc', 'turnoverTime_slow', 'turnoverTime_passive','struc_to_fast', 'struc_to_slow', 'meta_to_fast', 'fast_to_slow', 'meta_to_struc', 'meta_to_slow', 'meta_to_passive', 'struc_to_meta', 'struc_to_passive', 'fast_to_meta', 'slow_to_meta', 'passive_to_meta', 'fast_to_struc', 'slow_to_struc', 'passive_to_struc', 'fast_to_slow', 'fast_to_passive', 'slow_to_fast', 'slow_to_passive', 'passive_to_fast', 'passive_to_slow' ) %in% names(parms))){
    stop('You have a parameter missing that this function needs.')
  }
  
  u <- parms$inputs.fn(time=t, parms=parms)
  
  ###ODE###
  
  ans <- u * allocation_vector - transfer_matrix %*% decay_matrix %*% pools
  
  resp <- transfer_matrix %*% decay_matrix %*% pools
  total_resp <- sum(resp)
  
  #if we're not doing a zero-input run...
  if(u != 0){
    #check against a relative tolerance to ensure conservation of mass
    if(abs(u - (total_resp + sum(ans)))/u > rel_tol){
      #relative tolerance allowance
      stop('Conservation of mass does not hold')
    }
  }

  return(list(c(
    Cumulative_Respiration = total_resp,
    Lit_metabolic = ans[1],
    Lit_structural = ans[2],
    Soil_fast = ans[3],
    Soil_slow = ans[4],
    Soil_passive = ans[5]
  )
  )
  )
}