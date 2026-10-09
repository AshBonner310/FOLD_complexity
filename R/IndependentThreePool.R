#' Three-Pool First Order Soil Carbon Model with Linear Decay and Independent Pool Structure
#'
#' @param t 
#' A sequence of time values, at whatever resolution desired. This model is designed to run in a yearly timestep.
#' 
#' @param y 
#' A vector containing four named elements, CO2 and 3 soil pools (C_1, C_2, C_3), representing the initial value for an lsoda forward run/simulation. 
#' CO2 should be zero to start as this will track the Cumulative CO2 released.
#' 
#' @param parms 
#' Consists of a list containing 'ave_inputs', 'turnoverTime_1', 'turnoverTime_2', 'turnoverTime_3', 'decay_1', 'decay_2', 'decay_3', 'inputs.fn'. 
#' 'inputs' can be either a number for static carbon inputs or a vector representing inputs (mass concentration in kg per m$^2$) over time (years).
#' All 'turnoverTime_pool' are expressed in years.
#' All 'inputs_to_pool' are unitless fractions representing litter allocation.
#' Finally, any parms required for inputs.fn should also be included.
#' 
#' @param rel_tol 
#' default relative tolerance is set to 1e-10
#'
#' @return 
#' Returns cumulative CO2 release and the current value of the soil C pools.
#' 
#' @examples 
#' (none yet)
#' 
#' 
IndThreePool_ODE.fn <- function(t, #time of integration
                            y, #4-valued vector in the order of Cummulative Respiration, then pools in descending based on decay rate (faster-decaying pools first)
                            parms, #Parameter list containing input_type, inputs (yearly average), and all allocation of inputs into pools
                            rel_tol = 1e-10 #our default relative tolerance
){ 
  
  pools <- matrix(y[2:4], nrow=3)
  allocation_vector <- matrix(c(parms$input_1, 
                                parms$input_2, 
                                parms$input_3),
                              nrow = 3
  )
  
  decay_matrix <- diag(x= 1/c(parms$turnoverTime_1,
                            parms$turnoverTime_2,
                            parms$turnoverTime_3
                            )
                       )
  
  
  if(! all(c('ave_inputs', 'input_1', 'input_2', 'input_3', 
             'turnoverTime_1', 'turnoverTime_2', 'turnoverTime_3', 
             'inputs.fn') %in% names(parms))){
    stop('You have a parameter missing that this function needs.')
  }
  
  u <- parms$inputs.fn(time=t, parms=parms)
  
  ###ODE###
  
  ans <- u * allocation_vector - decay_matrix %*% pools
  
  resp <- decay_matrix %*% pools
  total_resp <- sum(resp)
  
  # #if we're not doing a zero-input run...
  # if(u != 0){
  #   #check against a relative tolerance to ensure conservation of mass
  #   if(abs(u - (total_resp + sum(ans)))/u > rel_tol){
  #     #relative tolerance allowance
  #     stop('Conservation of mass does not hold')
  #   }
  # }
  
  return(list(c(
    Cumulative_Respiration = total_resp,
    Soil_1 = ans[1],
    Soil_2 = ans[2],
    Soil_3 = ans[3]
  )
  )
  )
}