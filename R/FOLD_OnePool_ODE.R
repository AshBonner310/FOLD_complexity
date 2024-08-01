#' One-Pool First Order Soil Carbon Model with Linear Decay
#' 
#' Describes the model:
#' 
#' $$ \frac{dC}{dt} = u(t) - kC$$
#' 
#' where $u(t)$ represents carbon inputs (which can change over time due to seasonal values, regime shifts, experiments, etc.), and the decay rate $k$ represents the inverse of turnover time $\frac{1}{\tau}$ (turnover time being the average time carbon is expected to remain in the soil before cycling out).
#' This model 
#' 
#' @param t A sequence of time values, at whatever resolution desired. This model is designed to run in a monthly timestep.
#' 
#' @param y A vector containing two named elements, CO2 and soil, in that order, representing the initial value for an lsoda forward run/simulation. CO2 should be zero to start as this will track the Cumulative CO2 released.
#' 
#' @param parms Consists of a list containing 'inputs' and 'turnoverTime', where turnoverTime is in months and inputs can be either a number for static carbon inputs or a vector representing inputs over time.
#' 
#' @param rel_tol default relative tolerance is set to 1e-8
#'
#' @return Returns cumulative CO2 release and the current value of the C pool.
#' 
#' @examples (none yet)
#' 
OnePool_ODE.fn <- function(t, 
                           y, 
                           parms, 
                           rel_tol = 1e-8){ 
  
  y <- unname(y)
  CO2 <- y[1] 
  soil <- y[2]
  
  if(! all(c('turnoverTime', 'inputs') %in% names(parms))){
    stop('You have a parameter missing that this function needs.')
  }
  
  if(typeof(parms$inputs) == 'list'){
    
    u <- Scenario_Inputs.df[t, "Carbon"]
    
  }else{
    
    u <- parms$inputs
    
  }
  
  dCO2 <- soil / parms$turnoverTime
  dSoil <- u - soil / parms$turnoverTime
  
  #if we're not doing a zero-input run...
  if(u != 0){ 
    
    # print(paste("dC02", dCO2))
    # print(paste("dSoil", dSoil))
    # print(paste("Relative Tolerance", rel_tol))
    
    #make a relative tolerance to ensure conservation of mass
    if(abs(u - (dCO2 + dSoil))/u > rel_tol){ 
      #relative tolerance allowance
      stop('Conservation of mass does not hold')
    }
  }
  return(list(c(Cumulative_Respiration = dCO2, C_total = dSoil)))
}