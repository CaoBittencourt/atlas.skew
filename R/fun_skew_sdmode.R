# # [SETUP] ----------------------------------------------------------------
# # - Packages ---------------------------------------------------------------
# pkg <- c(
#   'modeest' #Mode
#   , 'Hmisc' #Weighted variance
# )
# 
# # Activate / install packages
# lapply(pkg, function(x)
#   if(!require(x, character.only = T))
#   {install.packages(x); require(x)})
# 
# # Package citation
# # lapply(pkg, function(x)
# #   {citation(package = x)})

# [FUNCTIONS] -------------------------------------------------------------
# - Sd-adjusted mode function --------------------------------------------
fun_skew_sdmode <- function(
    dbl_var
    , dbl_weights = NULL
    , dbl_scale_lb
    , dbl_scale_ub
    , dbl_discount = 0.25
    , lgc_sample_variance = F
){
  
  # Arguments validation
  stopifnot(
    "'dbl_var' must be numeric." = 
      is.numeric(dbl_var)
  )
  
  stopifnot(
    "'dbl_weights' must be either NULL or a numeric vector the same length as 'dbl_var'." = 
      any(
        all(
          is.numeric(dbl_weights)
          , length(dbl_weights) ==
            length(dbl_var)
        )
        , is.null(dbl_weights)
      )
  )
  
  stopifnot(
    "'dbl_scale_lb' must be numeric." = 
      is.numeric(dbl_scale_lb)
  )
  
  stopifnot(
    "'dbl_scale_ub' must be numeric." = 
      is.numeric(dbl_scale_ub)
  )
  
  stopifnot(
    "'dbl_discount' must be a numeric value between 0 and 1." = 
      all(
        is.numeric(dbl_discount)
        , dbl_discount >= 0
        , dbl_discount <= 1
      )
  )
  
  stopifnot(
    "'lgc_sample_variance' must be either TRUE or FALSE." = 
      all(
        is.logical(lgc_sample_variance)
        , !is.na(lgc_sample_variance)
      )
  )
  
  # Coerce bounds to one element
  dbl_scale_lb[[1]] -> dbl_scale_lb
  dbl_scale_ub[[1]] -> dbl_scale_ub
  
  # Prevent division by zero by increasing the scale by 1
  if(dbl_scale_ub == 0){

    dbl_scale_lb + 1 -> dbl_scale_lb
    dbl_scale_ub + 1 -> dbl_scale_ub
    dbl_var + 1 -> dbl_var

  }
  
  # Coerce the data into the bounds
  pmax(dbl_var, dbl_scale_lb) -> dbl_var
  pmin(dbl_var, dbl_scale_ub) -> dbl_var
  
  # Calculate upper limit for standard deviation
  sd(c(
    dbl_scale_lb
    , dbl_scale_ub
  )) -> dbl_sd_ub
  
  if(!lgc_sample_variance){
    
    # Population standard deviation
    dbl_sd_ub / sqrt(2) -> dbl_sd_ub
    
  }
  
  # Minimal weights
  if(length(dbl_weights)){
    
    dbl_weights / 
      min(
        dbl_weights
        , na.rm = T
      ) -> dbl_weights
    
    round(
      dbl_weights
    ) -> dbl_weights
    
  }
  
  # Calculate standard deviation
  # Sample standard deviation
  sqrt(wtd.var(
    dbl_var
    , dbl_weights
  )) -> dbl_sd
  
  if(is.na(dbl_sd)){
    
    0 -> dbl_sd
    
  }
  
  # Population standard deviation
  if(!lgc_sample_variance){
    
    dbl_sd * 
      sqrt((
        length(dbl_var) - 1
      ) / length(dbl_var)
      ) -> dbl_sd
    
  }
  
  # If weights are provided, repeat each element 'dbl_weights' times
  if(length(dbl_weights)){
    
    rep(
      dbl_var
      , times = 
        dbl_weights
    ) -> dbl_var
    
  }
  
  # Calculate mode
  unique(mlv(
    x = dbl_var
    , method = 'shorth'
  )) -> dbl_mode
  
  # Calculate sd-adjusted mode
  dbl_discount +
    (1 - dbl_discount) * (dbl_mode / dbl_scale_ub) +
    -dbl_discount * (dbl_sd / dbl_sd_ub) +
    -dbl_discount * (1 - (dbl_mode / dbl_scale_ub)) *
    (1 - 2 * (dbl_sd / dbl_sd_ub)) -> 
    dbl_skew
  
  # Output
  return(dbl_skew)
  
}

# # [TEST] ------------------------------------------------------------------
# # - Sd-adjusted mode ------------------------------------------------------
# fun_skew_sdmode(
#   dbl_var = pmax(rnorm(1000, 50, sd = 15), 0)
#   , dbl_weights = runif(1000, 25000, 250000)
#   , dbl_scale_lb = 0
#   , dbl_scale_ub = 100
#   , dbl_discount = 0.25
# )
