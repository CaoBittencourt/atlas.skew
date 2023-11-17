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
# - Shorth's mode function --------------------------------------------
fun_skew_mode <- function(
    dbl_var
    , dbl_weights = NULL
){

  # Arguments validation
  stopifnot(
    "'dbl_var' must be numeric." =
      is.numeric(dbl_var)
  )

  stopifnot(
    "'dbl_weights' must be either NULL or a numeric vector the same length as 'dbl_var'." =
      any(
        is.null(dbl_weights)
        , all(
          is.numeric(dbl_weights)
          , length(dbl_weights) ==
            nrow(cbind(dbl_var))
        )
      )
  )

  # Check if 'dbl_var' is named
  lgc_named <- F

  if(any(
    length(colnames(dbl_var)),
    length(names(dbl_var))
  )){

    lgc_named <- T

  }

  cbind(
    dbl_var
  ) -> dbl_var

  # If weights are provided, repeat each element 'dbl_weights' times
  if(length(dbl_weights)){

    # Minimal weights
    dbl_weights
    min(
      dbl_weights
      , na.rm = T
    ) -> dbl_weights

    round(
      dbl_weights
    ) -> dbl_weights

    dbl_var[rep(
      1:nrow(dbl_var)
      , times = dbl_weights
    ), ] -> dbl_var

  }

  # Calculate mode
  vapply(
    as.data.frame(dbl_var)
    , function(x){mlv(x, method = 'shorth')}
    , FUN.VALUE = numeric(1)
  ) -> dbl_mode

  rm(dbl_var)
  rm(dbl_weights)

  # If named, keep names
  if(!lgc_named){

    as.numeric(
      dbl_mode
    ) -> dbl_mode

  }

  rm(lgc_named)

  # Output
  return(dbl_mode)

}

# - Sd-adjusted mode function --------------------------------------------
fun_skew_sdmode <- function(
    dbl_var
    , dbl_weights = NULL
    , dbl_scale_lb = 0
    , dbl_scale_ub = 100
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
        is.null(dbl_weights)
        , all(
          is.numeric(dbl_weights)
          , length(dbl_weights) ==
            nrow(cbind(dbl_var))
        )
      )
  )

  stopifnot(
    "'dbl_scale_lb' must be numeric." =
      is.numeric(dbl_scale_lb)
  )

  stopifnot(
    "'dbl_scale_ub' must be numeric and greater than 'dbl_scale_lb'." =
      all(
        is.numeric(dbl_scale_ub)
        , dbl_scale_ub > 
          dbl_scale_lb
      )
    
  )
  
  stopifnot(
    "'lgc_sample_variance' must be either TRUE or FALSE." =
      all(
        is.logical(lgc_sample_variance)
        , !is.na(lgc_sample_variance)
      )
  )

  # Data wrangling
  dbl_scale_lb[[1]] -> dbl_scale_lb
  dbl_scale_ub[[1]] -> dbl_scale_ub

  # Check if 'dbl_var' is named
  lgc_named <- F

  if(any(
    length(colnames(dbl_var)),
    length(names(dbl_var))
  )){

    lgc_named <- T

  }

  cbind(
    dbl_var
  ) -> dbl_var

  dbl_var / (
    dbl_scale_ub -
      dbl_scale_lb
  ) -
    dbl_scale_lb / (
      dbl_scale_ub -
        dbl_scale_lb
    ) -> dbl_var

  rm(dbl_scale_lb)
  rm(dbl_scale_ub)

  # Calculate standard deviation
  # Sample standard deviation
  vapply(
    as.data.frame(dbl_var)
    , function(x){sqrt(wtd.var(
      x, dbl_weights, na.rm = T
    ))}
    , FUN.VALUE = numeric(1)
  ) -> dbl_sd

  # Normalize sd
  # Calculate upper limit for standard deviation
  sd(c(
    rep(
      0, (1 + floor((nrow(dbl_var) - 1) / 2))
    ),
    rep(
      1, nrow(dbl_var) - (1 + floor((nrow(dbl_var) - 1) / 2))
    )
  )) -> dbl_sd_ub

  if(!lgc_sample_variance){

    # Population standard deviation
    dbl_sd_ub *
      sqrt((
        nrow(dbl_var) - 1
      ) / nrow(dbl_var)
      ) -> dbl_sd_ub

    dbl_sd *
      sqrt((
        nrow(dbl_var) - 1
      ) / nrow(dbl_var)
      ) -> dbl_sd

  }

  dbl_sd / dbl_sd_ub -> dbl_sd

  rm(dbl_sd_ub)
  rm(lgc_sample_variance)

  # If weights are provided, repeat each element 'dbl_weights' times
  if(length(dbl_weights)){

    # Minimal weights
    dbl_weights /
    min(
      dbl_weights
      , na.rm = T
    ) -> dbl_weights

    round(
      dbl_weights
    ) -> dbl_weights

    dbl_var[rep(
      1:nrow(dbl_var)
      , times = dbl_weights
    ), ] -> dbl_var

  }

  # Calculate mode
  vapply(
    as.data.frame(dbl_var)
    , function(x){mlv(x, method = 'shorth')}
    , FUN.VALUE = numeric(1)
  ) -> dbl_mode

  rm(dbl_var)
  rm(dbl_weights)

  # Calculate sd-adjusted mode
  0.5 * dbl_sd +
    (1 - 0.5 * dbl_sd) * dbl_mode +
    -0.5 * dbl_sd * dbl_sd +
    -0.5 * dbl_sd * (1 - dbl_mode) *
    (1 - 2 * dbl_sd) ->
    dbl_skew

  rm(dbl_mode)
  rm(dbl_sd)

  # If named, keep names
  if(!lgc_named){

    as.numeric(
      dbl_skew
    ) -> dbl_skew

  }

  rm(lgc_named)

  # Output
  return(dbl_skew)

}

# # - Demode data (skewness recenter) ----------------------------------------------
# fun_skew_desdmode <- function(
    #     df_data
#     , dbl_weights = NULL
#     , dbl_scale_lb = 0
#     , dbl_scale_ub = 100
#     , lgc_sample_variance = F
#     , dbl_pct_remove = 1
# ){
#
#   # Arguments validation
#   stopifnot(
#     "'df_data' must be a data frame or a numeric matrix." =
#       any(
#         is.data.frame(df_data)
#         , all(
#           is.matrix(df_data),
#           is.numeric(df_data)
#         )
#       )
#   )
#
#   stopifnot(
#     "'dbl_scale_lb' must be numeric." =
#       is.numeric(dbl_scale_lb)
#   )
#
#   stopifnot(
#     "'dbl_scale_ub' must be numeric." =
#       is.numeric(dbl_scale_ub)
#   )
#
#   stopifnot(
#     "'dbl_pct_remove' must be a number between zero and one." =
#       all(
#         dbl_pct_remove >= 0,
#         dbl_pct_remove <= 1
#       )
#   )
#
#   # Data wrangling
#   dbl_scale_lb[[1]] -> dbl_scale_lb
#   dbl_scale_ub[[1]] -> dbl_scale_ub
#
#   as.matrix(
#     df_data
#   ) -> mtx_data
#
#   rm(df_data)
#
#   # Conditional output
#   if(dbl_pct_remove == 0){
#
#     # Output
#     return(mtx_data)
#
#   }
#
#   # Calculate bounded variable skewness
#   fun_skew_sdmode(
#     dbl_var =
#       mtx_data
#     , dbl_weights =
#       dbl_weights
#     , dbl_scale_lb =
#       dbl_scale_lb
#     , dbl_scale_ub =
#       dbl_scale_ub
#     , lgc_sample_variance =
#       lgc_sample_variance
#   ) -> mtx_skew
#
#   rm(dbl_weights)
#   rm(lgc_sample_variance)
#
#   # Subtract sdmode from data
#   mtx_data -
#     dbl_scale_ub *
#     dbl_pct_remove *
#     rbind(mtx_skew)[
#       rep(1, nrow(mtx_data))
#     ] -> mtx_centered
#
#   rm(mtx_skew)
#   rm(dbl_scale_ub)
#
#   # Truncate data
#   pmax(
#     mtx_centered
#     , dbl_scale_lb
#   ) -> mtx_centered
#
#   rm(dbl_scale_lb)
#
#   # Normalize by col max
#   mtx_centered /
#     apply(
#       # mtx_centered, 1
#       mtx_centered, 2
#       , max
#     ) -> mtx_centered
#
#   pmax(mtx_centered, 0) ->
#     mtx_centered
#
#   pmin(mtx_centered, 1) ->
#     mtx_centered
#
#   # Adjust original data
#   mtx_data *
#     mtx_centered ->
#     mtx_data
#
#   rm(mtx_centered)
#
#   as.data.frame(
#     mtx_data
#   ) -> mtx_data
#
#   # Output
#   return(mtx_data)
#
# }

# - Demode data (skewness recenter) ----------------------------------------------
fun_skew_desdmode <- function(
    df_data
    , dbl_weights = NULL
    , dbl_scale_lb = 0
    , dbl_scale_ub = 100
    , lgc_sample_variance = F
    , dbl_pct_remove = 1
){

  # Arguments validation
  stopifnot(
    "'df_data' must be a data frame or a numeric matrix." =
      any(
        is.data.frame(df_data)
        , all(
          is.matrix(df_data),
          is.numeric(df_data)
        )
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
    "'dbl_pct_remove' must be a number between zero and one." =
      all(
        dbl_pct_remove >= 0,
        dbl_pct_remove <= 1
      )
  )

  # Data wrangling
  dbl_scale_lb[[1]] -> dbl_scale_lb
  dbl_scale_ub[[1]] -> dbl_scale_ub

  as.matrix(
    df_data
  ) -> mtx_data

  rm(df_data)

  # Conditional output
  if(dbl_pct_remove == 0){

    # Output
    return(mtx_data)

  }

  # Calculate bounded variable skewness
  fun_skew_sdmode(
    dbl_var =
      mtx_data
    , dbl_weights =
      dbl_weights
    , dbl_scale_lb =
      dbl_scale_lb
    , dbl_scale_ub =
      dbl_scale_ub
    , lgc_sample_variance =
      lgc_sample_variance
  ) -> mtx_skew

  rm(dbl_weights)
  rm(lgc_sample_variance)

  # Subtract sdmode from data,
  # adjusting for sdmode
  mtx_data -
    dbl_scale_ub *
    dbl_pct_remove *
    (1 - mtx_skew) *
    rbind(mtx_skew)[
      rep(1, nrow(mtx_data))
    ] -> mtx_centered

  rm(mtx_skew)
  rm(dbl_scale_ub)
  rm(dbl_pct_remove)

  # Truncate data
  pmax(
    mtx_centered
    , dbl_scale_lb
  ) -> mtx_centered

  rm(dbl_scale_lb)

  # Normalize by col max
  mtx_centered /
    apply(
      mtx_centered, 2
      , max
    ) -> mtx_centered

  pmax(mtx_centered, 0) ->
    mtx_centered

  pmin(mtx_centered, 1) ->
    mtx_centered

  # Adjust original data
  mtx_data *
    mtx_centered ->
    mtx_data

  rm(mtx_centered)

  as.data.frame(
    mtx_data
  ) -> mtx_data

  # Output
  return(mtx_data)

}

# - Demode data (mode recenter) ----------------------------------------------
fun_skew_demode <- function(
    df_data
    , dbl_weights = NULL
    , dbl_scale_lb = 0
    , dbl_scale_ub = 100
    , dbl_pct_remove = 1
){

  # Arguments validation
  stopifnot(
    "'df_data' must be a data frame or a numeric matrix." =
      any(
        is.data.frame(df_data)
        , all(
          is.matrix(df_data),
          is.numeric(df_data)
        )
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
    "'dbl_pct_remove' must be a number between zero and one." =
      all(
        dbl_pct_remove >= 0,
        dbl_pct_remove <= 1
      )
  )


  # Data wrangling
  dbl_scale_lb[[1]] -> dbl_scale_lb
  dbl_scale_ub[[1]] -> dbl_scale_ub

  as.matrix(
    df_data
  ) -> mtx_data

  rm(df_data)

  # Conditional output
  if(dbl_pct_remove == 0){

    # Output
    return(mtx_data)

  }

  # Calculate bounded variable skewness
  fun_skew_mode(
    dbl_var =
      mtx_data
    , dbl_weights =
      dbl_weights
  ) -> mtx_mode

  rm(dbl_weights)

  # Subtract mode from data,
  # adjusting for mode
  mtx_data -
    dbl_pct_remove *
    (
      1 -
        mtx_mode / (
          dbl_scale_ub -
            dbl_scale_lb
        ) -
        dbl_scale_lb / (
          dbl_scale_ub -
            dbl_scale_lb
        )
    ) * rbind(mtx_mode)[
      rep(1, nrow(mtx_data))
    ] -> mtx_centered

  rm(mtx_mode)
  rm(dbl_scale_ub)
  rm(dbl_pct_remove)

  # Truncate data
  pmax(
    mtx_centered
    , dbl_scale_lb
  ) -> mtx_centered

  rm(dbl_scale_lb)

  # Normalize by col max
  mtx_centered /
    apply(
      mtx_centered, 2
      , max
    ) -> mtx_centered

  pmax(mtx_centered, 0) ->
    mtx_centered

  pmin(mtx_centered, 1) ->
    mtx_centered

  # Adjust original data
  mtx_data *
    mtx_centered ->
    mtx_data

  rm(mtx_centered)

  as.data.frame(
    mtx_data
  ) -> mtx_data

  # Output
  return(mtx_data)

}

# # [TEST] ------------------------------------------------------------------
# # - Mode ------------------------------------------------------------------
# fun_skew_mode(c(1,1,0))
# fun_skew_mode(c(1,0,0))
# fun_skew_mode(diag(5))

# # - Sd-adjusted mode 1 ------------------------------------------------------
# fun_skew_sdmode(
#   dbl_var = pmax(rnorm(1000, 50, sd = 15), 0)
#   , dbl_weights = runif(1000, 25000, 250000)
#   , dbl_scale_lb = 0
#   , dbl_scale_ub = 100
# )
#
# # - Sd-adjusted mode 2 ------------------------------------------------------
# fun_skew_sdmode(
#   dbl_var =
#     pmax(
#       cbind(
#         rnorm(1000, 50, sd = 15),
#         rnorm(1000, 50, sd = 15),
#         rnorm(1000, 50, sd = 15),
#         rnorm(1000, 50, sd = 15),
#         rnorm(1000, 50, sd = 15)
#       ), 0
#     )
#   , dbl_weights = runif(1000, 25000, 250000)
#   , dbl_scale_lb = 0
#   , dbl_scale_ub = 100
#   , dbl_discount = 0.25
# )
#
# # - Sd-adjusted mode 3 ------------------------------------------------------
# pmax(
#   cbind(
#     rnorm(1000, 50, sd = 15),
#     rnorm(1000, 50, sd = 15),
#     rnorm(1000, 50, sd = 15),
#     rnorm(1000, 50, sd = 15),
#     rnorm(1000, 50, sd = 15)
#   ), 0
# ) -> dsds
#
# colnames(dsds) <- letters[1:ncol(dsds)]
#
# fun_skew_sdmode(
#   dbl_var = dsds
#   , dbl_weights = runif(1000, 25000, 250000)
#   , dbl_scale_lb = 0
#   , dbl_scale_ub = 100
#   , dbl_discount = 0.25
# )
#
# # - Demode data -----------------------------------------------------------
# fun_skew_desdmode(
#   df_data = matrix(1, 10, 10) * runif(100, min = 0, max = 100)
#   , dbl_weights = runif(10, 1000, 250000)
#   , dbl_scale_lb = 0
#   , dbl_scale_ub = 100
#   , lgc_sample_variance = F
#   , dbl_pct_remove = 1
# )
