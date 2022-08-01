block_lnlp <- function (block, lib = NULL, pred = NULL, norm = 2, method = c("simplex", 
                                                               "s-map"), tp = 1, num_neighbors = switch(match.arg(method), 
                                                                                                        simplex = "e+1", `s-map` = 0), columns = NULL, 
          target_column = 1, stats_only = TRUE, first_column_time = FALSE, 
          exclusion_radius = NULL, epsilon = NULL, theta = NULL, silent = TRUE, 
          save_smap_coefficients = FALSE) 
{
  verbose = !silent
  if (norm != 2) {
    stop("block_lnlp(): L2-norm is the only metric currently available.")
  }
  if (!is.null(epsilon)) {
    stop("block_lnlp(): epsilon exlcusion not available.")
  }
  if (is.null(dim(block))) {
    dataFrame = data.frame(Index = seq(1:length(block)), 
                           Data = block)
    columns = "Data"
    target = "Data"
  }
  else if (ncol(block) >= 2) {
    if (first_column_time) {
      dataFrame = block
    }
    else {
      Index = seq(1:nrow(block))
      dataFrame = data.frame(Index = Index, cbind(block))
      first_column_time = TRUE
    }
    if (is.numeric(target_column)) {
      target_column = target_column + 1
      target = names(dataFrame)[target_column]
    }
    else {
      target = target_column
    }
    if (is.null(columns)) {
      columns = names(dataFrame)[3:ncol(dataFrame)]
    }
  }
  if (!is.character(columns) || length(columns) > 1) {
    E = length(columns)
  }
  else {
    E = length(strsplit(trimws(columns), "\\s+")[[1]])
  }
  if (verbose) {
    print(paste("block_lnlp(): Using target", target, 
                "columns", FlattenToString(columns), "E =", 
                E))
  }
  if (is.null(lib)) {
    lib = c(1, nrow(dataFrame))
  }
  if (is.null(pred)) {
    pred = c(1, nrow(dataFrame))
  }
  if (is.null(exclusion_radius)) {
    exclusionRadius = 0
  }
  else {
    exclusionRadius = exclusion_radius
  }
  if ("simplex" %in% method) {
    if ("e+1" %in% num_neighbors || "E+1" %in% 
        num_neighbors || "e + 1" %in% num_neighbors || 
        "E + 1" %in% num_neighbors) {
      knn = 0
    }
    else {
      knn = num_neighbors
    }
    smplx = Simplex(dataFrame = dataFrame, pathOut = "./", 
                    predictFile = "", lib = lib, pred = pred, E = E, 
                    Tp = tp, knn = knn, tau = -1, exclusionRadius = exclusionRadius, 
                    columns = columns, target = target, embedded = TRUE, 
                    const_pred = TRUE, verbose = verbose, validLib = vector(), 
                    generateSteps = 0, parameterList = FALSE, showPlot = FALSE)
    if (knn == 0) {
      knn = E + 1
    }
    stats = data.frame(cols = FlattenToString(columns))
    stats = cbind(stats, ComputeStats(list(smplx), E, 0, 
                                      tp, knn, NULL))
    if (stats_only) {
      return.object = stats
    }
    else {
      return.object = list(stats = stats, model_output = smplx)
    }
  }
  else if ("s-map" %in% method) {
    if (is.character(num_neighbors)) {
      knn = 0
    }
    else {
      knn = num_neighbors
    }
    if (is.null(theta)) {
      theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 
                0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
    }
    smapList = list()
    for (i in 1:length(theta)) {
      theta.i = theta[i]
      smapList[[i]] = SMap(pathIn = "./", dataFile = "", 
                           dataFrame = dataFrame, lib = lib, pred = pred, 
                           E = E, Tp = tp, knn = knn, tau = -1, theta = theta.i, 
                           exclusionRadius = exclusionRadius, columns = columns, 
                           target = target, smapFile = "", embedded = TRUE, 
                           const_pred = TRUE, verbose = verbose, validLib = vector(), 
                           generateSteps = 0, parameterList = FALSE, showPlot = FALSE)
    }
    names(smapList) = paste0("theta", theta)
    smapListPred = lapply(smapList, function(L) {
      L$predictions
    })
    stats = data.frame(cols = FlattenToString(columns))
    stats = cbind(stats, ComputeStats(smapListPred, E, 0, 
                                      tp, knn, theta)) %>% 
      mutate(across(where(is.list),as.numeric))
    if (stats_only) {
      return.object = stats
    }
    else {
      return.object = list(stats = stats, model_output = smapListPred)
    }
    if (save_smap_coefficients) {
      smapListCoef = lapply(smapList, function(L) {
        L$coefficients
      })
      smapListCov = lapply(smapList, function(L) {
        cols = ncol(L$coefficients)
        cov(L$coefficients[, 2:cols], use = "complete.obs")
      })
      return.object[["smap_coefficients"]] = smapListCoef
      return.object[["smap_coefficient_covariances"]] = smapListCov
    }
  }
  else {
    stop(paste("block_lnlp(): Invalid method:", method))
  }
  return(return.object)
}


ComputeStats = function( PredictList, E, tau, tp, knn.E, theta ) {
  #----------------------------------------------------------------------
  # rEDM 0.7 simplex stats_only = TRUE : data.frame E rows x 16 columns
  # rEDM 0.7 s_map   stats_only = TRUE : data.frame theta rows x 17 cols
  #----------------------------------------------------------------------
  #  "E"                "tau"                 "tp"
  #  "nn"   ("theta")   "num_pred"            "rho"
  #  "mae"              "rmse"                "perc"
  #  "p_val"            "const_pred_num_pred" "const_pred_rho"
  #  "const_pred_mae"   "const_pred_rmse"     "const_pred_perc"
  #  "const_p_val"
  #---------------------------------------------------------------------
  N = length( PredictList )
  
  # Here's the redundant part...
  stats = data.frame( E = E, tau = rep( tau, N ),
                      tp = rep( tp, N ), nn = knn.E )
  
  if ( ! is.null( theta ) ) {
    stats $ theta = theta
  }
  
  numPred      = sapply( PredictList, PredictN,             simplify = TRUE )
  errors       = sapply( PredictList, PredictError,         simplify = TRUE )
  constErrors  = sapply( PredictList, PredictConstError,    simplify = TRUE )
  percent      = sapply( PredictList, PercentSameSign,      simplify = TRUE ) 
  numConstPred = sapply( PredictList, PredictConstN,        simplify = TRUE )
  constPercent = sapply( PredictList, PercentConstSameSign, simplify = TRUE )
  pvals        = sapply( PredictList, PValue,               simplify = TRUE )
  constpvals   = sapply( PredictList, PValueConst,          simplify = TRUE )
  
  errors      = as.data.frame( t( errors ) )
  constErrors = as.data.frame( t( constErrors ) )
  
  stats $ num_pred            = numPred
  stats $ rho                 = errors $ rho
  stats $ mae                 = errors $ MAE
  stats $ rmse                = errors $ RMSE
  stats $ perc                = percent
  stats $ p_val               = pvals
  stats $ const_pred_num_pred = numConstPred
  stats $ const_pred_rho      = constErrors $ rho
  stats $ const_pred_mae      = constErrors $ MAE
  stats $ const_pred_rmse     = constErrors $ RMSE
  stats $ const_pred_perc     = constPercent
  stats $ const_p_val         = constpvals
  
  return( stats )
}

#------------------------------------------------------------------------
# sapply functions for ComputeStats
#------------------------------------------------------------------------
PValue = function( df ) {
  N.pred = length( which( !is.na( df $ Predictions ) ) )
  rho    = ComputeError( df $ Observations, df $ Predictions ) $ rho
  p_val  = NA
  if ( N.pred > 3 ) {
    # pnorm( q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE )
    p_val = max( 1E-10,
                 pnorm( atanh(rho), 0.0, 1 / sqrt(N.pred - 3), FALSE, FALSE ) )
  }
  return( p_val )
}
PValueConst = function( df ) {
  N.pred = length( which( !is.na( df $ Const_Predictions ) ) )
  rho    = ComputeError( df $ Observations, df $ Const_Predictions ) $ rho
  max( 1E-10, pnorm( atanh(rho), 0.0, 1 / sqrt(N.pred), FALSE, FALSE ) )
}
PercentSameSign = function( df ) {
  # Ratio of observations to predictions with same sign
  o = df $ Observations
  p = df $ Predictions
  N = length( which( ! is.na(o) ) )
  sum( abs( sign(o) + sign(p) ), na.rm = TRUE ) / 2 / N
}
PercentConstSameSign = function( df ) {
  o = df $ Observations
  p = df $ Const_Predictions
  N = length( which( ! is.na(o) ) )
  sum( abs( sign(o) + sign(p) ), na.rm = TRUE ) / 2 / N
}
PredictError = function( df ) {
  ComputeError( df $ Observations, df $ Predictions )
}
PredictConstError = function( df ) {
  ComputeError( df $ Observations, df $ Const_Predictions )
}
PredictN = function( df ) {
  length( which( ! is.na( df $ Predictions ) ) )
}
PredictConstN = function( df ) {
  length( which( ! is.na( df $ Const_Predictions ) ) )
}
FlattenToString = function( x ) {
  # R is wonderful... is.vector( list() ) is TRUE is.list( data.frame ) TRUE
  # Test for data.frame or matrix first, then list, then vector
  # or, use class string as selector
  if ( is.data.frame( x ) || is.matrix( x ) ) {
    s = ""
    for( row in 1:nrow( x ) ) {
      s = paste( s, paste( x[row,], collapse = " " ), collapse = " " )
    }
  }
  else if ( is.list( x ) ) {
    s = paste( unlist( x ),  collapse = " " )
  }
  else if ( is.vector( x ) ) {
    s = paste( x,  collapse = " " )
  }
  else {
    s = x
  }
  
  return ( s )
}
