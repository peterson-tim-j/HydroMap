get.variogramModel <- function(
  params, 
  param.names, 
  model=c('Mat'))  {

  # DEBUGGING
  op <- options("warn")
  on.exit(options(op))
  options(warn=1)  
  
  # Check if any parameters are nan or inf. If so, then later return obj function value inf
  hasImplauibleParams = is.nan(params) | is.infinite(params)

  # Check the format of the model input.
  if (class(model)[1] == "variogramModel") {
    
    # Get model types
    variogram.type = as.character(model$model)
    
    # Remove 'nug' from list of model types and add to 'model' variable
    filt = variogram.type != 'Nug'      
    modelType = variogram.type[filt]            
  } else if (is.character(model)) {
    modelType = model    
  } else {
    stop('When calling get.variogramModel(), fit.variogram.type must equal and the model should be a gstat variogram model or a vector of model types, eg c(\'Mat\')')
  }
  
  # Initialise variogram settings
  variogram.nugg = NULL
  variogram.psill = c()
  variogram.range = c()
  variogram.kappa = c()
  variogram.ang1 = c()
  variogram.anis1 = c()
  nModel = 0
  
  # Extract variogram settings.
  for (i in 1:length(param.names)) {
    if (!is.na(pmatch('psill', param.names[i])) ||
        !is.na(pmatch('range', param.names[i])) ||
        !is.na(pmatch('kappa', param.names[i])) ||
        !is.na(pmatch('nug', param.names[i]))  ||
        !is.na(pmatch('ang1', param.names[i])) ||
        !is.na(pmatch('anis1', param.names[i]))) {
      
      if (is.na(pmatch('nug', param.names[i]))) {
        iModel = as.numeric(substr(param.names[i], nchar(param.names[i]), nchar(param.names[i])))
        
        # Initialise some extra variogram parameters in case they're not used.        
        if (iModel>nModel) {
          variogram.kappa[iModel] = 0.5
          variogram.ang1[iModel] = 0.0
          variogram.anis1[iModel] = 1.0
        }
      }
      
      if (!is.na(pmatch('nug', param.names[i]))) {
        variogram.nug = max(0, params[i])
        
      } else if (!is.na(pmatch('psill', param.names[i]))) {
        variogram.psill[iModel] = max(0, params[i])
        
        if (iModel > nModel)
          nModel = iModel
        
      } else if (!is.na(pmatch('range', param.names[i]))) {
        variogram.range[iModel] = max(1e-6, params[i])
        
        if (iModel > nModel)
          nModel = iModel
        
      } else if (!is.na(pmatch('kappa', param.names[i]))) {
        variogram.kappa[iModel] = max(0, params[i])
        
        if (iModel > nModel)
          nModel = iModel
        
      } else if (!is.na(pmatch('ang1', param.names[i]))) {
        variogram.ang1[iModel] = params[i]
        
        if (iModel > nModel)
          nModel = iModel
        
      } else if (!is.na(pmatch('anis1', param.names[i]))) {
        variogram.anis1[iModel] = params[i]
        
        if (iModel > nModel)
          nModel = iModel
        
      }
    }
  }
  
  # Build variogram
  for (i in 1:nModel) {
    
    #message(paste('... DEBGUGGING: i, model type, nug, psill, range, kappa ',i, model[i], variogram.nug, variogram.psill[i], variogram.range[i], variogram.kappa[i]))
    
    if (i == 1 && !is.null(variogram.nug)) {
      vgm.model = vgm(
        nugget = variogram.nug,
        psill = variogram.psill[i],
        range = variogram.range[i],
        model = modelType[i],
        kappa = variogram.kappa[i],
        anis = c(variogram.ang1, variogram.anis1)
      )
    } else if (i == 1) {
      vgm.model = vgm(
        psill = variogram.psill[i],
        range = variogram.range[i],
        model = modelType[i],
        kappa = variogram.kappa[i],
        anis = c(variogram.ang1, variogram.anis1)
      )
    } else {
      vgm.model = vgm(
        psill = variogram.psill[i],
        range = variogram.range[i],
        model = modelType[i],
        add.to = vgm.model,
        kappa = variogram.kappa[i],
        anis = c(variogram.ang1, variogram.anis1)
      )
    }
  }
  return(vgm.model)
}
