get.objFunc <- function(
  params, 
  param.names, 
  formula, 
  grid, 
  grid.landtype.colname=NULL, 
  data, 
  data.fixedHead=NULL,
  newdata=NULL, 
  data.errvar.colname = NULL, 
  model=c('Mat'),
  nmin=0.2, 
  omax=NULL,
  fit.variogram.type=1, 
  objFunc.type=1, 
  lookup.table=NULL, 
  return.transformed.params = T, 
  return.predictions = F, 
  use.cluster=TRUE, 
  debug.level=0)  {

  # DEBUGGING
  op <- options("warn")
  on.exit(options(op))
  options(warn=1)  
  
  # Check if any parameters are nan or inf. If so, then later return obj function value inf
  hasImplauibleParams = is.nan(params) | is.infinite(params)

  # Loop through the input parameter names to update the values within the complete list of param names
  # If the lookup table is not empty, then back transform the parameters from an integer to true parameter value.
  paramsList = list(mrvbf.pslope=NULL,mrvbf.ppctl=NULL,smooth.std=NULL,nmax=NULL,nmax.fixedHead=NULL,maxdist=NULL,trendMaxDistFrac=NULL)
  params2Print = c()
  param.beyond.bounds = F;
  if (is.null(lookup.table)) {
    for (i in 1:length(param.names)) {
      # Skip if variogram variable and the variogram parameters are continuous and to be fitted.
      if (fit.variogram.type == 1 &&
          (
            !is.na(pmatch('psill', param.names[i])) ||
            !is.na(pmatch('range', param.names[i])) ||
            !is.na(pmatch('kappa', param.names[i])) ||
            !is.na(pmatch('nug', param.names[i])) ||
            !is.na(pmatch('ang1', param.names[i])) ||
            !is.na(pmatch('anis1', param.names[i]))
          )) {
        next
      }
      paramsList[[param.names[i]]] = params[i]
      params2Print[i] = paramsList[[param.names[i]]]
      
    }
  } else {
    lookup.table.names = names(lookup.table)
    # Loop through the loopup table and add default (ie non-calibrated) parameters to the list.
    # This is required to allow calibration of a partial subset of the parameters.
    for (i in 1:length(lookup.table.names)) {
      if (nrow(lookup.table[[lookup.table.names[i]]]) == 1) {
        paramsList[[lookup.table.names[i]]] = lookup.table[[lookup.table.names[i]]][1, 2];  
      }
    }
    
    # Handle the parameters to be calibrated.
    for (i in 1:length(param.names)) {
      # Skip if param is implausible. Note, this loop is still undertaken for the
      # plausible parameters to return the parameter transformation for integer parameters
      if (hasImplauibleParams[i]) {
        params2Print[i] = params[i]
        
        next
        
      }
      
      # Skip if variogram variable and the variogram parameters are continuous and to be fitted.
      if (fit.variogram.type == 1 &&
          (
            !is.na(pmatch('psill', param.names[i])) ||
            !is.na(pmatch('range', param.names[i])) ||
            !is.na(pmatch('kappa', param.names[i])) ||
            !is.na(pmatch('nug', param.names[i])) ||
            !is.na(pmatch('ang1', param.names[i])) ||
            !is.na(pmatch('anis1', param.names[i]))
          )) {
        next
      }
      
      # Assess if the parameter name is listed in the lookup table. If so, round the input
      # parameter to an integer and add the rounded value to the vector of transformed parameter values
      # to be returned to the optimisation routine
      isIntegerParam =  param.names[i] %in% lookup.table.names
      
      # Skip if the parameter is not in the lookup table (eg is a real parameter)
      if (!isIntegerParam) {
        paramsList[[param.names[i]]] = params[i]
        
        params2Print[i] = params[i]
        
        next
      }
      
      # Round parameter to an integer
      params[i] = round(params[i], 0)
      
      if (params[i] < 1 || params[i] > dim(lookup.table[[param.names[i]]])[1])
        param.beyond.bounds = T
      else {
          # Transform parameter from integer lookup value
          paramsList[[param.names[i]]] = lookup.table[[param.names[i]]][params[i], 2]
          
          params2Print[i] = paramsList[[param.names[i]]]
      }
    }
  }
  
  # Return Inf if any integer paramater is out of bounds or implausible.
  if (param.beyond.bounds || any(hasImplauibleParams)) {
    message('... Integer parameters are out of bounds. Returning Inf.');
    if (return.transformed.params) { 
      return(c(Inf, params))
    } else {
      return(Inf)
    } 
  }
    
  # Build variogram model
  if (fit.variogram.type == 1) {
  
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
      stop('When fit.variogram.type==1, the model should be a gstat variogram model or a vector of model types, eg c(\'Mat\')')
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
        
        params2Print[i] = params[i]
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
    model = vgm.model
    
    
    # Turn off fitting within get.variogram()
    fit.variogram.type = 3
    
  }

  # Call krige.head with the input parameter names.
  heads = do.call(krige.head, c(paramsList, list(formula=formula, grid=grid, grid.landtype.colname=grid.landtype.colname, data=data, data.fixedHead=data.fixedHead, newdata=newdata, data.errvar.colname=data.errvar.colname, model=model,
                nmin=nmin, omax=omax, fit.variogram.type=fit.variogram.type, use.cluster=use.cluster, debug.level=debug.level)))  

  if (return.predictions) {
    heads.to.return =  heads;
  }
  
  # Create filters for erroneous estimates
  notNA = !is.na(heads$head.pred) & !is.na(heads$head.var)
  notNAN = !is.nan(heads$head.pred) & !is.nan(heads$head.var)
  notNegVar = heads$head.var>0
  
  # Report on NAs filter
  if (any(!notNA))
    message(paste('... The following number of newdata estimates had a NA predictions or variances and have been removed from objective function calculation: ',sum(!notNA)),immediate.=T)

  # Report on NANs filter
  if (any(!notNAN))
    message(paste('... The following number of newdata estimates had a NAN predictions or variances and have been removed from objective function calculation: ',sum(!notNAN)),immediate.=T)

  # Report on variance filter
  if (any(!notNegVar))
    message(paste('... The following number of newdata estimates had a kriging variance <0 and have been removed from objective function calculation: ',sum(!notNegVar)),immediate.=T)

  # Apply filters
  heads = heads[notNA & notNAN & notNegVar,]
  newdata = newdata[notNA & notNAN & notNegVar,]
  
  # Add filter to returned data
  if (return.predictions) {
    heads.to.return$is.valid.est =  notNA & notNAN & notNegVar;
  }
  
  
  # Filter residuals to remove 1%ile and 99%ile results
  threshold_percentile = 1/100;
  filt = heads$resid/sqrt(heads$head.var) > quantile(heads$resid/sqrt(heads$head.var),threshold_percentile) & heads$resid/sqrt(heads$head.var) < quantile(heads$resid/sqrt(heads$head.var),1-threshold_percentile);
  newdata.at.tails = newdata[!filt,]
  message('... Number of outlier newdata estimates removed from objective function: ',sum(!filt),'.')
  message('    The mean and variance of the removed obs. newdata are ',mean(newdata.at.tails$head),' and ',var(newdata.at.tails$head) )  
  heads = heads[filt,];
  newdata = newdata[filt,];
  nObs = length(heads); 

  # Report summmary of est  
  if (debug.level>0) {
      message(paste('    The mean, min & max meters depth to water level is :',
                    mean(heads$DEM-heads$head.pred), min(heads$DEM-heads$head.pred), max(heads$DEM-heads$head.pred) ))    
    message(paste('    The mean, min & max meters head is :',
                  mean(heads$head.pred), min(heads$head.pred), max(heads$head.pred) ))    
    message(paste('    The mean, min & max meters DEM elevation is :',
                  mean(heads$DEM), min(heads$DEM), max(heads$DEM) ))    
    
  }              
  
  # Report the number of artesian points.
  is.artesian = heads$head.pred > heads$DEM & newdata$head < heads$DEM
  message(paste('... The number of artesian newdata estimates = ',sum(is.artesian)))
  
  # Calculate the requested objective functions
  objEst = vector(mode='numeric',length(objFunc.type))
  for (i in 1:length(objFunc.type)) {
      if (objFunc.type[i]==1) {
        # Calculate liklihood est divided by the number of obs (MLKV:mean log kriging variance). 
        # Taken from Samper and Neumen 1989 and Pardo-Iguzquiza & Dowd 2013.
        NSSE = sum((heads$resid^2)/heads$head.var);
        
        MLKV = sum(log(heads$head.var))
        
        objEst[i]= (log(2*pi) + MLKV + NSSE)/nObs
        
      } else if (objFunc.type[i]==2) {
        # As for type 1 BUT with an added penality for violoating a physical constraint.
        # For simplicity, below the constraint is simply the head being below the land surface.
        # This could be easily extended to include multiple spatially varying constraints.
        resid = heads$resid;
        var  = heads$head.var;    
        if (any(is.artesian)) {
          penalty.resid = heads$head.pred[is.artesian] - heads$DEM[is.artesian];
          resid[is.artesian] = abs(resid[is.artesian]) + abs(penalty.resid);
          message(paste('    The mean, min & max meters above the land surface from the artesian points is :',
                        sum(is.artesian), mean(penalty.resid), min(penalty.resid) ))
        }
        NSSE = sum((resid^2)/var);
        MLKV = sum(log(var))
        objEst[i]= (log(2*pi) + MLKV + NSSE)/nObs;
        
        
      } else if (objFunc.type[i]==3) {
        # RMSE
        objEst[i]= sqrt(mean(heads$resid^2));    
        
      } else if (objFunc.type[i]==4) {
        # As for type 3 BUT with an added penality for violoating a physical constraint.
        # For simplicity, below the constraint is simple the head being below the land surface.
        # This could be easily extended to include multiple spatially varying constraints.
        resid = heads$resid;
        if (any(is.artesian)) {
          penalty.resid = heads$head.pred[is.artesian] - heads$DEM[is.artesian];
          resid[is.artesian] = abs(resid[is.artesian]) + abs(penalty.resid);
          message(paste('    The mean, min & max meters above the land surface from the artesian points is :',
                        sum(is.artesian), mean(penalty.resid), min(penalty.resid) ))
        }    
        objEst[i]= sqrt(sum(resid^2)/nObs);    
        
      } else
        stop('objective function type unknown.')    
  }

    
  # Print paramaters to be analysed
  message('... Objective function results:')
  for (i in 1:length(objFunc.type)) {
    message(cat(format(c('Obj. Type',param.names,'Obj. func. val.'),width=16,justify='right')))
    message(cat(format(c(objFunc.type[i],params2Print,objEst[i]),width=16,justify='right')))  
  }
  message('                               ')
  
  if (return.transformed.params && return.predictions) { 
    return(list(obj.est =objEst, params.transformed = params, est = heads.to.return))
  } else if (!return.transformed.params && return.predictions) { 
    return(list(obj.est =objEst, est = heads.to.return))
  } else if (return.transformed.params && !return.predictions) { 
      return(c(objEst, params))
  } else {
    return(objEst)
  } 
  
    
}
