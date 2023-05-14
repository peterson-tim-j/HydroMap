get.variogram <- function(formula, data, cutoff=120000, width=500, model=c('Mat'), 
                          nmin=0, nmax=Inf, maxdist=Inf,  tol.RMSE = 1e-4, max.Its=10, 
                          fit.variogram.starts = 10, debug.level=0 ) {

  if (debug.level>0) 
    message('Building variogram model:')

  # Check model type is acceptable
  if (is.character(model)) {
    # get list of acceptable model 
    model.types = vgm();
    for (i in 1:length(model)) { 
      if (!model[i] %in% model)
	      stop(paste('The following model type is not an accepted by gstat vgm():',model[i]))
    }
  }
  
  # Check if head of depth to be est.
  if (debug.level>0) 
    message('... Checking formula terms.');
  var.names <- all.vars(formula);
  do.head.est <- any(match(var.names, 'head'),na.rm=TRUE);
  do.depth.est <- any(match(var.names, 'depth'),na.rm=TRUE);	   
  if (do.depth.est == FALSE && do.head.est ==FALSE)
    stop('   Variogram formula LHS must be "depth" or "head".')
     
  
  # Estimate initial resdiuals
  if (debug.level>0) 
    message('... Building initial regression model for residuals.');
  fit <- lm(formula, data=data);
  resid = residuals(fit);
  data$resid = resid;
  

  # Filter out NAs  
  filt_na = is.na(data$resid);
  if (any(filt_na)) {
    if (debug.level>0) 
      message('... Filtering out ', sum(filt_na), ' NA values from point data residuals.');
    dataFilt =  data[!filt_na,];    
  } else
    dataFilt =  data
  
  # Define/get initial variogram
  if (debug.level>0) 
    message('... Extracting input variogram settings.');
  if (class(model)[1] =="variogramModel") {
    filt = model$model!='Nug';
    variogram.type = model$model[filt] 
    variogram.range = model$range[filt] 
    variogram.psill = model$psill[filt]
    variogram.nugget = model$psill[!filt]
    variogram.kappa= model$kappa[filt]
    variogram.ang1= model$ang1[filt]
    variogram.anis1= model$anis1[filt]
  } else {
    if (is.character(model)) {
      variogram.type = model;
    } else
      variogram.type ='Mat';
    
    variogram.range = 25000;
    variogram.psill = 20;
    variogram.nugget = 1;
    variogram.kappa= 0.5;
  }
  
  # Fit variogam to initial residuals
  if (debug.level>0) 
    message('... Building initial variogram.');
  variogram_exp = variogram(resid ~ 1, dataFilt, cressie=FALSE, cutoff=cutoff, width=width);
  model = vgm(nugget=variogram.nugget,psill=variogram.psill, model=variogram.type, range=variogram.range, kappa=variogram.kappa);
  fit.kappa=F
  if (variogram.type=='Mat')
	  fit.kappa=T	
  if (debug.level>0) 
    message('... Fitting initial variogram.');
  model  = fit.variogram(variogram_exp, fit.sills=TRUE,  fit.ranges=TRUE, fit.kappa = fit.kappa, model = model);       

  # Update initial variogram
  if (debug.level>0) 
    message('... Extracting fitted initial variogram values.');
  variogram.type = model$model[2] 
  variogram.nugget = if(is.finite(model$psill[1])){model$psill[1]}else{1}
  variogram.range = if(is.finite(model$range[2])){model$range[2]}else{25000}
  variogram.psill = if(is.finite(model$psill[2])){model$psill[2]}else{20}
  variogram.kappa= if(is.finite(model$kappa[2])){model$kappa[2]}else{0.5}
  variogram.ang1= model$ang1[2]
  variogram.anis1= model$anis1[2]
  if (debug.level>0){
    message(paste('    type =',variogram.type))
    message(paste('    nugget =',variogram.nugget))
    message(paste('    psill =',variogram.psill))
    message(paste('    range =',variogram.range))
    message(paste('    kappa =',variogram.kappa))
  }
  
  # Return if the model is not to be refined.
  if (max.Its==0) 
    return(model)
                     
  # Filter out NAs  
  filt_na = is.na(data$resid);
  if (any(filt_na)) {
     dataFilt =  data[!filt_na,];    
  } else
     dataFilt =  data;	
    
    
  # Fit the model variogram using many initial starts and select the
  # model producing the lowest weights sum of squared error
  # --------------------------------------------------------------           

    # Initialisng data
    SSE_varg = matrix(Inf,nrow=fit.variogram.starts, ncol=5)

    # Cycle through each start. 
    if (debug.level>0) 
      message('... Doing multistart variogram fitting...');  
    for (j in 1:fit.variogram.starts) {

    	# Randomly sample all model parameters. The samples are 
    	# are +- one order from the parameter input.
    	sill = 0.1*variogram.psill + (10.0*variogram.psill - 0.1*variogram.psill) * runif(1);
    	nug = 0.1*variogram.nugget + (10.0*variogram.nugget - 0.1*variogram.nugget) * runif(1);
    	range = 0.1*variogram.range + (10.0*variogram.range - 0.1*variogram.range) * runif(1);
    	if (fit.kappa)
    	   kappa = 0.1*variogram.kappa + (10.0*variogram.kappa - 0.1*variogram.kappa) * runif(1);
    
    	# Build experimental and fit model variogram;	
    	if (fit.kappa) {
     	   model_sample = vgm(nugget=nug,psill=sill, model=variogram.type, range=range, kappa=kappa);
    	} else {
     	   model_sample = vgm(nugget=nug,psill=sill, model=variogram.type, range=range);
    	}
    	model_sample  = fit.variogram(variogram_exp, fit.sills=TRUE,  fit.ranges=TRUE, fit.kappa = fit.kappa, model = model_sample,  warn.if.neg = FALSE, debug.level=0);           
    
    	# Get the weighted sum of squares error and report to the user.
    	SSEerr = attr(model_sample, 'SSErr')
    	isSingular = attr(model_sample, 'singular')
    
    	if (!isSingular) {
    	  SSE_varg[j,1] = SSEerr;
    	  SSE_varg[j,2] = model_sample$psill[1];
    	  SSE_varg[j,3] =  model_sample$psill[2];
    	  SSE_varg[j,4] = model_sample$range[2];
    	  SSE_varg[j,5] = model_sample$kappa[2];
    	  
    	  if (debug.level>0) 
    	    message(paste('    Fit number ',j,' SSE =',SSEerr));
    	  
    	} else {
    	  if (debug.level>0) 
    	    message(paste('    Fit number ',j,' failed because of singularity.'));
    	  
    	}
  }
      
  # Find variogram with lowest SSE
  if (debug.level>0) 
    message('... Finding lowest SSE variogram fit.');
  filt = SSE_varg[,1] == min(SSE_varg[,1]);
  SSE_best = SSE_varg[filt,1]
  variogram.type = model$model[2] 
  variogram.nugget = SSE_varg[filt,2]
  variogram.psill = SSE_varg[filt,3]
  variogram.range = SSE_varg[filt,4]
  variogram.kappa = SSE_varg[filt,5]
  if (debug.level>0) {
    message(paste('    type =',variogram.type))
    message(paste('    nugget =',variogram.nugget))
    message(paste('    psill =',variogram.psill))
    message(paste('    range =',variogram.range))
    message(paste('    kappa =',variogram.kappa))
    message(paste('    SSE =',SSE_best))
    message('... Building lowest SSE variogram.');
  }
  
  model = vgm(nugget=variogram.nugget,psill=variogram.psill, model=variogram.type, range=variogram.range, kappa=variogram.kappa);

  return(model)
}
