get.elev <- function(data, grid, debug.level=0) {

  if (debug.level>0) 
    message('Getting DEM elevation point data:');
  
  # Check enviro variables are setup
  if (!exists('pkg.env') || is.null(pkg.env))
    stop('    The environment variable are not setup. Call set.env() with paths to SAGA.');	      
   
  # Check if 'data', 'DEM' or 'newdata' are strings. If so, import the data.
  if (is.character(data)) {
    if (debug.level>0) 
      message('... Reading in point data.');
    
    fname.data <- paste(pkg.env$working.path,data,sep=.Platform$file.sep);  
    data <- read.table(fname.data, sep=',',header=T);
    sp::coordinates(data) = ~Easting + Northing;
  }
  if (is.character(grid)) {
    if (debug.level>0) 
      message('... Reading in grid and writing DEM.asc');
    
    # Build path to DEM asc file
    fname.dem <- paste(pkg.env$working.path,grid,sep=.Platform$file.sep);  

    # Read in grid
    grid = sp::read.asciigrid(fname.dem,colname='DEM');
    sp::gridded(grid)=FALSE;
    grid = data.frame(grid);
    grid= grid[,1:3];
    names(grid) = c("DEM","Easting","Northing");
    filt=!is.na(grid[,"DEM"]);
    sp::coordinates(grid) = ~Easting + Northing;
    sp::gridded(grid)=TRUE;
    sp::fullgrid(grid)=TRUE;
    
    fname.dem <- paste(pkg.env$working.path,'DEM.asc',sep=.Platform$file.sep);  
  } else {
    # Build path to DEM asc file
    fname.dem <- paste(pkg.env$working.path,'DEM.asc',sep=.Platform$file.sep);  
  }      
          
  # Interpolate DEM to the required points. 
  if (!is.null(data) && (is.null(pkg.env$DEM.data) || length(pkg.env$DEM.data) != length(data))) {
    
    # Infill NA values of grid by taking the local average. This was essential to ensure 
    # DEM values at fixed head points beyond the mappng area (eg coastal points
    # with a fixed head of zero just beyond the DEM extent)
    dem.asRaster = raster::raster(grid,layer='DEM');
    data$DEM = NULL;
    if (debug.level>0) 
      message('... Interpolating grid to point locations.');    
    
    tmp = raster::extract(dem.asRaster, data, method='bilinear');    

    data$tmpName =  tmp;
    ncols = length(names(data))
    names(data)[ncols] = 'DEM'    
    
    # Append interpolated points to enviro variable
    if (is.null(pkg.env$DEM.data) || length(data) != length(pkg.env$DEM.data)) {      
  
      Easting = sp::coordinates(data)[1];
      Northing = sp::coordinates(data)[2];
      tmp = data.frame(Easting, Northing, tmp);
      names(tmp) = c('Easting','Northing','elev')
      sp::coordinates(tmp) = ~Easting + Northing
      pkg.env$DEM.data = tmp;	  	
	  
    } else {
      pkg.env$DEM.data$tmp = tmp;
    }
    ncols = length(names(pkg.env$DEM.data))
    names(pkg.env$DEM.data)[ncols] = 'DEM'
    rm(tmp);    

    return(data);
    
  } else {  
  
    data$DEM = pkg.env$DEM.data[['DEM']]	
    return(data);    
    
  }
}
