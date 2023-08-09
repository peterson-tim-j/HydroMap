#'
#' @export
set.env <- function(working.path=getwd(), saga.path=NULL, saga.modules=NULL, saga.lib.prefix='') {

  # Check RSAGA is loaded
  if (!require('RSAGA',quietly = T))
    stop('Please install the RSAGA package. It is required by HydroMap.')
  
  pkg.env <<- new.env()

  pkg.env$DEM <<- NULL
  pkg.env$DEM.data <<- NULL
  pkg.env$MrVBF.data <<- NULL
  pkg.env$MrVBF.grid <<- NULL
  pkg.env$MrVBF.grid.asraster <<- NULL
  pkg.env$MrRTF.data <<- NULL
  pkg.env$MrRTF.grid <<- NULL
  pkg.env$MrRTF.grid.asraster <<- NULL
  pkg.env$smoothDEM.data <<- NULL
  pkg.env$smoothDEM.grid <<- NULL
  pkg.env$working.path <<- NULL
  pkg.env$saga.settings<<- NULL
  pkg.env$saga.has.MrVFB <<- FALSE
  pkg.env$ind.subset  <<- NULL

  # Set working folder
  pkg.env$working.path <<- working.path;

  if (!is.null(saga.path) && !is.null(saga.modules)) {
    if ('RSAGA' %in% installed.packages()) {
      # Setup paths to saga_cmd and saga modules
      pkg.env$saga.settings <<- rsaga.env(path=saga.path,modules=saga.modules,lib.prefix=saga.lib.prefix,workspace=pkg.env$working.path)
      
      # Check Mr VBF module is available
      pkg.env$saga.has.MrVFB <<- rsaga.module.exists(lib = 'ta_morphometry', module='Multiresolution Index of Valley Bottom Flatness (MRVBF)', env = pkg.env$saga.settings);
      
      
      if (!pkg.env$saga.has.MrVFB) {
        message('RSAGA package does not appear to setup correctly. MrVBF and MrRTF based mapping will be unavailable.')
      }
    }
  } else {
      message('RSAGA package is not installed. Calculation of MrVBF and MrRTF grids will be unavailable.')
  }
}
