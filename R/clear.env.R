#'
#' @export
clear.env <- function() {

  if (!exists('pkg.env'))
    pkg.env <<- new.env()

  pkg.env$DEM.data <<- NULL
  pkg.env$MrVBF.data <<- NULL
  pkg.env$MrVBF.grid <<- NULL
  pkg.env$MrVBF.grid.asraster <<- NULL
  pkg.env$MrRTF.data <<- NULL
  pkg.env$MrRTF.grid <<- NULL
  pkg.env$MrRTF.grid.asraster <<- NULL
  pkg.env$smoothDEM.data <<- NULL
  pkg.env$smoothDEM.grid <<- NULL
}
