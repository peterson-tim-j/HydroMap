#' Reads in the point groundwater level observations file.
#'
#' \code{import.pointData} reads in a .csv file of the point groundwater level data.
#'
#' @details This function reads in .csv file of the point groundwater level data and
#' removing duplicates points. The coordinates must be projected.
#'
#' @param \code{data} character string for the point data .csv filename. The file must contain columns with the following names "easting", "northing" and "head".
#'
#' @param \code{min.sepDist} scalar value >0 defining the maximum distance between two observations that will result in a duplicate being detected. When a duplicate is detected then when
#' \code{orderby.head==TRUE}, the point with the greatest water level elevation will be kept ALso, the lowest elevation point will be kept. The default is 5m.
#'
#' ' @param \code{orderby.head} logical scalar defining if the highest elevation water level should be selected when a duplicate is found. The default is \code{TRUE}.
#'
#' @param \code{debug.level} Control the user messages. A value >0 outputs progress.
#'
#' @return
#' A point Spatial object  will be returned.
#'
#' @export
import.pointData <- function(data, min.sepDist = 5, orderby.head = TRUE, debug.level=0) {

	if (is.character(data)) {
	  if (debug.level>0)
	    message('... Reading in point obs. file');
	  if (dirname(data) == '.') {
	    fname.data <- paste(pkg.env$working.path,data,sep=.Platform$file.sep);
	  } else {
	    fname.data <- data;
	  }
	  data <- read.table(fname.data, sep=',',header=T);

	  # Chnage the column names to the lower casing
	  names(data) = tolower(names(data));

	  # Check the data has head or depth data
	  colnames = names(data)
	  if (orderby.head) {
	    if (!('head' %in% colnames))
	      stop('The input point data must contain a column named "head" when mapping heads.')
	  } else {
	    if (!('depth' %in% colnames))
	      stop('The input point data must contain a column named "depth" when mapping depth to water table.')
	  }

	  # Check minimum speration distance
	  if (min.sepDist<=0)
	    stop('The minimum seperation distance between points must be >0.')

	  # Removing duplicates
	  if (orderby.head) {
	    data = data[with(data, order(easting, northing, head, decreasing=TRUE)), ]
	  } else {
	    data = data[with(data, order(easting, northing, depth, decreasing=FALSE)), ]
	  }
	  coordinates(data) = ~easting + northing
	  ndata = length(data);
	  data = remove.duplicates(data, zero = min.sepDist)

	  if (debug.level>0)
	    message(paste('... Number duplicate location points removed:',ndata-length(data)))

	} else if (class(data)!='SpatialPointsDataFrame') {
	  stop('The input grid must be a file name or a SpatialPointsDataFrame')
	}

	return(data)
}
