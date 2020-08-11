##########################################
####  Calculate hypervolume overlap
##########################################
#### | Project name: Bioclimatic envelopes
#### | Creator: Mirza Cengic
#### | Contact: mirzaceng@gmail.com
##########################################

#' Calculate overlap percentage between hypervolumes
#'
#' @param x - Hypervolume set
#'
#' @return Numeric
#' @export
#'
#' @examples
hypervolume_overlap_percent <- function(x)
{
  stopifnot(inherits(x, "HypervolumeList"))
  # calculate overlap percentage
  union_vol <- x@HVList$Union@Volume 
  intersect_vol <- x@HVList$Intersection@Volume
  
  pct_overlap <- (intersect_vol /union_vol) * 100
  
  #
  return(pct_overlap)
}

