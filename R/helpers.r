#' finds if  estimated locations are on or off range
#' @param x a data frame of locations object
#' @param range a SpatialPolygonDataFrame of the range
#' @return a logical vector indicating if successive locations were on or off range
#' @export

setGeneric("onrange",
           function(x, range){
               standardGeneric("onrange")
           })


setMethod("onrange",
          c(x = "data.frame",range = "SpatialPolygonsDataFrame"),
          function(x,range){
              locs <- data.frame(lon = x$longitude, lat = x$latitude)
              on = spatstat::inside.owin(locs$lon,locs$lat,
                                         spatstat::as.owin(range))
              on
          }
          )

setMethod("onrange",
          c(x = "data.frame",range = "SpatialPolygons"),
          function(x,range){
              locs <- data.frame(x = x$x, y = x$y)
              on = spatstat::inside.owin(locs$x,locs$y,
                                         spatstat::as.owin(range))
              on
          }
          )
#' Function to get transition matrix
#' @param x a matrix of states
#' @param prob logical, should probabilities be returned and not transitions.
#' @export
setGeneric("trans.matrix",
           function(x, prob){
               standardGeneric("trans.matrix")
           })
setMethod("trans.matrix",
          c(x = "matrix", prob = "logical"),
          function(x, prob){
               tt <- table( c(x[,-ncol(x)]), c(x[,-1]) )
               if(prob) tt <- tt / rowSums(tt)
               res <- matrix(tt,ncol = ncol(tt),byrow = FALSE)
               colnames(res) <- rownames(res) <- attributes(tt)$dimnames[[1]]
               return(res)
          }
          )
