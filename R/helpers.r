#' Function to find if  estimated locations are on or off range
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
              locs = data.frame(lon = x$longitude, lat = x$latitude)
              on = spatstat::inside.owin(locs$lon,locs$lat,maptools::as.owin.SpatialPolygons(range))
              on
          }
          )

setMethod("onrange",
          c(x = "data.frame",range = "SpatialPolygons"),
          function(x,range){
              locs = data.frame(x = x$x, y = x$y)
              on = spatstat::inside.owin(locs$x,locs$y,maptools::as.owin.SpatialPolygons(range))
              on
          }
          )
#' Function to change a SpatialLinesDataFrame to a SpatialPolygonsDataFrame
#' @details This function is based on code found in the book Applied Spaatial Data
#' Analysis with R, 2nd edition, Roger S. Bivand, Edzer J. Pebesma and V. Gomez-Rubio
#' UseR! Series, Springer New York (section 2.6 pp 41)
#' @param x An object of class SpatialLinesDataFrame
#' @return A SpatialPolygonsDataFrame

setGeneric("sl2sp",
           function(x){
               standardGeneric("sl2sp")
           })
setMethod("sl2sp",
          c(x = "SpatialLinesDataFrame"),
            function(x) {
                lns <- slot(x,"lines")
                if(sum(table(sapply(lns,function(y)length(slot(y,"Lines")))))> 1) {
                    i <- sapply(lns,function(y) {
                        crds <- slot(slot(y,"Lines")[[1]],"coords")
                        identical(crds[1,],crds[nrow(crds),])
                    })
                    i2 <- x[i]
                    list_of_Lines <- slot(i2,"lines")
                    sp <- SpatialPolygons(lapply(list_of_Lines, function(y) {
                        Polygons(list(Polygon(slot(slot(y,"Lines")[[1]],
                                                   "coords"))),ID = slot(y,"ID"))
                    }), proj4string = CRS(proj4string(x)))
                }else{
                    sp <- SpatialPolygons(lapply(lns, function(y) {
                        Polygons(list(Polygon(slot(slot(y,"Lines")[[1]],
                                                   "coords"))),ID = slot(y,"ID"))
                    }), proj4string = CRS(proj4string(x)))
                }
                spdf <- as(sp, "SpatialPolygonsDataFrame")
                return(spdf)
            })
#' Function to make a SpatialPolygonsDataFrame from a matrix of points
#' @param points a 2 column matrix of points. The last point should be equal to the first
#' @param name a character name for the SPDF

setGeneric("coords2spdf",
           function(points, name){
               standardGeneric("coords2spdf")
           })

setMethod("coords2spdf",
          c(points = "matrix", name = "character"),
          function(points,name){
              p <- Polygon(points)
              polys <- SpatialPolygons(list(
                  Polygons(list(p), name)))
              res <- as(polys, "SpatialPolygonsDataFrame")
              return(res)
          }
          )

#' calculates a sequence matrix from consecutive TRUE/FALSE vector
#' @param x a logical vector of consecutive occurances
#' @return a matrix with transition probabilities 
setGeneric("seqMat",
           function(x){
               standardGeneric("seqMat")
           })
setMethod("seqMat",
          c(x = "logical"),
          function(x){
               TT <- length(which(x[1:(length(x)-1)]=="TRUE"&x[2:(length(x))]=="TRUE"))
               TF <- length(which(x[1:(length(x)-1)]=="TRUE"&x[2:(length(x))]=="FALSE"))
               FF <- length(which(x[1:(length(x)-1)]=="FALSE"&x[2:(length(x))]=="FALSE"))
               FT <- length(which(x[1:(length(x)-1)]=="FALSE"&x[2:(length(x))]=="TRUE"))
               mat <- matrix(c(FF,FT,TF,TT),ncol = 2,byrow = TRUE)
               colnames(mat) <- rownames(mat) <- c("FALSE","TRUE")
               return(prop.table(mat,margin = 1))
          }
          )
#' Function to turn off diagonal 2 state transition matrix elements into a transition probability matrix
#' @param q a numeric vector of length two of the off diagonal transition matrix (2 state)
#' @param t the time at which to calculate, by default is 1
#' @export
setGeneric("q2p",
           function(q, t){
               standardGeneric("q2p")
           })
setMethod("q2p",
          c(q = "numeric", t = "numeric"),
          function(q, t){
              p11 <- (q[1] * exp(t*(-q[2]-q[1])) + q[2])/sum(q)
              p22 <- (q[2] * exp(t * (-q[2] - q[1])) + q[1])/sum(q)
              p12 <- 1- p11
              p21 <- 1- p22
              Pt <- matrix(c(p11,p12,p21,p22),ncol = 2, byrow = TRUE)
              return(Pt)
          }
          )
setMethod("q2p",
          c(q = "numeric", t = "missing"),
          function(q, t ){
              p11 <- (q[1] * exp((-q[2]-q[1])) + q[2])/sum(q)
              p22 <- (q[2] * exp((-q[2] - q[1])) + q[1])/sum(q)
              p12 <- 1- p11
              p21 <- 1- p22
              Pt <- matrix(c(p11,p12,p21,p22),ncol = 2, byrow = TRUE)
              return(Pt)
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
               rownames(res) <- attributes(tt)$dimnames[[1]]
               colnames(res) <- attributes(tt)$dimnames[[2]]
               return(res)
          }
          )
#' Function to extract the log likelihood from the fitted models
#' @param fit object of class \code{mmre} returned by \code{fit.mmre}
#' @export
setGeneric("ll",
           function(fit){
               standardGeneric("ll")
           })
setMethod("ll",
          c(fit = "mmre"),
          function(fit){
              pars <- rbind(get.params(fit, FALSE),get.params(fit, TRUE))[,1]
              nll <- fit@fit$fn(pars)
              -nll
          }
          )

#' Function to get the transition propabilities from a fitted model
#' @param fit object of class \code{mmre} returned by \code{fit.mmre}
#' @param t the time index at which to calculate the probability matrix
#' @export
setGeneric("get.probs",
           function(fit, t){
               standardGeneric("get.probs")
           }
           )
setMethod("get.probs",
          c(fit = "mmre", t = "numeric"),
          function(fit,t){
              idx <- which(rownames(fit@sdreport) == "log_q")
              if(length(idx)==0){stop("no yet implemented for model with covariates")}
              qs <- fit@sdreport[idx,1]
              q2p(exp(qs),t)
          }
          )

#' Function to extract either the fixed or random parameters from model
#' @param fit object of class \code{mmre} returned by \code{fit.mmre}
#' @param random logical if TRUE will return random paramter estimates
#' @export

setGeneric("get.params",
           function(fit, random){
               standardGeneric("get.params")
           }
           )
setMethod("get.params",
          c(fit = "mmre", random = "logical"),
          function(fit,random){
              if(!random){
                  idx <- which(!(names(fit@parameters)%in%c("u","log_sig_u", "cov_par")))
              }else{
                  idx <- which(names(fit@parameters)%in%c("u","log_sig_u", "cov_par"))
              }
              nam <- names(fit@parameters)[idx]
              idx <- which(rownames(fit@sdreport)%in%nam)
              return(fit@sdreport[idx,])
          }
          )
#' Function to calculate AIC of mmre fitted model
#' @param fit object of class \code{mmre} returned by \code{fit.mmre}
#' @export
setGeneric("get.AIC",
           function(fit){
               standardGeneric("get.AIC")
           })
setMethod("get.AIC",
          c(fit = "mmre"),
          function(fit){
              ll <- ll(fit)
              k <- length(fit@fit$par)
              aic <- 2*(k - ll)
              aic
          }
          )
#' Function to do a likelihood ratio test between two models, null hypothesis
#' that there is no affect of alternative model \code{fit.alt}.
#' @param fit.full object of class \code{mmre} returned by \code{fit.mmre} for full model
#' @param fit.alt object of class \code{mmre} returned by \code{fit.mmre} for alternative model
#' @export
setGeneric("lr.test",
           function(fit.full,fit.alt){
               standardGeneric("lr.test")
           }
           )
setMethod("lr.test",
         c(fit.full = "mmre", fit.alt = "mmre"),
         function(fit.full, fit.alt){
             lrts <- 2*(ll(fit.full) - ll(fit.alt))
             k.full <- length(fit.full@fit$par)
             k.alt <- length(fit.alt@fit$par)
             k <- k.full - k.alt
             pval <- 1 - pchisq(lrts,k)
             cat(paste("H0: no altervnative effect. Pvalue = ",pval),"\n")
         }
         )


#' Imports
#' @importFrom sp SpatialPolygons Polygons proj4string as SpatialPolygonsDataFrame
#' @importFrom spatstat inside.owin
#' @importFrom msm msm
#' @importFrom ggmap get_map ggmap
#' @importFrom ggplo2 fortify geom_polygon aes geom_path geom_point geom_histogram ylab ggplot xlab xlim ylim annotation_custom layer_scales annotate scale_colour_manual geom_polygon geom_raster geom_contour coord_map geom_line geom_rect scale_x_discrete theme element_text
#' @importFrom gridExtra tableGrob grid.arrange
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom stats optim rlnorm rnorm
#' @importFrom expm expm
#' @importFrom mvtnorm mvnorm
#' @importFrom maptools as.owin.SpatialPolygons
#' @import TMB Rcpp 
#' @useDynLib mmre
