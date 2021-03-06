#' Function to location data

#' Method to map a SpatialPolygonsDataFrame onto google maps 
#' @param x A SpatialPolygonsDataFrame with latitude and longitude projection or a data.frame with longitude and latitude columns
#' @param zoom numeric 'zoom' of google map
#' @docType methods
#' @rdname map
#' @export

setGeneric("map",
           function(x,zoom){
               standardGeneric("map")
           })


#' @rdname map
#' @aliases map,SpatialPolygonsDataFrame,SpatialPolygonsDataFrame-method
#' @export
setMethod("map",
          c(x = "SpatialPolygonsDataFrame",zoom = "numeric"),
          function(x,zoom) {
              centre = rowMeans(x@bbox)
              map = ggmap::get_map(centre,zoom = zoom, maptype = "hybrid", source = "google")
              map = ggmap::ggmap(map)
              p = ggplot2::fortify(x)
              map = map + ggplot2::geom_polygon(ggplot2::aes(x=long, y=lat, group = group),
                                                fill='grey', size=2,color='black', data=p, alpha=0)
              map
          }
          )
#' @rdname map
#' @aliases map,data.frame,data.frame-method
#' @export
setMethod("map",
          c(x = "data.frame",zoom = "numeric"),
          function(x,zoom) {
              lonlat <- data.frame(longitude = x$longitude, latitude = x$latitude)
              centre = colMeans(lonlat)
              map = ggmap::get_map(centre,zoom = zoom, maptype = "hybrid", source = "google")
              map = ggmap::ggmap(map)
              map = map + ggplot2::geom_point(ggplot2::aes(x=longitude, y=latitude), data=lonlat)
              map
          }
          )


#' Function to plot a histogram of dates
#' @param t A vector of class POSIXct
#' @docType methods
#' @rdname hist.time
#' @name hist.time
#' @export

setGeneric("hist.time",
           function(t){
               standardGeneric("hist.time")
           })

setMethod("hist.time",
          c(t = "POSIXct"),
          function(t){
             dates = data.frame(Date = as.Date(t))
             ggplot2::ggplot(dates ,ggplot2::aes(x=Date)) +
                 ggplot2::geom_histogram(binwidth=1, colour="white") +
                     ggplot2::ylab("Frequency of transmissions") +
                     ggplot2::xlab("Month and Day")
          }
          )


#' Function to plot all animal tracks in one dataset along with the associated spatial ploygon range. Colours of tracks are created to
#' be as distinct as possible.
#' @param x a data frame with named elements \code{longitude}, \code{latitude}, and \code{TagID} of animal
#' @param range a SpatialPolygonsDataFrame of the associated spatial range
#' @export
#'
setGeneric("show.all",
           function(x,range){
               standardGeneric("show.all")
           })
setMethod("show.all",
          c(x = "data.frame",range = "SpatialPolygonsDataFrame"),
          function(x,range){
              p <- ggplot2::fortify(range)
              n1 <- grep("lat",names(x))[1]
              n2 <- grep("lon",names(x))[1]
              l <- data.frame(long = x[,n2],lat = x[,n1],TagID = x$TagID)
              n <- length(table(l$TagID))
              qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
              col_vector <- unlist(mapply(RColorBrewer::brewer.pal,
                                          qual_col_pals$maxcolors, rownames(qual_col_pals)))
              ggplot2::ggplot(ggplot2::aes(x = long,y = lat,group = TagID, colour = TagID),
                              data = l) +
                  ggplot2::geom_path(ggplot2::aes(x = long,y = lat),size = 0.2, data = l) +
                  ggplot2::geom_point(ggplot2::aes(x = long,y = lat), data = l)  +
                  ggplot2::scale_colour_manual(values = col_vector) +
                  ggplot2::geom_polygon(ggplot2::aes(x=long, y=lat, group = group),
                                        fill='grey', size=2,color='black', data=p, alpha=0)
          }
          )
#' Function to show duration of exposure data
#' @param x a data frame with the following named elements: \code{Start}, start time as a date
#' object of exposure start times; \code{End}, end time as a date object of exposure end times;
#' \code{Type}, a factor vector of exposure types.
#' @export
setGeneric("show.duration",
           function(x){
               standardGeneric("show.duration")
           })
setMethod("show.duration",
          c(x = "data.frame"),
          function(x){
              x$date <- strptime(x$Start, "%Y-%m-%d")
              x$start.new <- format(x$Start, format = "%H:%M:%S")
              x$end.new <- format(x$End, format = "%H:%M:%S")
              x$day <- factor(as.POSIXct(x$date))
              levels(x$day) <- 1:length(table(x$day))
              x$day <- as.numeric(as.character(x$day))
              x$event <- x$Type
              levels(x$event) <- 1:length(table(x$event))
              x$event <- as.numeric(as.character(x$event))
              
              plot <- ggplot2::ggplot(x, ggplot2::aes(day, Start)) +
                  ggplot2::geom_rect(ggplot2::aes(ymin = Start, ymax = End,
                                         xmin = (day - 0.45) + event/10,
                                         xmax = (day - 0.35) + event/10,
                                         fill = Type))+
                  ggplot2::scale_x_discrete("Date",labels = as.POSIXct(unique(x$date)),
                                   limits = 1:length(table(x$day))) +
                  ggplot2::ylab("Duration") +
                  ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,hjust=1))
              return(plot)
          }
          )
              

#' Function to show the individual normal random effects of fitted model retturned by \code{fit.mmre}
#' @param fit object of class \code{mmre} returned by \code{fit.mmre}
#' @export
#' 
setGeneric("show.random",
           function(fit){
               standardGeneric("show.random")
           })
setMethod("show.random",
          c(fit = "mmre"),
          function(fit){
              idx <- which(rownames(fit@sdreport)=="u")
              if(length(idx)==0){stop("no random effects")}
              n <- length(idx)/2
              us <- fit@sdreport[idx,]
              ran <- c(floor(range(us[,1])[1]),ceiling(range(us[,1])[2]))
              plot(us[1:n,1],pch = 20, type = "p",ylab = "",xlab = "",xaxt = "n",
                   ylim = ran,main = "",cex = 2)
              lines(us[1:n,1],col = "black",cex = 2)
              lines(us[(n+1):(2*n),1],col = "grey",cex = 2)
              points(us[(n+1):(2*n),1],pch = 20,col = "grey",cex = 2)
              legend("topright",col = c("black","grey"),pch = 20,lty = 1,
                     legend = c("transition 1->2","transition 2->1"),bty = "n")
              axis(1,at = 1:n,labels = as.factor(fit@summary$IDs),cex.axis = 0.7,las = 2)
              abline(h=0,lwd = 2,lty = 2)
          }
         )

#' Function to plot the (off diagonal) transition probabilties over a specified time period
#' @inheritParams show.random
#' @param time a numeric vector of times
#' @export
setGeneric("show.probs",
           function(fit, time){
               standardGeneric("show.probs")
           })
setMethod("show.probs",
          c(fit = "mmre", time = "numeric"),
          function(fit, time){
              off.on = on.off = numeric(length(time))
              for(i in 1:length(time)){
                  gp = get.probs(fit, t = time[i])
                  off.on[i] = gp[1,2]
                  on.off[i] = gp[2,1]
              }
              plot(1, type = "n", ylim = c(0,1), xlim = range(time))
              lines(time, off.on, lty = 2)
              lines(time, on.off, lty = 3)
              legend("topleft", lty = c(2,3), legend = c("1--2","2--1"), bty = "n")
          }
          )

