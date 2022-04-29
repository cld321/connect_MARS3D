# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
## FUNCTIONS TO COMPUTE CONNECTIVITY AND DISTANCE MATRIX   ##
## FROM EULERIAN SIMULATION OF TRACER DISPERSION-DIFFUSION ##
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #



###################### FUNCTION TO COMPUTE CONNECITIVTY MATRIX ######################

#extract connectivity matrix from subInfo 
 ## AGGREGATE CELLS IN ZONES BY SUM ##
  ## LABEL THE ROWS AND COLUMNS OF THE CONNECTIVITY MATRIX ##
  ## SOURCE ZONES AS COLUMNS / SINKS ZONES AS ROWS
  ## options: percentage and threshold
computeConnectivityMatrix <- function (subInfo, zones, threshold, percentage = T) {
    mat.connect <- aggregate(subInfo[,colnames(subInfo)[grep ('Traceur',colnames(subInfo))]], subInfo[,'name',drop=F], sum, na.rm=T)
    #create row.names in good order
    new.rownames<-paste("Traceur_", c(1:60),sep="")
    #arrange the rows in good order
    mat.connect <- mat.connect [match(new.rownames,mat.connect[,'name']),]
    mat.connect <- mat.connect [, 2:ncol(mat.connect) ] [zones,zones] #eliminate first column containing the names of sinkzones
    #mat.connect = true connectivity matrix
    #mat.percent = connectivity matrix expressed in percentage
    #subInfo extended, includes on last row total concentration of each traceur at release, extracted in masseTotal
    
    if (percentage) {
        #convert in percentage
        masseTotal <-as.numeric( subInfo[nrow(subInfo), colnames(subInfo)[grep ('Traceur',colnames(subInfo))]]) [zones]
        mat.percent <- do.call(cbind, lapply(1:ncol(mat.connect), function(x){ mat.connect[,x] * 100/ masseTotal[x] }))
        mat.percent [mat.percent < threshold] = 0 #apply threshold to data
        mat.final <- mat.percent
    }else{
        mat.connect [mat.connect < threshold] = 0 #apply threshold to data
        mat.final <- mat.connect}
    # asign the names of sinkzones
    colnames(mat.final) <- zones
    rownames(mat.final) <- zones
    return(mat.final)
}





###################### FUNCTIONS TO CALCULATE MATRIX OF DISTANCES #######################

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
   # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
   # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.

   if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
   if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
   else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
   else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
   m[tri] <- t(m)[tri]
   return(m)
}

GeoDistanceInKmMatrix <- function(df.geopoints){
   # Returns a matrix (M) of distances between geographic points.
   # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
   # (df.geopoints$lat[j], df.geopoints$lon[j]).
   # The row and column names are given by df.geopoints$name

   GeoDistanceInKMetres <- function(g1, g2){
      # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
      # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
      # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
      # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
      # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
     DistM <- function(g1, g2){
       require("Imap")
       #return()
       #g0<-data.frame(lat=48.40117,lon=-5.06079627) = zone 27, used for all 60 zones
       g0<-data.frame(lat=48.4614,lon=-5.0914)
       ifelse( g1$index > g2$index, 0, 
               ifelse( g1$index<=12 & g2$index<=12 | g1$index>12 & g2$index>12 , gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="km"),
                       gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g0$lat, lon.2=g0$lon, units="km") + gdist(lat.1=g0$lat, lon.1=g0$lon, lat.2=g2$lat, lon.2=g2$lon, units="km")))
       #return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="km")))
     }
     return(mapply(DistM, g1, g2))
   }

   n.geopoints <- nrow(df.geopoints)

   # The index column is used to ensure we only do calculations for the upper triangle of points
   df.geopoints$index <- 1:n.geopoints

   # Create a list of lists
   list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})

   # Get a matrix of distances (in metres)
   mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInKMetres), "lower")

   # Set the row and column names
   rownames(mat.distances) <- df.geopoints$zone
   colnames(mat.distances) <- df.geopoints$zone

   return(mat.distances)
}




#FUNCTION TO EXTRACT LAST CHARCATERS IN A STRING
substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
}
