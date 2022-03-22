# GRAPH THEORY #

source("FunctionsForConnectivity.R")
source("FunctionsForGraphs.R")
library("igraph")

# DISTANCE MATRIX ##################################
# dist.mat needed in node removal scenarios

# read any subInfo to extract geographical coordonates of zones, them remove file
subInfo <- read.table("Output60zones/subInfo2012_May_4w.txt", header=T)
df.points <- aggregate(subInfo [,c('lat','lon')],subInfo[,'zone',drop=F], mean,na.rm=T)[1:60,]
# assign the right order 
new.rownames<- c(1:60) 
df.points <- df.points [match(new.rownames,df.points$zone),]
dist.mat <- GeoDistanceInKmMatrix(df.points)
rm(subInfo, new.rownames)

# NODES ############################################
# select only zones with reefs, scale != 0
# add lat/lon and ID to them
habqual<- read.table("Data/scale_habitat_quality.txt", header=T)
habqual <- cbind(habqual, df.points[,2:3])
row.names(habqual) <- c(1:60)
habqual$zone<-as.factor(habqual$zone)

habqual <- habqual[habqual$scale != 0,]
habqual$ID <- 1:31

# select only ID of zones with reefs #
reefs <- which(read.table("Data/scale_habitat_quality.txt", header=T)$scale != 0)


# GRAPH$local $global ######################################
# read all subInfo file with dispersal values
# creates graph based on dispersal matrix and distance between connected zones
# calculate graphs metrics local=for each zone, global=for each grap/simulation
# threshold = ( 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)
filenames <- Sys.glob( path = "Output60zones/subInfo*", dirmark=T)

for (i in 1:length(filenames)){
  graph.ind <- computeGraphIndices(subInfo = read.table(filenames[[i]], header=T), 
                                   new.threshold = 0.01, zone = "all", g.strenght = "weak")
  filename = substr(filenames[[i]], 22, 22)
  if (i ==1){ 
    all.local <- cbind(zone = row.names(graph.ind$local), graph.ind$local, year=substr(filenames[[i]], 22, 25), 
                       month=substr(filenames[[i]], 27, 29), pld=substr(filenames[[i]], 31,32))
    all.global <- c(graph.ind$global, year=substr(filenames[[i]], 22, 25), 
                    month=substr(filenames[[i]], 27, 29), pld=substr(filenames[[i]], 31,32)) }
  if (i > 1){
    all.local <- rbind( all.local, cbind(zone = row.names(graph.ind$local), graph.ind$local, year=substr(filenames[[i]], 22, 25), 
                                         month=substr(filenames[[i]], 27, 29), pld=substr(filenames[[i]], 31,32)) )
    all.global <- rbind( all.global, c(graph.ind$global, year=substr(filenames[[i]], 22, 25), 
                                       month=substr(filenames[[i]], 27, 29), pld=substr(filenames[[i]], 31,32)))
  }
}


# PLD and threshold ##########################################################
# calculate graphs metrics local=for each zone, global=for each grap/simulation
# threshold = ( 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)
filenames <- list.files(path="Output60zones/", pattern="subInfo+.*_4w*")
filenames <- paste( "Output60zones/", filenames, sep="")

# the entire range tested
all.threshold = seq( from = 7, to = 0.001, by = -0.01)
# zoom on 0:1 range
all.threshold = c( 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)

# calculate graph's metrics over all simulations
# select the zone = c("all","manche","atlantic"); default is "all"
for (t in 1:length(all.threshold)) {
  for (i in 1:length(filenames)){
    graph.ind <- computeGraphIndices(subInfo = read.table(filenames[[i]], header=T), 
                                   new.threshold = all.threshold[t], zone = "all", g.strenght = "weak")
    filename = substr(filenames[[i]], 22, 32)
    if (i ==1){ 
      all.global <- data.frame(graph.ind$global, year=substr(filenames[[i]], 22, 25), 
                    month=substr(filenames[[i]], 27, 29), pld=substr(filenames[[i]], 31,32), threshold = all.threshold[t]) }
    if (i > 1){
       all.global <- rbind( all.global, data.frame(graph.ind$global, year=substr(filenames[[i]], 22, 25), 
                                       month=substr(filenames[[i]], 27, 29), pld=substr(filenames[[i]], 31,32), threshold = all.threshold[t]))
    }
  }
  ifelse (t == 1, data.threshold <- all.global, data.threshold <- rbind(data.threshold, all.global))
}

# repeat the above for each pld
# bind data from the three plds one by one
th3 <- pld.threshold # 3w
th4 <- pld.threshold # 4w
th6 <- pld.threshold # 6w




# LINKS FREQUENCY ######################################
filenames <- list.files(path="Output60zones/", pattern="subInfo+.*_3w*")
filenames <- paste( "Output60zones/", filenames, sep="")

list.mat <- lapply (1:length(filenames), function (i) {
  computeConnectivityMatrix (subInfo = read.table(filenames[[i]]), 
                             zones = reefs, 
                             percentage = T, 
                             threshold = 0.01) })

list.adj <- lapply(1:length(list.mat), function(i){ ifelse(list.mat[[i]] == 0, 0, 1) })
sd.mat <- apply(simplify2array(list.adj), 1:2, sum ) # check this again
row.names(sd.mat)<-c(1:31)
colnames(sd.mat)<-c(1:31)

sd.mat <- sd.mat / 27 # 27 simulations, hence 27 times max one connection can apear

sd.graph <-graph_from_adjacency_matrix(t(sd.mat), mode = c("directed"), weighted = T, diag = F)

lo<-cbind(habqual[,"lon"],habqual[,"lat"]) 
vertex_attr(sd.graph) <- list(index=row.names(sd.mat) , value=habqual[reefs,"scale"])

# use a high threshold (= 0.5)
sd.high <- graph_from_adjacency_matrix(t(ifelse(sd.mat == 1, 1, 0)), mode = c("directed"), weighted = T, diag = F)
plot(sd.high, layout=lo, edge.color = "black", edge.arrow.size=.2, edge.arrow.mode=0)
# use a high threshold (= 0.01)
sd.low <- graph_from_adjacency_matrix(t(ifelse(sd.mat < 0.1 & sd.mat != 0, 1, 0)), mode = c("directed"), weighted = T, diag = F)
plot(sd.low, layout=lo, edge.color = "black", edge.arrow.size=.2, edge.arrow.mode=0)




# Graphs ~ threshold ##############################
filenames <- list.files(path="Output60zones/", pattern="subInfo+.*_4w*")
filenames <- paste( "Output60zones/", filenames, sep="")

list.mat <- lapply (1:length(filenames), function (i) {
  computeConnectivityMatrix (subInfo = read.table(filenames[[i]]), zones=reefs, percentage = T, threshold = 0.01) })

th.mat <- apply(simplify2array(list.mat), 1:2, mean)
row.names(th.mat) <- c(1:31)
colnames(th.mat) <- c(1:31)
th.mat[th.mat < 0.01] <- 0





# SCENARIOS OF NODE REMOVAL########################
# requires reefs, dist.mat, habqual
filenames <- list.files(path="Output60zones/", pattern="subInfo+.*_4w*")
# add path to the file names subInfo
filenames <- paste( "Output60zones/", filenames, sep="")

library(plyr)
# requires dist.mat, habqual #
dist.mat <- read.table("distanceMatrix.txt") 
habqual<- read.table("scale_habitat_quality.txt", header=T)
new.threshold <- 0.01
reefs <- which(read.table("Data/scale_habitat_quality.txt", header=T)$scale != 0)

# SCENARIOS #
# for average scores per node
nodes.rm.random = computeRandomScenario(filenames, new.threshold, zones = reefs, habqual, dist.mat)
nodes.rm.weak = computeWeakScenario(filenames, new.threshold, zones = reefs, habqual, dist.mat)
nodes.rm.strong = computeStrongScenario(filenames, new.threshold, zones = reefs, habqual, dist.mat)
nodes.rm.weakHab = computeWeakHabScenario(filenames, new.threshold, zones = reefs, habqual, dist.mat)
nodes.rm.strongHab = computeStrongHabScenario(filenames, new.threshold, zones = reefs, habqual, dist.mat)

df <- rbind(nodes.rm.random[1:2], nodes.rm.weak, nodes.rm.strong, nodes.rm.weakHab, nodes.rm.strongHab)
df$scenario <- c(rep("Random", 31), 
                  rep("lowBC", 31), rep("highBC", 31), 
                  rep("smallReef", 31), rep("largeReef", 31))
df$node <- rep(c(1:31), 5)


# for random scenario, all data points
# to plot boxplots
random.all <- computeRandomScenarioAll(filenames, new.threshold, zones=reefs, habqual, dist.mat)
