# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# FUNCTIONS TO BUILT GRAPH/ESTIMATE GRAPHS METRICS #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# calculate graph indices local & global
# functions for node removal scenarios: 
# random, weakest/strongest links first and smallest/largest reefs first

# Note on constructing the graphs:
# g.disp = graph edges weighted with connectivity fluxes
# g.dist = graph edges weighted with distance matrix
# g.adj = graph edges are 0/1 (presence/absence of connection) from adjacent matrix

######################## Graph indices ########################

computeGraphIndices <- function (subInfo, new.threshold , zone = "all", g.strenght = "weak"){
  # subInfo = file with traceur concentrations
  # new.threshold = minimum concentration treshold for connectivity
  
  # compute connectivity matrix
  source("FunctionsForConnectivity.R")
  disp.mat <- computeConnectivityMatrix (subInfo = subInfo, zones = 1:60, percentage = T, threshold = new.threshold)
  
  # requires file with habitat scores in current directory
  # select only zones with reefs
  habqual<-read.table("Data/scale_habitat_quality.txt", header=T)
  disp.mat <- disp.mat[which(habqual$scale!=0),which(habqual$scale!=0)]
  
  # compute adjacent matrix, when link present=1 and link absent=0
  adj.mat <- ifelse(disp.mat == 0, 0, 1)
  
  # compute distance matrix between connected points
  df.points <- aggregate(subInfo [,c('lat','lon')],subInfo[,'zone',drop=F], mean,na.rm=T)[1:60,]
  new.rownames<- c(1:60) 
  df.points <- df.points [match(new.rownames,df.points$zone),]
  # function included in "FunctionsForConnectivity.R"
  dist.mat <- GeoDistanceInKmMatrix(df.points)
  dist.mat <- dist.mat[which(habqual$scale!=0),which(habqual$scale!=0)] * adj.mat

  # by default zone = "all" then all zones will be included
  if(zone == "all"){
    dist.mat <- dist.mat
    disp.mat <- disp.mat
    adj.mat <- adj.mat
    new.habqual <- habqual[habqual$scale!=0,"scale"]
    nb.zones <- 31
  }
  # select only zones on the ATlantic coast
  if(zone == "atlantic"){
    dist.mat <- dist.mat[12:31,12:31]
    disp.mat <- disp.mat[12:31,12:31]
    adj.mat <- adj.mat[12:31,12:31]
    new.habqual <- habqual[habqual$scale!=0,"scale"][12:31]
    nb.zones <- 20
  }
  # select only zones in the English Channel
  if(zone == "manche"){
    dist.mat <- dist.mat[1:11,1:11]
    disp.mat <- disp.mat[1:11,1:11]
    adj.mat <- adj.mat[1:11,1:11]
    new.habqual <- habqual[habqual$scale!=0,"scale"][1:11]
    nb.zones <- 11
  }

  require('igraph')
  # build directed graph (direction of links considered)
  # graph weighted, edges = transport fluxes from connectivity matrix (disp.mat)
  g.disp <- graph_from_adjacency_matrix(t(disp.mat), mode = c("directed"), weighted = T, diag = T)

  # build directed and weighted graph, edges = distances in km
  g.dist<-graph_from_adjacency_matrix(t(dist.mat), mode = c("directed"), weighted = T, diag = F)
  
  # build directed and weighted graph, edges = 1/0 (presence/absence of connection)
  g.adj<-graph_from_adjacency_matrix(t(adj.mat), mode = c("directed"), weighted = T, diag = F)
  # indices are based on either dispersal or distance matrix
  indice.local <- data.frame( degree = degree(g.disp),         # number of outgoing connections (links)
                              b_dist = betweenness(g.dist),    # betweeness centrality metric based on distances
                              b_adj = betweenness(g.adj) )     # betweeness centrality metric based on nb of steps
  #compute shortest path between nodes (in distance or number of steps)
  # results counts number of links in the path, so for number of nodes +1
  min.path <- shortest.paths(g.adj, algorithm = "dijkstra", mode = "out")+1
  # "Inf" is returned when there is no path connecting two nodes, replace it for mean calculation
  min.path[min.path=="Inf"]<-NA 
  
  # global indices are independent from links weight, only presence/absence of links counts
  indice.global <- data.frame(comp.no = components(g.disp, mode = g.strenght)$no,
                              comp.size = max(components(g.disp, mode = g.strenght)$csize)/nb.zones )
  
  return(list(local=indice.local, global=indice.global))
}




######################## COMPUTE GRAPH SCENARIOS ########################
# uses the function 'components' from the 'igraph' library
# all graphs edges are built on distance matrix

######################## Random Scenario ########################
# =================================================== #
### random node removal complete ###
computeRandomScenarioAll <- function (filenames, new.threshold, zones, habqual, dist.mat){   
  require('igraph')
  require('plyr')
  
  all.files <- lapply(1:length(filenames), function(s){
    
    s.mat <- computeConnectivityMatrix (subInfo = read.table(filenames[[s]]), zones, percentage = T, threshold = new.threshold)
    
    # add the attributes = distance (km)
    s.mat <- ifelse(s.mat == 0, 0, 1) * dist.mat[which(habqual$scale!=0),which(habqual$scale!=0)]
    nb.zones <- length(zones)
    row.names(s.mat) <- c(1:nb.zones)
    colnames(s.mat) <- c(1:nb.zones)
    
    random.all <- lapply(1:125, function(x){
      
      # start with complete graph
      g <- graph_from_adjacency_matrix(t(s.mat), mode = c("directed"), weighted = T, diag = F)
      nodes.id <- row.names( as.data.frame( components(g, mode = "weak")$membership))
      remove.order <- sample(nodes.id, length(nodes.id), replace = F)
      
      # compute metrics for baseline graph
      nodes.rm <- data.frame(#dens = round(edge_density(g), digits = 2), 
        #reciproc = round(reciprocity(g), digits = 2),
        #trans = round(transitivity(g), digits = 2),
        #diam = as.integer(diameter(g, directed = T)),
        comp.no = as.integer(components(g, mode = "weak")$no),
        comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
      for (i in 1:length(remove.order)){
        # node to remove
        g <- delete.vertices(g, remove.order[1] )
        new.nodes.rm <- data.frame( #dens = round(edge_density(g), digits = 2), 
          #reciproc = round(reciprocity(g), digits = 2),
          #trans = round(transitivity(g), digits = 2),
          #diam = as.integer(diameter(g, directed = T)),
          comp.no = as.integer(components(g, mode = "weak")$no),
          comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
        nodes.rm <- rbind(nodes.rm, new.nodes.rm)
        nodes.id <- row.names( as.data.frame( components(g, mode = "weak")$membership))
        remove.order <- sample(nodes.id, length(nodes.id), replace = F)
      }
      nodes.rm[nodes.rm=="NaN"] <- 1
      nodes.rm.random <- nodes.rm[1:(nb.zones),]
      x <- nodes.rm.random
    })
    
    nodes.rm.random <- data.frame(aaply(laply(random.all, as.matrix), c(2, 3), mean))
    #nodes.rm.random <- do.call(rbind, random.all)
    s <- nodes.rm.random
    
  }) 
  #nodes.rm.random <- data.frame(aaply(laply(all.files, as.matrix), c(2, 3), mean))
  nodes.rm.random <- do.call(rbind, all.files)
  return(nodes.rm.random)
}

# =================================================== #
# returns averages over all simulated events (N = 27) #
computeRandomScenario <- function (filenames, new.threshold, zones, habqual, dist.mat){   
  # filesnames: vector containing all path and name of files with raw dispersal data
  # new.threshols: connectivity threshold, as minimum flux of particle, below which = 0
  # zones: labels of coastal zones with reefs (31 out of 60)
  # habqual: scale of reef status (0, 0.25, 0.5, 0.75 and 1 the largest reefs)
  # dist.mat: full matrix with linear distances between all zones
  require('igraph')
  require('plyr')
  
all.files <- lapply(1:length(filenames), function(s){
  
  s.mat <- computeConnectivityMatrix (subInfo = read.table(filenames[[s]]), zones, percentage = T, threshold = new.threshold)
  
  # add the attributes = distance (km)
  s.mat <- ifelse(s.mat == 0, 0, 1) * dist.mat[which(habqual$scale!=0),which(habqual$scale!=0)]
  nb.zones <- length(zones)
  row.names(s.mat) <- c(1:nb.zones)
  colnames(s.mat) <- c(1:nb.zones)
  
  random.all <- lapply(1:125, function(x){
    
    # start with complete graph
    g <- graph_from_adjacency_matrix(t(s.mat), mode = c("directed"), weighted = T, diag = F)
    nodes.id <- row.names( as.data.frame( components(g, mode = "weak")$membership))
    remove.order <- sample(nodes.id, length(nodes.id), replace = F)
    
    # compute metrics for baseline graph
    nodes.rm <- data.frame(comp.no = as.integer(components(g, mode = "weak")$no),
                           comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
    for (i in 1:length(remove.order)){
      # node to remove
      g <- delete.vertices(g, remove.order[1] )
      new.nodes.rm <- data.frame(comp.no = as.integer(components(g, mode = "weak")$no),
                                comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
      
      nodes.rm <- rbind(nodes.rm, new.nodes.rm)
      nodes.id <- row.names( as.data.frame( components(g, mode = "weak")$membership))
      remove.order <- sample(nodes.id, length(nodes.id), replace = F)
    }
    nodes.rm[nodes.rm=="NaN"] <- 1
    nodes.rm.random <- nodes.rm[1:(nb.zones),]
    x <- nodes.rm.random
  })
  
  nodes.rm.random <- data.frame(aaply(laply(random.all, as.matrix), c(2, 3), mean))
  s <- nodes.rm.random
}) 
nodes.rm.random <- data.frame(aaply(laply(all.files, as.matrix), c(2, 3), mean))
return(nodes.rm.random)
}


######################## Weakest Links Scenario ########################
# remove the weakest links first
# betweeness metric estimated in graphs built on distance matrix

computeWeakScenario <- function (filenames, new.threshold, zones, habqual, dist.mat){   
  # filesnames: vector containing all path and name of files with raw dispersal data
  # new.threshols: connectivity threshold, as minimum flux of particle, below which = 0
  # zones: labels of coastal zones with reefs (31 out of 60)
  # habqual: scale of reef status (0, 0.25, 0.5, 0.75 and 1 the largest reefs)
  # dist.mat: full matrix with linear distances between all zones
  require('igraph')
  require('plyr')
  
all.files <- lapply(1:length(filenames), function(s){
  
  s.mat <- computeConnectivityMatrix (subInfo = read.table(filenames[[s]]), zones, percentage = T, threshold = new.threshold)
  # add the attributes = distance (km)
  s.mat <- dist.mat[which(habqual$scale!=0),which(habqual$scale!=0)] * ifelse(s.mat == 0, 0, 1)
  nb.zones <- length(zones)
  row.names(s.mat) <- c(1:nb.zones)
  colnames(s.mat) <- c(1:nb.zones)
  
  # start with complete graph
  g <- graph_from_adjacency_matrix(t(s.mat), mode = c("directed"), weighted = T, diag = F)
  # create the order with the weakest node
  nodes.id <- row.names( as.data.frame( components(g, mode = "weak")$membership))
  remove.order <- nodes.id[order(betweenness(g), decreasing = F)]
  
  # compute metrics for baseline graph
  nodes.rm <- data.frame(node = c("none"),
                         comp.no = as.integer(components(g, mode = "weak")$no),
                         comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
  for (i in 1:length(remove.order)){
    # node to remove
    g <- delete.vertices(g, remove.order[1] )
    new.nodes.rm <- data.frame(node = remove.order[1], 
                               comp.no = as.integer(components(g, mode = "weak")$no),
                               comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
    nodes.rm <- rbind(nodes.rm, new.nodes.rm)
    nodes.id <- row.names( as.data.frame( components(g, mode = "weak")$membership))
    remove.order <- nodes.id[order(betweenness(g), decreasing = F)]
  }
  nodes.rm[nodes.rm=="NaN"] <- 0
  nodes.rm[nodes.rm=="Na"] <- 0
  nodes.rm.weak <- nodes.rm[1:(nb.zones),2:4]
  s <- nodes.rm.weak
}) 
# node = remove.order is valid for last simulation only, as the order change between diff simu
nodes.rm.weak <- data.frame(aaply(laply(all.files, as.matrix), c(2, 3), mean))
return(nodes.rm.weak)
}


######################## Strongest Links Scenario ########################
# remove the strongest links
# betweeness metric estimated in graphs built on distance matrix

computeStrongScenario <- function (filenames, new.threshold, zones, habqual, dist.mat){   
  # filesnames: vector containing all path and name of files with raw dispersal data
  # new.threshols: connectivity threshold, as minimum flux of particle, below which = 0
  # zones: labels of coastal zones with reefs (31 out of 60)
  # habqual: scale of reef status (0, 0.25, 0.5, 0.75 and 1 the largest reefs)
  # dist.mat: full matrix with linear distances between all zones
  require('igraph')
  require('plyr')
  
all.files <- lapply(1:length(filenames), function(s){
  
  s.mat <- computeConnectivityMatrix (subInfo = read.table(filenames[[s]]), zones, percentage = T, threshold = new.threshold)
  # add the attributes = distance (km)
  s.mat <- dist.mat[which(habqual$scale!=0),which(habqual$scale!=0)] * ifelse(s.mat == 0, 0, 1)
  nb.zones <- length(zones)
  row.names(s.mat) <- c(1:nb.zones)
  colnames(s.mat) <- c(1:nb.zones)
  
  # start with complete graph
  g <- graph_from_adjacency_matrix(t(s.mat), mode = c("directed"), weighted = T, diag = F)
  # create the order with the weakest node
  nodes.id <- row.names( as.data.frame( components(g, mode = "weak")$membership))
  remove.order <- nodes.id[order(betweenness(g), decreasing = T)]
  
  # compute metrics for baseline graph
  nodes.rm <- data.frame(node = c("none"),
                         comp.no = as.integer(components(g, mode = "weak")$no),
                         comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
  for (i in 1:length(remove.order)){
    # node to remove
    g <- delete.vertices(g, remove.order[1] )
    new.nodes.rm <- data.frame(node = remove.order[1],
                               comp.no = as.integer(components(g, mode = "weak")$no),
                               comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
    nodes.rm <- rbind(nodes.rm, new.nodes.rm)
    nodes.id <- row.names( as.data.frame( components(g, mode = "weak")$membership))
    remove.order <- nodes.id[order(betweenness(g), decreasing = T)]
  }
  nodes.rm[nodes.rm=="NaN"] <- 0
  nodes.rm[nodes.rm=="Na"] <- 0
  nodes.rm.strong <- nodes.rm[1:(nb.zones),2:4]
  s <- nodes.rm.strong
}) 
# node = remove.order is valid for last simulation only, as the order change between diff simu
nodes.rm.strong <- data.frame(aaply(laply(all.files, as.matrix), c(2, 3), mean))
return(nodes.rm.strong)
}


######################## Smallest Habitat Scenario ########################
# random by smallest habitat

computeWeakHabScenario <- function (filenames, new.threshold, zones, habqual, dist.mat){   
  # filesnames: vector containing all path and name of files with raw dispersal data
  # new.threshols: connectivity threshold, as minimum flux of particle, below which = 0
  # zones: labels of coastal zones with reefs (31 out of 60)
  # habqual: scale of reef status (0, 0.25, 0.5, 0.75 and 1 the largest reefs)
  # dist.mat: full matrix with linear distances between all zones
  require('igraph')
  require('plyr')
  
all.files <- lapply(1:length(filenames), function(s){
  
  s.mat <- computeConnectivityMatrix (subInfo = read.table(filenames[[s]]), zones, percentage = T, threshold = new.threshold)
  # add the attributes = distance (km)
  s.mat <- dist.mat[which(habqual$scale!=0),which(habqual$scale!=0)] * ifelse(s.mat == 0, 0, 1)
  nb.zones <- length(zones)
  row.names(s.mat) <- c(1:nb.zones)
  colnames(s.mat) <- c(1:nb.zones)
  
  new.habqual <- habqual[habqual$scale != 0,]
  new.habqual$zone <- 1:nb.zones
  random.all <- lapply(1:100, function(x){
    
    # start with complete graph
    g <- graph_from_adjacency_matrix(t(s.mat), mode = c("directed"), weighted = T, diag = F)
    
    # for habqual score 0.25
    nodes.id <- as.character(new.habqual[new.habqual$scale==0.25, "zone"])
    remove.order <- sample(nodes.id, length(nodes.id), replace = F)
    
    # compute metrics for baseline graph
    nodes.rm <- data.frame(node = "none", 
                           comp.no = as.integer(components(g, mode = "weak")$no),
                           comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
    for (i in 1:length(remove.order)){
      # node to remove
      g <- delete.vertices(g, remove.order[i] )
      new.nodes.rm <- data.frame( node = remove.order[i],
                                  comp.no = as.integer(components(g, mode = "weak")$no),
                                  comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
      nodes.rm <- rbind(nodes.rm, new.nodes.rm)
    }
    # for habqual score 0.50
    nodes.id <- as.character(new.habqual[new.habqual$scale==0.50, "zone"])
    remove.order <- sample(nodes.id, length(nodes.id), replace = F)
    for (i in 1:length(remove.order)){
      # node to remove
      g <- delete.vertices(g, remove.order[i] )
      new.nodes.rm <- data.frame( node = remove.order[i],
                                  comp.no = as.integer(components(g, mode = "weak")$no),
                                  comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
      nodes.rm <- rbind(nodes.rm, new.nodes.rm)
    }
    # for habqual score 0.75
    nodes.id <- as.character(new.habqual[new.habqual$scale==0.75, "zone"])
    remove.order <- sample(nodes.id, length(nodes.id), replace = F)
    for (i in 1:length(remove.order)){
      # node to remove
      g <- delete.vertices(g, remove.order[i] )
      new.nodes.rm <- data.frame( node = remove.order[i],
                                  comp.no = as.integer(components(g, mode = "weak")$no),
                                  comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
      nodes.rm <- rbind(nodes.rm, new.nodes.rm)
    }
    # for habqual score 1
    nodes.id <- as.character(new.habqual[new.habqual$scale==1, "zone"])
    remove.order <- sample(nodes.id, length(nodes.id), replace = F)
    for (i in 1:length(remove.order)){
      # node to remove
      g <- delete.vertices(g, remove.order[i] )
      new.nodes.rm <- data.frame( node = remove.order[i],
                                  comp.no = as.integer(components(g, mode = "weak")$no),
                                  comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
      nodes.rm <- rbind(nodes.rm, new.nodes.rm)
    }
    nodes.rm[nodes.rm=="NaN"] <- 0
    x <- nodes.rm[1:(nb.zones),2:4]
  })
  
  nodes.rm.weakHab <- data.frame(aaply(laply(random.all, as.matrix), c(2, 3), mean))
  s <- nodes.rm.weakHab
}) 
nodes.rm.weakHab <- data.frame(aaply(laply(all.files, as.matrix), c(2, 3), mean))
return(nodes.rm.weakHab)
}


######################## Largest Reef Scenario ########################
# random by largest habitat

computeStrongHabScenario <- function (filenames, new.threshold, zones, habqual, dist.mat){   
  # filesnames: vector containing all path and name of files with raw dispersal data
  # new.threshols: connectivity threshold, as minimum flux of particle, below which = 0
  # zones: labels of coastal zones with reefs (31 out of 60)
  # habqual: scale of reef status (0, 0.25, 0.5, 0.75 and 1 the largest reefs)
  # dist.mat: full matrix with linear distances between all zones
  require('igraph')
  require('plyr')
  
all.files <- lapply(1:length(filenames), function(s){
  
  s.mat <- computeConnectivityMatrix (subInfo = read.table(filenames[[s]]), zones, percentage = T, threshold = new.threshold)
  # add the attributes = distance (km)
  s.mat <- dist.mat[which(habqual$scale!=0),which(habqual$scale!=0)] * ifelse(s.mat == 0, 0, 1)
  nb.zones <- length(zones)
  row.names(s.mat) <- c(1:nb.zones)
  colnames(s.mat) <- c(1:nb.zones)
  
  new.habqual <- habqual[habqual$scale != 0,]
  new.habqual$zone <- 1:nb.zones
  
  random.all <- lapply(1:100, function(x){
    
    # start with complete graph
    g <- graph_from_adjacency_matrix(t(s.mat), mode = c("directed"), weighted = T, diag = F)
    
    # for habqual score 1
    nodes.id <- as.character(new.habqual[new.habqual$scale==1, "zone"])
    remove.order <- sample(nodes.id, length(nodes.id), replace = F)
    
    # compute metrics for baseline graph
    nodes.rm <- data.frame(node = "none", 
                           comp.no = as.integer(components(g, mode = "weak")$no),
                           comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
    for (i in 1:length(remove.order)){
      # node to remove
      g <- delete.vertices(g, remove.order[i] )
      new.nodes.rm <- data.frame( node = remove.order[i],
                                  comp.no = as.integer(components(g, mode = "weak")$no),
                                  comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
      nodes.rm <- rbind(nodes.rm, new.nodes.rm)
    }
    # for habqual score 0.75
    nodes.id <- as.character(new.habqual[new.habqual$scale==0.75, "zone"])
    remove.order <- sample(nodes.id, length(nodes.id), replace = F)
    for (i in 1:length(remove.order)){
      # node to remove
      g <- delete.vertices(g, remove.order[i] )
      new.nodes.rm <- data.frame( node = remove.order[i],
                                  comp.no = as.integer(components(g, mode = "weak")$no),
                                  comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
      nodes.rm <- rbind(nodes.rm, new.nodes.rm)
    }
    # for habqual score 0.50
    nodes.id <- as.character(new.habqual[new.habqual$scale==0.50, "zone"])
    remove.order <- sample(nodes.id, length(nodes.id), replace = F)
    for (i in 1:length(remove.order)){
      # node to remove
      g <- delete.vertices(g, remove.order[i] )
      new.nodes.rm <- data.frame( node = remove.order[i],
                                  comp.no = as.integer(components(g, mode = "weak")$no),
                                  comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
      nodes.rm <- rbind(nodes.rm, new.nodes.rm)
    }
    # for habqual score 0.25
    nodes.id <- as.character(new.habqual[new.habqual$scale==0.25, "zone"])
    remove.order <- sample(nodes.id, length(nodes.id), replace = F)
    for (i in 1:length(remove.order)){
      # node to remove
      g <- delete.vertices(g, remove.order[i] )
      new.nodes.rm <- data.frame( node = remove.order[i],
                                  comp.no = as.integer(components(g, mode = "weak")$no),
                                  comp.size = round(max(components(g, mode = "weak")$csize)/nb.zones, digits = 2))
      nodes.rm <- rbind(nodes.rm, new.nodes.rm)
    }
    nodes.rm[nodes.rm=="NaN"] <- 0
    x <- nodes.rm[1:(nb.zones),2:4]
  })
  
  nodes.rm.strongHab <- data.frame(aaply(laply(random.all, as.matrix), c(2, 3), mean))
  s <- nodes.rm.strongHab
}) 
nodes.rm.strongHab <- data.frame(aaply(laply(all.files, as.matrix), c(2, 3), mean))
return(nodes.rm.strongHab)
}


