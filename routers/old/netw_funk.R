#generates the network "router2" from Cao et al. as a network object (requires library(network))
#Vertices have an attribute "rout" that is 1 if they are a router and 0 if they are an origin/destination.
#Edges has an attribute "dest" which is a vector of allowed destinations.
#One can plots the network with "plot.network(n2r)"
make_tworout <- function() {
    vnames <- c("a", "b", "c", "d", "e", "f", "g", "h", "r1", "r2")
    enames <- c("a_r1", "b_r1", "c_r1", "d_r1", "e_r2", "f_r2", "g_r2", "h_r2", "r1_a", "r1_b", "r1_c", "r1_d", "r2_e", "r2_e", "r2_f", "r2_g", "r1_r2", "r2_r1")
    r1r2 <- list("e", "f", "g", "h")
    r2r1 <- list("a", "b", "c", "d")
    elist <- matrix(0, nrow=18, ncol=2)
    for (j in 9:10){
        for (i in 1:4) {
            elist[((j-9)*4+i),1] <- (j-9)*4+i
            elist[((j-9)*4+i),2] <- j
        }
    }
    elist[9:16,1] <- elist[1:8,2]
    elist[9:16,2] <- elist[1:8,1]
    elist[17,1] <- 9
    elist[17,2] <- 10
    elist[18,1] <- 10
    elist[18,2] <- 9

    n2r <- network.initialize(10)
    n2r <- network.edgelist(elist, n2r, names.eval=enames)
    network.vertex.names(n2r) <- vnames
    set.vertex.attribute(n2r, "rout", c(rep(0,8), rep(1,2)))
    alldest <- which(get.vertex.attribute(n2r, "rout", unlist=F) == 0)
    vatt <- get.vertex.attribute(n2r, "rout", unlist=F)
    for (e in 1:nrow(elist)) {
        if (vatt[elist[e,1]] == 0) {
            set.edge.attribute(n2r, "dest", list(alldest), e)
        } else if (vatt[elist[e,2]] == 0) {
            set.edge.attribute(n2r, "dest", list(elist[e,2]), e)
        }
    }
    set.edge.attribute(n2r, "dest", list(c(5,6,7,8)), 17)
    set.edge.attribute(n2r, "dest", list(c(1,2,3,4)), 18)
    return(n2r)
}

#Computes the neighbors of a vertex where a route to a destination dest can pass.
neighb <- function(nw, vert, dest) {
    elist <- as.matrix.network(nw, "edgelist")
    idx <- which(elist[,1] == vert)
    nop <- NULL
    for (i in idx) {
        if (any(grepl(dest, get.edge.attribute(nw$mel[i], "dest"))) == F) {
            nop <- c(nop, i)
        }
    }
    idx <- setdiff(idx, nop)
    return(elist[idx,2])
}

#gets a list of paths from orig to dest (possibly more than one)
get_odpaths <- function(nw, orig, dest, lpaths) {
    onei <- neighb(nw, orig, dest)
    if (length(onei) == 0) {
        return(lpaths)
    }
    for (k in onei) {
        if (k == dest) {
            lpaths <- c(lpaths, list(c(orig, dest)))
        } else {
            kpaths <- get_odpaths(nw, k, dest, lpaths)
            for (p in seq(from = 1, by=1, length.out=length(kpaths))) {
                kpaths[[p]] <- c(orig, kpaths[[p]])
            }
            lpaths <- c(lpaths, kpaths)
        }
    }
    return(lpaths)
}

#gets all the paths from all the origins to all destinations
#prints a warning message if more than one path for a given origin/destination pair has been found
get_allpaths <- function(nw) {
    odlist <- which(get.vertex.attribute(nw, "rout", unlist=F) == 0)
    allpaths <- NULL
    for (i in odlist) {
        for (j in odlist) {
            odpaths <- get_odpaths(nw, i, j, NULL)
            if (length(odpaths) > 1) {
                cat("More than one path found from origin ", i, " to destination ", j,"\n")
            }
            allpaths <- c(allpaths, odpaths)
        }
    }
    return(allpaths)
}

#makes the matrix A for the network
make_A <- function(nw) {
    allpaths <- get_allpaths(nw)
    elist <- as.matrix.network(nw, "edgelist")
    n <- length(allpaths)
    m <- nrow(elist)
    A <- matrix(0, nrow=m, ncol=n)
    for (j in 1:n) {
        vidx <- allpaths[[j]]
        for (i in 2:length(vidx)) {
            e2 <- which(elist[,2] == vidx[i])
            e1 <- which(elist[e2,1] == vidx[i-1])
            A[e2[e1],j] <- 1
        }
    }
    return(A)
}









