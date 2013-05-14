resample <- function(x, ...) {
    return(x[sample.int(length(x), ...)])
}

swapel <- function(v, i, j) {
    temp <- v[i]
    v[i] <- v[j]
    v[j] <- temp
    return(v)
}

swapro <- function(M, i, j) {
    temp <- M[i,]
    M[i,] <- M[j,]
    M[j,] <- temp
    return(M)
}

getdata <- function(arith){
    #library(gmp)
    A <- as.matrix(read.table("A.txt"))
    m <<- nrow(A)
    n <<- ncol(A)
    qrA <- qr(A)
    if (qrA$rank < m) return("Error: A is not full rank")
    perm <<- qrA$pivot
    if (arith == "exact") {
        A <<- as.bigz(A)
        y <<- as.bigz(as.matrix(read.table("y.txt")))
        y1 <<- as.bigz(as.matrix(read.table("y1.txt")))
        y2 <<- as.bigz(as.matrix(read.table("y2.txt")))
        y3 <<- as.bigz(as.matrix(read.table("y3.txt")))
    } else {
        A <<- A
        y <<- as.matrix(read.table("y.txt"))
        y1 <<- as.matrix(read.table("y1.txt"))
        y2 <<- as.matrix(read.table("y2.txt"))
        y3 <<- as.matrix(read.table("y3.txt"))
    }
}

xpivot <- function(x1,Q12,i,j) {
  if (Q12[i,j] == 0) return("Error: Q12[i,j] == 0")
    t <- - x1[i] / Q12[i,j]
    x1 <- x1 + t * Q12[,j]
    x1[i] <- t
    return(x1)
}

geti <- function(x1, Q12, j) {
    I <- which(x1 * Q12[,j] < 0)
    if (length(I) == 0) return(-1)
    minel <- -max(x1[I]/Q12[I,j])
    i <- I[resample(which(x1[I]/Q12[I,j] == -minel),1)]
    return(i)
}

Qpivot <- function(Q12, perm, i, j) {
    m <- nrow(Q12)
    n <- m + ncol(Q12)
    Q12 <- rbind(Q12,matrix(0,ncol=n-m,nrow=1))
    Q12[m+1,j] <- 1
    for (k in setdiff(1:(n-m), j)) {
        Q12[,k] <- Q12[,k] - matrix(Q12[i,k],nrow=m+1,ncol=1) * Q12[,j] / matrix(Q12[i,j],nrow=m+1,ncol=1)
    }
    Q12[,j] <- Q12[,j]/matrix(Q12[i,j],nrow=m+1,ncol=1)
    Q12 <- swapro(Q12,i, m+1)
    perm <- swapel(perm, i, m+j)
    return(list(Q12=Q12[1:m,],perm=perm))
}

trygetij <- function(x1, Q12, I0, i0, sig) {
    m <- nrow(Q12)
    if (i0 > m) {
        if (sig == 1) {
            J = i0-m
        } else return(list(j=-1,J=NULL))
    } else {
        J <- which(sig * Q12[i0,] > 0)
    }
    if (length(J) == 0) return(list(i=-1,j=-1,J=J))
    if (length(I0) > 0) {
        K <- c()
        for (k in J) {
            if (any(Q12[I0,k] < 0)) K <- c(K,k)
        }
        #K <- J[unique(which(Q12[I0,J,drop=F] < 0, arr.ind=T)[,2])]
        JmK <- setdiff(J,K)
        if (length(JmK) == 0) return(list(i=-1,j=-1,J=J))
    } else {
        JmK <- J
    }
    idx <- sample.int(length(JmK), length(JmK))
    for (j in JmK[idx]) {
        i <- geti(x1, Q12, j)
        if (i >= 1) return(list(i=i, j=j, J=J))
    }
    return(list(i=-1, j=-1, J=J))
}

getij <- function(x1, Q12, perm, I0, i0, sig) {
    Qp <- list(Q12=Q12, perm=perm)
    ijJ <- trygetij(x1, Q12, I0, i0, sig)
    id <- perm[i0]
    while (ijJ$j < 0) {
        if (length(ijJ$J) == 0) return(list(i=-1, j=-1, Q12=Qp$Q12, perm=Qp$perm))
        j <- resample(ijJ$J,1)
        i <- I0[which(Qp$Q12[I0,j] < 0)[1]]
        Qp <- Qpivot(Qp$Q12,Qp$perm,i,j)
        i0 <- which(Qp$perm == id)
        ijJ <- trygetij(x1, Qp$Q12, I0, i0, sig)
    }
    return(list(i=ijJ$i, j=ijJ$j, Q12=Qp$Q12, perm=Qp$perm))
}

nextvert <- function(x1, Q12, perm, i0, sig) {
    I0 <- which(x1 == 0)
    L <- getij(x1, Q12, perm, I0, i0, sig)
    if (L$i == -1) return(list(x1=x1, Q12=Q12, perm=perm))
    x1 <- xpivot(x1, L$Q12, L$i, L$j)
    L <- Qpivot(L$Q12, L$perm, L$i, L$j)
    return(list(x1=x1, Q12=L$Q12, perm=L$perm))
}

getvert <- function(x1, Q12, perm) {
    L <- list(x1=x1, Q12=Q12, perm=perm)
    while (any(L$x1 < 0)) {
        Ineg <- which(L$x1 < 0)
        minel <- min(-L$x1[Ineg])
        i0 <- resample(which(L$x1 == -minel),1)
        L <- nextvert(L$x1, L$Q12, L$perm, i0, 1)
    }
    return(list(x1=L$x1, Q12=L$Q12, perm=L$perm))
}

getpv <- function(A, y, perm) {
    m <- nrow(A)
    n <- ncol(A)
    A1 <- A[,1:m]
    if (qr(A1)$rank < m) {
        print("Error:  A1 not full rank")
    }
    A2 <- A[,(m+1):n]
    Q11 <- solve(A1)
    x1 <- Q11 %*% y
    Q12 <- matrix(-Q11, dim(Q11)) %*% A2
    return(list(x1=x1, Q12=Q12, perm=perm))
}

getbound <- function(x1, Q12, perm, i0,sig) {
    L <- list(x1=x1, Q12=Q12, perm=perm)
    id <- perm[i0]
    #bound = -(max(abs(L$x1)) + 1)
    x <- rbind(L$x1, matrix(0,nrow=ncol(Q12),ncol=1))
    #while (sig*x[i0] > bound) {
    repeat {
        bound <- x[i0]
        L <- nextvert(L$x1, L$Q12, L$perm, i0, sig)
        i0 <- which(L$perm == id)
        x <- rbind(L$x1, matrix(0,nrow=ncol(Q12),ncol=1))
        if (bound == x[i0]) break
    }
    return(list(bound=bound, x1=L$x1, Q12=L$Q12, perm=L$perm))
}

get1bd <- function(A, y, perm, id) {
    n <- ncol(A)
    m <- nrow(A)
    L <- getpv(A,y,perm)
    L <- getvert(L$x1, L$Q12, L$perm)
    if (class(A) == "bigz") {
        xbd <- matrix.bigz(0, nrow=2, ncol=n)
        Vmin <- matrix.bigz(0, nrow=n, ncol=n)
        Vmax <- matrix.bigz(0, nrow=n, ncol=n)
    } else {
        xbd <- matrix(0, nrow=2, ncol=n)
        Vmin <- matrix(0, nrow=n, ncol=n)
        Vmax <- matrix(0, nrow=n, ncol=n)
    }
    i0 <- which(L$perm == id)
    B <- getbound(L$x1, L$Q12, L$perm, i0, -1)
    vmin <- rbind(B$x1, matrix(0,nrow=n-m,ncol=1))
    vmin[B$perm] <- matrix(vmin,ncol=1)
    B <- getbound(L$x1, L$Q12, L$perm, i0, +1)
    vmax <- rbind(B$x1, matrix(0,nrow=n-m,ncol=1))
    vmax[B$perm] <- matrix(vmax, ncol=1)
    xmin <- vmin[id]
    xmax <- vmax[id]
    return(list(xmin=xmin, xmax=xmax, vmin=vmin, vmax=vmax))
}

reduceA <- function(A,y,perm, id) {
    qrA <- qr(t(A))
    A <- A[qrA$pivot[1:qrA$rank],]
    y <- y[qrA$pivot[1:qrA$rank]]
    qrA <- qr(A)
    A <- A[,qrA$pivot]
    perm <- perm[qrA$pivot]
    nzcol <- 1:ncol(A)
    for (j in 1:ncol(A)) {
        if (all(A[,j] == 0)) nzcol <- setdiff(nzcol, j)
    }
    A <- A[, nzcol]
    perm <- perm[nzcol]
    i0 <- which(perm == id)
    perm <- sort.int(sort.int(perm, index.return=T)$ix, index.return=T)$ix
    if (length(i0) == 0) {
        print("Warning id = -1")
        id <- -1
    } else id <- perm[i0]
    return(list(A=A,y=y,perm=perm,id=id))
}


get1bdlp <- function(A, y, perm0, id) {
    m <- nrow(A)
    n <- ncol(A)
    of <- rep(0,n)
    of[id] <- 1
    B <- lp(direction="min", objective.in=of, const.mat=A, const.dir="==", const.rhs=y)
    xmin <- B$objval
    vmin <- B$solution
    B <- lp(direction="max", objective.in=of, const.mat=A, const.dir="==", const.rhs=y)
    xmax <- B$objval
    vmax <- B$solution
    return(list(xmin=xmin, xmax=xmax, vmin=vmin, vmax=vmax))
}


getbox <- function(A, y, perm) {
    n <- ncol(A)
    m <- nrow(A)
    L <- getpv(A,y,perm)
    L <- getvert(L$x1, L$Q12, L$perm)
    if (class(A) == "bigz") {
        xbd <- matrix.bigz(0, nrow=2, ncol=n)
        Vmin <- matrix.bigz(0, nrow=n, ncol=n)
        Vmax <- matrix.bigz(0, nrow=n, ncol=n)
    } else {
        xbd <- matrix(0, nrow=2, ncol=n)
        Vmin <- matrix(0, nrow=n, ncol=n)
        Vmax <- matrix(0, nrow=n, ncol=n)
    }
    for (i in 1:n) {
        id <- perm[i]
        i0 <- which(L$perm == id)
        B <- getbound(L$x1, L$Q12, L$perm, i0, -1)
        x <- rbind(B$x1,matrix(0,nrow=n-m,ncol=1))
        x[B$perm] <- x
        Vmin[,i] <- x
        B <- getbound(L$x1, L$Q12, L$perm, i0, +1)
        x <- rbind(B$x1,matrix(0,nrow=n-m,ncol=1))
        x[B$perm] <- x
        Vmax[,i] <- x
    }
    xbd[1,] <- diag(Vmin)
    xbd[2,] <- diag(Vmax)
    return(list(xbd=xbd, Vmin=Vmin, Vmax=Vmax))
}

#checkAxy <- function(x1, Q12, perm, A, y) {
#    m <- nrow(A)
#    n <- ncol(A)
#    cekv <- all(A[,perm[1:m]] %*% x1 == y)
#    if (cekv) {
#        cat("Ax = y\n")
#    } else {
#        cat("(Ax, y) = ")
#        print(cbind(A[,perm[1:m]] %*% x1, y))
#    }
#}


