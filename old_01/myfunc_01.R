swape <- function(vect,i,j) {
    temp <- vect[i]
    vect[i] <- vect[j]
    vect[j] <- temp
    return(vect)
}

swapc <- function(mat,i,j) {
    temp <- mat[,i]
    mat[,i] <- mat[,j]
    mat[,j] <- temp
    return(mat)
}

xchange <- function(A,B,a,b) {
    temp <- A[a]
    A[a] <- B[b]
    B[b] <- temp
    return(list(A=A,B=B))
}


first_vert <- function(A,y) {
    m <- nrow(A)
    n <- ncol(A)
    perm0 <- qr(A)$pivot
    A <- A[,perm0]
    A1 <- A[,1:m]
    A2 <- A[,(m+1):n]
    F <- solve(A1)
    W2 <- F %*% A2
    y_1 <- F %*% y
    J0 <- perm0[(m+1):n]
    cperm <- perm0
    Ineg <- which(y_1 < 0)
    while (length(Ineg) > 0) {
        ineg <- Ineg[1]
        Jneg <- which(W2[ineg,] < 0)
        J0 <- J0[Jneg]
        jneg <- Jneg[1]
        j0 <- J0[1]
        U <- diag(m)
        newcol <- -W2[,jneg]/W2[ineg,jneg]
        newcol[ineg] <- 1/W2[ineg,jneg]
        U[,ineg] <- newcol
        cperm <- swape(cperm,ineg,j0)
        W2 <- U %*% W2[,Jneg]
        W2[,1] <- newcol
        y_1 <- U %*% y_1
        Ineg <- which(y_1 < 0)
    }
    xperm0 <- sort(cperm[1:m])
    cperm <- c(xperm0, sort(cperm[(m+1):n]))
    x0 <- solve(A[,xperm0],y)
    v0 <- rbind(x0,matrix(0,nrow=n-m,ncol=1))
    v0[cperm] <- v0
    return(list(xperm0=xperm0, vperm0=cperm, x0=x0, v0=v0))
}

get_ineg <- function(W2) {
    ineg <- 0
    minn <- ncol(W2) + 1
    for (i in 1:m) {
        if (all(W2[i,] <= 0)) return(-i)
        Jneg <- which(W2[i,] < 0)
        lJ <- length(Jneg)
        if (lJ > 0 & lJ < minn) {
            minn <- lJ
            ineg <- i
        }
    }
    return(ineg)
}

red_pivot <- function(W2, ineg, perm1, perm2) {
    Jneg <- which(W2[ineg,] < 0)
    jneg <- Jneg[1]
    newcol <- -W2[,jneg]/W2[ineg,jneg]
    newcol[ineg] <- 1/W2[ineg,jneg]
    U <- diag(m)
    U[,ineg] <- newcol
    W2 <- U %*% W2[,Jneg]
    W2[,1] <- newcol
    upcol <- newcol - A[,perm1[ineg],drop=F]
    F <- F + 1/drop(1 + F[ineg,] %*% upcol) * F %*% upcol %*% F[ineg,,drop=F]
    temp <- perm1[ineg]
    perm2 <- perm2[Jneg]
    perm1[ineg] <- perm2[1]
    perm2[1] <- temp
    return(list(W2=W2, Jneg=Jneg, perm1=perm1, perm2=perm2, F=F))
}


