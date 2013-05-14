swape <- function(vect,i,j) {
    temp <- vect[i]
    vect[i] <- vect[j]
    vect[j] <- temp
    return(vect)
}

get_piv <- function(W2, y1) {
    #record the current minimum (in abs) of y negative coord.
    #initialize ipiv and jpiv to unrealistic values (for debug)
    miny0 <- -y1[which(y1 == max(y1[which(y1 < 0)]))]
    ipiv <- -1
    jpiv <- -1
    for (j in 1:ncol(W2)) {
        w <- W2[,j]
        Ineg <- which(y1 < 0 & w < 0);
        Izer <- which(y1 == 0 & w != 0);
        Ipos <- which(y1 > 0 & w != 0);
        #needs to pivot where y[i]<0 and w[i]<0
        #Moreover if y[i]==0 then need w[i]<0 and if y[i]>0 then need w[i]>0
        #This is to avoid creating more negative y coordinates out of non-neg ones
        #recall the formula:  y1 <- y - (y[i]/w[i])w and y1[i] <- y[i]/w[i]
        if (length(Ineg) > 0 & all(w[Izer] < 0) & all(w[Ipos] > 0)) {
            #test for all potential i-pivot
            for (i in Ineg) {
                #if we don't create any new negative coordinate by eliminating
                #the negative y[i] (recall that y1[i] <- y[i]/w[i] > 0)
                #the we got a pivot
                if (all(y[Ipos] * w[i] <= y[i] * w[Ipos])) {
                    return(list(i=i,j=j))
                #Otherwise we record the pivot that decreases the absolute value
                #of the smallest negative coordinate
                } else {
                    y2 <- y1 - y1[i]/w[i] * w
                    miny0 <- -y2[which(y2 == max(y2[which(y2 < 0)]))]
                    if (miny2 < miny0) {
                        miny0 <- miny2
                        ipiv <- i
                        jpiv <- j
                    }
                }
            }
        }
    }
    return(list(i=ipiv,j=jpiv))
}

rget_piv <- function(W2, y1) {
    Ineg <- which(y1 < 0)
    ineg <- Ineg[runif(1,min=1,max=length(Ineg))]
    Jneg <- which(W2[ineg,] <0)
    jneg <- Jneg[runif(1,min=1,max=length(Jneg))]
    return(list(i=ineg,j=jneg))
}

pivot <- function(W2, y1, i, j) {
    m <- nrow(W2)
    w <- W2[,j]
    U <- diag(m)
    U[,i] <- -w/w[i]
    U[i,i] <- 1/w[i]
    W2 <- U %*% W2
    W2[,j] <- U[,i]
    y1 <- U %*% y1
    return(list(W2=W2, y1=y1))
}

first_vert <- function(A,y) {
    m <- nrow(A)
    perm0 <- qr(A)$pivot
    A1 <- A[,perm0[1:m]]
    A2 <- A[,-perm0[1:m]]
    F <- solve(A1)
    W2 <- F %*% A2
    y1 <- F %*% y
    nit <- 0
    while(any(y1 < 0)) {
        ij <- get_piv(W2, y1)
        perm0 <- swape(perm0,ij$i, m + ij$j)
        piv <- pivot(W2, y1, ij$i, ij$j)
        W2 <- piv$W2
        y1 <- piv$y1
        nit <- nit + 1
    }
    return(list(perm=perm0,y=y1,W2=W2,nit=nit))
}

rfirst_vert <- function(A,y) {
    m <- nrow(A)
    perm0 <- qr(A)$pivot
    A1 <- A[,perm0[1:m]]
    A2 <- A[,-perm0[1:m]]
    F <- solve(A1)
    W2 <- F %*% A2
    y1 <- F %*% y
    nit <- 0
    while(any(y1 < 0)) {
        ij <- get_rpiv(W2, y1)
        perm0 <- swape(perm0,ij$i, m + ij$j)
        piv <- pivot(W2, y1, ij$i, ij$j)
        W2 <- piv$W2
        y1 <- piv$y1
        nit <- nit + 1
    }
    return(list(perm=perm0,y=y1,W2=W2,nit=nit))
}

rget_sepv <- function(A,y,perm,c) {
    perm <- swape(perm,1,c)
    perm0 <- qr(A[,perm])$pivot
    perm0 <- perm[perm0]
    F <- solve(A[,perm0[1:m]])
    W2 <- F %*% A[,-(perm0[1:m])]
    y1 <- F %*% y
    while(any(W2[1,] > 0)) {
print(y1)
        Jpos <- which(W2[1,] > 0)
        jpos <- Jpos[runif(1,min=1,max=length(Jpos))]
        if (any(W2[2:m,jpos] > 0)) {
            Ipos <- 1 + which(W2[2:m,jpos] > 0)
            i <- Ipos[runif(1,min=1,max=length(Ipos))]
            piv <- pivot(W2, y, i, jpos)
        } else {
            Ineg <- which(W2[,jpos] < 0)
            i <- Ineg[runif(1,min=1,max=length(Ineg))]
            piv <- pivot(W2, y, i, jpos)
        }
        W2 <- piv$W2
        y1 <- piv$y1
        perm0 <- swape(perm0,i,m+jpos)
    }
    return(list(W2=W2,y1=y1,perm=perm0))
}
