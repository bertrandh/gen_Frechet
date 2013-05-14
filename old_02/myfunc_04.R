swapel <- function(vect,i,j) {
    temp <- vect[i]
    vect[i] <- vect[j]
    vect[j] <- temp
    return(vect)
}

pivot <- function(A, i, j) {
    m <- nrow(A)
    w <- A[,j]
    U <- diag(m)
    U[,i] <- -w/w[i]
    U[i,i] <- 1/w[i]
    A <- U %*% A
    w <- A[,i]
    A[,i] <- round(A[,j])
    A[,j] <- w
    return(A)
}

normalize <- function(A,perm1) {
    n <- ncol(A)
    F <- solve(A[,perm1])
    perm <- c(perm1, sort(setdiff(1:n, perm1)))
    print(length(perm))
    A <- F %*% A[,perm]
    return(A)
}

get_sep_hyp <- function(A,y) {
    I <- which(y < 0)
    lI <- length(I)
    if (lI == 0) return(list(i=-1,A=A,y=y))
    for (i in I) {
        if (any(A[i,] < 0)) next
        return(list(i=i,A=A,y=y))
    }
    i <- I[ceiling(runif(1,0,lJ))]
    J <- which(A[i,] < 0)
    j <- J[ceiling(runif(1,0,length(J)))]
    n <- ncol(A)
    #browser()
    A <- pivot(cbind(A,y),i,j)
    y <- A[,n+1]
    A <- A[,1:n]
    return(get_sep_hyp(A,y))
}




