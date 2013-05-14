max0 <- function(a) {
    n <- length(a)
    if (n == 0) {
        print("Usage: minv(a,b) where a and b are two vectors of same positive length")
        return()
    }
    c <- rep(0,n)
    for (i in 1:n) {
        c[i] <- max(a[i], 0)
    }
    return(c)
}


minv <- function(a,b) {
    n <- length(a)
    if (n == 0 || length(b) != n) {
        print("Usage: minv(a,b) where a and b are two vectors of same positive length")
        return()
    }
    c <- rep(0,n)
    for (i in 1:n) {
        c[i] <- min(a[i], b[i])
    }
    return(c)
}

make_rr <- function(Y) {
    Aout <- colSums(Y[1:4,])
    Ain <- colSums(Y[9:12,])
    Bout <- colSums(Y[5:8,])
    Bin <- colSums(Y[13:16,])
    r1r2 <- (max0(Aout - Ain) + minv(Aout, Bin))/2
    r2r1 <- (max0(Bout - Bin) + minv(Bout, Ain))/2
    return(list(r1r2=r1r2, r2r1=r2r1))
}

complete_Y <- function(Y) {
    rr <- make_rr(Y)
    Y <- rbind(Y, rr$r1r2, rr$r2r1)
    return(Y)
}
