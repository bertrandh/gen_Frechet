library(lpSolve)

getlpbox <- function(A,y,perm) {
    n <- ncol(A)
    m <- nrow(A)
    xbd <- matrix(0, ncol=n, nrow=2)
    Vmin <- matrix(0,ncol=n, nrow=n)
    Vmax <- matrix(0,ncol=n, nrow=n)
    cmat <- diag(n)
    for (i in 1:n) {
        BV <- lp(direction="min", objective.in=cmat[i,], const.mat=A, const.dir="==", const.rhs=y)
        Vmin[,i] <- BV$solution
        BV <- lp(direction="max", objective.in=cmat[i,], const.mat=A, const.dir="==", const.rhs=y)
        Vmax[,i] <- BV$solution
    }
    xbd[1,] <- diag(Vmin)
    xbd[2,] <- diag(Vmax)
    return(list(xbd=xbd, Vmin=Vmin, Vmax=Vmax))
}

getlpbox_all <- function(Ainit, Yob) {
    nob <- ncol(Yob)
    perm <- qr(Ainit)$pivot
    Ap <- Ainit[,perm]
    for (i in 1:nob) {
        lpBV <- getlpbox(Ap, Yob[,i], perm)
    }
}

getbox_all <- function(Ainit, Yob) {
    nob <- ncol(Yob)
    perm <- qr(Ainit)$pivot
    Ap <- Ainit[,perm]
    for (i in 1:nob) {
        BV <- getbox(Ap, Yob[,i], perm)
    }
}


