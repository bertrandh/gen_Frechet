ftxt <- as.matrix(read.table("A.txt"))
ytxt <- as.matrix(read.table("y.txt"))
n <- ncol(Atxt)
m <- nrow(Atxt)
u <- 7
A <- Atxt[,-u]
y <- ytxt
n <- n-1
w <- Atxt[,u,drop=F]
source("myfunc.R")

cperm1 <- qr(A)$pivot[1:m]
cperm2 <- sort(setdiff(1:n, cperm1))
cperm <- c(cperm1, cperm2)
rperm <- 1:m

F <- solve(A[,cperm1])
A <- F %*% A[,c(cperm1, cperm2)]
y <- F %*% y
w <- F %*% w

t <- rep(0,m)
Y <- matrix(0,m,m)

finish <- 0
while (finish == 0) {
    for (i in 1:m) {
        t[i] <- y[i]/w[i]
        Y[,i] <- round(y - t[i] * w)
    }
    tperm <- sort.int(t,index.return=T)$ix
    t <- t[tperm]
    itperm <- rep(0,m)
    itperm[tperm] <- 1:m
    cperm1 <- cperm1[tperm]
    cperm <- c(cperm1, cperm2)
    rperm <- rperm[tperm]
    A <- A[tperm,]
    A[,1:m] <- A[,tperm]
    y <- y[tperm]
    w <- w[tperm]
    Y <- Y[tperm,tperm]
    print(y)

    for (i in 1:m) {
        if (any(Y[,i] < 0)) next
        if (any(A[i,] < 0)) next
        cat("tmin = ",t[i], "and i = ",i,"\n")
        finish <- 1
    }
    if (finish == 1) {
        print("finished!")
        break
    }
    i <- 0
    j <- 0
    while (i <= n & j == 0) {
        i <- i + 1
        J <- which(A[i,] < 0)
        lJ <- length(J)
        if (lJ > 0) {
            j <- J[ceiling(runif(1,0,lJ))]
            cat("i = ",i,"\n")
        }
    }
    if (i > n) print("Error: i > n")
    A <- pivot(cbind(A,y,w),i,j)
    cperm <- swapel(cperm,i,j)
    cperm1 <- cperm[1:m]
    cperm2 <- cperm[(m+1):n]
    w <- A[,n+2]
    y <- A[,n+1]
    A <- A[,1:n]
}


