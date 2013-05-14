Ainit <- as.matrix(read.table("A.txt"))
y <- as.matrix(read.table("y.txt"))
n <- ncol(Ainit)
m <- nrow(Ainit)
u <- 7
A <- Ainit[,-u]
n <- n-1
w <- Ainit[,u,drop=F]
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
    t.sorted <- sort.int(t,index.return=T)
    t <- t.sorted$x
    tperm <- t.sorted$ix
    itperm <- rep(0,m)
    itperm[tperm] <- 1:m
    cperm1 <- cperm1[tperm]
    cperm2 <- sort(setdiff(1:n, cperm1))
    cperm <- c(cperm1, cperm2)
    rperm <- rperm[tperm]
    A <- A[tperm,]
    A[,1:m] <- A[,tperm]
    y <- y[tperm]
    w <- w[tperm]
    Y <- Y[tperm,tperm]

    for (i in 1:m) {
        if (any(Y[,i] < 0)) next
        if (any(A[i,] < 0)) next
        print(t[i])
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
        }
    }
    if (i > n) print("Error: i > n")
    A <- pivot(cbind(A,y,w),i,j)
    cperm <- swapel(cperm,i,j)
    w <- A[,n+2]
    y <- A[,n+1]
    print(y)
    A <- A[,1:n]
}


