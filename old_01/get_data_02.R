A <- as.matrix(read.table("A.txt"))
y <- as.matrix(read.table("y.txt"))
m <- nrow(A)
n <- ncol(A)
d <- n-m
perm0 <- qr(A)$pivot
A <- A[,perm0]
A1 <- A[,1:m]
A2 <- A[,(m+1):n]

