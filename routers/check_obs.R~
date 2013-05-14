obs <- read.csv("obs2.txt")
testval <- rep(FALSE, 16)
nob <- nrow(obs)/16
Y <- matrix(0, nrow=16, ncol=nob)
for (i in 1:16) {
    test <- obs[seq(from=i, to=nrow(obs), by=16),3]
    testval[i] <- all(test == obs[i,3])
}
if (any(testval) == F) {
    print("Problem:  Not all y-coordinate correspond to obs.")
} else {
    for (j in 1:nob) {
        Y[,j] <- obs[(16*(j-1)+1):(16*j),2]
    }
}
