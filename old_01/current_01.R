source("get_data.R")
source("myfunc.R")
fv <- first_vert(A,y)

F <- solve(A1)
W2 <- F %*% A2
y1 <- F %*% y

get_piv <- function(W2, y1) {
    miny0 <- min(-y1 * (y1 < 0))
    ipiv <- -1
    jpiv <- -1
    for (j in 1:ncol(W2)) {
        w <- W2[,j]
        Ineg <- which(y1 < 0 & w != 0);
        Izer <- which(y1 == 0 & w != 0);
        Ipos <- which(y1 > 0 & w != 0);
        if (length(Ineg) > 0 & all(w[Izer] < 0) & all(w[Ipos] > 0)) {
            for (i in which(w[Ineg] < 0)) {
                if (all(y[Ipos] * w[i] <= y[i] * w[Ipos])) {
                    return(list(i=i,j=j))
                } else {
                    y2 <- y1 - y1[i]/w[i] * w
                    miny2 <- which(-y2 * (y2 < 0))
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
