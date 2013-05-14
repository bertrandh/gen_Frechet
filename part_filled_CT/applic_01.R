arr_ind <- function(k, m1, m2) {
    j <- ceiling(k/m1)
    i <- k - m1*(j-1)
    return(list(i=i,j=j))
}

make_ctex <- function(m1,m2,maxentry,distrib,nex) {
    source("../gfb_func.R")
    A <- matrix(0, nrow=m1+m2, ncol<-m1*m2)
    for (i in 1:m1) {
        A[i,((i-1)*m2+1):(i*m2)] <- 1
    }
    for (i in 1:m2) {
        A[m1+i, seq(from=i, to=m1*m2, by=m2)] <- 1
    }
    Ainit <- A
    Tindx <- matrix(1:(m1*m2), nrow=m1, ncol=m2)
    nrun <- 1
    while (1==1) {
        Tcell <- matrix(1, nrow=m1, ncol=m2)
        for (i in 1:floor(min(m1,m2)/2)) {
            Tcell[(2*i-1):(2*i),(2*i-1):(2*i)] <- 0
        }
        Tcell[(m1-1):m1,(m2-1):m2] <- 0
        if (m1 * m2 %% 2 == 1) {
            Tcell[m1-1, m2-1] <- 1
        }
        Tcell <- Tcell[sample.int(m1),sample.int(m2)]
        if (Tcell[1,1] == 1) {
            j <- which(Tcell[1,] == 0)[1]
            Tcell[,c(1,j)] <- Tcell[,c(j,1)]
        }
        idxfill <- which(Tcell == 1)
        nfill <- length(idxfill)

        if (distrib == "unif") {
            T0 <-  round(matrix(runif(m1*m2,0,maxentry), nrow=m1, ncol=m2))
        } else if (distrib == "beta") {
            T0 <-  round(matrix(maxentry*rbeta(m1*m2,0.1, 0.1), nrow=m1, ncol=m2))
        } else return()
        yinit <- A %*% as.vector(T0)

        ctab <- matrix("0", nrow=m1+1, ncol=m2+1)
        ctab[1:m1,m2+1] <- yinit[1:m1]
        ctab[1:m1, m2+1] <- yinit[(m1+1):(m1+m2)]
        ctab[m1+1,1:m2] <- yinit[1:m1]
        ctab[m1+1,m2+1] <- sum(as.numeric(ctab[1:m1,m2+1]))

        idxfill <- idxfill[sample.int(nfill, nfill)]
        idxfill0 <- idxfill
        A0 <- Ainit
        y0 <- yinit
        bds <- matrix(0, nrow=2, ncol=nfill+1)
        red <- reduceA(A0,y0,1:ncol(A0),1)
        B <- get1bd(red$A,red$y,red$perm,1)
        bds[1,1] <- B$xmin
        bds[2,1] <- B$xmax
        k <- 1
        while (length(idxfill > 0)) {
            id <- idxfill[1]
            red <- reduceA(A0,y0,1:ncol(A0),id)
            B <- get1bd(red$A,red$y,red$perm,red$id)
            t <- round((B$xmin + B$xmax)/2)
            ij <- arr_ind(idxfill0[k], m1, m2)
            ctab[ij$i,ij$j] <- t
            y0 <- y0 - t * A0[,id]
            if(any(y0 < 0)) browser()
            A0 <- A0[,-id]
            idx <- which(idxfill > id)
            idxfill[idx] <- idxfill[idx] - 1
            idxfill <- idxfill[-1]
            red <- reduceA(A0,y0,1:ncol(A0),1)
            B <- get1bd(red$A,red$y,red$perm,red$id)
            bds[1,k+1] <- B$xmin
            bds[2,k+1] <- B$xmax
            k <- k + 1
        }

        #bds[1,] <- bds[1,] + 100
        #bds[2,] <- bds[2,] + 100
        #barplot(bds, offset=-100, names.arg=c(1,idxfill0), beside=T)
        barplot(bds, names.arg=c(1,idxfill0), beside=T)
        print(Tindx * Tcell)
        yne <- readline("Good example? [y/n/e] ")
        if (substr(yne,1,1) == "y") {
            write.table(bds, paste("ct_", nex, "_bds.txt", sep=""), row.names=F, col.names=F)
            write.table(T0, paste("ct_", nex, "_T0.txt", sep=""), row.names=F, col.names=F)
            write.table(idxfill0, paste("ct_", nex, "_idxfill0.txt", sep=""), row.names=F, col.names=F)
            write.table(Tindx * Tcell, paste("ct_", nex, "_TindxTcell.txt", sep=""), row.names=F, col.names=F)
            write.table(ctab, paste("ct_", nex, "_ctab.txt", sep=""), row.names=F, col.names=F)
            pdf(paste("ct_", nex, "_barplot.pdf", sep=""))
            #barplot(bds, offset=-100, names.arg=c(1,idxfill0), beside=T)
            barplot(bds, names.arg=c(1,idxfill0), beside=T)
            dev.off()
            nex <- nex+1
        } else if (substr(yne,1,1) == "e") {
            cat("number of runs: ", nrun, "\n")
            return()
        }
        nrun <- nrun+1
    }
}



