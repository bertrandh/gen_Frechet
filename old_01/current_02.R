source("get_data.R")
source("myfunc.R")
#fv <- first_vert(A,y)

F <- solve(A1)
W2 <- F %*% A2
y1 <- F %*% y
