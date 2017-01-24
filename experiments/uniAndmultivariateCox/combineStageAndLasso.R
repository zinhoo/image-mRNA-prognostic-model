# combine stage and lasso-Cox
library("OIsurv")

mydata = read.table("rdata.txt", header = FALSE)
mydata = data.matrix(mydata)
s1 = nrow(mydata)
s2 = ncol(mydata)
surTime = mydata[, 1]
death = mydata[, 2]
stage = mydata[, 4]
lasso = mydata[, 5]

my.surv <- Surv(surTime, death)