# combine stage and lasso-Cox, 8 groups
library("OIsurv")

mydata = read.table("rdata.txt", header = FALSE)
mydata = data.matrix(mydata)
s1 = nrow(mydata)
s2 = ncol(mydata)
surTime = mydata[, 1]
death = mydata[, 2]
lasso = mydata[, 5]

mySurv <- Surv(surTime, death)

data2 = read.csv("D:/matlab/KIRC/integration_image_RNA_SM2/clinical/rdata.txt", header = FALSE);
stage = data2[[4]]; # stage 1, 2, 3, 4

types = c("Stage I", "Stage II", "Stage III", "Stage IV")
combineGroup = matrix(0, nrow=s1, ncol = 1)
for(i in 1:4){
	ind1 = stage==types[i] & lasso==1
	ind2 = stage==types[i] & lasso==2
	combineGroup[ind1] = i*2-1
	combineGroup[ind2] = i*2
}

# logrank
log1 = survdiff(mySurv ~ combineGroup)
p = pchisq(log1$chisq, 1, lower.tail=FALSE)
print(p)

# plot KM curve
fit = survfit(mySurv ~ combineGroup)
n1 = sum(combineGroup==1)
leg1 = paste("I & low(", n1, ")", sep = "")
n2 = sum(combineGroup==2)
leg2 = paste("I & high(", n2, ")", sep = "")
n3 = sum(combineGroup==3)
leg3 = paste("II & low(", n3, ")", sep = "")
n4 = sum(combineGroup==4)
leg4 = paste("II & high(", n4, ")", sep = "")
n5 = sum(combineGroup==5)
leg5 = paste("III & low(", n5, ")", sep = "")
n6 = sum(combineGroup==6)
leg6 = paste("III & high(", n6, ")", sep = "")
n7 = sum(combineGroup==7)
leg7 = paste("IV & low(", n7, ")", sep = "")
n8 = sum(combineGroup==8)
leg8 = paste("IV & high(", n8, ")", sep = "")

png(filename = "combineStageAndLasso_8.png", width = 5.5, height = 5.5,
	units = "cm", res = 300, pointsize = 7)
plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", lty = 1:8,
	col = 1:8, cex = 0.5)
grid()
legend(x = "topright", legend = c(leg1, leg2, leg3, leg4, leg5, leg6, leg7, leg8), lty = 1:8,
	col = 1:8, cex = 0.38)
text(10, 0.1, paste("p=", formatC(p, format="g", digits = 3), sep = ""),
	pos = 4, cex = 1)
dev.off()



