# Lasso-Cox within early stage (stage1,2)
library("OIsurv")

mydata = read.table("rdata.txt", header = FALSE)
mydata = data.matrix(mydata)
s1 = nrow(mydata)
s2 = ncol(mydata)
surTime = mydata[, 1]
death = mydata[, 2]
lasso = mydata[, 5]

mySurv <- Surv(surTime, death)
# focus on early stage
data2 = read.csv("D:/matlab/KIRC/integration_image_RNA_SM2/clinical/rdata.txt", header = FALSE);
stage = data2[[4]]; # stage 1, 2, 3, 4
ind = stage=="Stage I" | stage=="Stage II"
mySurv = mySurv[ind, ]
lasso = lasso[ind]

# logrank
log1 = survdiff(mySurv ~ lasso)
p = pchisq(log1$chisq, 1, lower.tail=FALSE)
print(p)

# plot KM curve
fit = survfit(mySurv ~ lasso)
n1 = sum(lasso==1)
leg1 = paste("Low risk(", n1, ")", sep = "")
n2 = sum(lasso==2)
leg2 = paste("High risk(", n2, ")", sep = "")

png(filename = "lassoInEarlyStage.png", width = 5.5, height = 5.5,
	units = "cm", res = 300, pointsize = 7)
plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", lty = 1:2,
	col = 1:2, cex = 0.5)
grid()
legend(x = "topright", legend = c(leg1, leg2), lty = 1:2,
	col = 1:2, cex = 0.55)
text(10, 0.1, paste("p=", formatC(p, format="g", digits = 3), sep = ""),
	pos = 4, cex = 1)
dev.off()



