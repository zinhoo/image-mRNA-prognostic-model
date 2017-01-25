# Lasso regularized Cox regression using 10-fold cross validation (CV).
# This program is only for examining which variables are selected by the 
# lasso Cox model. So the 10-fold CV is performed using the whole data set,
# and then the sparse regression coefficients are obtained from the selected 
# model.

library("OIsurv")
library("glmnet")

ptm <- proc.time()
set.seed(1)

mydata = read.table("rdata_imGene_related_norm.txt", header = FALSE)
mydata = data.matrix(mydata)
s1 = nrow(mydata)
s2 = ncol(mydata)
mySurv = Surv(mydata[, 1], mydata[, 2]);
x = mydata[, 3:s2]

# 10-fold CV
cvfit = cv.glmnet(x, mySurv, family = "cox")
plot(cvfit)
pred = predict(cvfit, newx = x, s = "lambda.min", type="response")

# write coefficients to file
coefMin = coef(cvfit, s = "lambda.min")
ind = which(coefMin!=0)
coefMin = cbind(ind, coefMin[ind])
write.table(coefMin, file = "coef_10CV_imGene_related_norm.txt",
            row.names = F, col.names = F, sep="\t")

group = cbind(numeric(s1))
mv = median(pred)
group[pred<mv] = 1
group[pred>=mv] = 2

# logrank
log1 = survdiff(mySurv ~ group)
p = pchisq(log1$chisq, 1, lower.tail=FALSE)
print(p)

# plot KM curve
fit = survfit(mySurv ~ group)
n1 = sum(group==1)
leg1 = paste("Low risk(", n1, ")", sep = "")
n2 = sum(group==2)
leg2 = paste("High risk(", n2, ")", sep = "")

png(filename = "KMCurve_10CV_imGene_related_norm.png", width = 5.5, height = 5.5,
	units = "cm", res = 300, pointsize = 7)
plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", lty = 1:2,
	col = 1:2, cex = 0.5)
grid()
legend(x = "topright", legend = c(leg1, leg2), lty = 1:2,
	col = 1:2, cex = 0.5)
text(10, 0.1, paste("p=", formatC(p, format="g", digits = 3), sep = ""),
	pos = 4, cex = 1)
dev.off()


print(proc.time() - ptm)