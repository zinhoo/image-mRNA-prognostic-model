# Lasso regularized Cox regression using leave-one-out (LOO) 
# cross validation (CV).
# During each run of LOO CV, 10-fold CV is performed on training set to 
# select the best model. Then the selected model is applied to the single
# held-out test sample to predict death risk.

library(OIsurv)
library(glmnet)
set.seed(1)

ptm <- proc.time()

mydata = read.table("rdata_imGene_related_norm.txt", header = FALSE)
mydata = data.matrix(mydata)
s1 = nrow(mydata)
s2 = ncol(mydata)
mySurv = Surv(mydata[, 1], mydata[, 2]);
x = mydata[, 3:s2]

# leave-one-out CV
group = cbind(numeric(s1))
ind = 1:s1
for(i in 1:s1){
  cvfit = cv.glmnet(x[ind!=i,], mySurv[ind!=i,], family = "cox", maxit=500000)
  preTrain = predict(cvfit, newx = x[ind!=i,], s = "lambda.min", type="response")
  mv = median(preTrain)
  preTest = predict(cvfit, newx = x[ind==i,], s = "lambda.min", type="response")
  if(preTest < mv){
    group[i] = 1
  }else{
    group[i] = 2
  }
  print(i)
}

write.table(group, file = "group_loo_imGene_related_norm.txt",
            row.names = F, col.names = F, sep="\t")

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

png(filename = "KMCurve_loo_imGene_related_norm.png", width = 5.5, height = 5.5,
	units = "cm", res = 300, pointsize = 7)
plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", lty = 1:2,
	col = 1:2, cex = 0.5)
grid()
legend(x = "topright", legend = c(leg1, leg2), lty = 1:2,
	col = 1:2, cex = 0.65)
text(10, 0.1, paste("p=", formatC(p, format="g", digits = 3), sep = ""),
	pos = 4, cex = 1)
dev.off()

print(proc.time() - ptm)