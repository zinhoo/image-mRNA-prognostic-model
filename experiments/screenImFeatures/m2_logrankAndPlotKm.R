# Do log rank test and plot KM curves for survival-associated image features
library(OIsurv)

imUni = read.table("rdata.txt");
surTime = imUni[[1]];
death = imUni[[2]];
mySurv = Surv(surTime, death);
s2 = ncol(imUni);

p = numeric(s2-2)
for(i in 3:s2){	
	group = imUni[[i]];
	if(length(unique(group)) == 1){
		p[i-2] = 3.14
	}else{
		log1 = survdiff(mySurv ~ group)
		p[i-2] = pchisq(log1$chisq, 1, lower.tail=FALSE)
	}		
}

logrankRes = cbind(which(p<0.05), p[p<0.05])
print(logrankRes)
write.table(logrankRes, file = "im_logrankRes.txt",
	row.names = F, col.names = F, sep="\t")

# plot KM curves for features with p value less than 0.05
s1res = nrow(logrankRes)
for(i in 1:s1res){
	feaInd = logrankRes[i, 1]
	group = imUni[[feaInd+2]]
	ng = length(unique(group))
	n1 = sum(group==1)
	leName1 = paste("Low group(", n1, ")", sep = "")
	n2 = sum(group==2)
	leName2 = paste("High group(", n2, ")", sep = "")
	
	fit = survfit(mySurv ~ group)
	fname = paste("KMCurve/", feaInd, ".png", sep = "")
  
	png(filename = fname, width = 5.5, height = 5.5,
		units = "cm", res = 300, pointsize = 7)
	plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", lty = 1:ng,
		col = 1:ng, cex = 0.5)
	grid()
	legend(x = "topright", legend = c(leName1, leName2), lty = 1:ng,
		col = 1:ng, cex = 0.65)
	text(10, 0.1, paste("p=", formatC(logrankRes[i, 2], format="g", digits = 3), sep = ""),
		pos = 4, cex = 1)
	dev.off()
}	
