# Log rank test of stage 1 vs stage 2, and plot KM curves

library(OIsurv)

cli = read.csv("rdata.txt", header = FALSE);
surTime = cli[[1]]
death = cli[[2]]
stage = cli[[4]]

# only keep stage 1 and stage 2
ind = stage == "Stage I" | stage == "Stage II"
surTime = surTime[ind]
death = death[ind]
stage = stage[ind]

mySurv = Surv(surTime, death);

log1 = survdiff(mySurv ~ stage) # stage
p = pchisq(log1$chisq, 1, lower.tail=FALSE)

logrankRes = data.frame(name = "stage 1 vs. 2", p = p) # write logrank res to txt file
write.table(logrankRes, file = "logrankRes_stage1vs2.txt",
            row.names = F, col.names = T, sep="\t") 

# save KM curve to png file for stage
types = c("Stage I", "Stage II")
ng = length(unique(stage))
fit = survfit(mySurv ~ stage)
fname = "KMCurve/stage1vs2.png"

# compute the number of patients in each group
num = numeric()
for(j in 1:ng){
  num[j] = sum(stage==types[j])
}

png(filename = fname, width = 5.5, height = 5.5,
    units = "cm", res = 300, pointsize = 7)
plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", lty = 1:ng,
	col = 1:ng, cex = 0.5)
grid()
legend(x = "topright", legend = paste(types, "(", num, ")", sep = ""), lty = 1:ng,
	col = 1:ng, cex = 0.5)
text(10, 0.1, paste("p=", formatC(p, format="g", digits = 3), sep = ""), pos = 4, cex = 1)
dev.off()

print(logrankRes)


