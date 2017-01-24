# Log rank test of grade and stage, and plot KM curves

library(OIsurv)

cli = read.csv("rdata.txt", header = FALSE);
surTime = cli[[1]]
death = cli[[2]]
grade = cli[[3]]
stage = cli[[4]]

mySurv = Surv(surTime, death);
p = numeric()
name = c("grade", "stage")


log1 = survdiff(mySurv ~ grade) # grade
p[1] = pchisq(log1$chisq, 1, lower.tail=FALSE)

log1 = survdiff(mySurv ~ stage) # stage
p[2] = pchisq(log1$chisq, 1, lower.tail=FALSE)

logrankRes = data.frame(name = name, p = p) # write logrank res to txt file
write.table(logrankRes, file = "logrankRes.txt",
            row.names = F, col.names = T, sep="\t") 


# save KM curve to png file for grade
types = c("G1", "G2", "G3", "G4")
ng = length(unique(grade))
fit = survfit(mySurv ~ grade)
fname = "KMCurve/grade.png"

# compute the number of patients in each group
num = numeric()
for(j in 1:ng){
  num[j] = sum(grade==types[j])
}

png(filename = fname, width = 5.5, height = 5.5,
    units = "cm", res = 300, pointsize = 7)
plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", lty = 1:ng,
	col = 1:ng, cex = 0.5)
grid()
legend(x = "topright", legend = paste(types, "(", num, ")", sep = ""), lty = 1:ng,
	col = 1:ng, cex = 0.5)
text(10, 0.1, paste("p=", formatC(p[1], format="g", digits = 3), sep = ""), pos = 4, cex = 1)
dev.off()


# save KM curve to png file for stage
types = c("Stage I", "Stage II", "Stage III", "Stage IV")
ng = length(unique(stage))
fit = survfit(mySurv ~ stage)
fname = "KMCurve/stage.png"

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
text(10, 0.1, paste("p=", formatC(p[2], format="g", digits = 3), sep = ""), pos = 4, cex = 1)
dev.off()

print(logrankRes)


