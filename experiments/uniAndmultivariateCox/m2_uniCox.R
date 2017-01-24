# univariate Cox regression

mydata = read.table("rdata.txt", header = FALSE)
mydata = data.matrix(mydata)
s1 = nrow(mydata)
s2 = ncol(mydata)
surTime = mydata[, 1]
death = mydata[, 2]
x = mydata[, 3:s2]

my.surv <- Surv(surTime, death)

s2_x = ncol(x)
rname = c("grade", "stage", "lasso-Cox", "CSNK2A1", "SPP1", "DEFB1",
	"PECAM1", "EDNRB", "TSPAN7", "VHL", "PBRM1", "BAP1", "SETD2", "TP53")
cname = c("exp(coef)", "exp(-coef)", "lower .95", "upper .95", "p")
resUni = matrix(nrow = s2_x, ncol = 5, dimnames = list(rname, cname))
for(i in 1:s2_x){
	coxph.fit <- coxph(my.surv ~ x[, i], method="breslow")
	mysum = summary(coxph.fit)
	resUni[i, 1:4] = mysum$conf.int
	resUni[i, 5] = mysum$coefficients[5]
}
print(resUni)

# write.table(resUni, file = "resUni.txt",
#            row.names = F, col.names = F, sep="\t")
save(resUni, file = "resUni.RData")