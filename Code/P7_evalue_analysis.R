
## SYNTAX FILE 7 - E-VALUE ANALYSIS                          ##

cloudstor <- "C:/Users/z3312911/Cloudstor/" # change to master file path
.libPaths(paste0(cloudstor,"R Library"))

libs <- c("EValue")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
library ("EValue")

# Load Cohen's d results created by primary analysis file
load(file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/cdresults.RData"))

e1 <- evalues.MD(est=cdresults[1,1], se=cdresults[1,2])
e2 <- evalues.MD(est=cdresults[1,3], se=cdresults[1,4])

e3 <- evalues.MD(est=cdresults[2,1], se=cdresults[2,2])
e4 <- evalues.MD(est=cdresults[2,3], se=cdresults[2,4])

e5 <- evalues.MD(est=cdresults[3,1], se=cdresults[3,2])
e6 <- evalues.MD(est=cdresults[3,3], se=cdresults[3,4])

evalue <- matrix(c(e1[2,1],e2[2,1],e3[2,1],e4[2,1],e5[2,1],e6[2,1]),byrow=TRUE,ncol=2,nrow=3)
rownames(evalue) <- c("Any alcohol consumption","Binge drinking","AUDIT-C Problematic drinking")
colnames(evalue) <- c("Pain Severity","Pain Interference")

save(evalue,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/evalue.RData"))
