
## SYNTAX FILE 2 - DESCRIPTIVE STATISTICS                    ##

cloudstor <- "C:/Users/z3312911/Cloudstor/" # change to master file path
.libPaths(paste0(cloudstor,"R Library"))

set.seed(710321,kind="L'Ecuyer-CMRG")

# Load imputed data
load(file=paste0(cloudstor,"PhD/Paper 6 - POINT Application/Data/Imputed data - long.RData"))

# Calculate proportions for 2nd (longitudinal) descriptives table
alc_use_prop1 <- matrix(unlist(lapply(impdatawide, function(x) {
  p2 <- mean(x$freq12)
  p3 <- mean(x$freq13)
  p4 <- mean(x$freq14)
  p5 <- mean(x$freq15)
  p6 <- mean(x$freq16)
  c(p2,p3,p4,p5,p6)
})),ncol=5,byrow=TRUE)
alc_use_prop2 <- matrix(unlist(lapply(impdatawide, function(x) {
  p2 <- mean(x$freq22)
  p3 <- mean(x$freq23)
  p4 <- mean(x$freq24)
  p5 <- mean(x$freq25)
  p6 <- mean(x$freq26)
  c(p2,p3,p4,p5,p6)
})),ncol=5,byrow=TRUE)
binge_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p2 <- mean(x$binge2)
  p3 <- mean(x$binge3)
  p4 <- mean(x$binge4)
  p5 <- mean(x$binge5)
  p6 <- mean(x$binge6)
  c(p2,p3,p4,p5,p6)
})),ncol=5,byrow=TRUE)
auditcprob_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p2 <- mean(x$auditcprob2)
  p3 <- mean(x$auditcprob3)
  p4 <- mean(x$auditcprob4)
  p5 <- mean(x$auditcprob5)
  p6 <- mean(x$auditcprob6)
  c(p2,p3,p4,p5,p6)
})),ncol=5,byrow=TRUE)
alc_pain_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p2 <- mean(as.numeric(x$alc_pain_12m2))
  p3 <- mean(as.numeric(x$alc_pain_12m3))
  p4 <- mean(as.numeric(x$alc_pain_12m4))
  p5 <- mean(as.numeric(x$alc_pain_12m5))
  p6 <- mean(as.numeric(x$alc_pain_12m6))
  c(p2,p3,p4,p5,p6)
})),ncol=5,byrow=TRUE)
opioid90_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p2 <- mean(x$opioid902)
  p3 <- mean(x$opioid903)
  p4 <- mean(x$opioid904)
  p5 <- mean(x$opioid905)
  p6 <- mean(x$opioid906)
  c(p2,p3,p4,p5,p6)
})),ncol=5,byrow=TRUE)

alc_use_prop_res1 <- colMeans(alc_use_prop1)
alc_use_prop_res2 <- colMeans(alc_use_prop2)
binge_prop_res <- colMeans(binge_prop)
auditcprob_prop_res <- colMeans(auditcprob_prop)
alc_pain_prop_res <- colMeans(alc_pain_prop)
opioid90_prop_res <- colMeans(opioid90_prop)

longprop <- rbind(alc_use_prop_res1,alc_use_prop_res2,binge_prop_res,auditcprob_prop_res,alc_pain_prop_res,opioid90_prop_res)
rm(list=c("alc_use_prop1","alc_use_prop2","binge_prop","auditcprob_prop","alc_pain_prop","opioid90_prop",
          "alc_use_prop_res1","alc_use_prop_res2","binge_prop_res","auditcprob_prop_res","alc_pain_prop_res","opioid90_prop_res"))

save(longprop,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/longprop.RData"))
