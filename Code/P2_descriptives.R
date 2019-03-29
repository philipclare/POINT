
## SYNTAX FILE 2 - DESCRIPTIVE STATISTICS                    ##

cloudstor <- "C:/Users/z3312911/Cloudstor/" # change to master file path
.libPaths(paste0(cloudstor,"R Library"))

set.seed(710321,kind="L'Ecuyer-CMRG")

# Load imputed data
load(file=paste0(cloudstor,"PhD/Paper 6 - POINT Application/Data/Imputed data - long.RData"))

# Calculate proportions for 2nd (longitudinal) descriptives table
alc_use_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p1 <- mean(x$alc_12m1)
  p2 <- mean(x$alc_12m2)
  p3 <- mean(x$alc_12m3)
  p4 <- mean(x$alc_12m4)
  p5 <- mean(x$alc_12m5)
  p6 <- mean(x$alc_12m6)
  c(p1,p2,p3,p4,p5,p6)
})),ncol=6,byrow=TRUE)
binge_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p2 <- mean(x$binge2)
  p3 <- mean(x$binge3)
  p4 <- mean(x$binge4)
  p5 <- mean(x$binge5)
  p6 <- mean(x$binge6)
  c(0,p2,p3,p4,p5,p6)
})),ncol=6,byrow=TRUE)
auditcprob_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p2 <- mean(x$auditcprob2)
  p3 <- mean(x$auditcprob3)
  p4 <- mean(x$auditcprob4)
  p5 <- mean(x$auditcprob5)
  p6 <- mean(x$auditcprob6)
  c(0,p2,p3,p4,p5,p6)
})),ncol=6,byrow=TRUE)
alc_pain_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p2 <- mean(as.numeric(x$alc_pain_12m2))
  p3 <- mean(as.numeric(x$alc_pain_12m3))
  p4 <- mean(as.numeric(x$alc_pain_12m4))
  p5 <- mean(as.numeric(x$alc_pain_12m5))
  p6 <- mean(as.numeric(x$alc_pain_12m6))
  c(0,p2,p3,p4,p5,p6)
})),ncol=6,byrow=TRUE)
opioid90_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p1 <- mean(x$opioid901)
  p2 <- mean(x$opioid902)
  p3 <- mean(x$opioid903)
  p4 <- mean(x$opioid904)
  p5 <- mean(x$opioid905)
  p6 <- mean(x$opioid906)
  c(p1,p2,p3,p4,p5,p6)
})),ncol=6,byrow=TRUE)
can_12m_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p1 <- mean(as.numeric(x$can_12m1)-1)
  p2 <- mean(as.numeric(x$can_12m2)-1)
  p3 <- mean(as.numeric(x$can_12m3)-1)
  p4 <- mean(as.numeric(x$can_12m4)-1)
  p5 <- mean(as.numeric(x$can_12m5)-1)
  p6 <- mean(as.numeric(x$can_12m6)-1)
  c(p1,p2,p3,p4,p5,p6)
})),ncol=6,byrow=TRUE)
cig_12m_prop <- matrix(unlist(lapply(impdatawide, function(x) {
  p1 <- mean(as.numeric(x$cig_12m1)-1)
  p2 <- mean(as.numeric(x$cig_12m2)-1)
  p3 <- mean(as.numeric(x$cig_12m3)-1)
  p4 <- mean(as.numeric(x$cig_12m4)-1)
  p5 <- mean(as.numeric(x$cig_12m5)-1)
  p6 <- mean(as.numeric(x$cig_12m6)-1)
  c(p1,p2,p3,p4,p5,p6)
})),ncol=6,byrow=TRUE)

alc_use_prop_res <- colMeans(alc_use_prop)
binge_prop_res <- colMeans(binge_prop)
auditcprob_prop_res <- colMeans(auditcprob_prop)
alc_pain_prop_res <- colMeans(alc_pain_prop)
opioid90_prop_res <- colMeans(opioid90_prop)
can_12m_prop_res <- colMeans(can_12m_prop)
cig_12m_prop_res <- colMeans(cig_12m_prop)

longprop <- cbind(alc_use_prop_res,binge_prop_res,auditcprob_prop_res,alc_pain_prop_res,opioid90_prop_res,can_12m_prop_res,cig_12m_prop_res)
rm(list=c("alc_use_prop","binge_prop","auditcprob_prop","alc_pain_prop","opioid90_prop","can_12m_prop","cig_12m_prop",
          "alc_use_prop_res","binge_prop_res","auditcprob_prop_res","alc_pain_prop_res","opioid90_prop_res","can_12m_prop_res","cig_12m_prop_res"))

save(longprop,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/longprop.RData"))
