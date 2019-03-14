## SYNTAX FILE 2 - FINAL DATA CREATION                       ##

cloudstor <- "C:/Users/z3312911/Cloudstor/" # change to master file path
.libPaths(paste0(cloudstor,"R Library"))

set.seed(710321,kind="L'Ecuyer-CMRG")

# Load imputed data created in Syntax File 2
load(file=paste0(cloudstor,"PhD/Paper 6 - POINT Application/Data/Imputed data.RData"))

# Define variable lists which determine variable order and varying variables for reshape
flist <- c("Participant_ID","time","b_Actual_age","b_employ3","b_sex","b_edu","b_MaritalStatus","b_Arth_12mR","b_Back_12mR","b_Head_12mR","b_Visc_12mR","b_Fibro_12mR","b_Cmplx_12mR","b_Shing_12mR","b_pain_duration_yrs3","PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","alc_12m","binge","auditcprob","BPI_PScore","BPI_interference")
varlist <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","alc_12m","binge","auditcprob","BPI_PScore","BPI_interference")

# Reshape data to wide format for use by LTMLE
impdatawide <- lapply(JMimpute_alc, function (x) {
  x <- as.data.frame(x[,flist])
  
  wide <- reshape(x,
                  v.names=varlist,
                  idvar="Participant_ID",
                  timevar="time",
                  sep="",
                  direction="wide")
  wide
})

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

dataalcuse1 <- lapply(impdatawide, function (x) {
  x <- x[,c(-27,-28,-30,-43,-44,-46,-59,-60,-62,-75,-76,-78,-91,-92,-94,-107,-108,-110)]
  x
})
dataalcuse2 <- lapply(impdatawide, function (x) {
  x <- x[,c(-27,-28,-29,-43,-44,-45,-59,-60,-61,-75,-76,-77,-91,-92,-93,-107,-108,-109)]
  x
})

databinge1 <- lapply(impdatawide, function (x) {
  x <- x[,c(-26,-28,-30,-42,-44,-46,-58,-60,-62,-74,-76,-78,-90,-92,-94,-106,-108,-110)]
  x
})
databinge2 <- lapply(impdatawide, function (x) {
  x <- x[,c(-26,-28,-29,-42,-44,-45,-58,-60,-61,-74,-76,-77,-90,-92,-93,-106,-108,-109)]
  x
})

dataauditprob1 <- lapply(impdatawide, function (x) {
  x <- x[,c(-26,-27,-30,-42,-43,-46,-58,-59,-62,-74,-75,-78,-90,-91,-94,-106,-107,-110)]
  x
})
dataauditprob2 <- lapply(impdatawide, function (x) {
  x <- x[,c(-26,-27,-29,-42,-43,-45,-58,-59,-61,-74,-75,-77,-90,-91,-93,-106,-107,-109)]
  x
})

save(dataalcuse1,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Data/alcuse1.RData"))
save(dataalcuse2,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Data/alcuse2.RData"))
save(databinge1,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Data/binge1.RData"))
save(databinge2,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Data/binge2.RData"))
save(dataauditprob1,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Data/auditprob1.RData"))
save(dataauditprob2,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Data/auditprob2.RData"))

rm(list=c("JMimpute_alc","impdatawide","dataalcuse1","dataalcuse2","databinge1","databinge2","dataauditprob1","dataauditprob2"))