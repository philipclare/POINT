
## SYNTAX FILE 6 - NAIVE ANALYSIS                            ##

cloudstor <- "C:/Users/z3312911/Cloudstor/" # change to master file path
.libPaths(paste0(cloudstor,"R Library"))

libs <- c("parallel","Amelia","lme4")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}

library ("parallel")
library("Amelia")

set.seed(546432,kind="L'Ecuyer-CMRG")

# Load imputed data
load(file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Data/Imputed data - wide.RData"))

# Set up parallel processing to reduce analysis time
numcores <- future::availableCores()-2
cl <- makeCluster(numcores)
parallel::clusterEvalQ(cl, cloudstor <- "C:/Users/z3312911/Cloudstor/")
parallel::clusterEvalQ(cl, .libPaths(paste0(cloudstor,"R Library")))
parallel::clusterEvalQ(cl, library("lme4"))
parallel::clusterEvalQ(cl, Lvars <- "time + b_Actual_age + b_employ3 + b_sex + b_edu + b_MaritalStatus + b_Arth_12mR + b_Back_12mR + b_Head_12mR + b_Visc_12mR + b_Fibro_12mR + b_Cmplx_12mR + b_Shing_12mR + b_pain_duration_yrs3 + PHQ9_Mod_sev + GADMod2Sev + Antidepressant_week + Antipsychotic_week + benzo_week + Nonopioid_analgesic_week + Pregabalin_week + can_12m + cig_12m + opioid90")
parallel::clusterEvalQ(cl, varlist_1 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","alc_12m","acum","BPI_PScore"))
parallel::clusterEvalQ(cl, varlist_2 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","alc_12m","acum","BPI_interference"))
parallel::clusterEvalQ(cl, varlist_3 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","binge","acum","BPI_PScore"))
parallel::clusterEvalQ(cl, varlist_4 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","binge","acum","BPI_interference"))
parallel::clusterEvalQ(cl, varlist_5 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","auditcprob","acum","BPI_PScore"))
parallel::clusterEvalQ(cl, varlist_6 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","auditcprob","acum","BPI_interference"))

## ANY ALCOHOL CONSUMPTION ANALYSIS ##
rialc_12m1 <- parLapply(cl=cl, impdatawide, function (x) {
  x <- x[,-15]
  x <- x[,c(-27,-28,-30,-43,-44,-46,-59,-60,-62,-75,-76,-78,-91,-92,-94,-107,-108,-110)]
  x$acum1 <- x$alc_12m1
  x$acum2 <- x$acum1 + x$alc_12m2
  x$acum3 <- x$acum2 + x$alc_12m3
  x$acum4 <- x$acum3 + x$alc_12m4
  x$acum5 <- x$acum4 + x$alc_12m5
  x$acum6 <- x$acum5 + x$alc_12m6
  x <- x[c(1:26,93,27:39,94,40:52,95,53:65,96,66:78,97,79:91,98,92)]
  x <- reshape(x,
               v.names=varlist_1,
               varying=list(c(15,29,43,57,71,85),
                            c(16,30,44,58,72,86),
                            c(17,31,45,59,73,87),
                            c(18,32,46,60,74,88),
                            c(19,33,47,61,75,89),
                            c(20,34,48,62,76,90),
                            c(21,35,49,63,77,91),
                            c(22,36,50,64,78,92),
                            c(23,37,51,65,79,93),
                            c(24,38,52,66,80,94),
                            c(25,39,53,67,81,95),
                            c(26,40,54,68,82,96),
                            c(27,41,55,69,83,97),
                            c(28,42,56,70,84,98)),
               idvar="Participant_ID",
               timevar="time",
               sep="",
               direction="long")
  x
  formula <- paste0(paste0("BPI_PScore ~ acum + ",Lvars)," + (1|Participant_ID)")
  lmer(data=x,
       formula=formula)
})

rialc_12mco1 <- matrix(unlist(parLapply(cl=cl, rialc_12m1, function (x) {coef(summary(x))[2]})),nrow=50,ncol=1)
rialc_12mse1 <- matrix(unlist(parLapply(cl=cl, rialc_12m1, function (x) {vcov(x)[2,2]})),nrow=50,ncol=1)
rialc_12mres1 <- matrix(unlist(mi.meld(q=rialc_12mco1, se=rialc_12mse1)),nrow=1,ncol=2)

rialc_12m2 <- parLapply(cl=cl, impdatawide, function (x) {
  x <- x[,-15]
  x <- x[,c(-27,-28,-29,-43,-44,-45,-59,-60,-61,-75,-76,-77,-91,-92,-93,-107,-108,-109)]
  x$acum1 <- x$alc_12m1
  x$acum2 <- x$acum1 + x$alc_12m2
  x$acum3 <- x$acum2 + x$alc_12m3
  x$acum4 <- x$acum3 + x$alc_12m4
  x$acum5 <- x$acum4 + x$alc_12m5
  x$acum6 <- x$acum5 + x$alc_12m6
  x <- x[c(1:26,93,27:39,94,40:52,95,53:65,96,66:78,97,79:91,98,92)]
  x <- reshape(x,
               v.names=varlist_2,
               varying=list(c(15,29,43,57,71,85),
                            c(16,30,44,58,72,86),
                            c(17,31,45,59,73,87),
                            c(18,32,46,60,74,88),
                            c(19,33,47,61,75,89),
                            c(20,34,48,62,76,90),
                            c(21,35,49,63,77,91),
                            c(22,36,50,64,78,92),
                            c(23,37,51,65,79,93),
                            c(24,38,52,66,80,94),
                            c(25,39,53,67,81,95),
                            c(26,40,54,68,82,96),
                            c(27,41,55,69,83,97),
                            c(28,42,56,70,84,98)),
               idvar="Participant_ID",
               timevar="time",
               sep="",
               direction="long")
  lmer(data=x,
       formula=paste0(paste0("BPI_interference ~ acum + ",Lvars)," + (1|Participant_ID)"))
})
rialc_12mco2 <- matrix(unlist(parLapply(cl=cl, rialc_12m2, function (x) {coef(summary(x))[2]})),nrow=50,ncol=1)
rialc_12mse2 <- matrix(unlist(parLapply(cl=cl, rialc_12m2, function (x) {vcov(x)[2,2]})),nrow=50,ncol=1)
rialc_12mres2 <- matrix(unlist(mi.meld(q=rialc_12mco2, se=rialc_12mse2)),nrow=1,ncol=2)

rm(list=c("rialc_12m1","rialc_12mco1","rialc_12mse1",
          "rialc_12m2","rialc_12mco2","rialc_12mse2"))

## BINGE DRINKING ANALYSIS ##
ribinge1 <- parLapply(cl=cl, impdatawide, function (x) {
  x <- x[,-15]
  x <- x[,c(-26,-28,-30,-42,-44,-46,-58,-60,-62,-74,-76,-78,-90,-92,-94,-106,-108,-110)]
  x$acum1 <- 0
  x$acum2 <- x$binge2
  x$acum3 <- x$acum2 + x$binge3
  x$acum4 <- x$acum3 + x$binge4
  x$acum5 <- x$acum4 + x$binge5
  x$acum6 <- x$acum5 + x$binge6
  x <- x[c(1:26,93,27:39,94,40:52,95,53:65,96,66:78,97,79:91,98,92)]
  x <- reshape(x,
               v.names=varlist_3,
               varying=list(c(15,29,43,57,71,85),
                            c(16,30,44,58,72,86),
                            c(17,31,45,59,73,87),
                            c(18,32,46,60,74,88),
                            c(19,33,47,61,75,89),
                            c(20,34,48,62,76,90),
                            c(21,35,49,63,77,91),
                            c(22,36,50,64,78,92),
                            c(23,37,51,65,79,93),
                            c(24,38,52,66,80,94),
                            c(25,39,53,67,81,95),
                            c(26,40,54,68,82,96),
                            c(27,41,55,69,83,97),
                            c(28,42,56,70,84,98)),
               idvar="Participant_ID",
               timevar="time",
               sep="",
               direction="long")
  x <- x[!x$time==1,]
  lmer(data=x,
       formula=paste0(paste0("BPI_PScore ~ acum + ",Lvars)," + (1|Participant_ID)"))
})
ribingeco1 <- matrix(unlist(parLapply(cl=cl, ribinge1, function (x) {coef(summary(x))[2]})),nrow=50,ncol=1)
ribingese1 <- matrix(unlist(parLapply(cl=cl, ribinge1, function (x) {vcov(x)[2,2]})),nrow=50,ncol=1)
ribingeres1 <- matrix(unlist(mi.meld(q=ribingeco1, se=ribingese1)),nrow=1,ncol=2)

ribinge2 <- parLapply(cl=cl, impdatawide, function (x) {
  x <- x[,-15]
  x <- x[,c(-26,-28,-29,-42,-44,-45,-58,-60,-61,-74,-76,-77,-90,-92,-93,-106,-108,-109)]
  x$acum1 <- 0
  x$acum2 <- x$binge2
  x$acum3 <- x$acum2 + x$binge3
  x$acum4 <- x$acum3 + x$binge4
  x$acum5 <- x$acum4 + x$binge5
  x$acum6 <- x$acum5 + x$binge6
  x <- x[c(1:26,93,27:39,94,40:52,95,53:65,96,66:78,97,79:91,98,92)]
  x <- reshape(x,
               v.names=varlist_4,
               varying=list(c(15,29,43,57,71,85),
                            c(16,30,44,58,72,86),
                            c(17,31,45,59,73,87),
                            c(18,32,46,60,74,88),
                            c(19,33,47,61,75,89),
                            c(20,34,48,62,76,90),
                            c(21,35,49,63,77,91),
                            c(22,36,50,64,78,92),
                            c(23,37,51,65,79,93),
                            c(24,38,52,66,80,94),
                            c(25,39,53,67,81,95),
                            c(26,40,54,68,82,96),
                            c(27,41,55,69,83,97),
                            c(28,42,56,70,84,98)),
               idvar="Participant_ID",
               timevar="time",
               sep="",
               direction="long")
  x <- x[!x$time==1,]
  lmer(data=x,
       formula=paste0(paste0("BPI_interference ~ acum + ",Lvars)," + (1|Participant_ID)"))
})
ribingeco2 <- matrix(unlist(parLapply(cl=cl, ribinge2, function (x) {coef(summary(x))[2]})),nrow=50,ncol=1)
ribingese2 <- matrix(unlist(parLapply(cl=cl, ribinge2, function (x) {vcov(x)[2,2]})),nrow=50,ncol=1)
ribingeres2 <- matrix(unlist(mi.meld(q=ribingeco2, se=ribingese2)),nrow=1,ncol=2)

rm(list=c("ribinge1","ribingeco1","ribingese1",
          "ribinge2","ribingeco2","ribingese2"))

## AUDIT-C HAZARDOUS DRINKING ANALYSIS ##
riaudit1 <- parLapply(cl=cl, impdatawide, function (x) {
  x <- x[,-15]
  x <- x[,c(-26,-27,-30,-42,-43,-46,-58,-59,-62,-74,-75,-78,-90,-91,-94,-106,-107,-110)]
  x$acum1 <- 0
  x$acum2 <- x$auditcprob2
  x$acum3 <- x$acum2 + x$auditcprob3
  x$acum4 <- x$acum3 + x$auditcprob4
  x$acum5 <- x$acum4 + x$auditcprob5
  x$acum6 <- x$acum5 + x$auditcprob6
  x <- x[c(1:26,93,27:39,94,40:52,95,53:65,96,66:78,97,79:91,98,92)]
  x <- reshape(x,
               v.names=varlist_5,
               varying=list(c(15,29,43,57,71,85),
                             c(16,30,44,58,72,86),
                             c(17,31,45,59,73,87),
                             c(18,32,46,60,74,88),
                             c(19,33,47,61,75,89),
                             c(20,34,48,62,76,90),
                             c(21,35,49,63,77,91),
                             c(22,36,50,64,78,92),
                             c(23,37,51,65,79,93),
                             c(24,38,52,66,80,94),
                             c(25,39,53,67,81,95),
                             c(26,40,54,68,82,96),
                             c(27,41,55,69,83,97),
                             c(28,42,56,70,84,98)),
               idvar="Participant_ID",
               timevar="time",
               sep="",
               direction="long")
  x <- x[!x$time==1,]
  lmer(data=x,
       formula=paste0(paste0("BPI_PScore ~ acum + ",Lvars)," + (1|Participant_ID)"))
})
riauditco1 <- matrix(unlist(parLapply(cl=cl, riaudit1, function (x) {coef(summary(x))[2]})),nrow=50,ncol=1)
riauditse1 <- matrix(unlist(parLapply(cl=cl, riaudit1, function (x) {vcov(x)[2,2]})),nrow=50,ncol=1)
riauditres1 <- matrix(unlist(mi.meld(q=riauditco1, se=riauditse1)),nrow=1,ncol=2)

riaudit2 <- parLapply(cl=cl, impdatawide, function (x) {
  x <- x[,-15]
  x <- x[,c(-26,-27,-29,-42,-43,-45,-58,-59,-61,-74,-75,-77,-90,-91,-93,-106,-107,-109)]
  x$acum1 <- 0
  x$acum2 <- x$auditcprob2
  x$acum3 <- x$acum2 + x$auditcprob3
  x$acum4 <- x$acum3 + x$auditcprob4
  x$acum5 <- x$acum4 + x$auditcprob5
  x$acum6 <- x$acum5 + x$auditcprob6
  x <- x[c(1:26,93,27:39,94,40:52,95,53:65,96,66:78,97,79:91,98,92)]
  x <- reshape(x,
               v.names=varlist_6,
               varying=list(c(15,29,43,57,71,85),
                            c(16,30,44,58,72,86),
                            c(17,31,45,59,73,87),
                            c(18,32,46,60,74,88),
                            c(19,33,47,61,75,89),
                            c(20,34,48,62,76,90),
                            c(21,35,49,63,77,91),
                            c(22,36,50,64,78,92),
                            c(23,37,51,65,79,93),
                            c(24,38,52,66,80,94),
                            c(25,39,53,67,81,95),
                            c(26,40,54,68,82,96),
                            c(27,41,55,69,83,97),
                            c(28,42,56,70,84,98)),
               idvar="Participant_ID",
               timevar="time",
               sep="",
               direction="long")
  x <- x[!x$time==1,]
  lmer(data=x,
       formula=paste0(paste0("BPI_interference ~ acum + ",Lvars)," + (1|Participant_ID)"))
})
riauditco2 <- matrix(unlist(parLapply(cl=cl, riaudit2, function (x) {coef(summary(x))[2]})),nrow=50,ncol=1)
riauditse2 <- matrix(unlist(parLapply(cl=cl, riaudit2, function (x) {vcov(x)[2,2]})),nrow=50,ncol=1)
riauditres2 <- matrix(unlist(mi.meld(q=riauditco2, se=riauditse2)),nrow=1,ncol=2)

rm(list=c("riaudit1","riauditco1","riauditse1",
          "riaudit2","riauditco2","riauditse2"))

riresults <- matrix(c(rialc_12mres1,rialc_12mres2,ribingeres1,ribingeres2,riauditres1,riauditres2),byrow=TRUE,ncol=4,nrow=3)
rownames(riresults) <- c("Any alcohol consumption","Binge drinking","AUDIT-C Problematic drinking")
colnames(riresults) <- c("BPI Pain Score log(OR)","BPI Pain Score SE","BPI Pain Int log(OR)","BPI Pain Int SE")

save(riresults,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/riresults.RData"))

