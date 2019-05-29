
## SYNTAX FILE 5 - NAIVE ANALYSIS                            ##

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
numcores <- future::availableCores()-1
cl <- makeCluster(numcores)
parallel::clusterEvalQ(cl, cloudstor <- "C:/Users/z3312911/Cloudstor/")
parallel::clusterEvalQ(cl, .libPaths(paste0(cloudstor,"R Library")))
parallel::clusterEvalQ(cl, library("lme4"))
parallel::clusterEvalQ(cl, Lvars <- "time + b_Actual_age + b_employ3 + b_sex + b_edu + b_MaritalStatus + b_Arth_12mR + b_Back_12mR + b_Head_12mR + b_Visc_12mR + b_Fibro_12mR + b_Cmplx_12mR + b_Shing_12mR + b_pain_duration_yrs3 + PHQ9_Mod_sev + GADMod2Sev + Antidepressant_week + Antipsychotic_week + benzo_week + Nonopioid_analgesic_week + Pregabalin_week + can_12m + cig_12m + opioid90 + PSEQ_Score")
parallel::clusterEvalQ(cl, varlist_1 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","PSEQ_Score","alc_pain_12m","alc_12m","acum","BPI_PScore"))
parallel::clusterEvalQ(cl, varlist_2 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","PSEQ_Score","alc_pain_12m","alc_12m","acum","BPI_interference"))
parallel::clusterEvalQ(cl, varlist_3 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","PSEQ_Score","alc_pain_12m","binge","acum","BPI_PScore"))
parallel::clusterEvalQ(cl, varlist_4 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","PSEQ_Score","alc_pain_12m","binge","acum","BPI_interference"))
parallel::clusterEvalQ(cl, varlist_5 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","PSEQ_Score","alc_pain_12m","auditcprob","acum","BPI_PScore"))
parallel::clusterEvalQ(cl, varlist_6 <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","PSEQ_Score","alc_pain_12m","auditcprob","acum","BPI_interference"))

## 2+ DRINKS PER WEEK ##
rialc_12m1a <- parLapply(cl=cl, impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6))
  x$acum2 <- x$freq12
  x$acum3 <- x$acum2 + x$freq13
  x$acum4 <- x$acum3 + x$freq14
  x$acum5 <- x$acum4 + x$freq15
  x$acum6 <- x$acum5 + x$freq16
  x <- x[c(1:27,99,28:41,100,42:55,101,56:69,102,70:83,103,84:97,104,98)]
  x <- reshape(x,
               v.names=varlist_1,
               varying=list(c(15,30,45,60,75,90),
                            c(16,31,46,61,76,91),
                            c(17,32,47,62,77,92),
                            c(18,33,48,63,78,93),
                            c(19,34,49,64,79,94),
                            c(20,35,50,65,80,95),
                            c(21,36,51,66,81,96),
                            c(22,37,52,67,82,97),
                            c(23,38,53,68,83,98),
                            c(24,39,54,69,84,99),
                            c(25,40,55,70,85,100),
                            c(26,41,56,71,86,101),
                            c(27,42,57,72,87,102),
                            c(28,43,58,73,88,103),
                            c(29,44,59,74,89,104)),
               idvar="Participant_ID",
               timevar="time",
               sep="",
               direction="long")
  formula <- paste0(paste0("BPI_PScore ~ acum + ",Lvars)," + (1|Participant_ID)")
  lmer(data=x,
       formula=formula)
})

rialc_12mco1a <- matrix(unlist(parLapply(cl=cl, rialc_12m1a, function (x) {coef(summary(x))[2]})),nrow=50,ncol=1)
rialc_12mse1a <- matrix(unlist(parLapply(cl=cl, rialc_12m1a, function (x) {vcov(x)[2,2]})),nrow=50,ncol=1)
rialc_12mres1a <- matrix(unlist(mi.meld(q=rialc_12mco1a, se=rialc_12mse1a)),nrow=1,ncol=2)

rialc_12m2a <- parLapply(cl=cl, impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6))
  x$acum2 <- x$freq12
  x$acum3 <- x$acum2 + x$freq13
  x$acum4 <- x$acum3 + x$freq14
  x$acum5 <- x$acum4 + x$freq15
  x$acum6 <- x$acum5 + x$freq16
  x <- x[c(1:27,99,28:41,100,42:55,101,56:69,102,70:83,103,84:97,104,98)]
  x <- reshape(x,
               v.names=varlist_2,
               varying=list(c(15,30,45,60,75,90),
                            c(16,31,46,61,76,91),
                            c(17,32,47,62,77,92),
                            c(18,33,48,63,78,93),
                            c(19,34,49,64,79,94),
                            c(20,35,50,65,80,95),
                            c(21,36,51,66,81,96),
                            c(22,37,52,67,82,97),
                            c(23,38,53,68,83,98),
                            c(24,39,54,69,84,99),
                            c(25,40,55,70,85,100),
                            c(26,41,56,71,86,101),
                            c(27,42,57,72,87,102),
                            c(28,43,58,73,88,103),
                            c(29,44,59,74,89,104)),
               idvar="Participant_ID",
               timevar="time",
               sep="",
               direction="long")
  lmer(data=x,
       formula=paste0(paste0("BPI_interference ~ acum + ",Lvars)," + (1|Participant_ID)"))
})
rialc_12mco2a <- matrix(unlist(parLapply(cl=cl, rialc_12m2a, function (x) {coef(summary(x))[2]})),nrow=50,ncol=1)
rialc_12mse2a <- matrix(unlist(parLapply(cl=cl, rialc_12m2a, function (x) {vcov(x)[2,2]})),nrow=50,ncol=1)
rialc_12mres2a <- matrix(unlist(mi.meld(q=rialc_12mco2a, se=rialc_12mse2a)),nrow=1,ncol=2)

rm(list=c("rialc_12m1a","rialc_12mco1a","rialc_12mse1a",
          "rialc_12m2a","rialc_12mco2a","rialc_12mse2a"))

## 4+ DRINKS PER WEEK ##
rialc_12m1b <- parLapply(cl=cl, impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6))
  x$acum2 <- x$freq22
  x$acum3 <- x$acum2 + x$freq23
  x$acum4 <- x$acum3 + x$freq24
  x$acum5 <- x$acum4 + x$freq25
  x$acum6 <- x$acum5 + x$freq26
  x <- x[c(1:27,99,28:41,100,42:55,101,56:69,102,70:83,103,84:97,104,98)]
  x <- reshape(x,
               v.names=varlist_1,
               varying=list(c(15,30,45,60,75,90),
                            c(16,31,46,61,76,91),
                            c(17,32,47,62,77,92),
                            c(18,33,48,63,78,93),
                            c(19,34,49,64,79,94),
                            c(20,35,50,65,80,95),
                            c(21,36,51,66,81,96),
                            c(22,37,52,67,82,97),
                            c(23,38,53,68,83,98),
                            c(24,39,54,69,84,99),
                            c(25,40,55,70,85,100),
                            c(26,41,56,71,86,101),
                            c(27,42,57,72,87,102),
                            c(28,43,58,73,88,103),
                            c(29,44,59,74,89,104)),
               idvar="Participant_ID",
               timevar="time",
               sep="",
               direction="long")
  formula <- paste0(paste0("BPI_PScore ~ acum + ",Lvars)," + (1|Participant_ID)")
  lmer(data=x,
       formula=formula)
})

rialc_12mco1b <- matrix(unlist(parLapply(cl=cl, rialc_12m1b, function (x) {coef(summary(x))[2]})),nrow=50,ncol=1)
rialc_12mse1b <- matrix(unlist(parLapply(cl=cl, rialc_12m1b, function (x) {vcov(x)[2,2]})),nrow=50,ncol=1)
rialc_12mres1b <- matrix(unlist(mi.meld(q=rialc_12mco1b, se=rialc_12mse1b)),nrow=1,ncol=2)

rialc_12m2b <- parLapply(cl=cl, impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6))
  x$acum2 <- x$freq22
  x$acum3 <- x$acum2 + x$freq23
  x$acum4 <- x$acum3 + x$freq24
  x$acum5 <- x$acum4 + x$freq25
  x$acum6 <- x$acum5 + x$freq26
  x <- x[c(1:27,99,28:41,100,42:55,101,56:69,102,70:83,103,84:97,104,98)]
  x <- reshape(x,
               v.names=varlist_2,
               varying=list(c(15,30,45,60,75,90),
                            c(16,31,46,61,76,91),
                            c(17,32,47,62,77,92),
                            c(18,33,48,63,78,93),
                            c(19,34,49,64,79,94),
                            c(20,35,50,65,80,95),
                            c(21,36,51,66,81,96),
                            c(22,37,52,67,82,97),
                            c(23,38,53,68,83,98),
                            c(24,39,54,69,84,99),
                            c(25,40,55,70,85,100),
                            c(26,41,56,71,86,101),
                            c(27,42,57,72,87,102),
                            c(28,43,58,73,88,103),
                            c(29,44,59,74,89,104)),
               idvar="Participant_ID",
               timevar="time",
               sep="",
               direction="long")
  lmer(data=x,
       formula=paste0(paste0("BPI_interference ~ acum + ",Lvars)," + (1|Participant_ID)"))
})
rialc_12mco2b <- matrix(unlist(parLapply(cl=cl, rialc_12m2b, function (x) {coef(summary(x))[2]})),nrow=50,ncol=1)
rialc_12mse2b <- matrix(unlist(parLapply(cl=cl, rialc_12m2b, function (x) {vcov(x)[2,2]})),nrow=50,ncol=1)
rialc_12mres2b <- matrix(unlist(mi.meld(q=rialc_12mco2b, se=rialc_12mse2b)),nrow=1,ncol=2)

rm(list=c("rialc_12m1b","rialc_12mco1b","rialc_12mse1b",
          "rialc_12m2b","rialc_12mco2b","rialc_12mse2b"))

## BINGE DRINKING ANALYSIS ##
ribinge1 <- parLapply(cl=cl, impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6))
  x$acum1 <- 0
  x$acum2 <- x$binge2
  x$acum3 <- x$acum2 + x$binge3
  x$acum4 <- x$acum3 + x$binge4
  x$acum5 <- x$acum4 + x$binge5
  x$acum6 <- x$acum5 + x$binge6
  x <- x[c(1:27,99,28:41,100,42:55,101,56:69,102,70:83,103,84:97,104,98)]
  x <- reshape(x,
               v.names=varlist_3,
               varying=list(c(15,30,45,60,75,90),
                            c(16,31,46,61,76,91),
                            c(17,32,47,62,77,92),
                            c(18,33,48,63,78,93),
                            c(19,34,49,64,79,94),
                            c(20,35,50,65,80,95),
                            c(21,36,51,66,81,96),
                            c(22,37,52,67,82,97),
                            c(23,38,53,68,83,98),
                            c(24,39,54,69,84,99),
                            c(25,40,55,70,85,100),
                            c(26,41,56,71,86,101),
                            c(27,42,57,72,87,102),
                            c(28,43,58,73,88,103),
                            c(29,44,59,74,89,104)),
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
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6))
  x$acum1 <- 0
  x$acum2 <- x$binge2
  x$acum3 <- x$acum2 + x$binge3
  x$acum4 <- x$acum3 + x$binge4
  x$acum5 <- x$acum4 + x$binge5
  x$acum6 <- x$acum5 + x$binge6
  x <- x[c(1:27,99,28:41,100,42:55,101,56:69,102,70:83,103,84:97,104,98)]
  x <- reshape(x,
               v.names=varlist_4,
               varying=list(c(15,30,45,60,75,90),
                            c(16,31,46,61,76,91),
                            c(17,32,47,62,77,92),
                            c(18,33,48,63,78,93),
                            c(19,34,49,64,79,94),
                            c(20,35,50,65,80,95),
                            c(21,36,51,66,81,96),
                            c(22,37,52,67,82,97),
                            c(23,38,53,68,83,98),
                            c(24,39,54,69,84,99),
                            c(25,40,55,70,85,100),
                            c(26,41,56,71,86,101),
                            c(27,42,57,72,87,102),
                            c(28,43,58,73,88,103),
                            c(29,44,59,74,89,104)),
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
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6))
  x$acum1 <- 0
  x$acum2 <- x$auditcprob2
  x$acum3 <- x$acum2 + x$auditcprob3
  x$acum4 <- x$acum3 + x$auditcprob4
  x$acum5 <- x$acum4 + x$auditcprob5
  x$acum6 <- x$acum5 + x$auditcprob6
  x <- x[c(1:27,99,28:41,100,42:55,101,56:69,102,70:83,103,84:97,104,98)]
  x <- reshape(x,
               v.names=varlist_5,
               varying=list(c(15,30,45,60,75,90),
                            c(16,31,46,61,76,91),
                            c(17,32,47,62,77,92),
                            c(18,33,48,63,78,93),
                            c(19,34,49,64,79,94),
                            c(20,35,50,65,80,95),
                            c(21,36,51,66,81,96),
                            c(22,37,52,67,82,97),
                            c(23,38,53,68,83,98),
                            c(24,39,54,69,84,99),
                            c(25,40,55,70,85,100),
                            c(26,41,56,71,86,101),
                            c(27,42,57,72,87,102),
                            c(28,43,58,73,88,103),
                            c(29,44,59,74,89,104)),
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
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6))
  x$acum1 <- 0
  x$acum2 <- x$auditcprob2
  x$acum3 <- x$acum2 + x$auditcprob3
  x$acum4 <- x$acum3 + x$auditcprob4
  x$acum5 <- x$acum4 + x$auditcprob5
  x$acum6 <- x$acum5 + x$auditcprob6
  x <- x[c(1:27,99,28:41,100,42:55,101,56:69,102,70:83,103,84:97,104,98)]
  x <- reshape(x,
               v.names=varlist_6,
               varying=list(c(15,30,45,60,75,90),
                            c(16,31,46,61,76,91),
                            c(17,32,47,62,77,92),
                            c(18,33,48,63,78,93),
                            c(19,34,49,64,79,94),
                            c(20,35,50,65,80,95),
                            c(21,36,51,66,81,96),
                            c(22,37,52,67,82,97),
                            c(23,38,53,68,83,98),
                            c(24,39,54,69,84,99),
                            c(25,40,55,70,85,100),
                            c(26,41,56,71,86,101),
                            c(27,42,57,72,87,102),
                            c(28,43,58,73,88,103),
                            c(29,44,59,74,89,104)),
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

riresults <- matrix(c(rialc_12mres1a,rialc_12mres2a,rialc_12mres1b,rialc_12mres2b,ribingeres1,ribingeres2,riauditres1,riauditres2),byrow=TRUE,ncol=4,nrow=4)
rownames(riresults) <- c("2+ drinks per week","4+ drinks per week","Binge drinking","AUDIT-C Problematic drinking")
colnames(riresults) <- c("BPI Pain Score log(OR)","BPI Pain Score SE","BPI Pain Int log(OR)","BPI Pain Int SE")

save(riresults,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/riresults.RData"))
