
## SYNTAX FILE 4 - LTMLE ANALYSIS EXCLUDING NON-DRINKERS     ##

cloudstor <- "C:/Users/z3312911/Cloudstor/" # change to master file path
.libPaths(paste0(cloudstor,"R Library"))

libs <- c("parallel","SuperLearner","ltmle","Amelia","ranger","gam")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}

library ("parallel")
library ("Amelia")

set.seed(857915,kind="L'Ecuyer-CMRG")

# Load imputed data
load(file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Data/Imputed data - wide.RData"))

# Set up parallel processing to reduce analysis time
numcores <- future::availableCores()
cl <- makeCluster(numcores)
parallel::clusterEvalQ(cl, cloudstor <- "C:/Users/z3312911/Cloudstor/")
parallel::clusterEvalQ(cl, .libPaths(paste0(cloudstor,"R Library")))
parallel::clusterEvalQ(cl, library("SuperLearner"))
parallel::clusterEvalQ(cl, library("ltmle"))
parallel::clusterEvalQ(cl, library("ranger"))
parallel::clusterEvalQ(cl, library("gam"))
parallel::clusterEvalQ(cl, Lvars1 <- c("PHQ9_Mod_sev1","GADMod2Sev1","Antidepressant_week1","Antipsychotic_week1","benzo_week1","Nonopioid_analgesic_week1","Pregabalin_week1","can_12m1","cig_12m1","opioid901","alc_pain_12m1",
                                       "PHQ9_Mod_sev2","GADMod2Sev2","Antidepressant_week2","Antipsychotic_week2","benzo_week2","Nonopioid_analgesic_week2","Pregabalin_week2","can_12m2","cig_12m2","opioid902","alc_pain_12m2",
                                       "PHQ9_Mod_sev3","GADMod2Sev3","Antidepressant_week3","Antipsychotic_week3","benzo_week3","Nonopioid_analgesic_week3","Pregabalin_week3","can_12m3","cig_12m3","opioid903","alc_pain_12m3",
                                       "PHQ9_Mod_sev4","GADMod2Sev4","Antidepressant_week4","Antipsychotic_week4","benzo_week4","Nonopioid_analgesic_week4","Pregabalin_week4","can_12m4","cig_12m4","opioid904","alc_pain_12m4",
                                       "PHQ9_Mod_sev5","GADMod2Sev5","Antidepressant_week5","Antipsychotic_week5","benzo_week5","Nonopioid_analgesic_week5","Pregabalin_week5","can_12m5","cig_12m5","opioid905","alc_pain_12m5",
                                       "PHQ9_Mod_sev6","GADMod2Sev6","Antidepressant_week6","Antipsychotic_week6","benzo_week6","Nonopioid_analgesic_week6","Pregabalin_week6","can_12m6","cig_12m6","opioid906","alc_pain_12m6"))
parallel::clusterEvalQ(cl, Lvars2 <- c("PHQ9_Mod_sev2","GADMod2Sev2","Antidepressant_week2","Antipsychotic_week2","benzo_week2","Nonopioid_analgesic_week2","Pregabalin_week2","can_12m2","cig_12m2","opioid902","alc_pain_12m2",
                                       "PHQ9_Mod_sev3","GADMod2Sev3","Antidepressant_week3","Antipsychotic_week3","benzo_week3","Nonopioid_analgesic_week3","Pregabalin_week3","can_12m3","cig_12m3","opioid903","alc_pain_12m3",
                                       "PHQ9_Mod_sev4","GADMod2Sev4","Antidepressant_week4","Antipsychotic_week4","benzo_week4","Nonopioid_analgesic_week4","Pregabalin_week4","can_12m4","cig_12m4","opioid904","alc_pain_12m4",
                                       "PHQ9_Mod_sev5","GADMod2Sev5","Antidepressant_week5","Antipsychotic_week5","benzo_week5","Nonopioid_analgesic_week5","Pregabalin_week5","can_12m5","cig_12m5","opioid905","alc_pain_12m5",
                                       "PHQ9_Mod_sev6","GADMod2Sev6","Antidepressant_week6","Antipsychotic_week6","benzo_week6","Nonopioid_analgesic_week6","Pregabalin_week6","can_12m6","cig_12m6","opioid906","alc_pain_12m6"))
parallel::clusterEvalQ(cl, regimeList1 <- list(function(row) c(1,1,1,1,1,1),
                                               function(row) c(0,1,1,1,1,1),
                                               function(row) c(0,0,1,1,1,1),
                                               function(row) c(0,0,0,1,1,1),
                                               function(row) c(0,0,0,0,1,1),
                                               function(row) c(0,0,0,0,0,1),
                                               function(row) c(0,0,0,0,0,0)))
parallel::clusterEvalQ(cl, sum.measures1 <- as.array(matrix(c(6,5,4,3,2,1,0),nrow=7,ncol=1)))
parallel::clusterEvalQ(cl, colnames(sum.measures1) <- "init")
parallel::clusterEvalQ(cl, regimeList2 <- list(function(row) c(1,1,1,1,1),
                                               function(row) c(0,1,1,1,1),
                                               function(row) c(0,0,1,1,1),
                                               function(row) c(0,0,0,1,1),
                                               function(row) c(0,0,0,0,1),
                                               function(row) c(0,0,0,0,0)))
parallel::clusterEvalQ(cl, sum.measures2 <- as.array(matrix(c(5,4,3,2,1,0),nrow=6,ncol=1)))
parallel::clusterEvalQ(cl, colnames(sum.measures2) <- "init")
parallel::clusterEvalQ(cl, create.Learner("SL.ranger", params = list(num.trees = 250)))
parallel::clusterEvalQ(cl, create.Learner("SL.randomForest", params = list(num.trees = 250)))
parallel::clusterEvalQ(cl, SLlib <- list(Q=c("SL.mean","SL.glm","SL.gam"),
                                         g=c("SL.mean","SL.glm","SL.gam","SL.ranger_1")))
parallel::clusterEvalQ(cl, flist <- c("Participant_ID","time","b_Actual_age","b_employ3","b_sex","b_edu","b_MaritalStatus","b_Arth_12mR","b_Back_12mR","b_Head_12mR","b_Visc_12mR","b_Fibro_12mR","b_Cmplx_12mR","b_Shing_12mR","b_pain_duration_yrs3","PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","alc_12m","binge","auditcprob","BPI_PScore","BPI_interference"))
parallel::clusterEvalQ(cl, varlist <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","alc_12m","binge","auditcprob","BPI_PScore","BPI_interference"))

## BINGE DRINKING ANALYSIS ##
databinge1_sub1 <- parLapply(cl,impdatawide, function (x) {
  x <- x[,-15]
  x <- x[,c(-28,-30,-44,-46,-60,-62,-76,-78,-92,-94,-108,-110)]
  x$acum <- x$alc_12m1 + x$alc_12m2 + x$alc_12m3 + x$alc_12m4 + x$alc_12m5 + x$alc_12m6
  x <- x[!x$acum==0,]
  x <- x[,c(-26,-40,-54,-68,-82,-96,-99,-100,-101,-102,-103,-104)]
  x
})
databinge2_sub1 <- parLapply(cl,impdatawide, function (x) {
  x <- x[,-15]
  x <- x[,c(-28,-29,-44,-45,-60,-61,-76,-77,-92,-93,-108,-109)]
  x$acum <- x$alc_12m1 + x$alc_12m2 + x$alc_12m3 + x$alc_12m4 + x$alc_12m5 + x$alc_12m6
  x <- x[!x$acum==0,]
  x <- x[,c(-26,-40,-54,-68,-82,-96,-99,-100,-101,-102,-103,-104)]
  x
})

jtbinge1 <- parLapply(cl=cl, databinge1_sub1, function (x) {ltmle(x,
                                                             Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeres1 <- parLapply(cl=cl, jtbinge1, function (x) {summary(x)})
jtbingeco1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$ATE$estimate})),nrow=50,ncol=1)
jtbingese1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$ATE$std.dev})),nrow=50,ncol=1)
jtbingeres1 <- matrix(unlist(mi.meld(q=jtbingeco1, se=jtbingese1)),nrow=1,ncol=2)

jtbinge2 <- parLapply(cl=cl, databinge2_sub1, function (x) {ltmle(x,
                                                             Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeres2 <- parLapply(cl=cl, jtbinge2, function (x) {summary(x)})
jtbingeco2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$ATE$estimate})),nrow=50,ncol=1)
jtbingese2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$ATE$std.dev})),nrow=50,ncol=1)
jtbingeres2 <- matrix(unlist(mi.meld(q=jtbingeco2, se=jtbingese2)),nrow=1,ncol=2)

rm(list=c("jtbinge1","jtbingeco1","jtbingese1","databinge1_sub1",
          "jtbinge2","jtbingeco2","jtbingese2","databinge2_sub1"))

## AUDIT PROBLEMATIC DRINKING ANALYSIS ##
dataauditprob1_sub1 <- parLapply(cl,impdatawide, function (x) {
  x <- x[,-15]
  x <- x[,c(-27,-30,-43,-46,-59,-62,-75,-78,-91,-94,-107,-110)]
  x$acum <- x$alc_12m1 + x$alc_12m2 + x$alc_12m3 + x$alc_12m4 + x$alc_12m5 + x$alc_12m6
  x <- x[!x$acum==0,]
  x <- x[,c(-26,-40,-54,-68,-82,-96,-99,-100,-101,-102,-103,-104)]
  x
})
dataauditprob2_sub1 <- parLapply(cl,impdatawide, function (x) {
  x <- x[,-15]
  x <- x[,c(-27,-29,-43,-45,-59,-61,-75,-77,-91,-93,-107,-109)]
  x$acum <- x$alc_12m1 + x$alc_12m2 + x$alc_12m3 + x$alc_12m4 + x$alc_12m5 + x$alc_12m6
  x <- x[!x$acum==0,]
  x <- x[,c(-26,-40,-54,-68,-82,-96,-99,-100,-101,-102,-103,-104)]
  x
})

jtaudit1 <- parLapply(cl=cl, dataauditprob1_sub1, function (x) {ltmle(x,
                                                                 Anodes=c("auditcprob2","auditcprob3","auditcprob4","auditcprob5","auditcprob6"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtauditres1 <- parLapply(cl=cl, jtaudit1, function (x) {summary(x)})
jtauditco1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$ATE$estimate})),nrow=50,ncol=1)
jtauditse1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$ATE$std.dev})),nrow=50,ncol=1)
jtauditres1 <- matrix(unlist(mi.meld(q=jtauditco1, se=jtauditse1)),nrow=1,ncol=2)

jtaudit2 <- parLapply(cl=cl, dataauditprob2_sub1, function (x) {ltmle(x,
                                                                 Anodes=c("auditcprob2","auditcprob3","auditcprob4","auditcprob5","auditcprob6"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtauditres2 <- parLapply(cl=cl, jtaudit2, function (x) {summary(x)})
jtauditco2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$ATE$estimate})),nrow=50,ncol=1)
jtauditse2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$ATE$std.dev})),nrow=50,ncol=1)
jtauditres2 <- matrix(unlist(mi.meld(q=jtauditco2, se=jtauditse2)),nrow=1,ncol=2)

rm(list=c("jtaudit1","jtauditco1","jtauditse1","dataauditprob1_sub1",
          "jtaudit2","jtauditco2","jtauditse2","dataauditprob2_sub1"))

## COMBINE RESULTS INTO MATRIX FOR EXCEL ##
jtresults_sub1 <- matrix(c(jtbingeres1,jtbingeres2,jtauditres1,jtauditres2),byrow=TRUE,ncol=4,nrow=2)
rownames(jtresults_sub1) <- c("Binge drinking","AUDIT-C Problematic drinking")
colnames(jtresults_sub1) <- c("BPI Pain Score log(OR)","BPI Pain Score SE","BPI Pain Int log(OR)","BPI Pain Int SE")

save(jtresults_sub1,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/jtresults_sub1.RData"))

