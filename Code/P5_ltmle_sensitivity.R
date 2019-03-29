
## SYNTAX FILE 5 - LTMLE ANALYSIS EXCLUDING EX-DRINKERS      ##

cloudstor <- "C:/Users/z3312911/Cloudstor/" # change to master file path
.libPaths(paste0(cloudstor,"R Library"))

libs <- c("parallel","SuperLearner","ltmle","Amelia","ranger","gam","MBESS")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}

library ("parallel")
library ("Amelia")

set.seed(725276,kind="L'Ecuyer-CMRG")

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
parallel::clusterEvalQ(cl, SLlib <- list(Q=c("SL.mean","SL.glm","SL.gam"),
                                         g=c("SL.mean","SL.glm","SL.gam","SL.ranger_1")))

## ANY ALCOHOL CONSUMPTION ANALYSIS ##
dataalcuse1_sub2 <- lapply(impdatawide, function (x) {
  x$alccum <- x$alc_12m1 + x$alc_12m2 + x$alc_12m3 + x$alc_12m4 + x$alc_12m5 + x$alc_12m6
  x <- x[,c(-28,-29,-31,-44,-45,-47,-60,-61,-63,-76,-77,-79,-92,-93,-95,-108,-109,-111)]
  x <- x[!(x$alccum==0 & x$b_alc_ever==1),]
  x <- x[,c(-15,-94)]
  x
})
dataalcuse2_sub2 <- lapply(impdatawide, function (x) {
  x$alccum <- x$alc_12m1 + x$alc_12m2 + x$alc_12m3 + x$alc_12m4 + x$alc_12m5 + x$alc_12m6
  x <- x[,c(-28,-29,-30,-44,-45,-46,-60,-61,-62,-76,-77,-78,-92,-93,-94,-108,-109,-110)]
  x <- x[!(x$alccum==0 & x$b_alc_ever==1),]
  x <- x[,c(-15,-94)]
  x
})

jtalcuse1 <- parLapply(cl=cl, dataalcuse1_sub2, function (x) {ltmle(x,
                                                               Anodes=c("alc_12m1","alc_12m2","alc_12m3","alc_12m4","alc_12m5","alc_12m6"),
                                                               Lnodes=Lvars1,
                                                               Ynodes=c("BPI_PScore1","BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                               survivalOutcome=FALSE,
                                                               SL.library=SLlib,
                                                               abar=list(c(1,1,1,1,1,1),c(0,0,0,0,0,0)))})
jtalcuseres1 <- parLapply(cl=cl, jtalcuse1, function (x) {summary(x)})
jtalcuseco1 <- matrix(unlist(parLapply(cl=cl, jtalcuseres1, function (x) {x$effect.measures$ATE$estimate})),nrow=50,ncol=1)
jtalcusese1 <- matrix(unlist(parLapply(cl=cl, jtalcuseres1, function (x) {x$effect.measures$ATE$std.dev})),nrow=50,ncol=1)
jtalcusetco1 <- matrix(unlist(parLapply(cl=cl, jtalcuseres1, function (x) {x$effect.measures$treatment$estimate})),nrow=50,ncol=1)
jtalcusetse1 <- matrix(unlist(parLapply(cl=cl, jtalcuseres1, function (x) {x$effect.measures$treatment$std.dev})),nrow=50,ncol=1)
jtalcusecco1 <- matrix(unlist(parLapply(cl=cl, jtalcuseres1, function (x) {x$effect.measures$control$estimate})),nrow=50,ncol=1)
jtalcusecse1 <- matrix(unlist(parLapply(cl=cl, jtalcuseres1, function (x) {x$effect.measures$control$std.dev})),nrow=50,ncol=1)
jtalcuseres1 <- matrix(c(unlist(mi.meld(q=jtalcuseco1,se=jtalcusese1)),
                         unlist(mi.meld(q=jtalcusetco1,se=jtalcusese1)),
                         unlist(mi.meld(q=jtalcusecco1,se=jtalcusecse1))),nrow=1,ncol=6)


jtalcuse2 <- parLapply(cl=cl, dataalcuse2_sub2, function (x) {ltmle(x,
                                                               Anodes=c("alc_12m1","alc_12m2","alc_12m3","alc_12m4","alc_12m5","alc_12m6"),
                                                               Lnodes=Lvars1,
                                                               Ynodes=c("BPI_interference1","BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                               survivalOutcome=FALSE,
                                                               SL.library=SLlib,
                                                               abar=list(c(1,1,1,1,1,1),c(0,0,0,0,0,0)))})
jtalcuseres2 <- parLapply(cl=cl, jtalcuse2, function (x) {summary(x)})
jtalcuseco2 <- matrix(unlist(parLapply(cl=cl, jtalcuseres2, function (x) {x$effect.measures$ATE$estimate})),nrow=50,ncol=1)
jtalcusese2 <- matrix(unlist(parLapply(cl=cl, jtalcuseres2, function (x) {x$effect.measures$ATE$std.dev})),nrow=50,ncol=1)
jtalcusetco2 <- matrix(unlist(parLapply(cl=cl, jtalcuseres2, function (x) {x$effect.measures$treatment$estimate})),nrow=50,ncol=1)
jtalcusetse2 <- matrix(unlist(parLapply(cl=cl, jtalcuseres2, function (x) {x$effect.measures$treatment$std.dev})),nrow=50,ncol=1)
jtalcusecco2 <- matrix(unlist(parLapply(cl=cl, jtalcuseres2, function (x) {x$effect.measures$control$estimate})),nrow=50,ncol=1)
jtalcusecse2 <- matrix(unlist(parLapply(cl=cl, jtalcuseres2, function (x) {x$effect.measures$control$std.dev})),nrow=50,ncol=1)
jtalcuseres2 <- matrix(c(unlist(mi.meld(q=jtalcuseco2,se=jtalcusese2)),
                         unlist(mi.meld(q=jtalcusetco2,se=jtalcusese2)),
                         unlist(mi.meld(q=jtalcusecco2,se=jtalcusecse2))),nrow=1,ncol=6)

rm(list=c("jtalcuse1","jtalcuseco1","jtalcusese1","jtalcusetco1","jtalcusetse1","jtalcusecco1","jtalcusecse1","dataalcuse1_sub2",
          "jtalcuse2","jtalcuseco2","jtalcusese2","jtalcusetco2","jtalcusetse2","jtalcusecco2","jtalcusecse2","dataalcuse2_sub2"))

## BINGE DRINKING ANALYSIS ##
databinge1_sub2 <- lapply(impdatawide, function (x) {
  x$alccum <- x$alc_12m1 + x$alc_12m2 + x$alc_12m3 + x$alc_12m4 + x$alc_12m5 + x$alc_12m6
  x <- x[,c(-27,-29,-31,-43,-45,-47,-59,-61,-63,-75,-77,-79,-91,-93,-95,-107,-109,-111)]
  x <- x[!(x$alccum==0 & x$b_alc_ever==1),]
  x <- x[,c(-15,-94)]
  x
})
databinge2_sub2 <- lapply(impdatawide, function (x) {
  x$alccum <- x$alc_12m1 + x$alc_12m2 + x$alc_12m3 + x$alc_12m4 + x$alc_12m5 + x$alc_12m6
  x <- x[,c(-27,-29,-30,-43,-45,-46,-59,-61,-62,-75,-77,-78,-91,-93,-94,-107,-109,-110)]
  x <- x[!(x$alccum==0 & x$b_alc_ever==1),]
  x <- x[,c(-15,-94)]
  x
})

jtbinge1 <- parLapply(cl=cl, databinge1_sub2, function (x) {ltmle(x,
                                                             Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeres1 <- parLapply(cl=cl, jtbinge1, function (x) {summary(x)})
jtbingeco1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$ATE$estimate})),nrow=50,ncol=1)
jtbingese1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$ATE$std.dev})),nrow=50,ncol=1)
jtbingetco1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$treatment$estimate})),nrow=50,ncol=1)
jtbingetse1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$treatment$std.dev})),nrow=50,ncol=1)
jtbingecco1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$control$estimate})),nrow=50,ncol=1)
jtbingecse1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$control$std.dev})),nrow=50,ncol=1)
jtbingeres1 <- matrix(c(unlist(mi.meld(q=jtbingeco1,se=jtbingese1)),
                        unlist(mi.meld(q=jtbingetco1,se=jtbingese1)),
                        unlist(mi.meld(q=jtbingecco1,se=jtbingecse1))),nrow=1,ncol=6)

jtbinge2 <- parLapply(cl=cl, databinge2_sub2, function (x) {ltmle(x,
                                                             Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeres2 <- parLapply(cl=cl, jtbinge2, function (x) {summary(x)})
jtbingeco2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$ATE$estimate})),nrow=50,ncol=1)
jtbingese2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$ATE$std.dev})),nrow=50,ncol=1)
jtbingetco2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$treatment$estimate})),nrow=50,ncol=1)
jtbingetse2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$treatment$std.dev})),nrow=50,ncol=1)
jtbingecco2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$control$estimate})),nrow=50,ncol=1)
jtbingecse2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$control$std.dev})),nrow=50,ncol=1)
jtbingeres2 <- matrix(c(unlist(mi.meld(q=jtbingeco2,se=jtbingese2)),
                        unlist(mi.meld(q=jtbingetco2,se=jtbingese2)),
                        unlist(mi.meld(q=jtbingecco2,se=jtbingecse2))),nrow=1,ncol=6)

rm(list=c("jtbinge1","jtbingeco1","jtbingese1","jtbingetco1","jtbingetse1","jtbingecco1","jtbingecse1","databinge1_sub2",
          "jtbinge2","jtbingeco2","jtbingese2","jtbingetco2","jtbingetse2","jtbingecco2","jtbingecse2","databinge2_sub2"))

## AUDIT PROBLEMATIC DRINKING ANALYSIS ##
dataauditprob1_sub2 <- lapply(impdatawide, function (x) {
  x$alccum <- x$alc_12m1 + x$alc_12m2 + x$alc_12m3 + x$alc_12m4 + x$alc_12m5 + x$alc_12m6
  x <- x[,c(-27,-28,-31,-43,-44,-47,-59,-60,-63,-75,-76,-79,-91,-92,-95,-107,-108,-111)]
  x <- x[!(x$alccum==0 & x$b_alc_ever==1),]
  x <- x[,c(-15,-94)]
  x
})
dataauditprob2_sub2 <- lapply(impdatawide, function (x) {
  x$alccum <- x$alc_12m1 + x$alc_12m2 + x$alc_12m3 + x$alc_12m4 + x$alc_12m5 + x$alc_12m6
  x <- x[,c(-27,-28,-30,-43,-44,-46,-59,-60,-62,-75,-76,-78,-91,-92,-94,-107,-108,-110)]
  x <- x[!(x$alccum==0 & x$b_alc_ever==1),]
  x <- x[,c(-15,-94)]
  x
})

jtaudit1 <- parLapply(cl=cl, dataauditprob1_sub2, function (x) {ltmle(x,
                                                                 Anodes=c("auditcprob2","auditcprob3","auditcprob4","auditcprob5","auditcprob6"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtauditres1 <- parLapply(cl=cl, jtaudit1, function (x) {summary(x)})
jtauditco1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$ATE$estimate})),nrow=50,ncol=1)
jtauditse1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$ATE$std.dev})),nrow=50,ncol=1)
jtaudittco1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$treatment$estimate})),nrow=50,ncol=1)
jtaudittse1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$treatment$std.dev})),nrow=50,ncol=1)
jtauditcco1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$control$estimate})),nrow=50,ncol=1)
jtauditcse1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$control$std.dev})),nrow=50,ncol=1)
jtauditres1 <- matrix(c(unlist(mi.meld(q=jtauditco1,se=jtauditse1)),
                        unlist(mi.meld(q=jtaudittco1,se=jtauditse1)),
                        unlist(mi.meld(q=jtauditcco1,se=jtauditcse1))),nrow=1,ncol=6)

jtaudit2 <- parLapply(cl=cl, dataauditprob2_sub2, function (x) {ltmle(x,
                                                                 Anodes=c("auditcprob2","auditcprob3","auditcprob4","auditcprob5","auditcprob6"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtauditres2 <- parLapply(cl=cl, jtaudit2, function (x) {summary(x)})
jtauditco2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$ATE$estimate})),nrow=50,ncol=1)
jtauditse2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$ATE$std.dev})),nrow=50,ncol=1)
jtaudittco2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$treatment$estimate})),nrow=50,ncol=1)
jtaudittse2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$treatment$std.dev})),nrow=50,ncol=1)
jtauditcco2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$control$estimate})),nrow=50,ncol=1)
jtauditcse2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$control$std.dev})),nrow=50,ncol=1)
jtauditres2 <- matrix(c(unlist(mi.meld(q=jtauditco2,se=jtauditse2)),
                        unlist(mi.meld(q=jtaudittco2,se=jtauditse2)),
                        unlist(mi.meld(q=jtauditcco2,se=jtauditcse2))),nrow=1,ncol=6)

rm(list=c("jtaudit1","jtauditco1","jtauditse1","jtaudittco1","jtaudittse1","jtauditcco1","jtauditcse1","dataauditprob1_sub2",
          "jtaudit2","jtauditco2","jtauditse2","jtaudittco2","jtaudittse2","jtauditcco2","jtauditcse2","dataauditprob2_sub2"))

## COMBINE RESULTS INTO MATRIX FOR EXCEL ##
jtresults_sub2 <- matrix(c(jtalcuseres1[1:2],jtalcuseres2[1:2],jtbingeres1[1:2],jtbingeres2[1:2],jtauditres1[1:2],jtauditres2[1:2]),byrow=TRUE,ncol=4,nrow=3)
rownames(jtresults_sub2) <- c("Any alcohol consumption","Binge drinking","AUDIT-C Problematic drinking")
colnames(jtresults_sub2) <- c("BPI Pain Score coef","BPI Pain Score SE","BPI Pain Int coef","BPI Pain Int SE")

save(jtresults_sub2,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/jtresults_sub2.RData"))



