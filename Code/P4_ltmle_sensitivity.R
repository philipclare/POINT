
## SYNTAX FILE 4 - SENSITIVITY ANALYSIS                      ##

cloudstor <- "C:/Users/z3312911/Cloudstor/" # change to master file path
.libPaths(paste0(cloudstor,"R Library"))

libs <- c("parallel","SuperLearner","ltmle","Amelia","ranger","gam","MBESS")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}

library ("parallel")
library ("Amelia")
library ("MBESS")

set.seed(599648,kind="L'Ecuyer-CMRG")

# Load imputed data
load(file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Data/Imputed data - wide.RData"))

# Obtain number of imputations from length of list for use in results matrices
nimpute <- length(impdatawide)

# Set up parallel processing to reduce analysis time
numcores <- future::availableCores()-1
cl <- makeCluster(numcores)
parallel::clusterEvalQ(cl, cloudstor <- "C:/Users/z3312911/Cloudstor/")
parallel::clusterEvalQ(cl, .libPaths(paste0(cloudstor,"R Library")))
parallel::clusterEvalQ(cl, library("SuperLearner"))
parallel::clusterEvalQ(cl, library("ltmle"))
parallel::clusterEvalQ(cl, library("ranger"))
parallel::clusterEvalQ(cl, library("gam"))
parallel::clusterEvalQ(cl, Lvars2 <- c("PHQ9_Mod_sev2","GADMod2Sev2","Antidepressant_week2","Antipsychotic_week2","benzo_week2","Nonopioid_analgesic_week2","Pregabalin_week2","can_12m2","cig_12m2","opioid902","PSEQ_Score2",
                                       "PHQ9_Mod_sev3","GADMod2Sev3","Antidepressant_week3","Antipsychotic_week3","benzo_week3","Nonopioid_analgesic_week3","Pregabalin_week3","can_12m3","cig_12m3","opioid903","PSEQ_Score3",
                                       "PHQ9_Mod_sev4","GADMod2Sev4","Antidepressant_week4","Antipsychotic_week4","benzo_week4","Nonopioid_analgesic_week4","Pregabalin_week4","can_12m4","cig_12m4","opioid904","PSEQ_Score4",
                                       "PHQ9_Mod_sev5","GADMod2Sev5","Antidepressant_week5","Antipsychotic_week5","benzo_week5","Nonopioid_analgesic_week5","Pregabalin_week5","can_12m5","cig_12m5","opioid905","PSEQ_Score5",
                                       "PHQ9_Mod_sev6","GADMod2Sev6","Antidepressant_week6","Antipsychotic_week6","benzo_week6","Nonopioid_analgesic_week6","Pregabalin_week6","can_12m6","cig_12m6","opioid906","PSEQ_Score6"))
parallel::clusterEvalQ(cl, create.Learner("SL.ranger", params = list(num.trees = 250)))
parallel::clusterEvalQ(cl, SLlib <- list(Q=c("SL.mean","SL.glm","SL.gam"),
                                         g=c("SL.mean","SL.glm","SL.gam","SL.ranger_1")))

## 2+ DRINKS PER WEEK ##
dataalcuse1a <- lapply(impdatawide, function (x) {
  x$alcpain <- alc_pain_12m1 + alc_pain_12m2 + alc_pain_12m3 + alc_pain_12m4 + alc_pain_12m5 + alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(b_alc_ever,alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
dataalcuse2a <- lapply(impdatawide, function (x) {
  x$alcpain <- alc_pain_12m1 + alc_pain_12m2 + alc_pain_12m3 + alc_pain_12m4 + alc_pain_12m5 + alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(b_alc_ever,alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jtalcuse1a <- parLapply(cl=cl, dataalcuse1a, function (x) {ltmle(x,
                                                                 Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtalcuseres1a <- parLapply(cl=cl, jtalcuse1a, function (x) {summary(x)})
jtalcuseco1a <- matrix(unlist(parLapply(cl=cl, jtalcuseres1a, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtalcusese1a <- matrix(unlist(parLapply(cl=cl, jtalcuseres1a, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtalcusetco1a <- matrix(unlist(parLapply(cl=cl, jtalcuseres1a, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtalcusetse1a <- matrix(unlist(parLapply(cl=cl, jtalcuseres1a, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtalcusecco1a <- matrix(unlist(parLapply(cl=cl, jtalcuseres1a, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtalcusecse1a <- matrix(unlist(parLapply(cl=cl, jtalcuseres1a, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtalcuseres1a <- matrix(c(unlist(mi.meld(q=jtalcuseco1a,se=jtalcusese1a)),
                          unlist(mi.meld(q=jtalcusetco1a,se=jtalcusese1a)),
                          unlist(mi.meld(q=jtalcusecco1a,se=jtalcusecse1a))),nrow=1,ncol=6)


jtalcuse2a <- parLapply(cl=cl, dataalcuse2a, function (x) {ltmle(x,
                                                                 Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtalcuseres2a <- parLapply(cl=cl, jtalcuse2a, function (x) {summary(x)})
jtalcuseco2a <- matrix(unlist(parLapply(cl=cl, jtalcuseres2a, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtalcusese2a <- matrix(unlist(parLapply(cl=cl, jtalcuseres2a, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtalcusetco2a <- matrix(unlist(parLapply(cl=cl, jtalcuseres2a, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtalcusetse2a <- matrix(unlist(parLapply(cl=cl, jtalcuseres2a, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtalcusecco2a <- matrix(unlist(parLapply(cl=cl, jtalcuseres2a, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtalcusecse2a <- matrix(unlist(parLapply(cl=cl, jtalcuseres2a, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtalcuseres2a <- matrix(c(unlist(mi.meld(q=jtalcuseco2a,se=jtalcusese2a)),
                          unlist(mi.meld(q=jtalcusetco2a,se=jtalcusese2a)),
                          unlist(mi.meld(q=jtalcusecco2a,se=jtalcusecse2a))),nrow=1,ncol=6)

rm(list=c("jtalcuse1a","jtalcuseco1a","jtalcusese1a","jtalcusetco1a","jtalcusetse1a","jtalcusecco1a","jtalcusecse1a","dataalcuse1a",
          "jtalcuse2a","jtalcuseco2a","jtalcusese2a","jtalcusetco2a","jtalcusetse2a","jtalcusecco2a","jtalcusecse2a","dataalcuse2a"))

## 4+ DRINKS PER WEEK ##
dataalcuse1b <- lapply(impdatawide, function (x) {
  x$alcpain <- alc_pain_12m1 + alc_pain_12m2 + alc_pain_12m3 + alc_pain_12m4 + alc_pain_12m5 + alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(b_alc_ever,alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
dataalcuse2b <- lapply(impdatawide, function (x) {
  x$alcpain <- alc_pain_12m1 + alc_pain_12m2 + alc_pain_12m3 + alc_pain_12m4 + alc_pain_12m5 + alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(b_alc_ever,alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jtalcuse1b <- parLapply(cl=cl, dataalcuse1b, function (x) {ltmle(x,
                                                                 Anodes=c("freq22","freq23","freq24","freq25","freq26"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtalcuseres1b <- parLapply(cl=cl, jtalcuse1b, function (x) {summary(x)})
jtalcuseco1b <- matrix(unlist(parLapply(cl=cl, jtalcuseres1b, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtalcusese1b <- matrix(unlist(parLapply(cl=cl, jtalcuseres1b, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtalcusetco1b <- matrix(unlist(parLapply(cl=cl, jtalcuseres1b, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtalcusetse1b <- matrix(unlist(parLapply(cl=cl, jtalcuseres1b, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtalcusecco1b <- matrix(unlist(parLapply(cl=cl, jtalcuseres1b, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtalcusecse1b <- matrix(unlist(parLapply(cl=cl, jtalcuseres1b, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtalcuseres1b <- matrix(c(unlist(mi.meld(q=jtalcuseco1b,se=jtalcusese1b)),
                          unlist(mi.meld(q=jtalcusetco1b,se=jtalcusese1b)),
                          unlist(mi.meld(q=jtalcusecco1b,se=jtalcusecse1b))),nrow=1,ncol=6)

jtalcuse2b <- parLapply(cl=cl, dataalcuse2b, function (x) {ltmle(x,
                                                                 Anodes=c("freq22","freq23","freq24","freq25","freq26"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtalcuseres2b <- parLapply(cl=cl, jtalcuse2b, function (x) {summary(x)})
jtalcuseco2b <- matrix(unlist(parLapply(cl=cl, jtalcuseres2b, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtalcusese2b <- matrix(unlist(parLapply(cl=cl, jtalcuseres2b, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtalcusetco2b <- matrix(unlist(parLapply(cl=cl, jtalcuseres2b, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtalcusetse2b <- matrix(unlist(parLapply(cl=cl, jtalcuseres2b, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtalcusecco2b <- matrix(unlist(parLapply(cl=cl, jtalcuseres2b, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtalcusecse2b <- matrix(unlist(parLapply(cl=cl, jtalcuseres2b, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtalcuseres2b <- matrix(c(unlist(mi.meld(q=jtalcuseco2b,se=jtalcusese2b)),
                          unlist(mi.meld(q=jtalcusetco2b,se=jtalcusese2b)),
                          unlist(mi.meld(q=jtalcusecco2b,se=jtalcusecse2b))),nrow=1,ncol=6)

rm(list=c("jtalcuse1b","jtalcuseco1b","jtalcusese1b","jtalcusetco1b","jtalcusetse1b","jtalcusecco1b","jtalcusecse1b","dataalcuse1b",
          "jtalcuse2b","jtalcuseco2b","jtalcusese2b","jtalcusetco2b","jtalcusetse2b","jtalcusecco2b","jtalcusecse2b","dataalcuse2b"))

## BINGE DRINKING ANALYSIS ##
databinge1 <- lapply(impdatawide, function (x) {
  x$alcpain <- alc_pain_12m1 + alc_pain_12m2 + alc_pain_12m3 + alc_pain_12m4 + alc_pain_12m5 + alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(b_alc_ever,alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
databinge2 <- lapply(impdatawide, function (x) {
  x$alcpain <- alc_pain_12m1 + alc_pain_12m2 + alc_pain_12m3 + alc_pain_12m4 + alc_pain_12m5 + alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(b_alc_ever,alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jtbinge1 <- parLapply(cl=cl, databinge1, function (x) {ltmle(x,
                                                             Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeres1 <- parLapply(cl=cl, jtbinge1, function (x) {summary(x)})
jtbingeco1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtbingese1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtbingetco1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtbingetse1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtbingecco1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtbingecse1 <- matrix(unlist(parLapply(cl=cl, jtbingeres1, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtbingeres1 <- matrix(c(unlist(mi.meld(q=jtbingeco1,se=jtbingese1)),
                        unlist(mi.meld(q=jtbingetco1,se=jtbingese1)),
                        unlist(mi.meld(q=jtbingecco1,se=jtbingecse1))),nrow=1,ncol=6)

jtbinge2 <- parLapply(cl=cl, databinge2, function (x) {ltmle(x,
                                                             Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeres2 <- parLapply(cl=cl, jtbinge2, function (x) {summary(x)})
jtbingeco2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtbingese2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtbingetco2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtbingetse2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtbingecco2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtbingecse2 <- matrix(unlist(parLapply(cl=cl, jtbingeres2, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtbingeres2 <- matrix(c(unlist(mi.meld(q=jtbingeco2,se=jtbingese2)),
                        unlist(mi.meld(q=jtbingetco2,se=jtbingese2)),
                        unlist(mi.meld(q=jtbingecco2,se=jtbingecse2))),nrow=1,ncol=6)

rm(list=c("jtbinge1","jtbingeco1","jtbingese1","jtbingetco1","jtbingetse1","jtbingecco1","jtbingecse1","databinge1",
          "jtbinge2","jtbingeco2","jtbingese2","jtbingetco2","jtbingetse2","jtbingecco2","jtbingecse2","databinge2"))

## AUDIT PROBLEMATIC DRINKING ANALYSIS ##
dataauditprob1 <- lapply(impdatawide, function (x) {
  x$alcpain <- alc_pain_12m1 + alc_pain_12m2 + alc_pain_12m3 + alc_pain_12m4 + alc_pain_12m5 + alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(b_alc_ever,alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
dataauditprob2 <- lapply(impdatawide, function (x) {
  x$alcpain <- alc_pain_12m1 + alc_pain_12m2 + alc_pain_12m3 + alc_pain_12m4 + alc_pain_12m5 + alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(b_alc_ever,alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jtaudit1 <- parLapply(cl=cl, dataauditprob1, function (x) {ltmle(x,
                                                                 Anodes=c("auditcprob2","auditcprob3","auditcprob4","auditcprob5","auditcprob6"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtauditres1 <- parLapply(cl=cl, jtaudit1, function (x) {summary(x)})
jtauditco1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtauditse1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtaudittco1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtaudittse1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtauditcco1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtauditcse1 <- matrix(unlist(parLapply(cl=cl, jtauditres1, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtauditres1 <- matrix(c(unlist(mi.meld(q=jtauditco1,se=jtauditse1)),
                        unlist(mi.meld(q=jtaudittco1,se=jtauditse1)),
                        unlist(mi.meld(q=jtauditcco1,se=jtauditcse1))),nrow=1,ncol=6)

jtaudit2 <- parLapply(cl=cl, dataauditprob2, function (x) {ltmle(x,
                                                                 Anodes=c("auditcprob2","auditcprob3","auditcprob4","auditcprob5","auditcprob6"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtauditres2 <- parLapply(cl=cl, jtaudit2, function (x) {summary(x)})
jtauditco2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtauditse2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtaudittco2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtaudittse2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtauditcco2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtauditcse2 <- matrix(unlist(parLapply(cl=cl, jtauditres2, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtauditres2 <- matrix(c(unlist(mi.meld(q=jtauditco2,se=jtauditse2)),
                        unlist(mi.meld(q=jtaudittco2,se=jtauditse2)),
                        unlist(mi.meld(q=jtauditcco2,se=jtauditcse2))),nrow=1,ncol=6)

rm(list=c("jtaudit1","jtauditco1","jtauditse1","jtaudittco1","jtaudittse1","jtauditcco1","jtauditcse1","dataauditprob1",
          "jtaudit2","jtauditco2","jtauditse2","jtaudittco2","jtaudittse2","jtauditcco2","jtauditcse2","dataauditprob2"))

## COMBINE RESULTS INTO MATRIX FOR EXCEL ##
sensresults <- matrix(c(jtalcuseres1a[1:2],jtalcuseres2a[1:2],jtalcuseres1b[1:2],jtalcuseres2b[1:2],jtbingeres1[1:2],jtbingeres2[1:2],jtauditres1[1:2],jtauditres2[1:2]),byrow=TRUE,ncol=4,nrow=3)
rownames(sensresults) <- c("Drinking 2+ per week","Drinking 4+ times per week","Binge drinking","AUDIT-C Problematic drinking")
colnames(sensresults) <- c("BPI Pain Score coef","BPI Pain Score SE","BPI Pain Int coef","BPI Pain Int SE")

save(sensresults,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/sensresults.RData"))



