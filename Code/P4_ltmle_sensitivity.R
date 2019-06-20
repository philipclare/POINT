
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

## 2+ TIMES PER WEEK ##
## ALCOHOL USE TO COPE WITH PAIN ##
datafreq2ag1 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
datafreq2bg1 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jtfreq2ag1 <- parLapply(cl=cl, datafreq2a, function (x) {ltmle(x,
                                                             Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq2resag1 <- parLapply(cl=cl, jtfreq2ag1, function (x) {summary(x)})
jtfreq2coag1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag1, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq2seag1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag1, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq2tcoag1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag1, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq2tseag1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag1, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq2ccoag1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag1, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq2cseag1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag1, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq2resag1 <- matrix(c(unlist(mi.meld(q=jtfreq2coag1,se=jtfreq2seag1)),
                        unlist(mi.meld(q=jtfreq2tcoag1,se=jtfreq2seag1)),
                        unlist(mi.meld(q=jtfreq2ccoag1,se=jtfreq2cseag1))),nrow=1,ncol=6)

jtfreq2bg1 <- parLapply(cl=cl, datafreq2b, function (x) {ltmle(x,
                                                             Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq2resbg1 <- parLapply(cl=cl, jtfreq2b, function (x) {summary(x)})
jtfreq2cobg1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg1, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq2sebg1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg1, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq2tcobg1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg1, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq2tsebg1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg1, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq2ccobg1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg1, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq2csebg1 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg1, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq2resbg1 <- matrix(c(unlist(mi.meld(q=jtfreq2cobg1,se=jtfreq2sebg1)),
                        unlist(mi.meld(q=jtfreq2tcobg1,se=jtfreq2sebg1)),
                        unlist(mi.meld(q=jtfreq2ccobg1,se=jtfreq2csebg1))),nrow=1,ncol=6)

rm(list=c("jtfreq2ag1","jtfreq2coag1","jtfreq2seag1","jtfreq2tcoag1","jtfreq2tseag1","jtfreq2ccoag1","jtfreq2cseag1","datafreq2ag1",
          "jtfreq2bg1","jtfreq2coag1","jtfreq2seag1","jtfreq2tcoag1","jtfreq2tseag1","jtfreq2ccoag1","jtfreq2cseag1","datafreq2bg1"))

## NO ALCOHOL USE TO COPE WITH PAIN ##
datafreq2ag2 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain==0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
datafreq2bg2 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain==0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jtfreq2ag2 <- parLapply(cl=cl, datafreq2a, function (x) {ltmle(x,
                                                               Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                               Lnodes=Lvars2,
                                                               Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                               survivalOutcome=FALSE,
                                                               SL.library=SLlib,
                                                               abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq2resag2 <- parLapply(cl=cl, jtfreq2ag2, function (x) {summary(x)})
jtfreq2coag2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag2, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq2seag2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag2, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq2tcoag2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag2, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq2tseag2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag2, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq2ccoag2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag2, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq2cseag2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resag2, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq2resag2 <- matrix(c(unlist(mi.meld(q=jtfreq2coag2,se=jtfreq2seag2)),
                          unlist(mi.meld(q=jtfreq2tcoag2,se=jtfreq2seag2)),
                          unlist(mi.meld(q=jtfreq2ccoag2,se=jtfreq2cseag2))),nrow=1,ncol=6)

jtfreq2bg2 <- parLapply(cl=cl, datafreq2b, function (x) {ltmle(x,
                                                               Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                               Lnodes=Lvars2,
                                                               Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                               survivalOutcome=FALSE,
                                                               SL.library=SLlib,
                                                               abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq2resbg2 <- parLapply(cl=cl, jtfreq2b, function (x) {summary(x)})
jtfreq2cobg2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg2, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq2sebg2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg2, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq2tcobg2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg2, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq2tsebg2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg2, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq2ccobg2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg2, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq2csebg2 <- matrix(unlist(parLapply(cl=cl, jtfreq2resbg2, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq2resbg2 <- matrix(c(unlist(mi.meld(q=jtfreq2cobg2,se=jtfreq2sebg2)),
                          unlist(mi.meld(q=jtfreq2tcobg2,se=jtfreq2sebg2)),
                          unlist(mi.meld(q=jtfreq2ccobg2,se=jtfreq2csebg2))),nrow=1,ncol=6)

rm(list=c("jtfreq2ag2","jtfreq2coag2","jtfreq2seag2","jtfreq2tcoag2","jtfreq2tseag2","jtfreq2ccoag2","jtfreq2cseag2","datafreq2ag2",
          "jtfreq2bg2","jtfreq2coag2","jtfreq2seag2","jtfreq2tcoag2","jtfreq2tseag2","jtfreq2ccoag2","jtfreq2cseag2","datafreq2bg2"))

## 4+ TIMES PER WEEK ##
## ALCOHOL USE TO COPE WITH PAIN ##
datafreq4ag1 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
datafreq4bg1 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jtfreq4ag1 <- parLapply(cl=cl, datafreq4ag1, function (x) {ltmle(x,
                                                             Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq4resag1 <- parLapply(cl=cl, jtfreq4ag1, function (x) {summary(x)})
jtfreq4coag1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag1, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq4seag1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag1, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq4tcoag1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag1, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq4tseag1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag1, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq4ccoag1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag1, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq4cseag1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag1, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq4resag1 <- matrix(c(unlist(mi.meld(q=jtfreq4coag1,se=jtfreq4seag1)),
                        unlist(mi.meld(q=jtfreq4tcoag1,se=jtfreq4seag1)),
                        unlist(mi.meld(q=jtfreq4ccoag1,se=jtfreq4cseag1))),nrow=1,ncol=6)

jtfreq4bg1 <- parLapply(cl=cl, datafreq4bg1, function (x) {ltmle(x,
                                                             Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq4resbg1 <- parLapply(cl=cl, jtfreq4bg1, function (x) {summary(x)})
jtfreq4cobg1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg1, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq4sebg1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg1, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq4tcobg1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg1, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq4tsebg1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg1, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq4ccobg1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg1, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq4csebg1 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg1, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq4resbg1 <- matrix(c(unlist(mi.meld(q=jtfreq4cobg1,se=jtfreq4sebg1)),
                        unlist(mi.meld(q=jtfreq4tcobg1,se=jtfreq4sebg1)),
                        unlist(mi.meld(q=jtfreq4ccobg1,se=jtfreq4csebg1))),nrow=1,ncol=6)

rm(list=c("jtfreq4ag1","jtfreq4coag1","jtfreq4seag1","jtfreq4tcoag1","jtfreq4tseag1","jtfreq4ccoag1","jtfreq4cseag1","datafreq4ag1",
          "jtfreq4bg1","jtfreq4cobg1","jtfreq4sebg1","jtfreq4tcobg1","jtfreq4tsebg1","jtfreq4ccobg1","jtfreq4csebg1","datafreq4bg1"))

## NO ALCOHOL USE TO COPE WITH PAIN ##
datafreq4ag2 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain==0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
datafreq4bg2 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain==0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jtfreq4ag2 <- parLapply(cl=cl, datafreq4ag2, function (x) {ltmle(x,
                                                                 Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq4resag2 <- parLapply(cl=cl, jtfreq4ag2, function (x) {summary(x)})
jtfreq4coag2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag2, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq4seag2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag2, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq4tcoag2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag2, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq4tseag2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag2, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq4ccoag2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag2, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq4cseag2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resag2, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq4resag2 <- matrix(c(unlist(mi.meld(q=jtfreq4coag2,se=jtfreq4seag2)),
                          unlist(mi.meld(q=jtfreq4tcoag2,se=jtfreq4seag2)),
                          unlist(mi.meld(q=jtfreq4ccoag2,se=jtfreq4cseag2))),nrow=1,ncol=6)

jtfreq4bg2 <- parLapply(cl=cl, datafreq4bg2, function (x) {ltmle(x,
                                                                 Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq4resbg2 <- parLapply(cl=cl, jtfreq4bg2, function (x) {summary(x)})
jtfreq4cobg2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg2, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq4sebg2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg2, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq4tcobg2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg2, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq4tsebg2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg2, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq4ccobg2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg2, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq4csebg2 <- matrix(unlist(parLapply(cl=cl, jtfreq4resbg2, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq4resbg2 <- matrix(c(unlist(mi.meld(q=jtfreq4cobg2,se=jtfreq4sebg2)),
                          unlist(mi.meld(q=jtfreq4tcobg2,se=jtfreq4sebg2)),
                          unlist(mi.meld(q=jtfreq4ccobg2,se=jtfreq4csebg2))),nrow=1,ncol=6)

rm(list=c("jtfreq4ag2","jtfreq4coag2","jtfreq4seag2","jtfreq4tcoag2","jtfreq4tseag2","jtfreq4ccoag2","jtfreq4cseag2","datafreq4ag2",
          "jtfreq4bg2","jtfreq4cobg2","jtfreq4sebg2","jtfreq4tcobg2","jtfreq4tsebg2","jtfreq4ccobg2","jtfreq4csebg2","datafreq4bg2"))

## BINGE DRINKING ANALYSIS ##
## ALCOHOL USE TO COPE WITH PAIN ##
databinge1g1 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
databinge2g1 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jtbingeag1 <- parLapply(cl=cl, databinge1g1, function (x) {ltmle(x,
                                                             Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeresag1 <- parLapply(cl=cl, jtbingeag1, function (x) {summary(x)})
jtbingecoag1 <- matrix(unlist(parLapply(cl=cl, jtbingeresag1, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtbingeseag1 <- matrix(unlist(parLapply(cl=cl, jtbingeresag1, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtbingetcoag1 <- matrix(unlist(parLapply(cl=cl, jtbingeresag1, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtbingetseag1 <- matrix(unlist(parLapply(cl=cl, jtbingeresag1, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtbingeccoag1 <- matrix(unlist(parLapply(cl=cl, jtbingeresag1, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtbingecseag1 <- matrix(unlist(parLapply(cl=cl, jtbingeresag1, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtbingeresag1 <- matrix(c(unlist(mi.meld(q=jtbingecoag1,se=jtbingeseag1)),
                        unlist(mi.meld(q=jtbingetcoag1,se=jtbingeseag1)),
                        unlist(mi.meld(q=jtbingeccoag1,se=jtbingecseag1))),nrow=1,ncol=6)

jtbingebg1 <- parLapply(cl=cl, databinge2g1, function (x) {ltmle(x,
                                                             Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeresbg1 <- parLapply(cl=cl, jtbinge2g1, function (x) {summary(x)})
jtbingecobg1 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg1, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtbingesebg1 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg1, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtbingetcobg1 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg1, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtbingetsebg1 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg1, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtbingeccobg1 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg1, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtbingecsebg1 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg1, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtbingeresbg1 <- matrix(c(unlist(mi.meld(q=jtbingecobg1,se=jtbingesebg1)),
                        unlist(mi.meld(q=jtbingetcobg1,se=jtbingesebg1)),
                        unlist(mi.meld(q=jtbingeccobg1,se=jtbingecsebg1))),nrow=1,ncol=6)

rm(list=c("jtbingeag1","jtbingecoag1","jtbingeseag1","jtbingetcoag1","jtbingetseag1","jtbingeccoag1","jtbingecseag1","databingeag1",
          "jtbingebg1","jtbingecobg1","jtbingesebg1","jtbingetcobg1","jtbingetsebg1","jtbingeccobg1","jtbingecsebg1","databingebg1"))

## NO ALCOHOL USE TO COPE WITH PAIN ##
databinge1g2 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain==0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
databinge2g2 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             auditcprob1,auditcprob2,auditcprob3,auditcprob4,auditcprob5,auditcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain==0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jtbingeag2 <- parLapply(cl=cl, databinge1g2, function (x) {ltmle(x,
                                                                 Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeresag2 <- parLapply(cl=cl, jtbingeag2, function (x) {summary(x)})
jtbingecoag2 <- matrix(unlist(parLapply(cl=cl, jtbingeresag2, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtbingeseag2 <- matrix(unlist(parLapply(cl=cl, jtbingeresag2, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtbingetcoag2 <- matrix(unlist(parLapply(cl=cl, jtbingeresag2, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtbingetseag2 <- matrix(unlist(parLapply(cl=cl, jtbingeresag2, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtbingeccoag2 <- matrix(unlist(parLapply(cl=cl, jtbingeresag2, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtbingecseag2 <- matrix(unlist(parLapply(cl=cl, jtbingeresag2, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtbingeresag2 <- matrix(c(unlist(mi.meld(q=jtbingecoag2,se=jtbingeseag2)),
                          unlist(mi.meld(q=jtbingetcoag2,se=jtbingeseag2)),
                          unlist(mi.meld(q=jtbingeccoag2,se=jtbingecseag2))),nrow=1,ncol=6)

jtbingebg2 <- parLapply(cl=cl, databinge2g2, function (x) {ltmle(x,
                                                                 Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                                 Lnodes=Lvars2,
                                                                 Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                                 survivalOutcome=FALSE,
                                                                 SL.library=SLlib,
                                                                 abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeresbg2 <- parLapply(cl=cl, jtbinge2g2, function (x) {summary(x)})
jtbingecobg2 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg2, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtbingesebg2 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg2, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtbingetcobg2 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg2, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtbingetsebg2 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg2, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtbingeccobg2 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg2, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtbingecsebg2 <- matrix(unlist(parLapply(cl=cl, jtbingeresbg2, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtbingeresbg2 <- matrix(c(unlist(mi.meld(q=jtbingecobg2,se=jtbingesebg2)),
                          unlist(mi.meld(q=jtbingetcobg2,se=jtbingesebg2)),
                          unlist(mi.meld(q=jtbingeccobg2,se=jtbingecsebg2))),nrow=1,ncol=6)

rm(list=c("jtbingeag2","jtbingecoag2","jtbingeseag2","jtbingetcoag2","jtbingetseag2","jtbingeccoag2","jtbingecseag2","databingeag2",
          "jtbingebg2","jtbingecobg2","jtbingesebg2","jtbingetcobg2","jtbingetsebg2","jtbingeccobg2","jtbingecsebg2","databingebg2"))

## AUDIT PROBLEMATIC DRINKING ANALYSIS ##
## ALCOHOL USE TO COPE WITH PAIN ##
datahazardag1 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
datahazardbg1 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain>0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jthazardag1 <- parLapply(cl=cl, datahazardag1, function (x) {ltmle(x,
                                                               Anodes=c("hazardcprob2","hazardcprob3","hazardcprob4","hazardcprob5","hazardcprob6"),
                                                               Lnodes=Lvars2,
                                                               Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                               survivalOutcome=FALSE,
                                                               SL.library=SLlib,
                                                               abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jthazardresag1 <- parLapply(cl=cl, jthazardag1, function (x) {summary(x)})
jthazardcoag1 <- matrix(unlist(parLapply(cl=cl, jthazardresag1, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jthazardseag1 <- matrix(unlist(parLapply(cl=cl, jthazardresag1, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jthazardtcoag1 <- matrix(unlist(parLapply(cl=cl, jthazardresag1, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jthazardtseag1 <- matrix(unlist(parLapply(cl=cl, jthazardresag1, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jthazardccoag1 <- matrix(unlist(parLapply(cl=cl, jthazardresag1, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jthazardcseag1 <- matrix(unlist(parLapply(cl=cl, jthazardresag1, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jthazardresag1 <- matrix(c(unlist(mi.meld(q=jthazardcoag1,se=jthazardseag1)),
                         unlist(mi.meld(q=jthazardtcoag1,se=jthazardseag1)),
                         unlist(mi.meld(q=jthazardccoag1,se=jthazardcseag1))),nrow=1,ncol=6)

jthazardbg1 <- parLapply(cl=cl, datahazardbg1, function (x) {ltmle(x,
                                                               Anodes=c("hazardcprob2","hazardcprob3","hazardcprob4","hazardcprob5","hazardcprob6"),
                                                               Lnodes=Lvars2,
                                                               Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                               survivalOutcome=FALSE,
                                                               SL.library=SLlib,
                                                               abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jthazardresbg1 <- parLapply(cl=cl, jthazardbg1, function (x) {summary(x)})
jthazardcobg1 <- matrix(unlist(parLapply(cl=cl, jthazardresbg1, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jthazardsebg1 <- matrix(unlist(parLapply(cl=cl, jthazardresbg1, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jthazardtcobg1 <- matrix(unlist(parLapply(cl=cl, jthazardresbg1, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jthazardtsebg1 <- matrix(unlist(parLapply(cl=cl, jthazardresbg1, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jthazardccobg1 <- matrix(unlist(parLapply(cl=cl, jthazardresbg1, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jthazardcsebg1 <- matrix(unlist(parLapply(cl=cl, jthazardresbg1, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jthazardresbg1 <- matrix(c(unlist(mi.meld(q=jthazardcobg1,se=jthazardsebg1)),
                         unlist(mi.meld(q=jthazardtcobg1,se=jthazardsebg1)),
                         unlist(mi.meld(q=jthazardccobg1,se=jthazardcsebg1))),nrow=1,ncol=6)

rm(list=c("jthazardag1","jthazardcoag1","jthazardseag1","jthazardtcoag1","jthazardtseag1","jthazardccoag1","jthazardcseag1","datahazardag1",
          "jthazardbg1","jthazardcobg1","jthazardsebg1","jthazardtcobg1","jthazardtsebg1","jthazardccobg1","jthazardcsebg1","datahazardbg1"))

## NO ALCOHOL USE TO COPE WITH PAIN ##
datahazardag2 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain==0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})
datahazardbg2 <- lapply(impdatawide, function (x) {
  x$alcpain <- x$alc_pain_12m2 + x$alc_pain_12m3 + x$alc_pain_12m4 + x$alc_pain_12m5 + x$alc_pain_12m6
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1))
  x <- x[!(x$alcpain==0),]
  x <- subset(x, select = -c(alcpain,alc_pain_12m1,alc_pain_12m2,alc_pain_12m3,alc_pain_12m4,alc_pain_12m5,alc_pain_12m6))
  x
})

jthazardag2 <- parLapply(cl=cl, datahazardag2, function (x) {ltmle(x,
                                                                   Anodes=c("hazardcprob2","hazardcprob3","hazardcprob4","hazardcprob5","hazardcprob6"),
                                                                   Lnodes=Lvars2,
                                                                   Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                                   survivalOutcome=FALSE,
                                                                   SL.library=SLlib,
                                                                   abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jthazardresag2 <- parLapply(cl=cl, jthazardag2, function (x) {summary(x)})
jthazardcoag2 <- matrix(unlist(parLapply(cl=cl, jthazardresag2, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jthazardseag2 <- matrix(unlist(parLapply(cl=cl, jthazardresag2, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jthazardtcoag2 <- matrix(unlist(parLapply(cl=cl, jthazardresag2, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jthazardtseag2 <- matrix(unlist(parLapply(cl=cl, jthazardresag2, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jthazardccoag2 <- matrix(unlist(parLapply(cl=cl, jthazardresag2, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jthazardcseag2 <- matrix(unlist(parLapply(cl=cl, jthazardresag2, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jthazardresag2 <- matrix(c(unlist(mi.meld(q=jthazardcoag2,se=jthazardseag2)),
                           unlist(mi.meld(q=jthazardtcoag2,se=jthazardseag2)),
                           unlist(mi.meld(q=jthazardccoag2,se=jthazardcseag2))),nrow=1,ncol=6)

jthazardbg2 <- parLapply(cl=cl, datahazardbg2, function (x) {ltmle(x,
                                                                   Anodes=c("hazardcprob2","hazardcprob3","hazardcprob4","hazardcprob5","hazardcprob6"),
                                                                   Lnodes=Lvars2,
                                                                   Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                                   survivalOutcome=FALSE,
                                                                   SL.library=SLlib,
                                                                   abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jthazardresbg2 <- parLapply(cl=cl, jthazardbg2, function (x) {summary(x)})
jthazardcobg2 <- matrix(unlist(parLapply(cl=cl, jthazardresbg2, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jthazardsebg2 <- matrix(unlist(parLapply(cl=cl, jthazardresbg2, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jthazardtcobg2 <- matrix(unlist(parLapply(cl=cl, jthazardresbg2, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jthazardtsebg2 <- matrix(unlist(parLapply(cl=cl, jthazardresbg2, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jthazardccobg2 <- matrix(unlist(parLapply(cl=cl, jthazardresbg2, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jthazardcsebg2 <- matrix(unlist(parLapply(cl=cl, jthazardresbg2, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jthazardresbg2 <- matrix(c(unlist(mi.meld(q=jthazardcobg2,se=jthazardsebg2)),
                           unlist(mi.meld(q=jthazardtcobg2,se=jthazardsebg2)),
                           unlist(mi.meld(q=jthazardccobg2,se=jthazardcsebg2))),nrow=1,ncol=6)

rm(list=c("jthazardag2","jthazardcoag2","jthazardseag2","jthazardtcoag2","jthazardtseag2","jthazardccoag2","jthazardcseag2","datahazardag2",
          "jthazardbg2","jthazardcobg2","jthazardsebg2","jthazardtcobg2","jthazardtsebg2","jthazardccobg2","jthazardcsebg2","datahazardbg2"))

## COMBINE RESULTS INTO MATRIX FOR EXCEL ##
## ALCOHOL USE TO COPE WITH PAIN ##
g1results <- matrix(c(jtfreq2resag1[1:2],jtfreq2resbg1[1:2],jtfreq4resag1[1:2],jtfreq4resbg1[1:2],
                      jtbingeresag1[1:2],jtbingeresbg1[1:2],jthazardresag1[1:2],jthazardresbg1[1:2]),byrow=TRUE,ncol=4,nrow=4)
rownames(g1results) <- c("Drinking 2+ per week","Drinking 4+ times per week","Binge drinking","AUDIT-C Problematic drinking")
colnames(g1results) <- c("BPI Pain Score coef","BPI Pain Score SE","BPI Pain Int coef","BPI Pain Int SE")

save(g1results,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/sensbresults.RData"))

## NO ALCOHOL USE TO COPE WITH PAIN ##
g2results <- matrix(c(jtfreq2resag2[1:2],jtfreq2resbg2[1:2],jtfreq4resag2[1:2],jtfreq4resbg2[1:2],
                      jtbingeresag2[1:2],jtbingeresbg2[1:2],jthazardresag2[1:2],jthazardresbg2[1:2]),byrow=TRUE,ncol=4,nrow=4)
rownames(g2results) <- c("Drinking 2+ per week","Drinking 4+ times per week","Binge drinking","AUDIT-C Problematic drinking")
colnames(g2results) <- c("BPI Pain Score coef","BPI Pain Score SE","BPI Pain Int coef","BPI Pain Int SE")

save(g2results,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/sensbresults.RData"))
