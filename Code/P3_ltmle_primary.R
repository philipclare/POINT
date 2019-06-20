
## SYNTAX FILE 3 - PRIMARY ANALYSIS                          ##

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
parallel::clusterEvalQ(cl, Lvars2 <- c("PHQ9_Mod_sev2","GADMod2Sev2","Antidepressant_week2","Antipsychotic_week2","benzo_week2","Nonopioid_analgesic_week2","Pregabalin_week2","can_12m2","cig_12m2","opioid902","PSEQ_Score2","alc_pain_12m2",
                                       "PHQ9_Mod_sev3","GADMod2Sev3","Antidepressant_week3","Antipsychotic_week3","benzo_week3","Nonopioid_analgesic_week3","Pregabalin_week3","can_12m3","cig_12m3","opioid903","PSEQ_Score3","alc_pain_12m3",
                                       "PHQ9_Mod_sev4","GADMod2Sev4","Antidepressant_week4","Antipsychotic_week4","benzo_week4","Nonopioid_analgesic_week4","Pregabalin_week4","can_12m4","cig_12m4","opioid904","PSEQ_Score4","alc_pain_12m4",
                                       "PHQ9_Mod_sev5","GADMod2Sev5","Antidepressant_week5","Antipsychotic_week5","benzo_week5","Nonopioid_analgesic_week5","Pregabalin_week5","can_12m5","cig_12m5","opioid905","PSEQ_Score5","alc_pain_12m5",
                                       "PHQ9_Mod_sev6","GADMod2Sev6","Antidepressant_week6","Antipsychotic_week6","benzo_week6","Nonopioid_analgesic_week6","Pregabalin_week6","can_12m6","cig_12m6","opioid906","PSEQ_Score6","alc_pain_12m6"))
parallel::clusterEvalQ(cl, create.Learner("SL.ranger", params = list(num.trees = 250)))
parallel::clusterEvalQ(cl, SLlib <- list(Q=c("SL.mean","SL.glm","SL.gam"),
                                         g=c("SL.mean","SL.glm","SL.gam","SL.ranger_1")))

## 2+ DRINKS PER WEEK ##
datafreq2a <- lapply(impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             hazardcprob1,hazardcprob2,hazardcprob3,hazardcprob4,hazardcprob5,hazardcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x
})
datafreq2b <- lapply(impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             hazardcprob1,hazardcprob2,hazardcprob3,hazardcprob4,hazardcprob5,hazardcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x
})

jtfreq2a <- parLapply(cl=cl, datafreq2a, function (x) {ltmle(x,
                                                             Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq2resa <- parLapply(cl=cl, jtfreq2a, function (x) {summary(x)})
jtfreq2coa <- matrix(unlist(parLapply(cl=cl, jtfreq2resa, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq2sea <- matrix(unlist(parLapply(cl=cl, jtfreq2resa, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq2tcoa <- matrix(unlist(parLapply(cl=cl, jtfreq2resa, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq2tsea <- matrix(unlist(parLapply(cl=cl, jtfreq2resa, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq2ccoa <- matrix(unlist(parLapply(cl=cl, jtfreq2resa, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq2csea <- matrix(unlist(parLapply(cl=cl, jtfreq2resa, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq2resa <- matrix(c(unlist(mi.meld(q=jtfreq2coa,se=jtfreq2sea)),
                  unlist(mi.meld(q=jtfreq2tcoa,se=jtfreq2sea)),
                  unlist(mi.meld(q=jtfreq2ccoa,se=jtfreq2csea))),nrow=1,ncol=6)

jtfreq2b <- parLapply(cl=cl, datafreq2b, function (x) {ltmle(x,
                                                             Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq2resb <- parLapply(cl=cl, jtfreq2b, function (x) {summary(x)})
jtfreq2cob <- matrix(unlist(parLapply(cl=cl, jtfreq2resb, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq2seb <- matrix(unlist(parLapply(cl=cl, jtfreq2resb, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq2tcob <- matrix(unlist(parLapply(cl=cl, jtfreq2resb, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq2tseb <- matrix(unlist(parLapply(cl=cl, jtfreq2resb, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq2ccob <- matrix(unlist(parLapply(cl=cl, jtfreq2resb, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq2cseb <- matrix(unlist(parLapply(cl=cl, jtfreq2resb, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq2resb <- matrix(c(unlist(mi.meld(q=jtfreq2cob,se=jtfreq2seb)),
                         unlist(mi.meld(q=jtfreq2tcob,se=jtfreq2seb)),
                         unlist(mi.meld(q=jtfreq2ccob,se=jtfreq2cseb))),nrow=1,ncol=6)

rm(list=c("jtfreq2a","jtfreq2coa","jtfreq2sea","jtfreq2tcoa","jtfreq2tsea","jtfreq2ccoa","jtfreq2csea","datafreq2a",
          "jtfreq2b","jtfreq2coa","jtfreq2sea","jtfreq2tcoa","jtfreq2tsea","jtfreq2ccoa","jtfreq2csea","datafreq2b"))

## 4+ DRINKS PER WEEK ##
datafreq4a <- lapply(impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             hazardcprob1,hazardcprob2,hazardcprob3,hazardcprob4,hazardcprob5,hazardcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x
})
datafreq4b <- lapply(impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             hazardcprob1,hazardcprob2,hazardcprob3,hazardcprob4,hazardcprob5,hazardcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x
})

jtfreq4a <- parLapply(cl=cl, datafreq4a, function (x) {ltmle(x,
                                                             Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq4resa <- parLapply(cl=cl, jtfreq4a, function (x) {summary(x)})
jtfreq4coa <- matrix(unlist(parLapply(cl=cl, jtfreq4resa, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq4sea <- matrix(unlist(parLapply(cl=cl, jtfreq4resa, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq4tcoa <- matrix(unlist(parLapply(cl=cl, jtfreq4resa, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq4tsea <- matrix(unlist(parLapply(cl=cl, jtfreq4resa, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq4ccoa <- matrix(unlist(parLapply(cl=cl, jtfreq4resa, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq4csea <- matrix(unlist(parLapply(cl=cl, jtfreq4resa, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq4resa <- matrix(c(unlist(mi.meld(q=jtfreq4coa,se=jtfreq4sea)),
                        unlist(mi.meld(q=jtfreq4tcoa,se=jtfreq4sea)),
                        unlist(mi.meld(q=jtfreq4ccoa,se=jtfreq4csea))),nrow=1,ncol=6)

jtfreq4b <- parLapply(cl=cl, datafreq4b, function (x) {ltmle(x,
                                                             Anodes=c("freq12","freq13","freq14","freq15","freq16"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtfreq4resb <- parLapply(cl=cl, jtfreq4b, function (x) {summary(x)})
jtfreq4cob <- matrix(unlist(parLapply(cl=cl, jtfreq4resb, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtfreq4seb <- matrix(unlist(parLapply(cl=cl, jtfreq4resb, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtfreq4tcob <- matrix(unlist(parLapply(cl=cl, jtfreq4resb, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtfreq4tseb <- matrix(unlist(parLapply(cl=cl, jtfreq4resb, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtfreq4ccob <- matrix(unlist(parLapply(cl=cl, jtfreq4resb, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtfreq4cseb <- matrix(unlist(parLapply(cl=cl, jtfreq4resb, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtfreq4resb <- matrix(c(unlist(mi.meld(q=jtfreq4cob,se=jtfreq4seb)),
                        unlist(mi.meld(q=jtfreq4tcob,se=jtfreq4seb)),
                        unlist(mi.meld(q=jtfreq4ccob,se=jtfreq4cseb))),nrow=1,ncol=6)

rm(list=c("jtfreq4a","jtfreq4coa","jtfreq4sea","jtfreq4tcoa","jtfreq4tsea","jtfreq4ccoa","jtfreq4csea","datafreq4a",
          "jtfreq4b","jtfreq4cob","jtfreq4seb","jtfreq4tcob","jtfreq4tseb","jtfreq4ccob","jtfreq4cseb","datafreq4b"))

## BINGE DRINKING ANALYSIS ##
databingea <- lapply(impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             hazardcprob1,hazardcprob2,hazardcprob3,hazardcprob4,hazardcprob5,hazardcprob6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x
})
databingeb <- lapply(impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             hazardcprob1,hazardcprob2,hazardcprob3,hazardcprob4,hazardcprob5,hazardcprob6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x
})

jtbingea <- parLapply(cl=cl, databinge1, function (x) {ltmle(x,
                                                             Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeresa <- parLapply(cl=cl, jtbingea, function (x) {summary(x)})
jtbingecoa <- matrix(unlist(parLapply(cl=cl, jtbingeresa, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtbingesea <- matrix(unlist(parLapply(cl=cl, jtbingeresa, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtbingetcoa <- matrix(unlist(parLapply(cl=cl, jtbingeresa, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtbingetsea <- matrix(unlist(parLapply(cl=cl, jtbingeresa, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtbingeccoa <- matrix(unlist(parLapply(cl=cl, jtbingeresa, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtbingecsea <- matrix(unlist(parLapply(cl=cl, jtbingeresa, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtbingeresa <- matrix(c(unlist(mi.meld(q=jtbingecoa,se=jtbingesea)),
                         unlist(mi.meld(q=jtbingetcoa,se=jtbingesea)),
                         unlist(mi.meld(q=jtbingeccoa,se=jtbingecsea))),nrow=1,ncol=6)

jtbingeb <- parLapply(cl=cl, databinge2, function (x) {ltmle(x,
                                                             Anodes=c("binge2","binge3","binge4","binge5","binge6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jtbingeresb <- parLapply(cl=cl, jtbinge2, function (x) {summary(x)})
jtbingecob <- matrix(unlist(parLapply(cl=cl, jtbingeresb, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jtbingeseb <- matrix(unlist(parLapply(cl=cl, jtbingeresb, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jtbingetcob <- matrix(unlist(parLapply(cl=cl, jtbingeresb, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jtbingetseb <- matrix(unlist(parLapply(cl=cl, jtbingeresb, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jtbingeccob <- matrix(unlist(parLapply(cl=cl, jtbingeresb, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jtbingecseb <- matrix(unlist(parLapply(cl=cl, jtbingeresb, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jtbingeresb <- matrix(c(unlist(mi.meld(q=jtbingecob,se=jtbingeseb)),
                        unlist(mi.meld(q=jtbingetcob,se=jtbingeseb)),
                        unlist(mi.meld(q=jtbingeccob,se=jtbingecseb))),nrow=1,ncol=6)

rm(list=c("jtbingea","jtbingecoa","jtbingesea","jtbingetcoa","jtbingetsea","jtbingeccoa","jtbingecsea","databingea",
          "jtbingeb","jtbingecob","jtbingeseb","jtbingetcob","jtbingetseb","jtbingeccob","jtbingecseb","databingeb"))

## hazard PROBLEMATIC DRINKING ANALYSIS ##
datahazarda <- lapply(impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             BPI_interference1,BPI_interference2,BPI_interference3,BPI_interference4,BPI_interference5,BPI_interference6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x
})
datahazardb <- lapply(impdatawide, function (x) {
  x <- subset(x, select = -c(b_alc_ever,
                             freq11,freq12,freq13,freq14,freq15,freq16,
                             freq21,freq22,freq23,freq24,freq25,freq26,
                             binge1,binge2,binge3,binge4,binge5,binge6,
                             BPI_PScore1,BPI_PScore2,BPI_PScore3,BPI_PScore4,BPI_PScore5,BPI_PScore6,
                             PHQ9_Mod_sev1,GADMod2Sev1,Antidepressant_week1,Antipsychotic_week1,benzo_week1,Nonopioid_analgesic_week1,
                             Pregabalin_week1,can_12m1,cig_12m1,opioid901,PSEQ_Score1,alc_pain_12m1))
  x
})

jthazarda <- parLapply(cl=cl, datahazarda, function (x) {ltmle(x,
                                                           Anodes=c("hazardcprob2","hazardcprob3","hazardcprob4","hazardcprob5","hazardcprob6"),
                                                           Lnodes=Lvars2,
                                                           Ynodes=c("BPI_PScore2","BPI_PScore3","BPI_PScore4","BPI_PScore5","BPI_PScore6"),
                                                           survivalOutcome=FALSE,
                                                           SL.library=SLlib,
                                                           abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jthazardresa <- parLapply(cl=cl, jthazarda, function (x) {summary(x)})
jthazardcoa <- matrix(unlist(parLapply(cl=cl, jthazardresa, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jthazardsea <- matrix(unlist(parLapply(cl=cl, jthazardresa, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jthazardtcoa <- matrix(unlist(parLapply(cl=cl, jthazardresa, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jthazardtsea <- matrix(unlist(parLapply(cl=cl, jthazardresa, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jthazardccoa <- matrix(unlist(parLapply(cl=cl, jthazardresa, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jthazardcsea <- matrix(unlist(parLapply(cl=cl, jthazardresa, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jthazardresa <- matrix(c(unlist(mi.meld(q=jthazardcoa,se=jthazardsea)),
                        unlist(mi.meld(q=jthazardtcoa,se=jthazardsea)),
                        unlist(mi.meld(q=jthazardccoa,se=jthazardcsea))),nrow=1,ncol=6)

jthazardb <- parLapply(cl=cl, datahazardb, function (x) {ltmle(x,
                                                             Anodes=c("hazardcprob2","hazardcprob3","hazardcprob4","hazardcprob5","hazardcprob6"),
                                                             Lnodes=Lvars2,
                                                             Ynodes=c("BPI_interference2","BPI_interference3","BPI_interference4","BPI_interference5","BPI_interference6"),
                                                             survivalOutcome=FALSE,
                                                             SL.library=SLlib,
                                                             abar=list(c(1,1,1,1,1),c(0,0,0,0,0)))})
jthazardresb <- parLapply(cl=cl, jthazardb, function (x) {summary(x)})
jthazardcob <- matrix(unlist(parLapply(cl=cl, jthazardresb, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
jthazardseb <- matrix(unlist(parLapply(cl=cl, jthazardresb, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
jthazardtcob <- matrix(unlist(parLapply(cl=cl, jthazardresb, function (x) {x$effect.measures$treatment$estimate})),nrow=nimpute,ncol=1)
jthazardtseb <- matrix(unlist(parLapply(cl=cl, jthazardresb, function (x) {x$effect.measures$treatment$std.dev})),nrow=nimpute,ncol=1)
jthazardccob <- matrix(unlist(parLapply(cl=cl, jthazardresb, function (x) {x$effect.measures$control$estimate})),nrow=nimpute,ncol=1)
jthazardcseb <- matrix(unlist(parLapply(cl=cl, jthazardresb, function (x) {x$effect.measures$control$std.dev})),nrow=nimpute,ncol=1)
jthazardresb <- matrix(c(unlist(mi.meld(q=jthazardcob,se=jthazardseb)),
                        unlist(mi.meld(q=jthazardtcob,se=jthazardseb)),
                        unlist(mi.meld(q=jthazardccob,se=jthazardcseb))),nrow=1,ncol=6)

rm(list=c("jthazarda","jthazardcoa","jthazardsea","jthazardtcoa","jthazardtsea","jthazardccoa","jthazardcsea","datahazarda",
          "jthazardb","jthazardcob","jthazardseb","jthazardtcob","jthazardtseb","jthazardccob","jthazardcseb","datahazardb"))

## COMBINE RESULTS INTO MATRIX FOR EXCEL ##
jtresults <- matrix(c(jtfreq2resa[1:2],jtfreq2resa[1:2],jtfreq4resb[1:2],jtfreq4resb[1:2],
                      jtbingeresa[1:2],jtbingeresb[1:2],jthazardresa[1:2],jthazardresb[1:2]),byrow=TRUE,ncol=4,nrow=4)
rownames(jtresults) <- c("Drinking 2+ per week","Drinking 4+ times per week","Binge drinking","hazard-C Problematic drinking")
colnames(jtresults) <- c("BPI Pain Score coef","BPI Pain Score SE","BPI Pain Int coef","BPI Pain Int SE")

# Cohen's d results
cdresults <- matrix(rep(NA,12),ncol=4,nrow=4)
cdresults[1,1] <- smd(Mean.1=jtfreq2resa[3],s.1=jtfreq2resa[4]*sqrt(1514),n.1=1514,Mean.2=jtfreq2resa[5],s.2=jtfreq2resa[6]*sqrt(1514),n.2=1514)
cdresults[1,2] <- sqrt(((3028/2292196)+(cdresults[1,1]^2/6052))*(3028/3026))
cdresults[1,3] <- smd(Mean.1=jtfreq2resb[3],s.1=jtfreq2resb[4]*sqrt(1514),n.1=1514,Mean.2=jtfreq2resb[5],s.2=jtfreq2resb[6]*sqrt(1514),n.2=1514)
cdresults[1,4] <- sqrt(((3028/2292196)+(cdresults[1,3]^2/6052))*(3028/3026))
cdresults[2,1] <- smd(Mean.1=jtfreq4resa[3],s.1=jtfreq4resa[4]*sqrt(1514),n.1=1514,Mean.2=jtfreq4resa[5],s.2=jtfreq4resa[6]*sqrt(1514),n.2=1514)
cdresults[2,2] <- sqrt(((3028/2292196)+(cdresults[1,1]^2/6052))*(3028/3026))
cdresults[2,3] <- smd(Mean.1=jtfreq4resb[3],s.1=jtfreq4resb[4]*sqrt(1514),n.1=1514,Mean.2=jtfreq4resb[5],s.2=jtfreq4resb[6]*sqrt(1514),n.2=1514)
cdresults[2,4] <- sqrt(((3028/2292196)+(cdresults[1,3]^2/6052))*(3028/3026))
cdresults[3,1] <- smd(Mean.1=jtbingeresa[3],s.1=jtbingeresa[4]*sqrt(1514),n.1=1514,Mean.2=jtbingeresa[5],s.2=jtbingeresa[6]*sqrt(1514),n.2=1514)
cdresults[3,2] <- sqrt(((3028/2292196)+(cdresults[2,1]^2/6052))*(3028/3026))
cdresults[3,3] <- smd(Mean.1=jtbingeresb[3],s.1=jtbingeresb[4]*sqrt(1514),n.1=1514,Mean.2=jtbingeresb[5],s.2=jtbingeresb[6]*sqrt(1514),n.2=1514)
cdresults[3,4] <- sqrt(((3028/2292196)+(cdresults[2,3]^2/6052))*(3028/3026))
cdresults[4,1] <- smd(Mean.1=jthazardresa[3],s.1=jthazardresa[4]*sqrt(1514),n.1=1514,Mean.2=jthazardresa[5],s.2=jthazardresa[6]*sqrt(1514),n.2=1514)
cdresults[4,2] <- sqrt(((3028/2292196)+(cdresults[3,1]^2/6052))*(3028/3026))
cdresults[4,3] <- smd(Mean.1=jthazardresb[3],s.1=jthazardresb[4]*sqrt(1514),n.1=1514,Mean.2=jthazardresb[5],s.2=jthazardresb[6]*sqrt(1514),n.2=1514)
cdresults[4,4] <- sqrt(((3028/2292196)+(cdresults[3,3]^2/6052))*(3028/3026))
rownames(cdresults) <- c("Drinking 2+ per week","Drinking 4+ times per week","Binge drinking","hazard-C Problematic drinking")
colnames(cdresults) <- c("BPI Pain Score d","BPI Pain Score d SE","BPI Pain Int d","BPI Pain Int d SE")

save(jtresults,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/jtresults.RData"))
save(cdresults,file=paste0(cloudstor,"PhD/Paper 6 - POINT application/Results/cdresults.RData"))
