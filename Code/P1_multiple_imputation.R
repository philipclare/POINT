
## SYNTAX FILE 1 - MULTIPLE IMPUTATION                       ##

cloudstor <- "C:/Users/z3312911/Cloudstor/" # change to master file path
.libPaths(paste0(cloudstor,"R Library"))

libs <- c("mitml","pan","haven","plyr","readr")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
library ("mitml")
library ("pan")
library ("haven")
library ("plyr")
library ("readr")

set.seed(368078)

# Load original dataset and drop daily drinking variable (not needed)
datalong_alc <- read_dta(file=paste0(cloudstor,"PhD/Paper 6 - POINT Application/Data/Point data long.dta"))
datalong_alc <- datalong_alc[,-28]

# Number of imputations to use
nimpute <- 50
burn <- 100
iter <- 10

# Define variable lists which determine variable order and varying variables for reshape to wide format
flist <- c("Participant_ID","time","b_Actual_age","b_employ3","b_sex","b_edu","b_MaritalStatus","b_Arth_12mR","b_Back_12mR","b_Head_12mR","b_Visc_12mR","b_Fibro_12mR","b_Cmplx_12mR","b_Shing_12mR","b_pain_duration_yrs3","b_alc_ever","PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","alc_12m","binge","auditcprob","BPI_PScore","BPI_interference")
varlist <- c("PHQ9_Mod_sev","GADMod2Sev","Antidepressant_week","Antipsychotic_week","benzo_week","Nonopioid_analgesic_week","Pregabalin_week","can_12m","cig_12m","opioid90","alc_pain_12m","alc_12m","binge","auditcprob","BPI_PScore","BPI_interference")

# Define factor variables to make the imputation work properly
datalong_alc$time <- factor(datalong_alc$time)
datalong_alc$b_employ3 <- factor(datalong_alc$b_employ3)
datalong_alc$b_sex <- factor(datalong_alc$b_sex)
datalong_alc$b_edu <- factor(datalong_alc$b_edu)
datalong_alc$b_MaritalStatus <- factor(datalong_alc$b_MaritalStatus)
datalong_alc$b_edu <- factor(datalong_alc$b_edu)
datalong_alc$b_Arth_12mR <- factor(datalong_alc$b_Arth_12mR)
datalong_alc$b_Back_12mR <- factor(datalong_alc$b_Back_12mR)
datalong_alc$b_Head_12mR <- factor(datalong_alc$b_Head_12mR)
datalong_alc$b_Visc_12mR <- factor(datalong_alc$b_Visc_12mR)
datalong_alc$b_Fibro_12mR <- factor(datalong_alc$b_Fibro_12mR)
datalong_alc$b_Cmplx_12mR <- factor(datalong_alc$b_Cmplx_12mR)
datalong_alc$b_Shing_12mR <- factor(datalong_alc$b_Shing_12mR)
datalong_alc$b_pain_duration_yrs3 <- factor(datalong_alc$b_pain_duration_yrs3)
datalong_alc$PHQ9_Mod_sev <- factor(datalong_alc$PHQ9_Mod_sev)
datalong_alc$GADMod2Sev <- factor(datalong_alc$GADMod2Sev)
datalong_alc$Antidepressant_week <- factor(datalong_alc$Antidepressant_week)
datalong_alc$Antipsychotic_week <- factor(datalong_alc$Antipsychotic_week)
datalong_alc$benzo_week <- factor(datalong_alc$benzo_week)
datalong_alc$Nonopioid_analgesic_week <- factor(datalong_alc$Nonopioid_analgesic_week)
datalong_alc$Pregabalin_week <- factor(datalong_alc$Pregabalin_week)
datalong_alc$can_12m <- factor(datalong_alc$can_12m)
datalong_alc$cig_12m <- factor(datalong_alc$cig_12m)
datalong_alc$alc_12m <- factor(datalong_alc$alc_12m)
datalong_alc$audit_1 <- factor(datalong_alc$audit_1)
datalong_alc$audit_1 <- factor(datalong_alc$audit_1)
datalong_alc$audit_2 <- factor(datalong_alc$audit_2)
datalong_alc$audit_3 <- factor(datalong_alc$audit_3)
datalong_alc$alc_pain_12m <- factor(datalong_alc$alc_pain_12m)
datalong_alc$b_alc_ever <- factor(datalong_alc$b_alc_ever)

# Define vector of variable types: -2=cluster var; 1=impute; 2=complete
type <- c(-2,2,2,1,2,2,2,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

# Run multiple imputation using Jomo via mitml and save as list of imputed data files
JMimpute_alc <- mitmlComplete(jomoImpute(datalong_alc,
                                         type=type,
                                         m=nimpute,
                                         n.burn=burn,
                                         n.iter=iter),print="all")

JMimpute_final <- lapply(JMimpute_alc, function (x) {
  x$alc_12m <- as.numeric(levels(x$alc_12m))[x$alc_12m]
  x$b_alc_ever <- as.numeric(levels(x$b_alc_ever))[x$b_alc_ever]
  x$alc_pain_12m <- as.numeric(levels(x$alc_pain_12m))[x$alc_pain_12m]
  x$binge <- ifelse(as.numeric(x$audit_3)==1,0,1)
  x$auditc <- as.numeric(x$audit_1) + 
    as.numeric(x$audit_2) + 
    as.numeric(x$audit_3)-3
  x$auditcprob <- ifelse(x$b_sex==1,
                         ifelse(x$auditc<4,0,1),
                         ifelse(x$auditc<3,0,1))
  x$opioid90 <- ifelse(x$totalopioiddose<90,0,1)
  x
})

# Reshape data to wide format for use by LTMLE
impdatawide <- lapply(JMimpute_final, function (x) {
  x <- as.data.frame(x[,flist])
  wide <- reshape(x,
                  v.names=varlist,
                  idvar="Participant_ID",
                  timevar="time",
                  sep="",
                  direction="wide")
})

# Save imputed datasets
save(JMimpute_final,file=paste0(cloudstor,"PhD/Paper 6 - POINT Application/Data/Imputed data - long.RData"))
save(impdatawide,file=paste0(cloudstor,"PhD/Paper 6 - POINT Application/Data/Imputed data - wide.RData"))