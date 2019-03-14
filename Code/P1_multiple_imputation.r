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

datalong_alc <- read_dta(file=paste0(cloudstor,"PhD/Paper 6 - POINT Application/Data/Point data long.dta"))
datalong_alc <- datalong_alc[,-27]

# Number of imputations to use
nimpute <- 50
burn <- 100
iter <- 10

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

# Define vector of variable types: -2=cluster var; 1=impute; 2=complete
type <- c(-2,2,2,1,2,2,2,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

# Run multiple imputation using Jomo via mitml and save as list of imputed data files
JMimpute_alc <- mitmlComplete(jomoImpute(datalong_alc,
                                         type=type,
                                         m=nimpute,
                                         n.burn=burn,
                                         n.iter=iter),print="all")

# For each imputed dataset, create the final exposure variables to be used in analysis
for (l in 1:nimpute){
  JMimpute_alc[[l]]$alc_12m <- as.numeric(levels(JMimpute_alc[[l]]$alc_12m))[JMimpute_alc[[l]]$alc_12m]
  JMimpute_alc[[l]]$binge <- ifelse(as.numeric(JMimpute_alc[[l]]$audit_1)==1,0,1)
  JMimpute_alc[[l]]$auditc <- as.numeric(JMimpute_alc[[l]]$audit_1) + 
    as.numeric(JMimpute_alc[[l]]$audit_2) + 
    as.numeric(JMimpute_alc[[l]]$audit_3)-3
  JMimpute_alc[[l]]$auditcprob <- ifelse(JMimpute_alc[[l]]$b_sex==1,
                                         ifelse(JMimpute_alc[[l]]$auditc<4,0,1),
                                         ifelse(JMimpute_alc[[l]]$auditc<3,0,1))
  JMimpute_alc[[l]]$opioid90 <- ifelse(JMimpute_alc[[l]]$totalopioiddose<90,0,1)
}

# Save imputed dataset
save(JMimpute_alc,file=paste0(cloudstor,"PhD/Paper 6 - POINT Application/Data/Imputed data.RData"))
