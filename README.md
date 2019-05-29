# The effect of prolonged alcohol use on chronic non-cancer pain in users of prescription opioids
## R Analysis Code

This repository contains R code used in the TMLE analysis of POINT data on the effect of alcohol exposure on chronic pain, by Clare et al. 2019

| Description | R-code |
| --- | --- |
| P1 - Multiple imputation | [Multiple imputation code](Code/P1_multiple_imputation.R) |
| P2 - Longitudinal sample descriptives (article Table 2) | [Longitudinal descriptives code](Code/P2_descriptives.R) |
| P3 - LTMLE primary analysis of alcohol consumption on pain using the package 'ltmle' (1). | [LTMLE analysis code](Code/P3_ltmle_primary.R) |
| P4 - LTMLE sensitivity analysis excluding those who used alcohol to cope with their pain | [LTMLE sensitivity analysis code](Code/P4_ltmle_sensitivity.R) |
| P5 - Naive analysis using random intercept models | [Naive analysis code](Code/P5_naive_analysis.R) |
| P6 - E-Value sensitivity analysis | [E-value analysis code](Code/P6_evalue_analysis.R) |

1. Lendle SD, Schwab J, Petersen ML, van der Laan MJ. ltmle: An R Package Implementing Targeted Minimum Loss-Based Estimation for Longitudinal Data. Journal of Statistical Software. 2017;81(1):1-21.
