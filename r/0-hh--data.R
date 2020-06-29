### Symptom combinations for detecting SARS-CoV-2 infection in nonhospitalized persons
### Read in source CSV data file; set up names, analytic data
# source data, test results, symptoms, preconstructed rules

# set working directory

## Source data
hh_195 <- read.csv("public_records.csv", stringsAsFactors = TRUE) # 195 x 31
hh <- subset(hh_195, !(pcr_pos==0 & sero_pos==1)) # 185 x 31
# hh_oth <- subset(hh_195, (pcr_pos==0 & sero_pos==1)) # 10 x 31
# dput(names(hh))

# Test results
test_cov2 <- as.integer(hh$any_pos)
test_cov2_adult <- test_cov2[hh$age_adult == 1]
test_cov2_child <- test_cov2[hh$age_adult == 0]
npos_cov2 <- sum(test_cov2); nneg_cov2 <- sum(1-test_cov2)
npos_cov2_adult <- sum(test_cov2_adult); nneg_cov2_adult <- sum(1-test_cov2_adult)
npos_cov2_child <- sum(test_cov2_child); nneg_cov2_child <- sum(1-test_cov2_child)

## Symptoms
# 15 column names for binary symptom values
symps <- c("wheeze", "throat", "sob", "nausea", "myalgia", "headache", 
   "fatigue", "discomf", "diarrhea", "cough", "chestpain", "abdpain",
   "fever_chills", "nasal_combo", "tastesmell_combo")
symp_labels <- c("wheeze", "throat", "sob", "nausea", "myalgia", "headache", 
   "fatigue", "discomfort", "diarrhea", "cough", "chest_pain", "abd_pain", 
   "fever_chills", "nasal", "taste_smell")

# Numeric matrices: 0=TN, 1=FP (rule=yes, test=no), 2=FN (rule=no, test=yes), 3=TP
sy_data <- data.matrix(hh[, symps]) # 185 x 15 matrix
sy_data_adult <- sy_data[hh$age_adult == 1,] # 122 x 15 matrix
sy_data_child <- sy_data[hh$age_adult == 0,] # 63 x 15 matrix
hh_cov2_symps <- 2*(hh$any_pos) + hh[, symps] # 195 x 15 data.frame
names(hh_cov2_symps) <- paste(symp_labels, "cov2", sep="_")

# symptom groupings for graphics; match enables indexing symps or symp_labels
cons_idx <- match(c("myalgia", "fatigue", "fever_chills"), symps)
upresp_idx <- match(c("throat", "nasal_combo"), symps)
loresp_idx <- match(c("wheeze", "sob", "discomf", "cough",
   "chestpain"), symps)
neuro_idx <- match(c("headache", "tastesmell_combo"), symps)
gastro_idx <- match(c("nausea", "diarrhea", "abdpain"), symps)

## Preconstructed rules
rules <- c("ili", "cdc", "ari", "cste", "cli", 
   "vaccine_a1", "vaccine_a2", "vaccine_a3", "vaccine_a_all")
rule_labels <- c("ILI", "CDC", "ARI", "CSTE", "CLI", 
   "Vaccine_A1", "Vaccine_A2", "Vaccine_A3", "Vaccine_A_all")

# 0=TN, 1=FP (rule=yes, test=no), 2=FN (rule=no, test=yes), 3=TP
rl_data <- data.matrix(hh[, rules]) # 185 x 9 matrix
hh_cov2_rules <- 2*(hh$any_pos) + hh[, rules] # 185 x 9 data.frame
names(hh_cov2_rules) <- paste(rule_labels, "cov2", sep="_")

# rule groupings for graphics
orule_idx <- match(c("ili", "cdc", "ari", "cste", "cli"), rules)
vaxa_idx <- match(c("vaccine_a1", "vaccine_a2", "vaccine_a3", "vaccine_a_all"), rules)

## Rule and symptom labels combined
rusy <- c(rules, symps)
rusy_labels <- c(rule_labels, symp_labels)
