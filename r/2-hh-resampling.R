### Symptom combinations for detecting SARS-CoV-2 infection in nonhospitalized persons
### Construct pseudosamples by resampling households

# resample households with replacement, 36 in State A, 25 in State B
# expand each pseudosample of households to pseudosample of cohabitants
# replicate original data into pseudosamples

# link cohabitant ID integer to household ID integer
hh_to_co <- with(hh, data.frame(hh_int = as.integer(factor(hh_id)), 
   co_int = as.integer(factor(study_id_merge))))
hh_rle <- rle(hh_to_co$hh_int)
co_by_hh <- with(hh_to_co, split(co_int, hh_int))

set.seed(24601)
nsamp <- 1e4
hh_idx_mat <- rbind(
   matrix(sample(1:36, 36*nsamp, replace=TRUE), nrow=36),
   matrix(sample(37:61, 25*nsamp, replace=TRUE), nrow=25))
hh_idx_mat <- apply(hh_idx_mat, 2, sort) # not necessary; makes later steps easier
# 61 x nsamp

# enumerate pseudosamples (psamp_num) and pseudohouseholds within them (psamp_id)
hh_idx_frame <- data.frame(psamp_num = c(col(hh_idx_mat)),
   psamp_hh_id = c(row(hh_idx_mat)),
   psamp_co_id = NA, # placeholder
   hh_idx = c(hh_idx_mat),
   hh_freq = hh_rle$lengths[c(hh_idx_mat)])
# 61*nsamp x 5
# table(hh_idx_frame$hh_idx)

# expand the pseudohousehold frame to a pseudo-individual (cohabiter) frame
#    repeat each household the number of times in the original frame
#    plug the list of household cohabitants into each occurrence of that household
co_idx_frame <- cbind(hh_idx_frame[rep(seq_len(61*nsamp), hh_idx_frame$hh_freq),], 
   co_id=unlist(co_by_hh[hh_idx_frame$hh_idx]))
# 1849909 x 6 (approx 185*nsamp)

# finally, expand the original data into pseudosamples
hh_psamp <- cbind(co_idx_frame,
   cbind(hh, hh_cov2_rules, hh_cov2_symps)[co_idx_frame$co_id,])
hh_psamp$psamp_co_id <- unlist(
   lapply(tabulate(hh_psamp$psamp_num, nbins=nsamp), seq_len))
# 1849909 x 61
# write.csv(hh_psamp, "hh_psamp_1e4.csv")
# write.csv(subset(hh_psamp, psamp_num < 401), "hh_psamp_400.csv")

# TN, FP, FN, TP for each rule or symptom ("rusy") in each pseudosample
# list of nsamp matrices, each 4 x 24
hh_psamp_rusy_adult <- lapply(
   split(
      x = subset(x=hh_psamp, subset=(age_adult==1), 
         select=c(names(hh_cov2_rules), names(hh_cov2_symps))),
      f = hh_psamp$psamp_num[hh_psamp$age_adult==1]),
   function(x) apply(x+1, 2, tabulate, nbins=4))

# Sensitivity and specificity for each rule in each pseudosample
# list of nsamp matrices, each 2 x 24
hh_psamp_rusy_adult_sesp <- lapply(hh_psamp_rusy_adult,
   function(x) rbind(sens = x[4,]/(x[4,]+x[3,]), spec = x[1,]/(x[1,]+x[2,]))
)

# TN, FP, FN, TP for each rule in each pseudosample
# list of nsamp matrices, each 4 x 24
hh_psamp_rusy_child <- lapply(
   split(
      x = subset(x=hh_psamp, subset=(age_adult==0), 
         select=c(names(hh_cov2_rules), names(hh_cov2_symps))),
      f = hh_psamp$psamp_num[hh_psamp$age_adult==0]),
   function(x) apply(x+1, 2, tabulate, nbins=4))

# Sensitivity and specificity for each rule in each pseudosample
# list of nsamp matrices, each 2 x 24
hh_psamp_rusy_child_sesp <- lapply(hh_psamp_rusy_child,
   function(x) rbind(sens = x[4,]/(x[4,]+x[3,]), spec = x[1,]/(x[1,]+x[2,]))
)

# rearrange lists of matrices as 4-way arrays
hh_psamp_rusy_adult_sesp_array <- hh_psamp_rusy_child_sesp_array <- 
   array(NA, dim=c(2, 24, nsamp), dimnames = list(c("sens", "spec"), 
      c(rule_labels, symp_labels), seq_len(nsamp)))
hh_psamp_rusy_adult_sesp_array[] <- unlist(hh_psamp_rusy_adult_sesp)
hh_psamp_rusy_child_sesp_array[] <- unlist(hh_psamp_rusy_child_sesp)
hh_psamp_rusy_adult_sesp_array <- aperm(hh_psamp_rusy_adult_sesp_array, c(3, 2, 1))
hh_psamp_rusy_child_sesp_array <- aperm(hh_psamp_rusy_child_sesp_array, c(3, 2, 1))
# hh_psamp_rusy_adult_sesp_array[m, , ] == hh_psamp_rusy_adult_sesp[[m]]
hh_psamp_rusy_adult_sesp_array[1, , 1] # sensitivity re: PCR result in pseudosample 1
hh_psamp_rusy_adult_sesp_array[1, , 2] # specificity re: PCR result in pseudosample 1
# adult - child difference; need to reverse sign to graph specificity
hh_psamp_rusy_diff_sesp_array <- 
   hh_psamp_rusy_adult_sesp_array - hh_psamp_rusy_child_sesp_array

# bootstrap mean sensitivity: adult, child, difference
sens_cov2_adult_bsmean <- colMeans(hh_psamp_rusy_adult_sesp_array[, , 1])
sens_cov2_child_bsmean <- colMeans(hh_psamp_rusy_child_sesp_array[, , 1])
sens_cov2_diff_bsmean <- colMeans(hh_psamp_rusy_diff_sesp_array[, , 1])
# sens_cov2_diff_bsmean <- sens_cov2_adult_bsmean - sens_cov2_child_bsmean
# bootstrap mean specificity: adult, child, difference
spec_cov2_adult_bsmean <- colMeans(hh_psamp_rusy_adult_sesp_array[, , 2])
spec_cov2_child_bsmean <- colMeans(hh_psamp_rusy_child_sesp_array[, , 2])
spec_cov2_diff_bsmean <- colMeans(hh_psamp_rusy_diff_sesp_array[, , 2])
# spec_cov2_diff_bsmean <- spec_cov2_adult_bsmean - spec_cov2_child_bsmean

# bootstrap covariance of sensitivity and specificity

# 2 x 2 (sens by spec) bootstrap covariance for each rule
sesp_cov2_adult_bscov <- lapply(rusy_labels, function(x) 
   cov(cbind(hh_psamp_rusy_adult_sesp_array[, x, 1], 
      hh_psamp_rusy_adult_sesp_array[, x, 2])))
sesp_cov2_child_bscov <- lapply(rusy_labels, function(x) 
   cov(cbind(hh_psamp_rusy_child_sesp_array[, x, 1], 
      hh_psamp_rusy_child_sesp_array[, x, 2])))
sesp_cov2_diff_bscov <- lapply(rusy_labels, function(x) 
   cov(cbind(hh_psamp_rusy_diff_sesp_array[, x, 1], 
      hh_psamp_rusy_diff_sesp_array[, x, 2])))
names(sesp_cov2_adult_bscov) <- names(sesp_cov2_child_bscov) <- 
   names(sesp_cov2_diff_bscov) <- rusy_labels

