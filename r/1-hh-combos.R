### Symptom combinations for detecting SARS-CoV-2 infection in nonhospitalized persons
### Calculate performance of preconstructed rules and of symptom combinations

# TN, FP, FN, and TP for rules; all 2 x 9 matrices
test_cov2_pos.rules <- apply(rl_data[test_cov2==1, rules]+1, 2, tabulate, nbins=2)/npos_cov2
test_cov2_neg.rules <- apply(rl_data[test_cov2==0, rules]+1, 2, tabulate, nbins=2)/nneg_cov2
test_cov2_adult_pos.rules <- apply(rl_data[test_cov2==1 & hh$age_adult==1, rules]+1, 
   2, tabulate, nbins=2)/npos_cov2_adult
test_cov2_adult_neg.rules <- apply(rl_data[test_cov2==0 & hh$age_adult==1, rules]+1, 
   2, tabulate, nbins=2)/nneg_cov2_adult
test_cov2_child_pos.rules <- apply(rl_data[test_cov2==1 & hh$age_adult==0, rules]+1, 
   2, tabulate, nbins=2)/npos_cov2_child
test_cov2_child_neg.rules <- apply(rl_data[test_cov2==0 & hh$age_adult==0, rules]+1, 
   2, tabulate, nbins=2)/nneg_cov2_child

# Performance measures for rules
sens_cov2.rl <- test_cov2_pos.rules[2,]; spec_cov2.rl <- test_cov2_neg.rules[1,]
sens_cov2_adult.rl <- test_cov2_adult_pos.rules[2,]
spec_cov2_adult.rl <- test_cov2_adult_neg.rules[1,]
sens_cov2_child.rl <- test_cov2_child_pos.rules[2,]
spec_cov2_child.rl <- test_cov2_child_neg.rules[1,]

# Utility function: compute # indicants in matrix sym for all combinations of size n
combn.tab <- function(n, sym) {
   combn(x = dimnames(sym)[[2]], m = n, 
      FUN = function(x) tabulate(rowSums(sym[, x], na.rm=TRUE)+1, nbins=n+1))
}

# Problem sizes
# choose(15, 1:15)
# c(15, 105, 455, 1365, 3003, 5005, 6435, 6435, 5005, 3003, 1365, 455, 105, 15, 1)

# Objects to store compact names, tables of # indicants for each combination
names.sx.list <- vector(mode="list", 15)
names.sx.list[[1]] <- LETTERS[1:15]
# initialize lists
test_cov2_neg.sx.list <- test_cov2_pos.sx.list <- 
   test_cov2_adult_neg.sx.list <- test_cov2_adult_pos.sx.list <- 
   test_cov2_child_neg.sx.list <- test_cov2_child_pos.sx.list <- vector(mode="list", 15)

test_cov2_pos.sx.list[[1]] <- apply(sy_data[test_cov2==1, symps]+1, 
   2, tabulate, nbins=2)/npos_cov2
test_cov2_neg.sx.list[[1]] <- apply(sy_data[test_cov2==0, symps]+1, 
   2, tabulate, nbins=2)/nneg_cov2
test_cov2_adult_pos.sx.list[[1]] <- apply(sy_data[test_cov2==1 & hh$age_adult==1, symps]+1, 
   2, tabulate, nbins=2)/npos_cov2_adult
test_cov2_adult_neg.sx.list[[1]] <- apply(sy_data[test_cov2==0 & hh$age_adult==1, symps]+1, 
   2, tabulate, nbins=2)/nneg_cov2_adult
test_cov2_child_pos.sx.list[[1]] <- apply(sy_data[test_cov2==1 & hh$age_adult==0, symps]+1, 
   2, tabulate, nbins=2)/npos_cov2_child
test_cov2_child_neg.sx.list[[1]] <- apply(sy_data[test_cov2==0 & hh$age_adult==0, symps]+1, 
   2, tabulate, nbins=2)/nneg_cov2_child

print(system.time({
for(.i. in 2:15) {
   print(.i.)
   names.sx.list[[.i.]] <- apply(combn(LETTERS[1:15], .i.), 2, paste, collapse="")
   test_cov2_pos.sx.list[[.i.]] <- combn.tab(.i., sy_data[test_cov2==1, symps])/npos_cov2
   test_cov2_neg.sx.list[[.i.]] <- combn.tab(.i., sy_data[test_cov2==0, symps])/nneg_cov2
   test_cov2_adult_pos.sx.list[[.i.]] <- 
      combn.tab(.i., sy_data[test_cov2==1 & hh$age_adult==1, symps])/npos_cov2_adult
   test_cov2_adult_neg.sx.list[[.i.]] <- 
      combn.tab(.i., sy_data[test_cov2==0 & hh$age_adult==1, symps])/nneg_cov2_adult
   test_cov2_child_pos.sx.list[[.i.]] <- 
      combn.tab(.i., sy_data[test_cov2==1 & hh$age_adult==0, symps])/npos_cov2_child
   test_cov2_child_neg.sx.list[[.i.]] <- 
      combn.tab(.i., sy_data[test_cov2==0 & hh$age_adult==0, symps])/nneg_cov2_child
}
})) # 3 seconds

# Sensitivity and specificity for each combination at each cutoff
sens_cov2.sx.list <- lapply(test_cov2_pos.sx.list[1:14], function(x)
   apply(x[nrow(x):1,], 2, cumsum)[(nrow(x)-1):1, , drop=FALSE])
spec_cov2.sx.list <- lapply(test_cov2_neg.sx.list[1:14], function(x) 
   apply(x, 2, cumsum)[-nrow(x), , drop=FALSE])

sens_cov2_adult.sx.list <- lapply(test_cov2_adult_pos.sx.list[1:14], function(x)
   apply(x[nrow(x):1,], 2, cumsum)[(nrow(x)-1):1, , drop=FALSE])
spec_cov2_adult.sx.list <- lapply(test_cov2_adult_neg.sx.list[1:14], function(x) 
   apply(x, 2, cumsum)[-nrow(x), , drop=FALSE])

sens_cov2_child.sx.list <- lapply(test_cov2_child_pos.sx.list[1:14], function(x)
   apply(x[nrow(x):1,], 2, cumsum)[(nrow(x)-1):1, , drop=FALSE])
spec_cov2_child.sx.list <- lapply(test_cov2_child_neg.sx.list[1:14], function(x) 
   apply(x, 2, cumsum)[-nrow(x), , drop=FALSE])

# Likelihood ratios positive and negative, diagnostic odds ratios
lrp_cov2.sx.list <- lrn_cov2.sx.list <- dor_cov2.sx.list <- 
   lrp_cov2_adult.sx.list <- lrn_cov2_adult.sx.list <- dor_cov2_adult.sx.list <- 
   lrp_cov2_child.sx.list <- lrn_cov2_child.sx.list <- dor_cov2_child.sx.list <- 
   vector(mode="list", 14)
for(.i. in 1:14) {
   lrp_cov2.sx.list[[.i.]] <- sens_cov2.sx.list[[.i.]]/(1-spec_cov2.sx.list[[.i.]])
   lrn_cov2.sx.list[[.i.]] <- (1-sens_cov2.sx.list[[.i.]])/spec_cov2.sx.list[[.i.]]
   dor_cov2.sx.list[[.i.]] <- lrp_cov2.sx.list[[.i.]]/lrn_cov2.sx.list[[.i.]]
   lrp_cov2_adult.sx.list[[.i.]] <- 
      sens_cov2_adult.sx.list[[.i.]]/(1-spec_cov2_adult.sx.list[[.i.]])
   lrn_cov2_adult.sx.list[[.i.]] <- 
      (1-sens_cov2_adult.sx.list[[.i.]])/spec_cov2_adult.sx.list[[.i.]]
   dor_cov2_adult.sx.list[[.i.]] <- 
      lrp_cov2_adult.sx.list[[.i.]]/lrn_cov2_adult.sx.list[[.i.]]
   lrp_cov2_child.sx.list[[.i.]] <- 
      sens_cov2_child.sx.list[[.i.]]/(1-spec_cov2_child.sx.list[[.i.]])
   lrn_cov2_child.sx.list[[.i.]] <- 
      (1-sens_cov2_child.sx.list[[.i.]])/spec_cov2_child.sx.list[[.i.]]
   dor_cov2_child.sx.list[[.i.]] <- 
      lrp_cov2_child.sx.list[[.i.]]/lrn_cov2_child.sx.list[[.i.]]
}

