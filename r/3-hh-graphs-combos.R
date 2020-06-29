### Symptom combinations for detecting SARS-CoV-2 infection in nonhospitalized persons
### Graphics: how symptom combinations perform

library("viridis")
pal <- viridis(6)

# Calculate convex hulls for all combos of 15 symptoms
hull_cov2.sx.list <-  hull_cov2_adult.sx.list <-  hull_cov2_child.sx.list <- 
   vector(mode="list", 105)
hull_names <- character(105)
k <- 1
print(system.time({
for(i in 1:14) {
   for(j in 1:i) {
      # print(c(j,i,k))
      hull_names[k] <- paste(j,i,sep="/")
      hull_cov2.sx.list[[k]] <- 
         chull(sens_cov2.sx.list[[i]][j,], spec_cov2.sx.list[[i]][j,])
      hull_cov2_adult.sx.list[[k]] <- 
         chull(sens_cov2_adult.sx.list[[i]][j,], spec_cov2_adult.sx.list[[i]][j,])
      hull_cov2_child.sx.list[[k]] <- 
         chull(sens_cov2_child.sx.list[[i]][j,], spec_cov2_child.sx.list[[i]][j,])
      k <- k+1
   }
}
})) # < 1 sec
rm(i,j,k)

# Standard graphic set-up
setup.sens.spec.plot <- function(base.algs=0, sens, spec) {
   plot(c(0, 1), c(0, 1), type="n", bty="n", xaxt="n", yaxt="n", 
      xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
   axis(1, at=(0:5)/5, labels=paste0((0:5)*20, "%"))
   axis(2, at=(0:5)/5, labels=paste0((0:5)*20, "%"), las=1)
   segments(x0=0, x1=1, y0=0, y1=1, col="gray", lty="dotted", lwd=2)
   if(base.algs) {
      points(1-spec, sens, col=c("black", "gray")[base.algs])
      text(1-spec, sens - 0.02, rule_labels[orule_idx], 
         adj=c(0.5, 1), col=c("black", "gray")[base.algs], 
         cex=c(1,0.8)[base.algs])
   }
   invisible()
}
#setup.sens.spec.plot(base.algs=0)
setup.sens.spec.plot(base.algs=1, sens=sens_cov2.rl[orule_idx], spec=spec_cov2.rl[orule_idx])
setup.sens.spec.plot(base.algs=2, sens=sens_cov2.rl[orule_idx], spec=spec_cov2.rl[orule_idx])

pdf("hh-combos_cov2.pdf")
# 1. Baseline rules
setup.sens.spec.plot(base.algs=1, sens=sens_cov2.rl[orule_idx], spec=spec_cov2.rl[orule_idx])
title(main="Performance of existing rules, PCR, all ages", adj=0)

# 2. Cloud of single symptoms
setup.sens.spec.plot(base.algs=2, sens=sens_cov2.rl[orule_idx], spec=spec_cov2.rl[orule_idx])
title(main="All single-symptom rules, PCR, all ages", adj=0)
polygon(1-spec_cov2.sx.list[[1]][1, hull_cov2.sx.list[[1]]],
   sens_cov2.sx.list[[1]][1, hull_cov2.sx.list[[1]]],
   lwd=1, border="gray")
text(1-spec_cov2.sx.list[[1]][1,], sens_cov2.sx.list[[1]][1,], 
   names.sx.list[[1]], cex=0.95, col=pal[2])
legend(x="bottomright", legend=c(symp_labels), title="Single symptoms",
   pch=LETTERS[1:15], cex=0.85,
   col="black", bty="n")

# 3. Cloud for 1+/2 + convex hull
setup.sens.spec.plot(base.algs=2, sens=sens_cov2.rl[orule_idx], spec=spec_cov2.rl[orule_idx])
title(main="All combinations, either of 2 symptoms, PCR, all ages", adj=0)
polygon(1-spec_cov2.sx.list[[2]][1, hull_cov2.sx.list[[2]]],
   sens_cov2.sx.list[[2]][1, hull_cov2.sx.list[[2]]], 
   lwd=1, border="gray")
text(1-spec_cov2.sx.list[[2]][1,], sens_cov2.sx.list[[2]][1,], 
   names.sx.list[[2]], cex=0.5, col=pal[2], srt=45)
legend(x="bottomright", legend=c(symp_labels), title="Single symptoms",
   pch=LETTERS[1:15], cex=0.85,
   col="black", bty="n")

# 4. Convex hulls for 1+/2, 2/2 + baseline rules
setup.sens.spec.plot(base.algs=2, sens=sens_cov2.rl[orule_idx], spec=spec_cov2.rl[orule_idx])
title(main="All combinations of 2 symptoms, PCR, all ages", adj=0)
# lines(1-rf_all_spec_bsmean, rf_all_sens_bsmean, col="gray")
polygon(1-spec_cov2.sx.list[[2]][1, hull_cov2.sx.list[[2]]], 
   sens_cov2.sx.list[[2]][1, hull_cov2.sx.list[[2]]], border=pal[1])
polygon(1-spec_cov2.sx.list[[2]][2, hull_cov2.sx.list[[3]]], 
   sens_cov2.sx.list[[2]][2, hull_cov2.sx.list[[3]]], lty=2, border=pal[3])
legend(x="bottomright", hull_names[2:3], lty=1:2, col=pal[c(1,3)], bty="n")

# 5. Convex hulls for 1+/3, 2+/3, 3/3 + baseline rules
setup.sens.spec.plot(base.algs=2, sens=sens_cov2.rl[orule_idx], spec=spec_cov2.rl[orule_idx])
title(main="All combinations of 3 symptoms, PCR, all ages", adj=0)
# lines(1-rf_all_spec_bsmean, rf_all_sens_bsmean, col="gray")
polygon(1-spec_cov2.sx.list[[3]][1, hull_cov2.sx.list[[4]]], 
   sens_cov2.sx.list[[3]][1, hull_cov2.sx.list[[4]]], border=pal[1])
polygon(1-spec_cov2.sx.list[[3]][2, hull_cov2.sx.list[[5]]], 
   sens_cov2.sx.list[[3]][2, hull_cov2.sx.list[[5]]], lty=2, border=pal[3])
polygon(1-spec_cov2.sx.list[[3]][3, hull_cov2.sx.list[[6]]], 
   sens_cov2.sx.list[[3]][3, hull_cov2.sx.list[[6]]], lty=3, lwd=2, border=pal[5])
legend(x="bottomright", hull_names[4:6], lty=1:3, lwd=2, col=pal[c(1, 3, 5)], bty="n")

# 6. Convex hulls for 1+/4, 2+/4, 3+/4, 4/4
setup.sens.spec.plot(base.algs=2, sens=sens_cov2.rl[orule_idx], 
   spec=spec_cov2.rl[orule_idx])
title(main="All combinations of 4 symptoms, PCR, all ages", adj=0)
# lines(1-rf_all_spec_bsmean, rf_all_sens_bsmean, col="gray")
polygon(1-spec_cov2.sx.list[[4]][1, hull_cov2.sx.list[[7]]], 
   sens_cov2.sx.list[[4]][1, hull_cov2.sx.list[[7]]], border=pal[1])
polygon(1-spec_cov2.sx.list[[4]][2, hull_cov2.sx.list[[8]]], 
   sens_cov2.sx.list[[4]][2, hull_cov2.sx.list[[8]]], lty=2, border=pal[2])
polygon(1-spec_cov2.sx.list[[4]][3, hull_cov2.sx.list[[9]]], 
   sens_cov2.sx.list[[4]][3, hull_cov2.sx.list[[9]]], lty=3, border=pal[3], lwd=2)
polygon(1-spec_cov2.sx.list[[4]][4, hull_cov2.sx.list[[10]]], 
   sens_cov2.sx.list[[4]][4, hull_cov2.sx.list[[10]]], lty=4, border=pal[4])
legend(x="bottomright", hull_names[7:10], lty=1:4, col=pal[1:4], lwd=2, bty="n")

# 7. Convex hulls for 1+/5, 2+/5, 3+/5, 4+/5, 5/5
setup.sens.spec.plot(base.algs=2, sens=sens_cov2.rl[orule_idx], 
   spec=spec_cov2.rl[orule_idx])
title(main="All combinations of 5 symptoms, PCR, all ages", adj=0)
# lines(1-rf_all_spec_bsmean, rf_all_sens_bsmean, col="gray")
polygon(1-spec_cov2.sx.list[[5]][1, hull_cov2.sx.list[[11]]], 
   sens_cov2.sx.list[[5]][1, hull_cov2.sx.list[[11]]], border=pal[1])
polygon(1-spec_cov2.sx.list[[5]][2, hull_cov2.sx.list[[12]]], 
   sens_cov2.sx.list[[5]][2, hull_cov2.sx.list[[12]]], lty=2, border=pal[2])
polygon(1-spec_cov2.sx.list[[5]][3, hull_cov2.sx.list[[13]]], 
   sens_cov2.sx.list[[5]][3, hull_cov2.sx.list[[13]]], lty=3, border=pal[3], lwd=2)
polygon(1-spec_cov2.sx.list[[5]][4, hull_cov2.sx.list[[14]]], 
   sens_cov2.sx.list[[5]][4, hull_cov2.sx.list[[14]]], lty=4, border=pal[4])
polygon(1-spec_cov2.sx.list[[5]][5, hull_cov2.sx.list[[15]]], 
   sens_cov2.sx.list[[5]][5, hull_cov2.sx.list[[15]]], lty=5, border=pal[5])
legend(x="bottomright", hull_names[11:15], lty=1:5, col=pal[1:5], lwd=2, bty="n")
dev.off()

pdf("hh-combos_cov2_adult.pdf")
# 1. Baseline rules
setup.sens.spec.plot(base.algs=1, sens=sens_cov2_adult.rl[orule_idx], 
   spec=spec_cov2_adult.rl[orule_idx])
title(main="Performance of existing rules, PCR, adults", adj=0)

# 2. Cloud of single symptoms
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_adult.rl[orule_idx], 
   spec=spec_cov2_adult.rl[orule_idx])
title(main="All single-symptom rules, PCR, adults", adj=0)
polygon(1-spec_cov2_adult.sx.list[[1]][1, hull_cov2_adult.sx.list[[1]]], 
   sens_cov2_adult.sx.list[[1]][1, hull_cov2_adult.sx.list[[1]]],
   lwd=1, border="gray")
text(1-spec_cov2_adult.sx.list[[1]][1,], sens_cov2_adult.sx.list[[1]][1,], 
   names.sx.list[[1]], cex=0.95, col=pal[2])
legend(x="bottomright", legend=c(symp_labels), title="Single symptoms",
   pch=LETTERS[1:15], cex=0.85,
   col="black", bty="n")

# 3. Cloud for 1+/2 + convex hull
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_adult.rl[orule_idx], 
   spec=spec_cov2_adult.rl[orule_idx])
title(main="All combinations, either of 2 symptoms, PCR, adults", adj=0)
polygon(1-spec_cov2_adult.sx.list[[2]][1, hull_cov2_adult.sx.list[[2]]], 
   sens_cov2_adult.sx.list[[2]][1, hull_cov2_adult.sx.list[[2]]],
   lwd=1, border="gray")
text(1-spec_cov2_adult.sx.list[[2]][1,], sens_cov2_adult.sx.list[[2]][1,], 
   names.sx.list[[2]], cex=0.5, col=pal[1], srt=45)
legend(x="bottomright", legend=c(symp_labels), title="Single symptoms",
   pch=LETTERS[1:15], cex=0.85,
   col="black", bty="n")

# 4. Convex hulls for 1+/2, 2/2 + baseline rules
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_adult.rl[orule_idx], 
   spec=spec_cov2_adult.rl[orule_idx])
title(main="All combinations of 2 symptoms, PCR, adults", adj=0)
polygon(1-spec_cov2_adult.sx.list[[2]][1, hull_cov2_adult.sx.list[[2]]], 
   sens_cov2_adult.sx.list[[2]][1, hull_cov2_adult.sx.list[[2]]], 
   border=pal[1])
polygon(1-spec_cov2_adult.sx.list[[2]][2, hull_cov2_adult.sx.list[[3]]], 
   sens_cov2_adult.sx.list[[2]][2, hull_cov2_adult.sx.list[[3]]], 
   lty=2, border=pal[3])
legend(x="bottomright", hull_names[2:3], lty=1:2, col=pal[c(1,3)], lwd=2, bty="n")

# 5. Convex hulls for 1+/3, 2+/3, 3/3 + baseline rules
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_adult.rl[orule_idx], 
   spec=spec_cov2_adult.rl[orule_idx])
title(main="All combinations of 3 symptoms, PCR, adults", adj=0)
polygon(1-spec_cov2_adult.sx.list[[3]][1, hull_cov2_adult.sx.list[[4]]], 
   sens_cov2_adult.sx.list[[3]][1, hull_cov2_adult.sx.list[[4]]],
   border=pal[1])
polygon(1-spec_cov2_adult.sx.list[[3]][2, hull_cov2_adult.sx.list[[5]]], 
   sens_cov2_adult.sx.list[[3]][2, hull_cov2_adult.sx.list[[5]]], 
   lty=2, border=pal[3])
polygon(1-spec_cov2_adult.sx.list[[3]][3, hull_cov2_adult.sx.list[[6]]], 
   sens_cov2_adult.sx.list[[3]][3, hull_cov2_adult.sx.list[[6]]], 
   lty=3, lwd=2, border=pal[5])
legend(x="bottomright", hull_names[4:6], lty=1:3, col=pal[c(1,3,5)], lwd=2, bty="n")

# 6. Convex hulls for 1+/4, 2+/4, 3+/4, 4/4
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_adult.rl[orule_idx], 
   spec=spec_cov2_adult.rl[orule_idx])
title(main="All combinations of 4 symptoms, PCR, adults", adj=0)
polygon(1-spec_cov2_adult.sx.list[[4]][1, hull_cov2_adult.sx.list[[7]]], 
   sens_cov2_adult.sx.list[[4]][1, hull_cov2_adult.sx.list[[7]]],
   border=pal[1])
polygon(1-spec_cov2_adult.sx.list[[4]][2, hull_cov2_adult.sx.list[[8]]], 
   sens_cov2_adult.sx.list[[4]][2, hull_cov2_adult.sx.list[[8]]], 
   lty=2, border=pal[2])
polygon(1-spec_cov2_adult.sx.list[[4]][3, hull_cov2_adult.sx.list[[9]]], 
   sens_cov2_adult.sx.list[[4]][3, hull_cov2_adult.sx.list[[9]]], 
   lty=3, lwd=2, border=pal[3])
polygon(1-spec_cov2_adult.sx.list[[4]][4, hull_cov2_adult.sx.list[[10]]], 
   sens_cov2_adult.sx.list[[4]][4, hull_cov2_adult.sx.list[[10]]], 
   lty=4, border=pal[4])
legend(x="bottomright", hull_names[7:10], lty=1:4, col=pal[1:4], lwd=2, bty="n")

# 7. Convex hulls for 1+/5, 2+/5, 3+/5, 4+/5, 5/5
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_adult.rl[orule_idx], 
   spec=spec_cov2_adult.rl[orule_idx])
title(main="All combinations of 5 symptoms, PCR, adults", adj=0)
polygon(1-spec_cov2_adult.sx.list[[5]][1, hull_cov2_adult.sx.list[[11]]], 
   sens_cov2_adult.sx.list[[5]][1, hull_cov2_adult.sx.list[[11]]],
   border=pal[1])
polygon(1-spec_cov2_adult.sx.list[[5]][2, hull_cov2_adult.sx.list[[12]]], 
   sens_cov2_adult.sx.list[[5]][2, hull_cov2_adult.sx.list[[12]]], 
   lty=2, border=pal[2])
polygon(1-spec_cov2_adult.sx.list[[5]][3, hull_cov2_adult.sx.list[[13]]], 
   sens_cov2_adult.sx.list[[5]][3, hull_cov2_adult.sx.list[[13]]], 
   lty=3, lwd=2, border=pal[3])
polygon(1-spec_cov2_adult.sx.list[[5]][4, hull_cov2_adult.sx.list[[14]]], 
   sens_cov2_adult.sx.list[[5]][4, hull_cov2_adult.sx.list[[14]]], 
   lty=4, border=pal[4])
polygon(1-spec_cov2_adult.sx.list[[5]][5, hull_cov2_adult.sx.list[[15]]], 
   sens_cov2_adult.sx.list[[5]][5, hull_cov2_adult.sx.list[[15]]], 
   lty=5, border=pal[5])
legend(x="bottomright", hull_names[11:15], lty=1:5, col=pal[1:5], lwd=2, bty="n")
dev.off()

pdf("hh-combos_cov2_child.pdf")
# 1. Baseline rules
setup.sens.spec.plot(base.algs=1, sens=sens_cov2_child.rl[orule_idx], 
   spec=spec_cov2_child.rl[orule_idx])
title(main="Performance of existing rules, PCR, children", adj=0)

# 2. Cloud of single symptoms
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_child.rl[orule_idx], 
   spec=spec_cov2_child.rl[orule_idx])
title(main="All single-symptom rules, PCR, children", adj=0)
polygon(1-spec_cov2_child.sx.list[[1]][1, hull_cov2_child.sx.list[[1]]], 
   sens_cov2_child.sx.list[[1]][1, hull_cov2_child.sx.list[[1]]],
   lwd=1, border="gray")
text(1-spec_cov2_child.sx.list[[1]][1,], 
   sens_cov2_child.sx.list[[1]][1,], names.sx.list[[1]], cex=0.95, col=pal[2])
legend(x="bottomright", legend=c(symp_labels), title="Single symptoms",
   pch=LETTERS[1:15], cex=0.85,
   col="black", bty="n")

# 3. Cloud for 1+/2 + convex hull
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_child.rl[orule_idx], 
   spec=spec_cov2_child.rl[orule_idx])
title(main="All combinations, either of 2 symptoms, PCR, children", adj=0)
polygon(1-spec_cov2_child.sx.list[[2]][1, hull_cov2_child.sx.list[[2]]], 
   sens_cov2_child.sx.list[[2]][1, hull_cov2_child.sx.list[[2]]],
   lwd=1, border="gray")
text(1-spec_cov2_child.sx.list[[2]][1,], 
   sens_cov2_child.sx.list[[2]][1,], names.sx.list[[2]], cex=0.5, col=pal[1], srt=45)
legend(x="bottomright", legend=c(symp_labels), title="Single symptoms",
   pch=LETTERS[1:15], cex=0.85,
   col="black", bty="n")

# 4. Convex hulls for 1+/2, 2/2 + baseline rules
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_child.rl[orule_idx], 
   spec=spec_cov2_child.rl[orule_idx])
title(main="All combinations of 2 symptoms, PCR, children", adj=0)
polygon(1-spec_cov2_child.sx.list[[2]][1, hull_cov2_child.sx.list[[2]]], 
   sens_cov2_child.sx.list[[2]][1, hull_cov2_child.sx.list[[2]]],
   border=pal[1])
polygon(1-spec_cov2_child.sx.list[[2]][2, hull_cov2_child.sx.list[[3]]], 
   sens_cov2_child.sx.list[[2]][2, hull_cov2_child.sx.list[[3]]], 
   lty=2, border=pal[3])
legend(x="bottomright", hull_names[2:3], lty=1:2, col=pal[c(1,3)], lwd=2, bty="n")

# 5. Convex hulls for 1+/3, 2+/3, 3/3 + baseline rules
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_child.rl[orule_idx], 
   spec=spec_cov2_child.rl[orule_idx])
title(main="All combinations of 3 symptoms, PCR, children", adj=0)
polygon(1-spec_cov2_child.sx.list[[3]][1, hull_cov2_child.sx.list[[4]]], 
   sens_cov2_child.sx.list[[3]][1, hull_cov2_child.sx.list[[4]]],
   border=pal[1])
polygon(1-spec_cov2_child.sx.list[[3]][2, hull_cov2_child.sx.list[[5]]], 
   sens_cov2_child.sx.list[[3]][2, hull_cov2_child.sx.list[[5]]], 
   lty=2, border=pal[3])
polygon(1-spec_cov2_child.sx.list[[3]][3, hull_cov2_child.sx.list[[6]]], 
   sens_cov2_child.sx.list[[3]][3, hull_cov2_child.sx.list[[6]]], 
   lty=3, lwd=2, border=pal[5])
legend(x="bottomright", hull_names[4:6], lty=1:3, col=pal[c(1,3,5)], lwd=2, bty="n")

# 6. Convex hulls for 1+/4, 2+/4, 3+/4, 4/4
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_child.rl[orule_idx], 
   spec=spec_cov2_child.rl[orule_idx])
title(main="All combinations of 4 symptoms, PCR, children", adj=0)
polygon(1-spec_cov2_child.sx.list[[4]][1, hull_cov2_child.sx.list[[7]]], 
   sens_cov2_child.sx.list[[4]][1, hull_cov2_child.sx.list[[7]]],
   border=pal[1])
polygon(1-spec_cov2_child.sx.list[[4]][2, hull_cov2_child.sx.list[[8]]], 
   sens_cov2_child.sx.list[[4]][2, hull_cov2_child.sx.list[[8]]], 
   lty=2, border=pal[2])
polygon(1-spec_cov2_child.sx.list[[4]][3, hull_cov2_child.sx.list[[9]]], 
   sens_cov2_child.sx.list[[4]][3, hull_cov2_child.sx.list[[9]]], 
   lty=3, lwd=2, border=pal[3])
polygon(1-spec_cov2_child.sx.list[[4]][4, hull_cov2_child.sx.list[[10]]], 
   sens_cov2_child.sx.list[[4]][4, hull_cov2_child.sx.list[[10]]], 
   lty=4, border=pal[4])
legend(x="bottomright", hull_names[7:10], lty=1:4, col=pal[1:4], lwd=2, bty="n")

# 7. Convex hulls for 1+/5, 2+/5, 3+/5, 4+/5, 5/5
setup.sens.spec.plot(base.algs=2, sens=sens_cov2_child.rl[orule_idx], 
   spec=spec_cov2_child.rl[orule_idx])
title(main="All combinations of 5 symptoms, PCR, children", adj=0)
polygon(1-spec_cov2_child.sx.list[[5]][1, hull_cov2_child.sx.list[[11]]], 
   sens_cov2_child.sx.list[[5]][1, hull_cov2_child.sx.list[[11]]],
   border=pal[1])
polygon(1-spec_cov2_child.sx.list[[5]][2, hull_cov2_child.sx.list[[12]]], 
   sens_cov2_child.sx.list[[5]][2, hull_cov2_child.sx.list[[12]]], 
   lty=2, border=pal[2])
polygon(1-spec_cov2_child.sx.list[[5]][3, hull_cov2_child.sx.list[[13]]], 
   sens_cov2_child.sx.list[[5]][3, hull_cov2_child.sx.list[[13]]], 
   lty=3, lwd=2, border=pal[3])
polygon(1-spec_cov2_child.sx.list[[5]][4, hull_cov2_child.sx.list[[14]]], 
   sens_cov2_child.sx.list[[5]][4, hull_cov2_child.sx.list[[14]]], 
   lty=4, border=pal[4])
polygon(1-spec_cov2_child.sx.list[[5]][5, hull_cov2_child.sx.list[[15]]], 
   sens_cov2_child.sx.list[[5]][5, hull_cov2_child.sx.list[[15]]], 
   lty=5, border=pal[5])
legend(x="bottomright", hull_names[11:15], lty=1:5, col=pal[1:5], lwd=2, bty="n")
dev.off()
