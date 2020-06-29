### Symptom combinations for detecting SARS-CoV-2 infection in nonhospitalized persons
### Graphics: how selected rules perform in adults and children

# Simple graphics: 5 symptom groupings, 3 rule groupings

library("viridis")

# viridis(2) # "#440154FF" "#FDE725FF"

# Individual symptoms, superposed
plot(c(0, 1), c(0, 1), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
axis(1, at=(0:5)/5, labels=paste0((0:5)*20, "%"))
axis(2, at=(0:5)/5, labels=paste0((0:5)*20, "%"), las=1)
segments(x0=0, x1=1, y0=0, y1=1, col="gray", lty="dotted", lwd=2)

segments(x0=1-spec_cov2_adult.sx.list[[1]], y0=sens_cov2_adult.sx.list[[1]],
   x1=1-spec_cov2_child.sx.list[[1]], y1=sens_cov2_child.sx.list[[1]],
   col="gray")
points(1-spec_cov2_adult.sx.list[[1]], sens_cov2_adult.sx.list[[1]], 
   pch=21, bg="#440154FF", col="black")
points(1-spec_cov2_child.sx.list[[1]], sens_cov2_child.sx.list[[1]], 
   pch=21, bg="#FDE725FF", col="#440154FF")
text(labels=LETTERS[1:15], 1-spec_cov2_adult.sx.list[[1]], 
   sens_cov2_adult.sx.list[[1]]+0.025)
text(labels=letters[1:15], 1-spec_cov2_child.sx.list[[1]], 
   sens_cov2_child.sx.list[[1]]-0.025)
legend(x="topright", legend=c("Adults", "Children"), 
   pch=21, col=c("black", "#440154FF"), pt.bg=c("#440154FF", "#FDE725FF"),
   bty="n")
legend(x="bottomright", legend=c(symp_labels), title="Single symptoms",
   pch=LETTERS[1:15], cex=0.85,
   col="black", bty="n")

# Rules, superposed
plot(c(0, 1), c(0, 1), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
axis(1, at=(0:5)/5, labels=paste0((0:5)*20, "%"))
axis(2, at=(0:5)/5, labels=paste0((0:5)*20, "%"), las=1)
segments(x0=0, x1=1, y0=0, y1=1, col="gray", lty="dotted", lwd=2)

segments(x0=1-spec_cov2_adult.rl[1:5], y0=sens_cov2_adult.rl[1:5],
   x1=1-spec_cov2_child.rl[1:5], y1=sens_cov2_child.rl[1:5],
   col="gray")
points(1-spec_cov2_adult.rl[1:5], sens_cov2_adult.rl[1:5], 
   pch=21, bg="#440154FF", col="black")
points(1-spec_cov2_child.rl[1:5], sens_cov2_child.rl[1:5], 
   pch=21, bg="#FDE725FF", col="#440154FF")
text(labels=rule_labels[1:5], 1-spec_cov2_adult.rl[1:5], 
   sens_cov2_adult.rl[1:5]+0.025)
text(labels=tolower(rule_labels[1:5]), 1-spec_cov2_child.rl[1:5], 
   sens_cov2_child.rl[1:5]-0.025)
legend(x="bottomright", legend=c("Adults", "Children"), 
   pch=21, col=c("black", "#440154FF"), pt.bg=c("#440154FF", "#FDE725FF"),
   bty="n")

### Rules and symptoms in groups

## Rules in groups

pdf("rules-adult-child.pdf")
# Existing rules via orule_idx ILI, CDC, ARI, CSTE, CLI
plot(c(0, 0.79), c(0.25, 1.04), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
title(main="Existing symptom-based rules", adj=0)
axis(1, seq(0, 0.75, 0.05), labels=FALSE); axis(2, seq(0.25, 1, 0.05), labels=FALSE)
axis(1, at=seq(0, 0.7, 0.1), labels=paste0(seq(0, 70, 10), "%"))
axis(2, at=seq(0.3, 1, 0.1), labels=paste0(seq(30, 100, 10), "%"), las=1)
segments(x0=0, x1=1, y0=0, y1=1, col="gray", lty="dotted", lwd=2)
legend(x="bottomright", legend=c("Adults", "Children"), 
   pch=21, col=c("black", "#440154FF"), pt.bg=c("#440154FF", "#FDE725FF"),
   bty="n")
segments(x0=1-spec_cov2_adult.rl[orule_idx], y0=sens_cov2_adult.rl[orule_idx],
   x1=1-spec_cov2_child.rl[orule_idx], y1=sens_cov2_child.rl[orule_idx],
   col="gray")
points(1-spec_cov2_adult.rl[orule_idx], sens_cov2_adult.rl[orule_idx], 
   pch=21, bg="#440154FF", col="black")
points(1-spec_cov2_child.rl[orule_idx], sens_cov2_child.rl[orule_idx], 
   pch=21, bg="#FDE725FF", col="#440154FF")
# blend1 <- c(1/2, 3/4, 1/2, 1/3, 1/2)
# text(labels=rule_labels[1:5],
   # 1-(blend1*spec_cov2_adult.rl[1:5]+(1-blend1)*spec_cov2_child.rl[1:5]),
   # blend1*sens_cov2_adult.rl[1:5]+(1-blend1)*sens_cov2_child.rl[1:5])
text(labels=paste0(
   rule_labels[orule_idx], "\n(",
   round(100*spec_cov2_adult.rl[orule_idx]), ", ", 
   round(100*sens_cov2_adult.rl[orule_idx]), ")"),
   1-spec_cov2_adult.rl[orule_idx], sens_cov2_adult.rl[orule_idx]+c(1,1,1,-1,1)*0.025,
   cex = 0.5)
text(labels=paste0(
   rule_labels[orule_idx], "\n(",
   round(100*spec_cov2_child.rl[orule_idx]), ", ", 
   round(100*sens_cov2_child.rl[orule_idx]), "}"), 
   1-spec_cov2_child.rl[orule_idx], sens_cov2_child.rl[orule_idx]+c(-1,1,-1,-1,-1)*0.025,
   cex = 0.5)

# "Vaccine A" rules
plot(c(0, 0.75), c(0.25, 1), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
title(main="\"Vaccine A\" symptom-based rules", adj=0)
axis(1, seq(0, 0.75, 0.05), labels=FALSE); axis(2, seq(0.25, 1, 0.05), labels=FALSE)
axis(1, seq(0, 0.7, 0.1)); axis(2, seq(0.3, 1, 0.1), las=1)
segments(x0=0, x1=1, y0=0, y1=1, col="gray", lty="dotted", lwd=2)
legend(x="bottomright", legend=c("Adults", "Children"), 
   pch=21, col=c("black", "#440154FF"), pt.bg=c("#440154FF", "#FDE725FF"),
   bty="n")
segments(x0=1-spec_cov2_adult.rl[vaxa_idx], y0=sens_cov2_adult.rl[vaxa_idx],
   x1=1-spec_cov2_child.rl[vaxa_idx], y1=sens_cov2_child.rl[vaxa_idx],
   col="gray")
points(1-spec_cov2_adult.rl[vaxa_idx], sens_cov2_adult.rl[vaxa_idx], 
   pch=21, bg="#440154FF", col="black")
points(1-spec_cov2_child.rl[vaxa_idx], sens_cov2_child.rl[vaxa_idx], 
   pch=21, bg="#FDE725FF", col="#440154FF")
# blend1 <- c(1/2, 3/4, 1/2, 1/3, 1/2)
# text(labels=rule_labels[1:5],
# 1-(blend1*spec_cov2_adult.rl[1:5]+(1-blend1)*spec_cov2_child.rl[1:5]),
# blend1*sens_cov2_adult.rl[1:5]+(1-blend1)*sens_cov2_child.rl[1:5])
text(labels=paste0(
   rule_labels[vaxa_idx], "\n(",
   round(100*spec_cov2_adult.rl[vaxa_idx]), ", ", 
   round(100*sens_cov2_adult.rl[vaxa_idx]), ")"),
   1-spec_cov2_adult.rl[vaxa_idx], sens_cov2_adult.rl[vaxa_idx]+c(1,1,1,1)*0.025,
   cex = 0.5)
text(labels=paste0(
   rule_labels[vaxa_idx], "\n(",
   round(100*spec_cov2_child.rl[vaxa_idx]), ", ", 
   round(100*sens_cov2_child.rl[vaxa_idx]), "}"), 
   1-spec_cov2_child.rl[vaxa_idx], sens_cov2_child.rl[vaxa_idx]+c(-1,-1,-1,-1)*0.025,
   cex = 0.5)
dev.off()

## Symptoms in groups

# pdf("symptoms-adult-child.pdf")
# constitutional: cons_idx
# myalgia, fatigue, fever_chills
plot(c(0, 0.6), c(0.4, 1), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
title(main="Constitutional symptoms", adj=0)
axis(1, seq(0, 0.6, 0.05), labels=FALSE); axis(2, seq(0.4, 1, 0.05), labels=FALSE)
axis(1, seq(0, 0.6, 0.1)); axis(2, seq(0.4, 1, 0.1), las=1)
segments(x0=0, x1=1, y0=0, y1=1, col="gray", lty="dotted", lwd=2)
legend(x="bottomright", legend=c("Adults", "Children"), 
   pch=21, col=c("black", "#440154FF"), pt.bg=c("#440154FF", "#FDE725FF"),
   bty="n")
segments(x0=1-spec_cov2_adult.sx.list[[1]][cons_idx], 
   y0=sens_cov2_adult.sx.list[[1]][cons_idx],
   x1=1-spec_cov2_child.sx.list[[1]][cons_idx], 
   y1=sens_cov2_child.sx.list[[1]][cons_idx],
   col="gray")
points(1-spec_cov2_adult.sx.list[[1]][cons_idx], 
   sens_cov2_adult.sx.list[[1]][cons_idx], 
   pch=21, bg="#440154FF", col="black")
points(1-spec_cov2_child.sx.list[[1]][cons_idx], 
   sens_cov2_child.sx.list[[1]][cons_idx], 
   pch=21, bg="#FDE725FF", col="#440154FF")
text(labels=paste0(
   symp_labels[cons_idx], "\n(",
   round(100*spec_cov2_adult.sx.list[[1]][cons_idx]), ", ", 
   round(100*sens_cov2_adult.sx.list[[1]][cons_idx]), ")"),
   1-spec_cov2_adult.sx.list[[1]][cons_idx]+c(1,0,0)*0.03, 
   sens_cov2_adult.sx.list[[1]][cons_idx]+c(0,1,1)*0.025,
   cex = 0.5)
text(labels=paste0(
   symp_labels[cons_idx], "\n(",
   round(100*spec_cov2_child.sx.list[[1]][cons_idx]), ", ", 
   round(100*sens_cov2_child.sx.list[[1]][cons_idx]), "}"), 
   1-spec_cov2_child.sx.list[[1]][cons_idx]+c(0,0,-1)*0.03, 
   sens_cov2_child.sx.list[[1]][cons_idx]+c(-1,-1,0)*0.025,
   cex = 0.5)

# Upper respiratory
# sore throat, nasal
plot(c(0, 1), c(0, 1), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
title(main="Upper respiratory symptoms", adj=0)
axis(1, seq(0, 1, 0.05), labels=FALSE); axis(2, seq(0, 1, 0.05), labels=FALSE)
axis(1, seq(0, 1, 0.1)); axis(2, seq(0, 1, 0.1), las=1)
segments(x0=0, x1=1, y0=0, y1=1, col="gray", lty="dotted", lwd=2)
legend(x="bottomright", legend=c("Adults", "Children"), 
   pch=21, col=c("black", "#440154FF"), pt.bg=c("#440154FF", "#FDE725FF"),
   bty="n")
segments(x0=1-spec_cov2_adult.sx.list[[1]][upresp_idx], 
   y0=sens_cov2_adult.sx.list[[1]][upresp_idx],
   x1=1-spec_cov2_child.sx.list[[1]][upresp_idx], 
   y1=sens_cov2_child.sx.list[[1]][upresp_idx],
   col="gray")
points(1-spec_cov2_adult.sx.list[[1]][upresp_idx], 
   sens_cov2_adult.sx.list[[1]][upresp_idx], 
   pch=21, bg="#440154FF", col="black")
points(1-spec_cov2_child.sx.list[[1]][upresp_idx], 
   sens_cov2_child.sx.list[[1]][upresp_idx], 
   pch=21, bg="#FDE725FF", col="#440154FF")
text(labels=paste0(
   symp_labels[upresp_idx], "\n(",
   round(100*spec_cov2_adult.sx.list[[1]][upresp_idx]), ", ", 
   round(100*sens_cov2_adult.sx.list[[1]][upresp_idx]), ")"),
   1-spec_cov2_adult.sx.list[[1]][upresp_idx]+c(1,1)*0.06, 
   sens_cov2_adult.sx.list[[1]][upresp_idx]+c(0,0)*0.03,
   cex = 0.5)
text(labels=paste0(
   symp_labels[upresp_idx], "\n(",
   round(100*spec_cov2_child.sx.list[[1]][upresp_idx]), ", ", 
   round(100*sens_cov2_child.sx.list[[1]][upresp_idx]), "}"), 
   1-spec_cov2_child.sx.list[[1]][upresp_idx]+c(-1,-1)*0.06, 
   sens_cov2_child.sx.list[[1]][upresp_idx]+c(0,0)*0.03,
   cex = 0.5)

# Lower respiratory
# wheezing, shortness of breath, pain breathing, cough, chest pain
plot(c(0, 1), c(0, 1), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
title(main="Lower respiratory symptoms", adj=0)
axis(1, seq(0, 1, 0.05), labels=FALSE); axis(2, seq(0, 1, 0.05), labels=FALSE)
axis(1, seq(0, 1, 0.1)); axis(2, seq(0, 1, 0.1), las=1)
segments(x0=0, x1=1, y0=0, y1=1, col="gray", lty="dotted", lwd=2)
legend(x="bottomright", legend=c("Adults", "Children"), 
   pch=21, col=c("black", "#440154FF"), pt.bg=c("#440154FF", "#FDE725FF"),
   bty="n")
segments(x0=1-spec_cov2_adult.sx.list[[1]][loresp_idx], 
   y0=sens_cov2_adult.sx.list[[1]][loresp_idx],
   x1=1-spec_cov2_child.sx.list[[1]][loresp_idx], 
   y1=sens_cov2_child.sx.list[[1]][loresp_idx],
   col="gray")
points(1-spec_cov2_adult.sx.list[[1]][loresp_idx], 
   sens_cov2_adult.sx.list[[1]][loresp_idx], 
   pch=21, bg="#440154FF", col="black")
points(1-spec_cov2_child.sx.list[[1]][loresp_idx], 
   sens_cov2_child.sx.list[[1]][loresp_idx], 
   pch=21, bg="#FDE725FF", col="#440154FF")
text(labels=paste0(
   symp_labels[loresp_idx], "\n(",
   round(100*spec_cov2_adult.sx.list[[1]][loresp_idx]), ", ", 
   round(100*sens_cov2_adult.sx.list[[1]][loresp_idx]), ")"),
   1-spec_cov2_adult.sx.list[[1]][loresp_idx]+c(0,1,0,0,0)*0.06, 
   sens_cov2_adult.sx.list[[1]][loresp_idx]+c(-1,0,1,1,1)*0.03,
   cex = 0.5)
text(labels=paste0(
   symp_labels[loresp_idx], "\n(",
   round(100*spec_cov2_child.sx.list[[1]][loresp_idx]), ", ", 
   round(100*sens_cov2_child.sx.list[[1]][loresp_idx]), "}"), 
   1-spec_cov2_child.sx.list[[1]][loresp_idx]+c(0,0,0,0,1)*0.06, 
   sens_cov2_child.sx.list[[1]][loresp_idx]+c(-1,1,1,1,0)*0.03,
   cex = 0.5)

# Neurological
# headache, taste_smell
# pdf("symptoms-othr-adult-child_20200529.pdf")
plot(c(0, 0.6), c(0.4, 1), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
title(main="Neurological symptoms", adj=0)
axis(1, seq(0, 0.6, 0.05), labels=FALSE); axis(2, seq(0.4, 1, 0.05), labels=FALSE)
axis(1, seq(0, 0.6, 0.1)); axis(2, seq(0.4, 1, 0.1), las=1)
segments(x0=0, x1=1, y0=0, y1=1, col="gray", lty="dotted", lwd=2)
legend(x="bottomright", legend=c("Adults", "Children"), 
   pch=21, col=c("black", "#440154FF"), pt.bg=c("#440154FF", "#FDE725FF"),
   bty="n")
segments(x0=1-spec_cov2_adult.sx.list[[1]][neuro_idx], 
   y0=sens_cov2_adult.sx.list[[1]][neuro_idx],
   x1=1-spec_cov2_child.sx.list[[1]][neuro_idx], 
   y1=sens_cov2_child.sx.list[[1]][neuro_idx],
   col="gray")
points(1-spec_cov2_adult.sx.list[[1]][neuro_idx], 
   sens_cov2_adult.sx.list[[1]][neuro_idx], 
   pch=21, bg="#440154FF", col="black")
points(1-spec_cov2_child.sx.list[[1]][neuro_idx], 
   sens_cov2_child.sx.list[[1]][neuro_idx], 
   pch=21, bg="#FDE725FF", col="#440154FF")
text(labels=paste0(
   symp_labels[neuro_idx], "\n(",
   round(100*spec_cov2_adult.sx.list[[1]][neuro_idx]), ", ", 
   round(100*sens_cov2_adult.sx.list[[1]][neuro_idx]), ")"),
   1-spec_cov2_adult.sx.list[[1]][neuro_idx]+c(0,0)*0.03, 
   sens_cov2_adult.sx.list[[1]][neuro_idx]+c(1,1)*0.025,
   cex = 0.5)
text(labels=paste0(
   symp_labels[neuro_idx], "\n(",
   round(100*spec_cov2_child.sx.list[[1]][neuro_idx]), ", ", 
   round(100*sens_cov2_child.sx.list[[1]][neuro_idx]), "}"), 
   1-spec_cov2_child.sx.list[[1]][neuro_idx]+c(0,0)*0.03, 
   sens_cov2_child.sx.list[[1]][neuro_idx]+c(-1,-1)*0.025,
   cex = 0.5)

# Gastrointestinal
# nausea, diarrhea, abdominal pain
# pdf("symptoms-othr-adult-child_20200529.pdf")
plot(c(0.0, 0.4), c(0.1, 0.5), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
title(main="Gastrointestinal symptoms", adj=0)
axis(1, seq(0, 0.4, 0.05), labels=FALSE); axis(2, seq(0.1, 0.5, 0.05), labels=FALSE)
axis(1, seq(0, 0.4, 0.1)); axis(2, seq(0.1, 0.5, 0.1), las=1)
segments(x0=0, x1=1, y0=0, y1=1, col="gray", lty="dotted", lwd=2)
legend(x="bottomright", legend=c("Adults", "Children"), 
   pch=21, col=c("black", "#440154FF"), pt.bg=c("#440154FF", "#FDE725FF"),
   bty="n")
segments(x0=1-spec_cov2_adult.sx.list[[1]][gastro_idx], 
   y0=sens_cov2_adult.sx.list[[1]][gastro_idx],
   x1=1-spec_cov2_child.sx.list[[1]][gastro_idx], 
   y1=sens_cov2_child.sx.list[[1]][gastro_idx],
   col="gray")
points(1-spec_cov2_adult.sx.list[[1]][gastro_idx], 
   sens_cov2_adult.sx.list[[1]][gastro_idx], 
   pch=21, bg="#440154FF", col="black")
points(1-spec_cov2_child.sx.list[[1]][gastro_idx], 
   sens_cov2_child.sx.list[[1]][gastro_idx], 
   pch=21, bg="#FDE725FF", col="#440154FF")
text(labels=paste0(
   symp_labels[gastro_idx], "\n(",
   round(100*spec_cov2_adult.sx.list[[1]][gastro_idx]), ", ", 
   round(100*sens_cov2_adult.sx.list[[1]][gastro_idx]), ")"),
   1-spec_cov2_adult.sx.list[[1]][gastro_idx]+c(-1,0,0)*0.02, 
   sens_cov2_adult.sx.list[[1]][gastro_idx]+c(0,1,1)*0.015,
   cex = 0.5)
text(labels=paste0(
   symp_labels[gastro_idx], "\n(",
   round(100*spec_cov2_child.sx.list[[1]][gastro_idx]), ", ", 
   round(100*sens_cov2_child.sx.list[[1]][gastro_idx]), "}"), 
   1-spec_cov2_child.sx.list[[1]][gastro_idx]+c(1,0,0)*0.02, 
   sens_cov2_child.sx.list[[1]][gastro_idx]+c(0,-1,-1)*0.015,
   cex = 0.5)
# dev.off()

##### graphics
pi_seq <- seq(0, 2*pi, length=181)
crit <- sqrt(2*qf(c(0.95, 0.975, 0.99, 0.995), df1=2, df2=Inf))

# pdf("hh-resampling_20200529.pdf")
setup.sens.spec.plot(base.algs=1, sens=sens_cov2_adult.rl[orule_idx], 
   spec=spec_cov2_adult.rl[orule_idx])
title(main="Rule performance, PCR-confirmed, adults\nwith 95% bootstrap confidence ellipse", adj=0)
for(j in orule_idx) {
   points(1-spec_cov2_adult_bsmean[j], sens_cov2_adult_bsmean[j], cex=0.5)
   polygon(1-(spec_cov2_adult_bsmean[j] + 
         sqrt(sesp_cov2_adult_bscov[[j]][2, 2])*crit[1]*
            cos(pi_seq+acos(sesp_cov2_adult_bscov[[j]][1,2]))),
      sens_cov2_adult_bsmean[j] + 
         sqrt(sesp_cov2_adult_bscov[[j]][1, 1])*crit[1]*cos(pi_seq))
}

setup.sens.spec.plot(base.algs=1, sens=sens_cov2_child.rl[orule_idx], 
   spec=spec_cov2_child.rl[orule_idx])
title(main="Rule performance, PCR-confirmed, children\nwith 95% bootstrap confidence ellipse", adj=0)
for(j in orule_idx) {
   points(1-spec_cov2_child_bsmean[j], sens_cov2_child_bsmean[j], cex=0.5)
   polygon(1-(spec_cov2_child_bsmean[j] + 
         sqrt(sesp_cov2_child_bscov[[j]][2, 2])*crit[1]*
         cos(pi_seq+acos(sesp_cov2_child_bscov[[j]][1,2]))),
      sens_cov2_child_bsmean[j] + 
         sqrt(sesp_cov2_child_bscov[[j]][1, 1])*crit[1]*cos(pi_seq))
}
# dev.off()


# Bivariate differences between adults and children
# apply(apply(hh_psamp_rusy_diff_sesp_array[,,,"pcr"], c(2,3), mean), 2, range)
# apply(apply(hh_psamp_rusy_diff_sesp_array[,,,"pcr"], c(2,3), quantile, pr=0.01), 2, range)
# bounding box for means of differences: c(-.17, -.13), c(.32, .45)
# bounding box for bivariate bootstrap distributions of differences:
#    c(-.5, -.5), c(.5, .5) or c(-.6, -.6), c(.6, .6)

# for each rule or symptom, graph a subset of differences (101?),
#    95% confidence ellipse, bootstrap mean difference, observed difference

pal <- viridis(6)
bb <- 0.5 # half-width of bounding box

# single example -- ILI
plot(bb*c(-1, 1), bb*c(-1, 1), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
axis(1, at=seq(-bb, bb, by=0.1), labels=paste0(100*seq(-bb, bb, by=0.1), "%"))
axis(2, at=seq(-bb, bb, by=0.1), labels=paste0(100*seq(-bb, bb, by=0.1), "%"), las=1)
segments(x0=-bb, x1=bb, y0=0, col="gray", lty="dashed")
segments(x0=0, y0=-bb, y1=bb, col="gray", lty="dashed")
title(main=paste("Adult/child differences:", rusy_labels[1]), adj=0)
points(-hh_psamp_rusy_diff_sesp_array[1:401, rusy_labels[1], "spec"],
   hh_psamp_rusy_diff_sesp_array[1:401, rusy_labels[1], "sens"],
   cex=0.5, pch=1, col="gray") 
polygon(
   -(spec_cov2_diff_bsmean[1] + sqrt(sesp_cov2_diff_bscov[[1]][2, 2])*crit[1]*
      cos(pi_seq+acos(sesp_cov2_diff_bscov[[1]][1,2]))),
   sens_cov2_diff_bsmean[1] + sqrt(sesp_cov2_diff_bscov[[1]][1, 1])*crit[1]*cos(pi_seq))
points(-spec_cov2_diff_bsmean[1], sens_cov2_diff_bsmean[1], cex=2,
   pch=21, bg=pal[3], col=pal[1])

# Sequence of comparisons
# pdf("hh-adult-child-diffs-5.pdf")
for(.i. in 1:24) {
plot(c(-0.5, 0.5), c(-0.5, 0.5), type="n", bty="n", xaxt="n", yaxt="n", 
   xlab=expression("100%" - "Specificity"), ylab="Sensitivity")
axis(1, at=seq(-5, 5)/10, labels=paste0(seq(-5, 5)*10, "%"))
axis(2, at=seq(-5, 5)/10, labels=paste0(seq(-5, 5)*10, "%"), las=1)
segments(x0=-0.5, x1=0.5, y0=0, col="gray", lty="dashed")
segments(x0=0, y0=-0.5, y1=0.56, col="gray", lty="dashed")
title(main=paste("Adult/child differences:", rusy_labels[.i.]), adj=0)
points(-hh_psamp_rusy_diff_sesp_array[1:401, rusy_labels[.i.], "spec"],
   hh_psamp_rusy_diff_sesp_array[1:401, rusy_labels[.i.], "sens"],
   cex=0.5, pch=1, col="gray") 
polygon(
   -(spec_cov2_diff_bsmean[.i.] + sqrt(sesp_cov2_diff_bscov[[.i.]][2, 2])*crit[1]*
         cos(pi_seq+acos(sesp_cov2_diff_bscov[[.i.]][1,2]))),
   sens_cov2_diff_bsmean[.i.] + sqrt(sesp_cov2_diff_bscov[[.i.]][1, 1])*crit[1]*cos(pi_seq))
points(-spec_cov2_diff_bsmean[.i.], sens_cov2_diff_bsmean[.i.], cex=2,
   pch=21, bg=pal[3], col=pal[1])
}
# dev.off()
