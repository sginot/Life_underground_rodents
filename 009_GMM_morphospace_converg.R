source("005_phylogeny_and_contrasts_generic.R")
source("007_import_landmarks.R")
library(scales)
library(viridis)
library(geomorph)
library(paleomorph)

# Opt to print plots or not:
print.plots <- F
if (print.plots & !dir.exists("plots")) {dir.create("plots")} #Create output folder if it does not exist

# Load project metadata
disparat <- read.csv("LISTE_IODINE_RODENT - Main_list.csv", h=T, dec=".", sep=",")

# Subset the metadata
d <- disparat[match(specimens$genus, disparat$Genre),]
dd <- disparat[match(spec$genus, disparat$Genre),]
ddd <- disparat[match(spem$genus, disparat$Genre),]

##############################################
#Morphospace occupation by main rodent clades#
##############################################

if (print.plots) {
pdf(file = paste("plots", "pca_craniumLM_genera.pdf", sep = "/"), width = 10, height = 10)
plot(pcac$x[,1:2], asp = 1, pch = 21, bg = c(1:3)[as.factor(d$Clade)])
text(pcac$x[,1:2], labels = specimens$genus, cex = 0.5, pos = 2)
dev.off()

pdf(file = paste("plots", "pca_allLMandcurves_genera.pdf", sep = "/"), width = 10, height = 10)
plot(pcaM$x[,1:2], asp = 1, pch = 21, bg = c(1:3)[as.factor(dd$Clade)])
text(pcaM$x[,1:2], labels = spec$genus, cex = 0.5, pos = 2)
dev.off()

pdf(file = paste("plots", "pca_mandibleLM_genera.pdf", sep = "/"), width = 10, height = 10)
plot(pcamd$x[,1:2], asp = 1, pch = 21, bg = c(1:3)[as.factor(ddd$Clade)])
text(pcamd$x[,1:2], labels = spem$genus, cex = 0.5, pos = 2)
dev.off()
}

#######################################
#Morphospace occupation clade by clade#
#######################################
if (print.plots) {
pdf(file = paste("plots", "clade_by_clade_space.pdf", sep = "/"), width = 15, height = 10)
layout(matrix(1:9, ncol = 3, byrow=F))
par(mar = c(1,3,1,0.2))

# Ctenohystrica
plot(pcac$x[which(d$Clade == "CTENO"),1:2],pch = 21, bg = 1, ylab="", xlab="", main = "Ctenohystrica", xaxt = "n", yaxt = "n", xlim = c(min(pcac$x[,1]), max(pcac$x[,1])), ylim = c(min(pcac$x[,2]), max(pcac$x[,2])))
text(pcac$x[which(d$Clade == "CTENO"),1:2], labels = d$Genre[which(d$Clade == "CTENO")], cex = 0.5, pos = 2)
title(ylab="Cranium landmarks", line=1, cex.lab=1.5)
abline(h = 0, v = 0, col ="grey", lty = 2)

plot(pcamd$x[which(ddd$Clade == "CTENO"),1:2],pch = 21, bg = 1, ylab="", xlab="", main = "Ctenohystrica", xaxt = "n", yaxt = "n", xlim = c(min(pcamd$x[,1]), max(pcamd$x[,1])), ylim = c(min(pcamd$x[,2]), max(pcamd$x[,2])))
text(pcamd$x[which(ddd$Clade == "CTENO"),1:2], labels = ddd$Genre[which(ddd$Clade == "CTENO")], cex = 0.5, pos = 2)
title(ylab="Mandible landmarks", line=1, cex.lab=1.5)
abline(h = 0, v = 0, col ="grey", lty = 2)

plot(pcaM$x[which(dd$Clade == "CTENO"),1:2], pch = 21, bg = 1, ylab="", xlab="", main = "Ctenohystrica", xaxt = "n", yaxt = "n", xlim = c(min(pcaM$x[,1]), max(pcaM$x[,1])), ylim = c(min(pcaM$x[,2]), max(pcaM$x[,2])))
text(pcaM$x[which(dd$Clade == "CTENO"),1:2], labels = dd$Genre[which(dd$Clade == "CTENO")], cex = 0.5, pos = 2)
title(ylab="All landmarks and curves", line=1, cex.lab=1.5)
abline(h = 0, v = 0, col ="grey", lty = 2)

# Mouse related clade
plot(pcac$x[which(d$Clade == "MOUS"),1:2],pch = 21, bg = 2, ylab="", xlab="", main = "Mouse-related", xaxt = "n", yaxt = "n", xlim = c(min(pcac$x[,1]), max(pcac$x[,1])), ylim = c(min(pcac$x[,2]), max(pcac$x[,2])))
text(pcac$x[which(d$Clade == "MOUS"),1:2], labels = d$Genre[which(d$Clade == "MOUS")], cex = 0.5, pos = 2)
title(ylab="Cranium landmarks", line=1, cex.lab=1.5)
abline(h = 0, v = 0, col ="grey", lty = 2)

plot(pcamd$x[which(ddd$Clade == "MOUS"),1:2],pch = 21, bg = 2, ylab="", xlab="", main = "Mouse-related", xaxt = "n", yaxt = "n", xlim = c(min(pcamd$x[,1]), max(pcamd$x[,1])), ylim = c(min(pcamd$x[,2]), max(pcamd$x[,2])))
text(pcamd$x[which(ddd$Clade == "MOUS"),1:2], labels = ddd$Genre[which(ddd$Clade == "MOUS")], cex = 0.5, pos = 2)
title(ylab="Mandible landmarks", line=1, cex.lab=1.5)
abline(h = 0, v = 0, col ="grey", lty = 2)

plot(pcaM$x[which(dd$Clade == "MOUS"),1:2], pch = 21, bg = 2, ylab="", xlab="", main = "Mouse-related", xaxt = "n", yaxt = "n", xlim = c(min(pcaM$x[,1]), max(pcaM$x[,1])), ylim = c(min(pcaM$x[,2]), max(pcaM$x[,2])))
text(pcaM$x[which(dd$Clade == "MOUS"),1:2], labels = dd$Genre[which(dd$Clade == "MOUS")], cex = 0.5, pos = 2)
title(ylab="All landmarks and curves", line=1, cex.lab=1.5)
abline(h = 0, v = 0, col ="grey", lty = 2)

# Squirrel related clade
plot(pcac$x[which(d$Clade == "SQUIR"),1:2],pch = 21, bg = 3, ylab="", xlab="", main = "Squirrel-related", xaxt = "n", yaxt = "n", xlim = c(min(pcac$x[,1]), max(pcac$x[,1])), ylim = c(min(pcac$x[,2]), max(pcac$x[,2])))
text(pcac$x[which(d$Clade == "SQUIR"),1:2], labels = d$Genre[which(d$Clade == "SQUIR")], cex = 0.5, pos = 2)
title(ylab="Cranium landmarks", line=1, cex.lab=1.5)
abline(h = 0, v = 0, col ="grey", lty = 2)

plot(pcamd$x[which(ddd$Clade == "SQUIR"),1:2],pch = 21, bg = 3, ylab="", xlab="", main = "Squirrel-related", xaxt = "n", yaxt = "n", xlim = c(min(pcamd$x[,1]), max(pcamd$x[,1])), ylim = c(min(pcamd$x[,2]), max(pcamd$x[,2])))
text(pcamd$x[which(ddd$Clade == "SQUIR"),1:2], labels = ddd$Genre[which(ddd$Clade == "SQUIR")], cex = 0.5, pos = 2)
title(ylab="Mandible landmarks", line=1, cex.lab=1.5)
abline(h = 0, v = 0, col ="grey", lty = 2)

plot(pcaM$x[which(dd$Clade == "SQUIR"),1:2], pch = 21, bg = 3, ylab="", xlab="", main = "Squirrel-related", xaxt = "n", yaxt = "n", xlim = c(min(pcaM$x[,1]), max(pcaM$x[,1])), ylim = c(min(pcaM$x[,2]), max(pcaM$x[,2])))
text(pcaM$x[which(dd$Clade == "SQUIR"),1:2], labels = dd$Genre[which(dd$Clade == "SQUIR")], cex = 0.5, pos = 2)
title(ylab="All landmarks and curves", line=1, cex.lab=1.5)
abline(h = 0, v = 0, col ="grey", lty = 2)

dev.off()
}
########################################
#Make plot coloring PCAs group by group#
########################################

if (print.plots) {
pdf(file = paste("plots", "combin_PCA_bygroup.pdf", sep = "/"), width = 7, height = 7)
layout(matrix(1:9, ncol = 3, byrow=T))
par(mar = c(1,3,1,0.2))

#Cranium LM row
plot(pcac$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(ylab="Cranium landmarks", line=1, cex.lab=1.5)
points(pcac$x[which(d$Clade == "CTENO"),1:2], cex = 1.5, pch = 21, bg = alpha(1, 0.7))
plot(pcac$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
points(pcac$x[which(d$Clade == "MOUS"),1:2], cex = 1.5, pch = 21, bg = alpha(2, 0.7))
plot(pcac$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
points(pcac$x[which(d$Clade == "SQUIR"),1:2], cex = 1.5, pch = 21, bg = alpha(3, 0.7))

#Mandible LM row
plot(pcamd$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(ylab="Mandible landmarks", line=1, cex.lab=1.5)
points(pcamd$x[which(ddd$Clade == "CTENO"),1:2], cex = 1.5, pch = 21, bg = alpha(1, 0.7))
plot(pcamd$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
points(pcamd$x[which(ddd$Clade == "MOUS"),1:2], cex = 1.5, pch = 21, bg = alpha(2, 0.7))
plot(pcamd$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
points(pcamd$x[which(ddd$Clade == "SQUIR"),1:2], cex = 1.5, pch = 21, bg = alpha(3, 0.7))

#LM and curves row
plot(pcaM$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(ylab="All curves and landmarks", line=1, cex.lab=1.5)
points(pcaM$x[which(dd$Clade == "CTENO"),1:2], cex = 1.5, pch = 21, bg = alpha(1, 0.7))
plot(pcaM$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
points(pcaM$x[which(dd$Clade == "MOUS"),1:2], cex = 1.5, pch = 21, bg = alpha(2, 0.7))
plot(pcaM$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
points(pcaM$x[which(dd$Clade == "SQUIR"),1:2], cex = 1.5, pch = 21, bg = alpha(3, 0.7))
dev.off()
}

##############################################
#Morphospace occupation by rodent morphotypes#
##############################################

if (print.plots) {
pdf(file = paste("plots", "pca_craniumLM_morphotypes.pdf", sep = "/"), width = 10, height = 10)
plot(pcac$x[,1:2], asp = 1, pch = 21, bg = c(4:7)[as.factor(d$MORPHOTYPE)])
text(pcac$x[,1:2], labels = specimens$genus, cex = 0.5, pos = 2)
dev.off()

pdf(file = paste("plots", "pca_allLMandcurves_morphotypes.pdf", sep = "/"), width = 10, height = 10)
plot(pcaM$x[,1:2], asp = 1, pch = 21, bg = c(4:7)[as.factor(dd$MORPHOTYPE)])
text(pcaM$x[,1:2], labels = spec$genus, cex = 0.5, pos = 2)
dev.off()

pdf(file = paste("plots", "pca_mandibleLM_morphotypes.pdf", sep = "/"), width = 10, height = 10)
plot(pcamd$x[,1:2], asp = 1, pch = 21, bg = c(4:7)[as.factor(ddd$MORPHOTYPE)])
text(pcamd$x[,1:2], labels = spem$genus, cex = 0.5, pos = 2)
dev.off()
}

##################################################
#Make plot coloring PCAs morphotype by morphotype#
##################################################

if (print.plots) {
pdf(file = paste("plots", "combin_PCA_morphotype.pdf", sep = "/"), width = 10, height = 8)
layout(matrix(1:12, ncol = 4, byrow=T))
par(mar = c(1,3,3,0.2))

f <- as.factor(d$Clade)

#Cranium LM row
plot(pcac$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(main = "Hystricomorphous", ylab="Cranium landmarks", line=1, cex.lab=1.5)
hys <- which(d$MORPHOTYPE == "HYSTRI")
points(pcac$x[hys,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[hys]])
mis <- which(d$Clade != "CTENO" & d$MORPHOTYPE == "HYSTRI")
text(pcac$x[mis,1:2], labels = d$Genre[mis], pos = 2)

plot(pcac$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(main = "Myomorphous", line=1, cex.lab=1.5)
myo <- which(d$MORPHOTYPE == "MYOM")
points(pcac$x[myo,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[myo]])
mis <- which(d$Clade != "MOUS" & d$MORPHOTYPE == "MYOM")
text(pcac$x[mis,1:2], labels = d$Genre[mis], pos = 2)

plot(pcac$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(main = "Protrogomorphous", line=1, cex.lab=1.5)
pro <- which(d$MORPHOTYPE == "PROTRO")
points(pcac$x[pro,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[pro]])
text(pcac$x[pro,1:2], labels = d$Genre[pro], pos = 2)

plot(pcac$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(main = "Sciuromorphous", line=1, cex.lab=1.5)
sci <- which(d$MORPHOTYPE == "SCIURO")
points(pcac$x[sci,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[sci]])
mis <- which(d$Clade != "SQUIR" & d$MORPHOTYPE == "SCIURO")
text(pcac$x[mis,1:2], labels = d$Genre[mis], pos = 2)

#Mandible LM row
f <- as.factor(ddd$Clade)

plot(pcamd$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(ylab="Mandible landmarks", line=1, cex.lab=1.5)
hys <- which(ddd$MORPHOTYPE == "HYSTRI")
points(pcamd$x[hys,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[hys]])
mis <- which(ddd$Clade != "CTENO" & ddd$MORPHOTYPE == "HYSTRI")
text(pcamd$x[mis,1:2], labels = ddd$Genre[mis], pos = 2)

plot(pcamd$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
myo <- which(ddd$MORPHOTYPE == "MYOM")
points(pcamd$x[myo,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[myo]])
mis <- which(ddd$Clade != "MOUS" & ddd$MORPHOTYPE == "MYOM")
text(pcamd$x[mis,1:2], labels = ddd$Genre[mis], pos = 2)

plot(pcamd$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
pro <- which(ddd$MORPHOTYPE == "PROTRO")
points(pcamd$x[pro,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[pro]])
text(pcamd$x[pro,1:2], labels = ddd$Genre[pro], pos = 2)

plot(pcamd$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
sci <- which(ddd$MORPHOTYPE == "SCIURO")
points(pcamd$x[sci,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[sci]])
mis <- which(ddd$Clade != "SQUIR" & ddd$MORPHOTYPE == "SCIURO")
text(pcamd$x[mis,1:2], labels = ddd$Genre[mis], pos = 2)

#LM and curves row
f <- as.factor(dd$Clade)

plot(pcaM$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(ylab="All curves and landmarks", line=1, cex.lab=1.5)
hys <- which(dd$MORPHOTYPE == "HYSTRI")
points(pcaM$x[hys,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[hys]])
mis <- which(dd$Clade != "CTENO" & dd$MORPHOTYPE == "HYSTRI")
text(pcaM$x[mis,1:2], labels = dd$Genre[mis], pos = 2)

plot(pcaM$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
myo <- which(dd$MORPHOTYPE == "MYOM")
points(pcaM$x[myo,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[myo]])
mis <- which(dd$Clade != "MOUS" & dd$MORPHOTYPE == "MYOM")
text(pcaM$x[mis,1:2], labels = dd$Genre[mis], pos = 2)

plot(pcaM$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
pro <- which(dd$MORPHOTYPE == "PROTRO")
points(pcaM$x[pro,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[pro]])
text(pcaM$x[pro,1:2], labels = dd$Genre[pro], pos = 2)

plot(pcaM$x[,1:2], asp = 1, pch = 20, cex = 1.2, col = alpha("gray", 0.5), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
sci <- which(dd$MORPHOTYPE == "SCIURO")
points(pcaM$x[sci,1:2], cex = 2, pch = 21, bg = alpha(1:3, 0.7)[f[sci]])
mis <- which(dd$Clade != "SQUIR" & dd$MORPHOTYPE == "SCIURO")
text(pcaM$x[mis,1:2], labels = dd$Genre[mis], pos = 2)

dev.off()
}
#######################################
#Disparity calculation and comparisons#
#######################################

#Apply functions
SORC <- SOR.fac(pcac$x, fac = as.factor(d$Clade), bootstrap = T, it = 1000)
SOVC <- SOV.fac(pcac$x, fac = as.factor(d$Clade), bootstrap = T, it = 1000)

SORmd <- SOR.fac(pcamd$x, fac = as.factor(ddd$Clade), bootstrap = T, it = 1000)
SOVmd <- SOV.fac(pcamd$x, fac = as.factor(ddd$Clade), bootstrap = T, it = 1000)
    
SORM <- SOR.fac(pcaM$x, fac = as.factor(dd$Clade), bootstrap = T, it = 1000)
SOVM <- SOV.fac(pcaM$x, fac = as.factor(dd$Clade), bootstrap = T, it = 1000)

if (print.plots) {
pdf(file = paste("plots", "disparity_boxplots.pdf", sep = "/"), width = 6, height = 6)
layout(matrix(1:6, ncol=3, byrow=T))
par(mar=c(2,2,3,0.2))

boxplot(SOVC$bsSOV, col = alpha(1:3,0.7), names = names(SOVC$SOV), yaxt ="n")
title(main = "Cranium LM", ylab = "Sum of variances", line = 1)
points(1:3, SOVC$SOV, pch = 20, col = "black")

boxplot(SOVmd$bsSOV, col = alpha(1:3,0.7), names = names(SOVmd$SOV), yaxt ="n")
title(main = "Mandible LM", line = 1)
points(1:3, SOVmd$SOV, pch = 20, col = "black")

boxplot(SOVM$bsSOV, col = alpha(1:3,0.7), names = names(SOVM$SOV), yaxt ="n")
title(main = "All curves and LM", line = 1)
points(1:3, SOVM$SOV, pch = 20, col = "black")

boxplot(SORC$bsSOR, col = alpha(1:3,0.7), names = names(SORC$SOR), yaxt ="n")
title(main = "Cranium LM", ylab = "Sum of Ranges", line = 1)
points(1:3, SORC$SOR, pch = 20, col = "black")

boxplot(SORmd$bsSOR, col = alpha(1:3,0.7), names = names(SORmd$SOR), yaxt ="n")
title(main = "Mandible LM", line = 1)
points(1:3, SORmd$SOR, pch = 20, col = "black")

boxplot(SORM$bsSOR, col = alpha(1:3,0.7), names = names(SORM$SOR), yaxt ="n")
title(main = "All curves and LM", line = 1)
points(1:3, SORM$SOR, pch = 20, col = "black")

dev.off()
}

############################################################################
# Integration and modularity patterns through correlation/congruence matrix#
############################################################################

#Average at the generic level
Mc_gen <- apply(X = Mc, 2, FUN = tapply, INDEX = d$Genre, mean)
M_gen <- apply(X = M, 2, FUN = tapply, INDEX = dd$Genre, mean)
Mmd_gen <- apply(X = Mmd, 2, FUN = tapply, INDEX = ddd$Genre, mean)

#Compute correlation matrix
correl <- cor(M_gen)
#plotting optional (not as readable as congruence matrix)
#image(abs(correl))

#Compute congruence matrix
A_gen <- array(NA, dim = c(dim(M_gen)[2]/3, 3, dim(M_gen)[1]))
for (i in 1:dim(A_gen)[3]) {A_gen[,,i] <- matrix(M_gen[i,], nrow = dim(M_gen)[2]/3, ncol=3, byrow=T)}

congru <- dotcorr(A_gen)

if (print.plots) {
layout(1:2)
#Plot congruence matrix with all landmarks and curves
par(mar = c(6.5,6.5,1,1))
image(abs(congru), xaxt = "n", yaxt = "n", col = hcl.colors(n=12, palette = "viridis"))
limits <- c(20.5, 45.5, 70.5, 95.5, 120.5, 145.5)
abline(h = 1/dim(A_gen)[1]*limits, v = 1/dim(A_gen)[1]*limits)
tck <- c(10, 32.5, 57.5, 82.5, 107.5, 132.5, 151)
labs <- c("Cranium", "Vault", "L Dorsal Zyg", "L Ventral Zyg","R Dorsal Zyg", "R Ventral Zyg", "Mandible")
axis(side = 1, at = 1/dim(A_gen)[1]*tck, labels = labs, las=2)
axis(side = 2, at = 1/dim(A_gen)[1]*tck, labels = labs, las=2)

#Plot only cranium and mandible LMs
par(mar = c(6.5,6.5,1,1))
image(abs(congru[c(1:20,145:157), c(1:20,145:157)]), xaxt = "n", yaxt = "n", col = hcl.colors(n=12, palette = "viridis"))
limits <- c(20.5)
abline(h = 1/32*limits, v = 1/32*limits)
tck <- c(10, 26)
labs <- c("Cranium", "Mandible")
axis(side = 1, at = 1/32*tck, labels = labs, las=2)
axis(side = 2, at = 1/32*tck, labels = labs, las=2)
}

#Compute phylogenetic independent contrasts

prnd.M <- data.pruning(tree=gen_tre, dat=M_gen, split.t="_", split.d=" ", positions=1)

M_gen_pics <- pics.df(prnd.M)

Apic <- array(NA, dim = c(dim(M_gen_pics)[2]/3, 3, dim(M_gen_pics)[1]))
for (i in 1:dim(Apic)[3]) {Apic[,,i] <- matrix(M_gen_pics[i,], nrow = dim(M_gen_pics)[2]/3, ncol=3, byrow=T)}

congru_pic <- dotcorr(Apic)

#Matrix of DIFFERENCE between congruence and PICs congruence
difcong <- abs(congru)-abs(congru_pic)

#Plot congruence matrix with all landmarks and curves and its PIC equivalent, and the difference matrix
if (print.plots) {
pdf(file = paste("plots/", "congruence_matrices_integration.pdf"), height = 7, width = 20)

layout(matrix(1:3, ncol=3))
par(mar = c(6.5,6.5,3,0))

image(abs(congru), xaxt = "n", yaxt = "n", main = "Landmark correlation matrix", col = hcl.colors(n=12, palette = "viridis"))
limits <- c(19.5, 44.5, 70, 95, 120.5, 145.5)
abline(h = 1/dim(Apic)[1]*limits, v = 1/dim(Apic)[1]*limits)
tck <- c(10, 32.5, 57.5, 82.5, 107.5, 132.5, 151)
labs <- c("Cranium", "Vault", "L Dorsal Zyg", "L Ventral Zyg","R Dorsal Zyg", "R Ventral Zyg", "Mandible")
axis(side = 1, at = 1/dim(Apic)[1]*tck, labels = labs, las=2)
axis(side = 2, at = 1/dim(Apic)[1]*tck, labels = labs, las=2)

par(mar = c(6.5,5.5,3,0))
image(abs(congru_pic), xaxt = "n", yaxt = "n", main = "PICs correlation matrix", col = hcl.colors(n=12, palette = "viridis"))
limits <- c(19.5, 44.5, 70, 95, 120.5, 145.5)
abline(h = 1/dim(Apic)[1]*limits, v = 1/dim(Apic)[1]*limits)
tck <- c(10, 32.5, 57.5, 82.5, 107.5, 132.5, 151)
labs <- c("Cranium", "Vault", "L Dorsal Zyg", "L Ventral Zyg","R Dorsal Zyg", "R Ventral Zyg", "Mandible")
axis(side = 1, at = 1/dim(Apic)[1]*tck, labels = labs, las=2)

par(mar = c(6.5,5.5,3,1))
image(z=difcong, col=hcl.colors(n=12, "viridis"), xaxt = "n", yaxt = "n", main = "Correlation difference matrix")
limits <- c(19.5, 44.5, 70, 95, 120.5, 145.5)
abline(h = 1/dim(Apic)[1]*limits, v = 1/dim(Apic)[1]*limits)
tck <- c(10, 32.5, 57.5, 82.5, 107.5, 132.5, 151)
labs <- c("Cranium", "Vault", "L Dorsal Zyg", "L Ventral Zyg","R Dorsal Zyg", "R Ventral Zyg", "Mandible")
axis(side = 1, at = 1/dim(Apic)[1]*tck, labels = labs, las=2)

dev.off()
}

######################################
#(Phylo)morphospaces and deformations#
######################################
linmd <- c(1:8,11,1,NA, 7, 12, NA, 8, 9)
linmand <- c(146:153,156,146,NA, 152, 157, NA,153,154)
curcran <- c(21:45,NA,46:70,NA,71:95,NA,96:120,NA,121:145)
lincrlat <- c(1:3,NA,17:15,12,5,20,1)
lincrandors <- c(1,3,18,45,19,4,1, NA, 6,8,NA, 7,9, NA, 10,13,NA, 11,14)
lincranlat <- c(1:3,21:45,17:15,12,5,20,1, NA, 46:70,NA,71:95)

# 1. Cranium, landmarks only

o <- cbind(Mc,d)

prnd.cranium <- array.pruning(gen_tre, o, A = pAcc$rotated, split.t = "_", split.d = "_", positions = 1)

data <- prnd.cranium$data[,1:60]
metadata <- prnd.cranium$data[,61:65]
tree <- prnd.cranium$tree
fac <- "Clade"

A <- prnd.cranium$array
mshp <- mshape(A)

dimnames(A)[[3]] <- tree$tip.label
pca.geomorph <- gm.prcomp(A= A, phy = tree)
phypca.geomorph <- gm.prcomp(A= A, phy = tree, align.to.phy = T)

if (print.plots) {

pdf(file=paste("plots/","cranium_morphospace_and_phylospace.pdf", sep=""), height=10, width=10)
plotnames <- T
layout(matrix(1:4,ncol=2))
par(mar = c(4,4,0.5,0.5))

plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F,anc.states=F), cex=0.1)
points(pca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
pca.geomorph$x[,1] > 0.8 * max(pca.geomorph$x[,1]) | 
pca.geomorph$x[,1] < 0.5 * min(pca.geomorph$x[,1])) | 
(pca.geomorph$x[,2] > 0.5 * max(pca.geomorph$x[,2]) | 
pca.geomorph$x[,2] < 0.5 * min(pca.geomorph$x[,2])))
if (plotnames) {text(pca.geomorph$x[namtoplot,1:2], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

plot(phypca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F,anc.states=F), cex=0.1)
points(phypca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
phypca.geomorph$x[,1] > 0.5 * max(phypca.geomorph$x[,1]) | 
phypca.geomorph$x[,1] < 0.8 * min(phypca.geomorph$x[,1])) | 
(phypca.geomorph$x[,2] > 0.5 * max(phypca.geomorph$x[,2]) | 
phypca.geomorph$x[,2] < 0.8 * min(phypca.geomorph$x[,2])))
if (plotnames) {text(phypca.geomorph$x[namtoplot,1:2], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), axis1=3, axis2=4, cex=0.1)
points(pca.geomorph$x[,3:4], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
pca.geomorph$x[,3] > 0.5 * max(pca.geomorph$x[,3]) | 
pca.geomorph$x[,3] < 0.5 * min(pca.geomorph$x[,3])) | 
(pca.geomorph$x[,4] > 0.8 * max(pca.geomorph$x[,4]) | 
pca.geomorph$x[,4] < 0.4 * min(pca.geomorph$x[,4])))
if (plotnames) {text(pca.geomorph$x[namtoplot,3:4], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

plot(phypca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), axis1=3, axis2=4, cex = 0.1)
points(phypca.geomorph$x[,3:4], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
phypca.geomorph$x[,3] > 0.5 * max(phypca.geomorph$x[,3]) | 
phypca.geomorph$x[,3] < 0.5 * min(phypca.geomorph$x[,3])) | 
(phypca.geomorph$x[,4] > 0.5 * max(phypca.geomorph$x[,4]) | 
phypca.geomorph$x[,4] < 0.5 * min(phypca.geomorph$x[,4])))
if (plotnames) {text(phypca.geomorph$x[namtoplot,3:4], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

dev.off()
}

if (print.plots) {
pdf(paste("plots/", "cranium_space_with_deform.pdf"), width = 8, height =6)
layout(matrix(c(1,3,3,2,3,3,0,4,5), ncol=3, byrow=T))

par(mar = c(0,2,0,2))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 2)
plot(mshp[,c(3,1)], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = c(shp$max[i,1], mshp[i,1]), col=alpha("red", 0.5), lwd=3)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = c(shp$min[i,1], mshp[i,1]), col=alpha("blue", 0.5), lwd=3)}


plot(mshp[,3], -mshp[,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[lincrlat,3], -mshp[lincrlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = -c(shp$max[i,2], mshp[i,2]), col=alpha("red", 0.5), lwd=3)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = -c(shp$min[i,2], mshp[i,2]), col=alpha("blue", 0.5), lwd=3)}


par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), cex = 0.1)
points(pca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)

par(mar = c(0,1,0,1))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 1)
plot(mshp[,c(3,1)], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = c(shp$max[i,1], mshp[i,1]), col=alpha("red", 0.5), lwd=3)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = c(shp$min[i,1], mshp[i,1]), col=alpha("blue", 0.5), lwd=3)}

plot(mshp[,3], -mshp[,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[lincrlat,3], -mshp[lincrlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = -c(shp$max[i,2], mshp[i,2]), col=alpha("red", 0.5), lwd=3)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = -c(shp$min[i,2], mshp[i,2]), col=alpha("blue", 0.5), lwd=3)}

dev.off()
}

# 2. Cranium: landmarks and curves

o <- cbind(MC,dd)

prnd.C <- array.pruning(gen_tre, o, A = pAC$rotated, split.t = "_", split.d = "_", positions = 1)

data <- prnd.C$data[,1:435]
metadata <- prnd.C$data[,436:440]
tree <- prnd.C$tree
fac <- "Clade"

A <- prnd.C$array
mshp <- mshape(A)

dimnames(A)[[3]] <- tree$tip.label
pca.geomorph <- gm.prcomp(A= A, phy = tree)
phypca.geomorph <- gm.prcomp(A= A, phy = tree, align.to.phy = T)

if (print.plots) {

pdf(file=paste("plots/","cranium_curves_morphospace_and_phylospace.pdf", sep=""), height=10, width=10)
plotnames <- T
layout(matrix(1:4,ncol=2))

par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), cex = 0.1)
points(pca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
pca.geomorph$x[,1] > 0.5 * max(pca.geomorph$x[,1]) | 
pca.geomorph$x[,1] < 0.7 * min(pca.geomorph$x[,1])) | 
(pca.geomorph$x[,2] > 0.7 * max(pca.geomorph$x[,2]) | 
pca.geomorph$x[,2] < 0.3 * min(pca.geomorph$x[,2])))
if (plotnames) {text(pca.geomorph$x[namtoplot,1:2], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

plot(phypca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F,anc.states=F), cex=0.1)
points(phypca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
phypca.geomorph$x[,1] > 0.9 * max(phypca.geomorph$x[,1]) | 
phypca.geomorph$x[,1] < 0.1 * min(phypca.geomorph$x[,1])) | 
(phypca.geomorph$x[,2] > 0.3 * max(phypca.geomorph$x[,2]) | 
phypca.geomorph$x[,2] < 0.3 * min(phypca.geomorph$x[,2])))
if (plotnames) {text(phypca.geomorph$x[namtoplot,1:2], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), axis1=3, axis2=4, cex = 0.1)
points(pca.geomorph$x[,3:4], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
pca.geomorph$x[,3] > 0.3 * max(pca.geomorph$x[,3]) | 
pca.geomorph$x[,3] < 0.8 * min(pca.geomorph$x[,3])) | 
(pca.geomorph$x[,4] > 0.7 * max(pca.geomorph$x[,4]) | 
pca.geomorph$x[,4] < 0.5 * min(pca.geomorph$x[,4])))
if (plotnames) {text(pca.geomorph$x[namtoplot,3:4], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

plot(phypca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), axis1=3, axis2=4, cex = 0.1)
points(phypca.geomorph$x[,3:4], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
phypca.geomorph$x[,3] > 0.2 * max(phypca.geomorph$x[,3]) | 
phypca.geomorph$x[,3] < 0.6 * min(phypca.geomorph$x[,3])) | 
(phypca.geomorph$x[,4] > 0.3 * max(phypca.geomorph$x[,4]) | 
phypca.geomorph$x[,4] < 0.3 * min(phypca.geomorph$x[,4])))
if (plotnames) {text(phypca.geomorph$x[namtoplot,3:4], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

dev.off()
}

if (print.plots) {
pdf(paste("plots/", "cranium_curves_space_with_deform.pdf"), width = 8, height =6)
layout(matrix(c(1,3,3,2,3,3,0,4,5), ncol=3, byrow=T))

par(mar = c(0,2,0,2))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 2)
plot(mshp[,c(3,1)], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = c(shp$max[i,1], mshp[i,1]), col=alpha("red", 0.5), lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = c(shp$min[i,1], mshp[i,1]), col=alpha("blue", 0.5), lwd=2)}


plot(mshp[,3], -mshp[,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[lincranlat,3], -mshp[lincranlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = -c(shp$max[i,2], mshp[i,2]), col=alpha("red", 0.5), lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = -c(shp$min[i,2], mshp[i,2]), col=alpha("blue", 0.5), lwd=2)}


par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), cex = 0.1)
points(pca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)

par(mar = c(0,1,0,1))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 1)
plot(mshp[,c(3,1)], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = c(shp$max[i,1], mshp[i,1]), col=alpha("red", 0.5), lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = c(shp$min[i,1], mshp[i,1]), col=alpha("blue", 0.5), lwd=2)}

plot(mshp[,3], -mshp[,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[lincranlat,3], -mshp[lincranlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = -c(shp$max[i,2], mshp[i,2]), col=alpha("red", 0.5), lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = -c(shp$min[i,2], mshp[i,2]), col=alpha("blue", 0.5), lwd=2)}

dev.off()
}

# 3. Mandible space

o <- cbind(Mmd,ddd)

prnd.md <- array.pruning(gen_tre, o, A = pAmc$rotated, split.t = "_", split.d = "_", positions = 1)

data <- prnd.md$data[,1:36]
metadata <- prnd.md$data[,37:41]
tree <- prnd.md$tree
fac <- "Clade"

A <- prnd.md$array
mshp <- mshape(A)

dimnames(A)[[3]] <- tree$tip.label
pca.geomorph <- gm.prcomp(A= A, phy = tree)
phypca.geomorph <- gm.prcomp(A= A, phy = tree, align.to.phy = T)

if (print.plots) {
pdf(file=paste("plots/","mandible_morphospace_and_phylospace.pdf", sep=""), height=10, width=10)
plotnames <- T
layout(matrix(1:4,ncol=2))

par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), cex = 0.1)
points(pca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
pca.geomorph$x[,1] > 0.1 * max(pca.geomorph$x[,1]) | 
pca.geomorph$x[,1] < 0.8 * min(pca.geomorph$x[,1])) | 
(pca.geomorph$x[,2] > 0.3 * max(pca.geomorph$x[,2]) | 
pca.geomorph$x[,2] < 0.5 * min(pca.geomorph$x[,2])))
if (plotnames) {text(pca.geomorph$x[namtoplot,1:2], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

plot(phypca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F,anc.states=F), cex=0.1)
points(phypca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
phypca.geomorph$x[,1] > 0.9 * max(phypca.geomorph$x[,1]) | 
phypca.geomorph$x[,1] < 0.7 * min(phypca.geomorph$x[,1])) | 
(phypca.geomorph$x[,2] > 0.8 * max(phypca.geomorph$x[,2]) | 
phypca.geomorph$x[,2] < 0.8 * min(phypca.geomorph$x[,2])))
if (plotnames) {text(phypca.geomorph$x[namtoplot,1:2], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), axis1=3, axis2=4, cex = 0.1)
points(pca.geomorph$x[,3:4], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
pca.geomorph$x[,3] > 0.5 * max(pca.geomorph$x[,3]) | 
pca.geomorph$x[,3] < 0.6 * min(pca.geomorph$x[,3])) | 
(pca.geomorph$x[,4] > 0.9 * max(pca.geomorph$x[,4]) | 
pca.geomorph$x[,4] < 0.7 * min(pca.geomorph$x[,4])))
if (plotnames) {text(pca.geomorph$x[namtoplot,3:4], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

plot(phypca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), axis1=3, axis2=4, cex = 0.1)
points(phypca.geomorph$x[,3:4], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
phypca.geomorph$x[,3] > 0.4 * max(phypca.geomorph$x[,3]) | 
phypca.geomorph$x[,3] < 0.9 * min(phypca.geomorph$x[,3])) | 
(phypca.geomorph$x[,4] > 0.4 * max(phypca.geomorph$x[,4]) | 
phypca.geomorph$x[,4] < 0.5 * min(phypca.geomorph$x[,4])))
if (plotnames) {text(phypca.geomorph$x[namtoplot,3:4], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

dev.off()
}

if (print.plots) {
pdf(paste("plots/", "mandible_space_with_deform.pdf"), width = 8, height =6)
layout(matrix(c(1,3,3,2,3,3,0,4,5), ncol=3, byrow=T))

par(mar = c(0,2,0,2))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 2)
plot(mshp[,1], -mshp[,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[linmd,1], -mshp[linmd,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,1], mshp[i,1]), y = -c(shp$max[i,2], mshp[i,2]), col=alpha("red", 0.5), lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,1], mshp[i,1]), y = -c(shp$min[i,2], mshp[i,2]), col=alpha("blue", 0.5), lwd=2)}


plot(mshp[,1], mshp[,3], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[linmd,1], mshp[linmd,3], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,1], mshp[i,1]), y = c(shp$max[i,3], mshp[i,3]), col=alpha("red", 0.5), lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,1], mshp[i,1]), y = c(shp$min[i,3], mshp[i,3]), col=alpha("blue", 0.5), lwd=2)}


par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), cex = 0.1)
points(pca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)

par(mar = c(0,1,0,1))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 1)
plot(mshp[,1], -mshp[,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[linmd,1], -mshp[linmd,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,1], mshp[i,1]), y = -c(shp$max[i,2], mshp[i,2]), col=alpha("red", 0.5), lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,1], mshp[i,1]), y = -c(shp$min[i,2], mshp[i,2]), col=alpha("blue", 0.5), lwd=2)}

plot(mshp[,1], mshp[,3], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5))
lines(mshp[linmd,1], mshp[linmd,3], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,1], mshp[i,1]), y = c(shp$max[i,3], mshp[i,3]), col=alpha("red", 0.5), lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,1], mshp[i,1]), y = c(shp$min[i,3], mshp[i,3]), col=alpha("blue", 0.5), lwd=2)}

dev.off()
}


# 4. All landmarks and curves (combine cranium and mandible)

o <- cbind(M,dd)
a <- array(NA, dim = c(dim(M)[2]/3, 3, dim(M)[1]))
for (i in 1:dim(M)[1]) {a[,,i] <- matrix(M[i,], ncol=3, byrow=T)}
dimnames(a)[[3]] <- rownames(o)

prnd.all <- array.pruning(gen_tre, o, A = a, split.t = "_", split.d = "_", positions = 1)

data <- prnd.all$data[,1:471]
metadata <- prnd.all$data[,472:476]
tree <- prnd.all$tree
fac <- "Clade"

A <- prnd.all$array
mshp <- mshape(A)

dimnames(A)[[3]] <- tree$tip.label
pca.geomorph <- gm.prcomp(A= A, phy = tree)
phypca.geomorph <- gm.prcomp(A= A, phy = tree, align.to.phy = T)


if (print.plots) {
pdf(file=paste("plots/","all_skull_morphospace_and_phylospace.pdf", sep=""), height=10, width=10)
plotnames <- T
layout(matrix(1:4,ncol=2))

par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), cex = 0.1)
points(pca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
pca.geomorph$x[,1] > 0.7 * max(pca.geomorph$x[,1]) | 
pca.geomorph$x[,1] < 0.6 * min(pca.geomorph$x[,1])) | 
(pca.geomorph$x[,2] > 0.7 * max(pca.geomorph$x[,2]) | 
pca.geomorph$x[,2] < 0.6 * min(pca.geomorph$x[,2])))
if (plotnames) {text(pca.geomorph$x[namtoplot,1:2], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

plot(phypca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F,anc.states=F), cex=0.1)
points(phypca.geomorph$x[,1:2], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
phypca.geomorph$x[,1] > 0.85 * max(phypca.geomorph$x[,1]) | 
phypca.geomorph$x[,1] < 0.7 * min(phypca.geomorph$x[,1])) | 
(phypca.geomorph$x[,2] > 0.5 * max(phypca.geomorph$x[,2]) | 
phypca.geomorph$x[,2] < 0.5 * min(phypca.geomorph$x[,2])))
if (plotnames) {text(phypca.geomorph$x[namtoplot,1:2], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), axis1=3, axis2=4, cex = 0.1)
points(pca.geomorph$x[,3:4], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
pca.geomorph$x[,3] > 0.5 * max(pca.geomorph$x[,3]) | 
pca.geomorph$x[,3] < 0.7 * min(pca.geomorph$x[,3])) | 
(pca.geomorph$x[,4] > 0.5 * max(pca.geomorph$x[,4]) | 
pca.geomorph$x[,4] < 0.7 * min(pca.geomorph$x[,4])))
if (plotnames) {text(pca.geomorph$x[namtoplot,3:4], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

plot(phypca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F), axis1=3, axis2=4, cex = 0.1)
points(phypca.geomorph$x[,3:4], bg = c(1:3)[as.factor(metadata$Clade)], pch=21, cex = 2)
namtoplot <- which((
phypca.geomorph$x[,3] > 0.2 * max(phypca.geomorph$x[,3]) | 
phypca.geomorph$x[,3] < 0.7 * min(phypca.geomorph$x[,3])) | 
(phypca.geomorph$x[,4] > 0.4 * max(phypca.geomorph$x[,4]) | 
phypca.geomorph$x[,4] < 0.5 * min(phypca.geomorph$x[,4])))
if (plotnames) {text(phypca.geomorph$x[namtoplot,3:4], labels = names(namtoplot), pos=c(1,3), cex = 0.7)}

dev.off()
}

if (print.plots) {
pdf(file=paste("plots/","deformation_PCA_allskull.pdf", sep=""), height=10, width=20)

layout(matrix(1:18, nrow = 3))
par(mar = c(0,0,0,0))

#PC1
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 1)

plot(mshp[1:145,3], mshp[1:145,1], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = c(shp$max[1:145,1][i], mshp[1:145,1][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = c(shp$min[1:145,1][i], mshp[1:145,1][i]), col=alpha("blue", 0.5), lwd = 2)}
title(main = "PC1", line = -1)

plot(mshp[1:145,3], -mshp[1:145,2], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincranlat,3], -mshp[lincranlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$max[1:145,2][i], mshp[1:145,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$min[1:145,2][i], mshp[1:145,2][i]), col=alpha("blue", 0.5), lwd = 2)}


plot(mshp[146:157,1], -mshp[146:157,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[linmand,1], -mshp[linmand,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$max[146:157,2][i], mshp[146:157,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$min[146:157,2][i], mshp[146:157,2][i]), col=alpha("blue", 0.5), lwd = 2)}

#PC2
shp <- shp.deform(shp = mshape(A), PCA = pca.geomorph, axis = 2)

plot(mshp[1:145,3], mshp[1:145,1], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = c(shp$max[1:145,1][i], mshp[1:145,1][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = c(shp$min[1:145,1][i], mshp[1:145,1][i]), col=alpha("blue", 0.5), lwd = 2)}
title(main = "PC2", line = -1)

plot(mshp[1:145,3], -mshp[1:145,2], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincranlat,3], -mshp[lincranlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$max[1:145,2][i], mshp[1:145,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$min[1:145,2][i], mshp[1:145,2][i]), col=alpha("blue", 0.5), lwd = 2)}


plot(mshp[146:157,1], -mshp[146:157,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[linmand,1], -mshp[linmand,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$max[146:157,2][i], mshp[146:157,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$min[146:157,2][i], mshp[146:157,2][i]), col=alpha("blue", 0.5), lwd = 2)}

#PC3
shp <- shp.deform(shp = mshape(A), PCA = pca.geomorph, axis = 3)

plot(mshp[1:145,3], mshp[1:145,1], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = c(shp$max[1:145,1][i], mshp[1:145,1][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = c(shp$min[1:145,1][i], mshp[1:145,1][i]), col=alpha("blue", 0.5), lwd = 2)}
title(main = "PC3", line = -1)

plot(mshp[1:145,3], -mshp[1:145,2], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincranlat,3], -mshp[lincranlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$max[1:145,2][i], mshp[1:145,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$min[1:145,2][i], mshp[1:145,2][i]), col=alpha("blue", 0.5), lwd = 2)}


plot(mshp[146:157,1], -mshp[146:157,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[linmand,1], -mshp[linmand,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$max[146:157,2][i], mshp[146:157,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$min[146:157,2][i], mshp[146:157,2][i]), col=alpha("blue", 0.5), lwd = 2)}

#PC4
shp <- shp.deform(shp = mshape(A), PCA = pca.geomorph, axis = 4)

plot(mshp[1:145,3], mshp[1:145,1], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = c(shp$max[1:145,1][i], mshp[1:145,1][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = c(shp$min[1:145,1][i], mshp[1:145,1][i]), col=alpha("blue", 0.5), lwd = 2)}
title(main = "PC4", line = -1)

plot(mshp[1:145,3], -mshp[1:145,2], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincranlat,3], -mshp[lincranlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$max[1:145,2][i], mshp[1:145,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$min[1:145,2][i], mshp[1:145,2][i]), col=alpha("blue", 0.5), lwd = 2)}


plot(mshp[146:157,1], -mshp[146:157,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[linmand,1], -mshp[linmand,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$max[146:157,2][i], mshp[146:157,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$min[146:157,2][i], mshp[146:157,2][i]), col=alpha("blue", 0.5), lwd = 2)}

#PC5
shp <- shp.deform(shp = mshape(A), PCA = pca.geomorph, axis = 5)

plot(mshp[1:145,3], mshp[1:145,1], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = c(shp$max[1:145,1][i], mshp[1:145,1][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = c(shp$min[1:145,1][i], mshp[1:145,1][i]), col=alpha("blue", 0.5), lwd = 2)}
title(main = "PC5", line = -1)

plot(mshp[1:145,3], -mshp[1:145,2], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincranlat,3], -mshp[lincranlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$max[1:145,2][i], mshp[1:145,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$min[1:145,2][i], mshp[1:145,2][i]), col=alpha("blue", 0.5), lwd = 2)}


plot(mshp[146:157,1], -mshp[146:157,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[linmand,1], -mshp[linmand,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$max[146:157,2][i], mshp[146:157,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$min[146:157,2][i], mshp[146:157,2][i]), col=alpha("blue", 0.5), lwd = 2)}

#PC6
shp <- shp.deform(shp = mshape(A), PCA = pca.geomorph, axis = 6)

plot(mshp[1:145,3], mshp[1:145,1], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = c(shp$max[1:145,1][i], mshp[1:145,1][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = c(shp$min[1:145,1][i], mshp[1:145,1][i]), col=alpha("blue", 0.5), lwd = 2)}
title(main = "PC6", line = -1)

plot(mshp[1:145,3], -mshp[1:145,2], asp = 1, bty = "n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[lincranlat,3], -mshp[lincranlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$max[1:145,2][i], mshp[1:145,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[1:145,3][i], mshp[1:145,3][i]), y = -c(shp$min[1:145,2][i], mshp[1:145,2][i]), col=alpha("blue", 0.5), lwd = 2)}


plot(mshp[146:157,1], -mshp[146:157,2], asp = 1, bty ="n", xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1, col = alpha("black",0.5), xlim = minetmax(unlist(shp)[-c(436:471, 907:942)]), ylim = minetmax(unlist(shp)[-c(436:471, 907:942)]))
lines(mshp[linmand,1], -mshp[linmand,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$max[146:157,2][i], mshp[146:157,2][i]), col=alpha("red", 0.5), lwd = 2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[146:157,1][i], mshp[146:157,1][i]), y = -c(shp$min[146:157,2][i], mshp[146:157,2][i]), col=alpha("blue", 0.5), lwd = 2)}

dev.off()
}




