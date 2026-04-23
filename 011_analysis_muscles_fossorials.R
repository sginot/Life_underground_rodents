source("Rfunctions1.txt")
source("CUSTOM_R_FUNCTIONS.R")
source("004_phylogeny_and_contrasts_specific.R")
source("007_import_landmarks.R")
source("009_GMM_morphospace_converg.R")
library(scales)
library(nlme)
library(convevol)
library(viridis)

# Opt to print plots or not:
print.plots <- T
if (print.plots & !dir.exists("plots")) {dir.create("plots")} #Create output folder if it does not exist

# Opt to plot figures with species names 
plot.names <- F

# Opt to sink all output to log file
sinking <- F
if (sinking) {sink(file = "output_fossorial_analyses.R", append = F, type = "output")}

# Define number of iterations of simulations to compute Stayton's convergence indices
NSIM <- 4

##################################
# Analyze the bite force dataset #
##################################

bfdat <- read.csv("Data bite force rodents.csv", sep = "\t", dec = ",", stringsAsFactors = T, na.strings = "")
#str(bfdat)

# OLS linear models
modbf1 <- lm(bfdat$logBF ~ bfdat$logmass)
modbf2 <- aov(bfdat$logBF ~ bfdat$logmass * bfdat$life.style.disparat)
modbf3 <- aov(bfdat$logBF ~ bfdat$logmass + bfdat$life.style.disparat)
modbf4 <- lm(bfdat$logBF ~ bfdat$logmass + bfdat$life.style.disparat)
resmod1 <- residuals(modbf1)

print("ANOVA Bite Force againt mass * lifestyle")
summary(modbf2)
print("ANOVA Bite Force againt mass + lifestyle")
summary(modbf3)
print("Linear regression Bite Force againt mass + lifestyle")
summary(modbf4)

# Simple biplot
cols1 <- viridis(4)

if (print.plots) {
pdf(paste("plots/", if (plot.names) {"named_"}, "biplot_mass_BF.pdf", sep = ""))
if (plot.names) {par(mar = c(4,4,1,8))}
plot(bfdat$logmass, bfdat$logBF, bg = alpha(cols1, 0.8)[bfdat$life.style.disparat], pch=c(21:22)[bfdat$chisel], cex = 2, xlab = "Log10 Body mass (g)", ylab = "Log10 Bite Force (N)")

if (plot.names) {
labels <- abbreviate(paste(bfdat$genus, bfdat$species, sep=" "))
text(bfdat$logmass, bfdat$logBF, labels = labels, cex=0.6, pos=1)
axis(side = 4, at = seq(from = min(bfdat$logBF), to = max(bfdat$logBF), length.out = length(labels)), labels = sort(paste(labels, names(labels), sep="-")), las = 2, cex.axis = 0.5, tick = T)
}

legend("topleft", pch = c(rep(21,nlevels(bfdat$life.style.disparat)), 22), pt.bg = c(cols1, NA), legend = c(levels(bfdat$life.style.disparat), "Chisel tooth digger"), pt.cex = 2, bty="n")
plot.lm.CI(bfdat$logmass, bfdat$logBF, alpha=0.35, border = F)
abline(modbf1, lwd = 2, col = "black")

par(fig = c(0.7,0.94,0.15,0.6), new = T, mar=c(0,2,0,0.1))
if (plot.names) {par(fig = c(0.54,0.75,0.13,0.54), new = T, mar=c(0,2,0,0.1))}
boxplot(resmod1 ~ bfdat$life.style.disparat, col = cols1, las = 2, yaxt = "n", xaxt="n", cex.lab = 0.5)
title(ylab = "Residual bite force", line = 0.5)
dev.off()
}

# Phylogenetic GLS models
rownames(bfdat) <- paste(bfdat$genus, bfdat$species, sep = " ")
bf.prnd <- data.pruning(tree = sp_tre, dat = bfdat, split.t = "_", split.d = " ", positions = c(1,2))

lfstl <- relevel(bf.prnd$pruned.dat$life.style.disparat, ref = "subterranean")
pglsbf <- gls(logBF ~ logmass + lfstl, data = bf.prnd$pruned.data, correlation = corBrownian(phy = bf.prnd$pruned.tree))

print("PGLS Bite force against mass + lifestyle, with subterranean as reference group")
summary(pglsbf)

# Post-hoc Pairwise Tukey test (taken from phytools blog)

library(multcomp)
model.matrix.gls<-function(object,...){
        model.matrix(terms(object),data=getData(object),...)
}
model.frame.gls<-function(object,...){
        model.frame(formula(object),data=getData(object),...)
}
terms.gls<-function(object,...){
        terms(model.frame(object),...)
}
post.hoc<-glht(pglsbf,linfct=mcp(lfstl="Tukey"))

print("Tukey posthoc pgls test of pairwise bite force differences between groups")
summary(post.hoc)
 
# Phenogram aka traitgram

x <- residuals(lm(bf.prnd$pruned.data$logBF ~ bf.prnd$pruned.data$logmass))
tree <- bf.prnd$pruned.tree
isbur <- rep("NO", length(bf.prnd$pruned.data$life.style.disparat))
isbur[bf.prnd$pruned.data$life.style.disparat == "subterranean"] <- "subterranean"
isbur[bf.prnd$pruned.data$life.style.disparat == "burrower"] <- "burrower"
names(x) <- names(isbur) <- bf.prnd$pruned.tree$tip.label
mtree <- make.simmap(tree = tree, x = isbur, message = F)
cols<-c(alpha("black", 0.1), cols1[2], cols1[3]); names(cols)<-c( "NO","burrower","subterranean")
 
if (print.plots) {
pdf(paste("plots/", "traitgram_residual_BF.pdf", sep = ""), height=8, width = 5)
par(mar = c(4,4,1,1))
phenogram(x = x, tree = mtree, colors = cols, fsize = 0.7, ylab = "Body mass-residual bite force")
dev.off()
}

# Convergence indices for residuals of linear regression of bite force against mass

foct <- names(x)[which(isbur != "NO")]

#print("Stayton C1-C4 convergence indices applied to body mass residual bite force") 
#cvsg <- convSig(phy = tree, traits = matrix(x), focaltaxa = foct, nsim = NSIM)
#print(cvsg)

p <- phylosig(tree = tree, x = x, method = "lambda")
print("Pagel's lambda applied to body mass residual bite force") 
print(p)

####################################################################################
# Create species level data frame containing size, muscle, and lifestyle variables #
####################################################################################

disparat <- read.csv("LISTE_IODINE_RODENT - Main_list.csv", h=T, dec=".", sep=",", na.strings = c("NA", "", "/"))
lifestyles <- disparat[,c("Clade","Genre", "espèces", "LIFESTYLE", "LIFESTYLE_bis", "Life.mode..Alhajeri.et.al..2020.", "DIET", "DIET.BROAD")]
colnames(lifestyles) <- c("Clade","Genre", "espèces", "LIFESTYLE_old", "LIFESTYLE", "Life.mode", "DIET", "DIET.BROAD")
#lifestyles[grep("fossorial", lifestyles$LIFESTYLE), "LIFESTYLE"] <- "terricolous"
lifestyles[grep("arboreal", lifestyles$LIFESTYLE), c("LIFESTYLE_old", "LIFESTYLE")] <- "arboreal"
lifestyles[grep("omnivorous \\(\\A\\?\\)", lifestyles$DIET.BROAD), "DIET.BROAD"] <- "omnivorous (A)"

#plot(gmean_mass[-attr(df_sum_pcsa, "na.action")], total_pcsa)
#identify(gmean_mass[-attr(df_sum_pcsa, "na.action")], total_pcsa)

sp_centroid_size <- data.frame(species = paste(specimens$genus, specimens$sp, sep=" "), size = sizAcc)
sp_avg_csize <- tapply(X=sp_centroid_size$size, FUN = mean, INDEX = sp_centroid_size$species)

index_lifestyles <- as.factor(paste(lifestyles$Genre, lifestyles$espèces, sep=" "))

sp_clade <- as.factor(tapply(X = lifestyles$Clade, FUN = unique, INDEX = index_lifestyles))
sp_lifestyle <- as.factor(tapply(X = lifestyles$LIFESTYLE, FUN = unique, INDEX = index_lifestyles))
sp_lifemode <- as.factor(tapply(X = lifestyles$Life.mode, FUN = unique, INDEX = index_lifestyles))
sp_diet <- as.factor(tapply(X = lifestyles$DIET, FUN = unique, INDEX = index_lifestyles))
sp_dietbroad <- as.factor(tapply(X = lifestyles$DIET.BROAD, FUN = unique, INDEX = index_lifestyles))

splfst <- sp_lifestyle
splfst[which(sp_lifestyle == "subterranean (C)")] <- "subterranean"
splfst[which(sp_lifestyle == "burrower (C)")] <- "burrower"
splfst[which(sp_lifestyle == "semiaquatic/burrower")] <- "burrower"
splfst <- droplevels(splfst)

mlifestyle <- data.frame(clade = sp_clade, lifestyle = sp_lifestyle, lifestyle_simple = splfst, lifemode = sp_lifemode, diet = sp_diet, diet.broad = sp_dietbroad)
ms <- as.matrix(sp_avg_csize)

sp_mtch <- mat.mtch(mlifestyle, ms)

#boxplot(sp_mtch$s2 ~ sp_mtch$lifestyle, las= 2, xlab="")

mscl_mtch <- vec.mtch(total_mass, total_pcsa)
index_mscl <- shorten.labels(rownames(mscl_mtch), positions = c(1,2), split=" ")

sp_avg_total_mass <- tapply(X = mscl_mtch$total_mass, INDEX = index_mscl, FUN = mean)
sp_avg_total_pcsa <- tapply(X = mscl_mtch$total_pcsa, INDEX = index_mscl, FUN = mean)

sp_mscl_mtch <- vec.mtch(sp_avg_total_mass, sp_avg_total_pcsa)

mat_mscl_lifestyle <- mat.mtch(sp_mtch,sp_mscl_mtch)
colnames(mat_mscl_lifestyle) <- c("clade","lifestyle", "lifestyle_simple", "lifemode", "diet", "diet.broad", "size", "mass", "pcsa")
mat_mscl_lifestyle$size <- as.numeric(mat_mscl_lifestyle$size)
mat_mscl_lifestyle$lifestyle <- as.factor(mat_mscl_lifestyle$lifestyle)
mat_mscl_lifestyle$lifemode <- as.factor(mat_mscl_lifestyle$lifemode)
mat_mscl_lifestyle$diet.broad <- as.factor(mat_mscl_lifestyle$diet.broad)
#str(mat_mscl_lifestyle)

#plot(log(mat_mscl_lifestyle$size), log(mat_mscl_lifestyle$mass), col = c(1:nlevels(mat_mscl_lifestyle$lifestyle))[mat_mscl_lifestyle$lifestyle], pch=20, cex=2)

#mat_mscl_lifestyle[6,] 
#Ammospermohpilus cranium has a problem of scale...

# m is the final data matrix, NB : continuous variables are log values !
m <- mat_mscl_lifestyle[-which(rownames(mat_mscl_lifestyle) == "Ammospermophilus leucurus"),]
m[,c("size","mass","pcsa")] <- log(m[,c("size","mass","pcsa")])

################################################################
# Analyses and plots of size against total muscle mass or pcsa #
################################################################

#Make linear models
mod_mass <- lm(m$mass ~ m$size)
mod_pcsa <- lm(m$pcsa ~ m$size)

mod_mass2 <- lm(m$mass ~ m$size * m$lifestyle_simple)
mod_pcsa2 <- lm(m$pcsa ~ m$size * m$lifestyle_simple)

mod_mass3 <- lm(m$mass ~ m$size + m$lifestyle_simple)
mod_pcsa3 <- lm(m$pcsa ~ m$size + m$lifestyle_simple)

# Compute residuals from models
resid_mass <- residuals(mod_mass)
resid_pcsa <- residuals(mod_pcsa)

#Compute confidence intervals
ci_mass <- confint(mod_mass)
ci_pcsa <- confint(mod_pcsa)

# Print summaries
print("ANOVA of total muscle mass against size * lifestyle")
summary(aov(mod_mass2))
print("Linear regression of total muscle mass against size * lifestyle")
summary(mod_mass2)
print("Linear regression of total muscle mass against size + lifestyle")
summary(mod_mass3)
print("ANOVA of total muscle pcsa against size * lifestyle")
summary(aov(mod_pcsa2))
print("Linear regression of total muscle pcsa against size * lifestyle")
summary(mod_pcsa2)
print("Linear regression of total muscle pcsa against size * lifestyle")
summary(mod_pcsa3)

#Define taxa which are chisel tooth diggers
chisels <- unique(paste(lifestyles$Genre[grep(pattern="(C)", lifestyles$LIFESTYLE)], lifestyles$espèces[grep(pattern="(C)", lifestyles$LIFESTYLE)], sep = " "))

#Make logical to select chisel tooth diggers
chisel_dig <- rep(F, dim(m)[1])
chisel_dig[grep(paste(chisels, collapse="|"), x = rownames(m))] <- T

# Make simpler lifestyle factor
lfst <- as.factor(m$lifestyle_simple)

m <- data.frame(m, chisel = chisel_dig)

cols <- viridis(nlevels(lfst))

if (print.plots) {
# Plot size vs muscle mass and size vs pcsa 
pdf(paste("plots/", if (plot.names) {"named_"}, "muscle_size_regression.pdf", sep =""), height = 6, width = 10)

#palette(c("forestgreen", "brown", "green", "blue", "gray80", "gold", "red", "black"))

if (plot.names) {layout(matrix(c(1,2,1,2,3,3), ncol = 3)); par(mar = c(4,4,1,1))} else {layout(matrix(1:2, ncol=2))}

plot(m$size, m$mass, bg = alpha(cols, alpha=0.85)[lfst], pch=c(21,22)[as.factor(chisel_dig)], cex=1.5, xlab = "Log cranium centroid size", ylab = "Log total muscle mass")
#identify(m[,2:3], labels = rownames(m), cex=0.6)

if (plot.names) {
labels <- abbreviate(rownames(m))
text(m$size, m$mass, labels = labels, cex = 0.5, pos = 1)}

plot.lm.CI(m$size, m$mass, line = F, alpha=0.5, border = F)
abline(mod_mass, lwd = 2)

if (!plot.names) {text(m$size[which(lfst == "subterranean")], m$mass[which(lfst == "subterranean")], labels = shorten.labels(rownames(m),positions=1, split=" ")[which(lfst == "subterranean")], cex = 0.5, pos = 1, font = 2)}

legend("topleft", bty = "n", pt.bg = alpha(c(cols, NA), alpha=0.85), pch=c(rep(21, nlevels(lfst)), 22), cex = 0.8, pt.cex=2, legend = c(levels(lfst), "Chisel-tooth digger"))

plot(m$size, m$pcsa, bg = alpha(cols, alpha=0.85)[lfst], pch=c(21,22)[as.factor(chisel_dig)], cex=1.5, xlab = "Log cranium centroid size", ylab = "Log total muscle PCSA")
#identify(m[,c(2,4)], labels = rownames(m), cex=0.6)
plot.lm.CI(m$size, m$pcsa, line = F, alpha=0.5, border = F)
abline(mod_pcsa, lwd = 2)

if (plot.names) {
labels <- abbreviate(rownames(m))
text(m$size, m$pcsa, labels = labels, cex = 0.5, pos = 1)
par(mar = c(0,0,0,0))
plot(rep_len(c(1,2), length.out=length(rownames(m))), rep_len(1:(length(rownames(m))/2), length.out =length(rownames(m))), xlab = "", ylab = "", type = "n", axes = F, xlim = c(0,3))
text(sort(rep_len(c(1,2), length.out=length(rownames(m)))),rep_len(1:(length(rownames(m))/2), length.out =length(rownames(m))), labels = paste(labels, names(labels), sep = "-"), cex = 0.65)}

if (!plot.names) {text(m$size[which(lfst == "subterranean")], m$pcsa[which(lfst == "subterranean")], labels = shorten.labels(rownames(m),positions=1, split=" ")[which(lfst == "subterranean")], cex = 0.5, pos = 1, font = 2)}

if (!plot.names) {
par(fig = c(0.83,0.95,0.18,0.6), new = T, mar=c(0,2,0,0.1))
boxplot(resid_pcsa ~ lfst, col = cols, las = 2, yaxt = "n", xaxt="n", cex.lab = 0.5)
title(ylab = "Residual muscle PCSA", line = 0.5)
par(fig = c(0.33,0.45,0.18,0.6), new = T, mar=c(0,2,0,0.1))
boxplot(resid_mass ~ lfst, col = cols, las = 2, yaxt = "n", xaxt="n", cex.lab = 0.5)
title(ylab = "Residual muscle mass", line = 0.5)}

dev.off()
}

if (print.plots) {
pdf(paste("plots/", "combined_biplots.pdf", sep = ""), width = 8.5, height = 4)

layout(matrix(1:3, ncol=3))
par(mar = c(4.5,4.5,3,0.5))

plot(bfdat$logmass, bfdat$logBF, bg = alpha(cols[-3], 0.8)[bfdat$life.style.disparat], pch=c(21:22)[bfdat$chisel], cex = 1.5, xlab = "Log10 Body mass (g)", ylab = "Log10 Bite Force (N)", cex.lab = 1.5)
mtext("A.", line = 0.5, adj = 0, cex = 1.5, font = 2)
legend("topleft", pch = c(rep(21,nlevels(bfdat$life.style.disparat)), 22), pt.bg = c(cols[-3], NA), legend = c(levels(bfdat$life.style.disparat), "Chisel-tooth digger"), pt.cex = 2, bty="n")
plot.lm.CI(bfdat$logmass, bfdat$logBF, alpha=0.35, border = F)
abline(modbf1, lwd = 2, col = "black")

plot(m$size, m$mass, bg = alpha(cols, alpha=0.85)[lfst], pch=c(21,22)[as.factor(chisel_dig)], cex=1.5, xlab = "Log cranium centroid size", ylab = "Log total muscle mass", cex.lab = 1.5)
plot.lm.CI(m$size, m$mass, line = F, alpha=0.5, border = F)
mtext("B.", line = 0.5, adj = 0, cex = 1.5, font = 2)
abline(mod_mass, lwd = 2)
legend("topleft", bty = "n", pt.bg = alpha(c(cols, NA), alpha=0.85), pch=c(rep(21, nlevels(lfst)), 22), pt.cex=2, legend = c(levels(lfst), "Chisel-tooth digger"))

plot(m$size, m$pcsa, bg = alpha(cols, alpha=0.85)[lfst], pch=c(21,22)[as.factor(chisel_dig)], cex=1.5, xlab = "Log cranium centroid size", ylab = "Log total muscle PCSA", cex.lab = 1.5)
mtext("C.", line = 0.5, adj = 0, cex = 1.5, font = 2)
plot.lm.CI(m$size, m$pcsa, line = F, alpha=0.5, border = F)
abline(mod_pcsa, lwd = 2)


par(fig = c(0.235,0.325,0.16,0.6), new = T, mar=c(0,2,0,0.1))
boxplot(resmod1 ~ bfdat$life.style.disparat, col = cols[-3], las = 2, yaxt = "n", xaxt="n", cex.lab = 0.5)
title(ylab = "Residual bite force", line = 0.5)
par(fig = c(0.56,0.655,0.16,0.6), new = T, mar=c(0,2,0,0.1))
boxplot(resid_mass ~ lfst, col = cols, las = 2, yaxt = "n", xaxt="n", cex.lab = 0.5)
title(ylab = "Residual muscle mass", line = 0.5)
par(fig = c(0.89,0.99,0.16,0.6), new = T, mar=c(0,2,0,0.1))
boxplot(resid_pcsa ~ lfst, col = cols, las = 2, yaxt = "n", xaxt="n", cex.lab = 0.5)
title(ylab = "Residual muscle PCSA", line = 0.5)

dev.off()
}


# Boxplots to check differences in size-corrected total muscle mass or pcsa with different lifestyle classifications
if (print.plots) {
par(mar = c(7,5,1,1))
layout(matrix(1:4, ncol=2, byrow=2))
boxplot(resid_mass ~ lfst, col = cols, las = 2, xlab="", ylab = "Size-residual total muscle mass")
points(resid_mass[which(chisel_dig)] ~ lfst[which(chisel_dig)], pch = 22)

boxplot(resid_pcsa ~ lfst, col = cols, las=2, xlab="", ylab = "Size-residual total muscle PCSA")
points(resid_pcsa[which(chisel_dig)] ~ lfst[which(chisel_dig)], pch = 22)

boxplot(resid_mass ~ m$lifemode, col = c(1,3,2,5), las = 2, xlab="", ylab = "Size-residual total muscle mass")
points(resid_mass[which(chisel_dig)] ~ m$lifemode[which(chisel_dig)], pch = 22)

boxplot(resid_pcsa ~ m$lifemode, col = c(1,3,2,5), las=2, xlab="", ylab = "Size-residual total muscle PCSA")
points(resid_pcsa[which(chisel_dig)] ~ m$lifemode[which(chisel_dig)], pch = 22)
}
#plot(m$size, resid_pcsa,bg = alpha(c(1:8), alpha=0.85)[m$lifestyle], pch=c(21,24)[as.factor(chisel_dig)], cex=2)
#identify(x=m$sp_avg_csize, y=resid_pcsa, labels=rownames(m))

# Ordered barplots to illustrate which lifestyles have proportionally more muscle, with species names to check individual taxa.
if (print.plots) {
pdf(paste("plots/","barplots_ordered_resid_muscles.pdf", sep = ""), height = 16, width = 8)
layout(matrix(1:2, ncol=2))
par(mar = c(3,8,0,0.2))
barplot(sort(resid_mass), col = alpha(cols, alpha=0.85)[lfst[order(resid_mass)]], names.arg = rownames(m)[order(resid_mass)], las = 2, cex.names=0.6, horiz=T)
barplot(sort(resid_pcsa), col = alpha(cols, alpha=0.85)[lfst[order(resid_pcsa)]], names.arg = rownames(m)[order(resid_pcsa)], las = 2, cex.names=0.6, horiz=T)
dev.off()
}

# Ordered barplots to illustrate which DIET (simplified) have proportionally more muscle, with species names to check individual taxa.
dietsimple <- as.factor(unlist(lapply(X = strsplit(as.character(m$diet.broad), split = " "), FUN = fun <- function(x) {x[[1]]})))

if (print.plots) {
pdf(paste("plots/","barplots_ordered_resid_muscles_diet.pdf", sep = ""), height = 16, width = 8)
layout(matrix(1:2, ncol=2))
par(mar = c(3,8,0,0.2))
barplot(sort(resid_mass), col = cols[dietsimple[order(resid_mass)]], names.arg = rownames(m)[order(resid_mass)], las = 2, cex.names=0.6, horiz=T)
barplot(sort(resid_pcsa), col = cols[dietsimple[order(resid_pcsa)]], names.arg = rownames(m)[order(resid_pcsa)], las = 2, cex.names=0.6, horiz=T)
dev.off()
}

#################################################
# Make Traitgrams for total mass and total PCSA #
#################################################

mprun <- data.pruning(tree = sp_tre, dat = m[-which(is.na(m$lifestyle_simple)),], split.t = "_", split.d = " ", positions = c(1,2))
resmass <- residuals(lm(mprun$pruned.data$mass ~ mprun$pruned.data$size))
respcsa <- residuals(lm(mprun$pruned.data$pcsa ~ mprun$pruned.data$size))
isbur <- mprun$pruned.data$lifestyle_simple
isbur[which(!(mprun$pruned.data$lifestyle_simple == "subterranean" | mprun$pruned.data$lifestyle_simple == "burrower"))] <- "NO"
isbur <- as.factor(isbur)
names(resmass) <- names(respcsa) <- names(isbur) <- mprun$pruned.tree$tip.label
mtree <- make.simmap(tree = mprun$pruned.tree, x = isbur, message = F)
cols<-c(cols[2], alpha("black", 0.1), cols[4]); names(cols)<-levels(isbur)

if (print.plots) {
if (plot.names) {
pdf(paste("plots/", "named_traitgram_residual_muscle_mass.pdf", sep = ""), height=12, width = 5)
par(mar = c(4,4,1,1))
phenogram(x = resmass, tree = mtree, colors = cols, fsize = 0.5, ylab = "Size-residual muscle mass")
dev.off()

pdf(paste("plots/", "named_traitgram_residual_muscle_pcsa.pdf", sep = ""), height=12, width = 5)
par(mar = c(4,4,1,1))
phenogram(x = respcsa, tree = mtree, colors = cols, fsize = 0.5, ylab = "Size-residual muscle PCSA")
dev.off()
} else {

pdf(paste("plots/", "traitgram_residual_muscle_mass_pcsa.pdf", sep = ""), height= 5, width = 5)
layout(matrix(1:2, ncol=2))
par(mar = c(1,5,1,0))
phenogram(x = resmass, tree = mtree, colors = cols, fsize = 0.01, spread.labels =F,xlim = c(0,70), axes = F)
axis(side = 2)
mtext(side = 2, text = "Size-residual muscle mass", line = 3)
par(mar = c(1,0,1,5))
phenogram(x = respcsa, tree = mtree, colors = cols, fsize = 0.01, spread.labels =F,xlim = c(70,0), axes = F)
axis(side = 4)
mtext(side = 4, text = "Size-residual muscle PCSA", line = 3)
dev.off()
}
}

# Compute convergence signal

foct <- names(isbur[which(isbur == "burrower" | isbur == "subterranean")])

#print("Stayton's C1-C4 indices for size-residual total muscle mass and total PCSA")
#cvsg_mass <- convSig(phy = mtree, traits = matrix(resmass), focaltaxa = foct, nsim = NSIM)
#print(cvsg_mass)
#cvsg_pcsa <- convSig(phy = mtree, traits = matrix(respcsa), focaltaxa = foct, nsim = NSIM)
#print(cvsg_pcsa)
print("Pagel's Lambda for size-residual total muscle mass and total PCSA")
phylosig(tree = mtree, x = resmass, method = "lambda")
phylosig(tree = mtree, x = respcsa, method = "lambda")

##################################
# Make phylogenetic heatmap plot #
##################################

cols <- c(viridis(5)[2], "black", viridis(5)[4]); names(cols)<-levels(isbur)
X <- cbind(resmass, respcsa)

if (print.plots) {
pdf(paste("plots/", "phylo_heatmap_mass_PCSA.pdf", sep = ""), height=12, width = 5)
phylo.heatmap(tree = mtree, X = X, split=c(0.5,0.5), fsize=c(0.01,1,0.5), standardize=F, colors = viridis(5))
par(fig = c(0.05,0.612,0.096,0.904), new = T)
plot(mtree, add = T, fsize = 0.01, lwd = 3, colors = cols)
dev.off()
}

########################################
# Include individual muscle group data #
########################################

index_prop <- shorten.labels(rownames(muscle_prop), positions = c(1,2), split=" ")
sp_prop <- avg.matrix.fac(mat = muscle_prop, INDEX = index_prop)
prop_match <- mat.mtch(mlifestyle, sp_prop)

mscl_lsr_mass <- read.csv("tables/adductor_LSR_mass.csv", h = T, sep = ",", dec = ".", row.names = 1)
index_lsr_mass <- shorten.labels(rownames(mscl_lsr_mass), positions = c(1,2), split=" ")
sp_lsr_mass <- avg.matrix.fac(mat = mscl_lsr_mass, INDEX = index_lsr_mass )
lsr_match <- mat.mtch(mlifestyle, sp_lsr_mass)

dat_dig <- data.frame(dat_dig)
index_dig_mass <- shorten.labels(rownames(dat_dig), positions = c(1,2), split=" ")
sp_dig_mass <- na.omit(avg.matrix.fac(mat = as.matrix(dat_dig$lsr_mass), INDEX = index_dig_mass))

all_lsr_match <- mat.mtch(lsr_match, sp_dig_mass)
all_lsr_match[,7:13] <- apply(all_lsr_match[,7:13], 2, as.numeric)

##############################################
# Discriminant analyses of muscle group data #
##############################################

matlda <- as.matrix(prop_match[,c("DM","ePT", "iPT", "SM","T","ZM")])

pdf(paste("plots/", "LDA_lifemode.pdf"))
lda_lifemode <- ldawithplot(mat = matlda, grouping=as.factor(prop_match$lifemode), pos.inset = c(0.15,0.45,0.15,0.45))
dev.off()

pdf(paste("plots/", "LDA_lifestyle.pdf"))
lda_lifestyle <- ldawithplot(mat = matlda, grouping=as.factor(prop_match$lifestyle_simple), pos.inset = c(0.15,0.45,0.15,0.45))
dev.off()

pdf(paste("plots/", "LDA_diet.pdf"))
lda_diet <- ldawithplot(mat = matlda, grouping=as.factor(prop_match$diet.broad), pos.inset = c(0.15,0.45,0.15,0.45))
dev.off()

# Do the same for LSR data instead of proportions
pdf(paste("plots/", "LDA_lsr_lifemode.pdf"))
lda_lsr_lifemode <- ldawithplot(mat = all_lsr_match[,7:13], grouping=as.factor(all_lsr_match$lifemode), pos.inset = c(0.15,0.45,0.15,0.45))
dev.off()

pdf(paste("plots/", "LDA_lsr_lifestyle.pdf"))
lda_lsr_lifestyle <- ldawithplot(mat = all_lsr_match[,7:13], grouping=as.factor(all_lsr_match$lifestyle_simple), pos.inset = c(0.15,0.45,0.15,0.45))
dev.off()

pdf(paste("plots/", "LDA_lsr_diet.pdf"))
lda_lsr_diet <- ldawithplot(mat = all_lsr_match[,7:13], grouping=as.factor(all_lsr_match$diet.broad), pos.inset = c(0.15,0.45,0.15,0.45))
dev.off()

prop.correct.pred <- function(LDA, grouping) {
predictgrp <- LDA$predictions$class
correct <- na.omit(grouping)
df <- data.frame(correct = correct, predicted = predictgrp)
rownames(df) <- rownames(LDA$predictions$posterior)
tab <- table(df)
tab/apply(tab, 1, sum)*100
}

lifemode_pred <- prop.correct.pred(lda_lifemode, grouping = as.factor(prop_match$lifemode))
lifestyle_pred <- prop.correct.pred(lda_lifestyle, grouping = as.factor(prop_match$lifestyle_simple))
diet_pred <- prop.correct.pred(lda_diet, grouping = as.factor(prop_match$diet.broad))

lifemode_pred2 <- prop.correct.pred(lda_lsr_lifemode, grouping = as.factor(all_lsr_match$lifemode))
lifestyle_pred2 <- prop.correct.pred(lda_lsr_lifestyle, grouping = as.factor(all_lsr_match$lifestyle_simple))
diet_pred2 <- prop.correct.pred(lda_lsr_diet, grouping = as.factor(all_lsr_match$diet.broad))

#Do LDA with total mass / total PCSA only

lif <- as.factor(m$lifestyle_simple[-which(is.na(m$lifestyle))])
mat <- cbind(resid_mass, resid_pcsa)[-which(is.na(m$lifestyle)),]

pdf(paste("plots/", "LDA_totalmassandPCSA_lifestyle.pdf"))
lda_mass <- ldawithplot(mat = mat, grouping = lif)
dev.off()

mass_pred <- prop.correct.pred(lda_mass, grouping = lif)

# Print "confusion matrices" for the LDAs

print("Confusion matrix for lifestyle LDA using only residual total mass and PCSA")
print(mass_pred)
print("Confusion matrix for lifemode LDA using muscle propotions")
print(lifemode_pred)
print("Confusion matrix for lifestyle LDA using muscle propotions")
print(lifestyle_pred)

################################################
# Compare muscle proportions along phylogeny   #
################################################

DIGmatch <- sp_dig_mass[match(rownames(prop_match), rownames(sp_dig_mass))]
scaledDIG <- scale(DIGmatch)
#prop <- prop_match[, c("DM","ePT", "iPT", "SM","T","ZM")]
prop_match <- data.frame(prop_match, scaledDIG)

cols <- c(
	hsv(0.897, 0.728, 0.952),
	hsv(0.726, 0.440, 0.885),
	hsv(0.167, 0.714, 0.878),
	hsv(0.093, 1, 1),
	hsv(0.483, 0.714, 0.878),
	hsv(0.593, 1, 1),
	hsv(0, 0.690, 0.788))
names(cols) <- c("SM", "DM", "ZM", "T", "ePT", "iPT", "DIG")

prop_pruned <- data.pruning(tree = sp_tre, dat = prop_match, split.t = "_", split.d = " ", positions = c(1,2))

suborders <- as.factor(c(rep("Sciuromorpha", 25), 
		rep("Hystricomorpha", 41),
		rep("Supramyomorpha", 144)))	
			
edg <- find.edgecols(phy = prop_pruned$pruned.tree, factor = suborders, color=c("gray30","chartreuse3","red"))
burrowers <- prop_pruned$pruned.data$lifestyle_simple == "subterranean" | prop_pruned$pruned.data$lifestyle_simple == "burrower"
fnt <- rep(3, length(prop_pruned$pruned.tree$tip.label))
fnt[which(burrowers)] <- 4

colss <- matrix(rep(alpha(cols, 0.5), length(prop_pruned$pruned.tree$tip.label)), ncol=7, byrow=T)
for (i in 1:length(cols)) {colss[which(burrowers),i] <- cols[i]}

if (print.plots) {
pdf(paste("plots/", "all_phylo_muscle_prop.pdf"), width = 15, height = 22)
par(mar = c(0,0,0,0))
plot(prop_pruned$pruned.tree, show.tip.label=T, edge.color=edg[[2]], lwd=2, x.lim = 400, cex=0.8, font = fnt)
phydataplot(x = prop_pruned$pruned.data$SM, phy = prop_pruned$pruned.tree, offset = 55, style = "bars", beside = F, col = colss[,1], border ="white")
phydataplot(x = prop_pruned$pruned.data$DM, phy = prop_pruned$pruned.tree, offset =100, style = "bars", beside = F, col = colss[,2], border ="white")
phydataplot(x = prop_pruned$pruned.data$ZM, phy = prop_pruned$pruned.tree, offset = 170, style = "bars", beside = F, col = colss[,3], border ="white")
phydataplot(x = prop_pruned$pruned.data$T, phy = prop_pruned$pruned.tree, offset = 220, style = "bars", beside = F, col = colss[,4], border ="white")
phydataplot(x = prop_pruned$pruned.data$ePT, phy = prop_pruned$pruned.tree, offset = 270, style = "bars", beside = F, col = colss[,5], border ="white")
phydataplot(x = prop_pruned$pruned.data$iPT, phy = prop_pruned$pruned.tree, offset = 300, style = "bars", beside = F, col = colss[,6], border ="white")
#phydataplot(x = prop_pruned$pruned.data$scaledDIG*5, phy = prop_pruned$pruned.tree, offset = 320, style = "bars", beside = F, col = colss[,7], border ="white")
dev.off()


pdf(paste("plots/", "all_phylo_muscle_prop_scaled.pdf"), width = 15, height = 22)
par(mar = c(0,0,0,0))
scaling <- 6
scale <- T
plot(prop_pruned$pruned.tree, show.tip.label=T, edge.color=edg[[2]], lwd=2, x.lim = 400, cex=0.8, font = fnt)
phydataplot(x = scale(prop_pruned$pruned.data$SM, scale = scale), phy = prop_pruned$pruned.tree, offset = 65, style = "bars", beside = F, col = colss[,1], border ="white", scaling=scaling)
phydataplot(x = scale(prop_pruned$pruned.data$DM, scale = scale), phy = prop_pruned$pruned.tree, offset = 105, style = "bars", beside = F, col = colss[,2], border ="white", scaling=scaling)
phydataplot(x = scale(prop_pruned$pruned.data$ZM, scale = scale), phy = prop_pruned$pruned.tree, offset = 140, style = "bars", beside = F, col = colss[,3], border ="white", scaling=scaling)
phydataplot(x = scale(prop_pruned$pruned.data$T, scale = scale), phy = prop_pruned$pruned.tree, offset = 190, style = "bars", beside = F, col = colss[,4], border ="white", scaling=scaling)
phydataplot(x = scale(prop_pruned$pruned.data$ePT, scale = scale), phy = prop_pruned$pruned.tree, offset = 230, style = "bars", beside = F, col = colss[,5], border ="white", scaling=scaling)
phydataplot(x = scale(prop_pruned$pruned.data$iPT, scale = scale), phy = prop_pruned$pruned.tree, offset = 270, style = "bars", beside = F, col = colss[,6], border ="white", scaling=scaling)
phydataplot(x = prop_pruned$pruned.data$scaledDIG, phy = prop_pruned$pruned.tree, offset = 310, style = "bars", beside = F, col = colss[,7], border ="white", scaling=scaling)
dev.off()
}

##################################
# Make phylogenetic heatmap plot #
##################################

ppd <- prop_pruned$pruned.data[-which(is.na(prop_pruned$pruned.data$lifestyle_simple)),]
pptree <- drop.tip(phy = prop_pruned$pruned.tree, tip = rownames(prop_pruned$pruned.data)[which(is.na(prop_pruned$pruned.data$lifestyle_simple))])

isbur <- ppd$lifestyle_simple
isbur[which(!(ppd$lifestyle_simple == "subterranean" | ppd$lifestyle_simple == "burrower"))] <- "NO"
isbur <- as.factor(isbur)
names(isbur) <- pptree$tip.label
mtree <- make.simmap(tree = pptree, x = isbur, message = F)
colors<-c(viridis(5)[2], "black", viridis(5)[4]); names(colors)<-levels(isbur)

if (print.plots) {
pdf(paste("plots/", "phylo_heatmap.pdf", sep = ""), height=12, width = 5)
phylo.heatmap(tree = mtree, X = ppd[7:12], split=c(0.75,0.25), fsize=c(0.01,1,0.5), standardize=T, colors = viridis(8))
par(fig = c(0.05,0.61,0.13,0.905), new = T)
plot(mtree, add = T, fsize = 0.01, lwd = 3, colors = colors)
dev.off()
}

###################
# Make traitgrams #
###################

if (print.plots) {
if (plot.names) {
pdf(paste("plots/", "traitgram_propDM.pdf", sep = ""), height=12, width = 5)
par(mar = c(4,4,1,1))
x <- ppd$DM; names(x) <- mtree$tip.label
col <- c(alpha(cols[2], 0.6), alpha(cols[2], 0.15), cols[2]); names(col)<-c( "burrower","NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, ylab = "Proportion DM")
dev.off()

pdf(paste("plots/", "traitgram_propSM.pdf", sep = ""), height=12, width = 5)
par(mar = c(4,4,1,1))
x <- ppd$SM; names(x) <- mtree$tip.label
col <- c(alpha(cols[1], 0.6), alpha(cols[1], 0.15), cols[1]); names(col)<-c( "burrower","NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, ylab = "Proportion SM")
dev.off()

pdf(paste("plots/", "traitgram_propZM.pdf", sep = ""), height=12, width = 5)
par(mar = c(4,4,1,1))
x <- ppd$ZM; names(x) <- mtree$tip.label
col <- c(alpha(cols[3], 0.6), alpha(cols[3], 0.25), cols[3]); names(col)<-c( "burrower","NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, ylab = "Proportion ZM")
dev.off()

pdf(paste("plots/", "traitgram_propT.pdf", sep = ""), height=12, width = 5)
par(mar = c(4,4,1,1))
x <- ppd$T; names(x) <- mtree$tip.label
col <- c(alpha(cols[4], 0.6), alpha(cols[4], 0.15), cols[4]); names(col)<-c( "burrower","NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, ylab = "Proportion T")
dev.off()

pdf(paste("plots/", "traitgram_propePT.pdf", sep = ""), height=12, width = 5)
par(mar = c(4,4,1,1))
x <- ppd$ePT; names(x) <- mtree$tip.label
col <- c(alpha(cols[5], 0.6), alpha(cols[5], 0.15), cols[5]); names(col)<-c( "burrower","NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, ylab = "Proportion ePT")
dev.off()

pdf(paste("plots/", "traitgram_propiPT.pdf", sep = ""), height=12, width = 5)
par(mar = c(4,4,1,1))
x <- ppd$iPT; names(x) <- mtree$tip.label
col <- c(alpha(cols[6], 0.6), alpha(cols[6], 0.15), cols[6]); names(col)<-c( "burrower","NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, ylab = "Proportion iPT")
dev.off()

#pdf(paste("plots/", "traitgram_propDIG.pdf", sep = ""), height=12, width = 5)
#par(mar = c(4,4,1,1))
#x <- ppd$DIG; names(x) <- mtree$tip.label
#col <- c(alpha(cols[7], 0.2), cols[7]); names(col)<-c( "FALSE","TRUE")
#phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, ylab = "Proportion DIG")
#dev.off()
} else {
pdf(paste("plots/", "all_traitgram_prop_nonames_sub_only.pdf", sep = ""), height=10, width = 5)
layout(matrix(1:6, ncol=2))
par(mar = c(0,5,0,0))
x <- ppd$DM; names(x) <- mtree$tip.label
col <- c(alpha(cols[2], 0.15), alpha(cols[2], 0.15), cols[2]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.01, spread.labels =F,xlim = c(0,70), axes = F)
axis(side = 2)
mtext(side = 2, text = "Proportion DM", line = 3)

x <- ppd$ZM; names(x) <- mtree$tip.label
col <- c(alpha(cols[3], 0.15), alpha(cols[3], 0.15), cols[3]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.01, spread.labels =F,xlim = c(0,70), axes = F)
axis(side = 2)
mtext(side = 2, text = "Proportion ZM", line = 3)

x <- ppd$ePT; names(x) <- mtree$tip.label
col <- c(alpha(cols[5], 0.15), alpha(cols[5], 0.15), cols[5]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.01, spread.labels =F,xlim = c(0,70), axes = F)
axis(side = 2)
mtext(side = 2, text = "Proportion ePT", line = 3)

par(mar = c(0,0,0,5))
x <- ppd$SM; names(x) <- mtree$tip.label
col <- c(alpha(cols[1], 0.15), alpha(cols[1], 0.15), cols[1]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.01, spread.labels =F,xlim = c(70,0), axes = F)
axis(side = 4)
mtext(side = 4, text = "Proportion SM", line = 3)

x <- ppd$T; names(x) <- mtree$tip.label
col <- c(alpha(cols[4], 0.15), alpha(cols[4], 0.15), cols[4]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.01, spread.labels =F,xlim = c(70,0), axes = F)
axis(side = 4)
mtext(side = 4, text = "Proportion T", line = 3)

x <- ppd$iPT; names(x) <- mtree$tip.label
col <- c(alpha(cols[6], 0.15), alpha(cols[6], 0.15), cols[6]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.01, spread.labels =F,xlim = c(70,0), axes = F)
axis(side = 4)
mtext(side = 4, text = "Proportion iPT", line = 3)
dev.off()
}
}

if (print.plots) {
pdf(paste("plots/", "traitgram_propDIG.pdf", sep = ""), height=12, width = 5)
x <- ppd$scaledDIG; names(x) <- mtree$tip.label
nas <- which(is.na(x))
col <- c(alpha(cols[7], 0.15), alpha(cols[7], 0.15), cols[7]); names(col)<-c("burrower", "NO", "subterranean")
par(mar = c(4,4,1,1))
phenogram(x = x[-nas], tree = drop.tip(mtree, names(nas)), colors = col, fsize = 0.5, spread.labels = T, xlim = c(0,70))
dev.off()
}

if (print.plots & plot.names) {

pdf(paste("plots/", "all_traitgram_with_names.pdf", sep = ""), height=15, width = 8)
layout(matrix(1:6, ncol=2))
par(mar = c(0,5,0,1))
x <- ppd$DM; names(x) <- mtree$tip.label
col <- c(alpha(cols[2], 0.15), alpha(cols[2], 0.15), cols[2]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, spread.labels =F, axes = F)
axis(side = 2)
mtext(side = 2, text = "Proportion DM", line = 3)

x <- ppd$ZM; names(x) <- mtree$tip.label
col <- c(alpha(cols[3], 0.15), alpha(cols[3], 0.15), cols[3]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, spread.labels =F, axes = F)
axis(side = 2)
mtext(side = 2, text = "Proportion ZM", line = 3)

x <- ppd$ePT; names(x) <- mtree$tip.label
col <- c(alpha(cols[5], 0.15), alpha(cols[5], 0.15), cols[5]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, spread.labels =F, axes = F)
axis(side = 2)
mtext(side = 2, text = "Proportion ePT", line = 3)

par(mar = c(0,0,0,5))
x <- ppd$SM; names(x) <- mtree$tip.label
col <- c(alpha(cols[1], 0.15), alpha(cols[1], 0.15), cols[1]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, spread.labels =F, axes = F)
axis(side = 4)
mtext(side = 4, text = "Proportion SM", line = 3)

x <- ppd$T; names(x) <- mtree$tip.label
col <- c(alpha(cols[4], 0.15), alpha(cols[4], 0.15), cols[4]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, spread.labels =F, axes = F)
axis(side = 4)
mtext(side = 4, text = "Proportion T", line = 3)

x <- ppd$iPT; names(x) <- mtree$tip.label
col <- c(alpha(cols[6], 0.15), alpha(cols[6], 0.15), cols[6]); names(col)<-c("burrower", "NO", "subterranean")
phenogram(x = x, tree = mtree, colors = col, fsize = 0.7, spread.labels =F, axes = F)
axis(side = 4)
mtext(side = 4, text = "Proportion iPT", line = 3)
dev.off()
}

################################################
# Compute C1-C4 indices for individual muscles #
################################################

foct <- names(isbur[which(isbur == "burrower" | isbur == "subterranean")])

#print("Stayton C1-C4 indices for individual muscles proportions, for burrowers and subterranean")
#cvsg_SM <- convSig(phy = mtree, traits = matrix(ppd$SM), focaltaxa = foct, nsim = NSIM)
#print(cvsg_SM)
#cvsg_DM <- convSig(phy = mtree, traits = matrix(ppd$DM), focaltaxa = foct, nsim = NSIM)
#print(cvsg_DM)
#cvsg_ZM <- convSig(phy = mtree, traits = matrix(ppd$ZM), focaltaxa = foct, nsim = NSIM)
#print(cvsg_ZM)
#cvsg_ePT <- convSig(phy = mtree, traits = matrix(ppd$ePT), focaltaxa = foct, nsim = NSIM)
#print(cvsg_ePT)
#cvsg_iPT <- convSig(phy = mtree, traits = matrix(ppd$iPT), focaltaxa = foct, nsim = NSIM)
#print(cvsg_iPT)
#cvsg_T <- convSig(phy = mtree, traits = matrix(ppd$T), focaltaxa = foct, nsim = NSIM)
#print(cvsg_T)

print("Stayton C1-C4 indices for combined muscles proportions, for burrowers and subterranean")
cvsg_prop <- convSig(phy = mtree, traits = as.matrix(ppd[,7:12]), focaltaxa = foct, nsim = NSIM)
print(cvsg_prop)

foctsub <- names(isbur[which(isbur == "subterranean")])

#print("Stayton C1-C4 indices for individual muscles proportions, for subterranean only")
#cvsg_SM_sub <- convSig(phy = mtree, traits = matrix(ppd$SM), focaltaxa = foctsub, nsim = NSIM)
#print(cvsg_SM_sub)
#cvsg_DM_sub <- convSig(phy = mtree, traits = matrix(ppd$DM), focaltaxa = foctsub, nsim = NSIM)
#print(cvsg_DM_sub)
#cvsg_ZM_sub <- convSig(phy = mtree, traits = matrix(ppd$ZM), focaltaxa = foctsub, nsim = NSIM)
#print(cvsg_ZM_sub)
#cvsg_ePT_sub <- convSig(phy = mtree, traits = matrix(ppd$ePT), focaltaxa = foctsub, nsim = NSIM)
#print(cvsg_ePT_sub)
#cvsg_iPT_sub <- convSig(phy = mtree, traits = matrix(ppd$iPT), focaltaxa = foctsub, nsim = NSIM)
#print(cvsg_iPT_sub)
#cvsg_T_sub <- convSig(phy = mtree, traits = matrix(ppd$T), focaltaxa = foctsub, nsim = NSIM)
#print(cvsg_T_sub)

print("Stayton C1-C4 indices for combined muscles proportions, for subterranean only")
cvsg_prop_sub <- convSig(phy = mtree, traits = as.matrix(ppd[,7:12]), focaltaxa = foctsub, nsim = NSIM)
print(cvsg_prop_sub)

print("Pagel lambda indices for individual muscles proportions")
phylosig(tree = mtree, x = ppd$DM, method = "lambda")
phylosig(tree = mtree, x = ppd$ZM, method = "lambda")
phylosig(tree = mtree, x = ppd$ePT, method = "lambda")
phylosig(tree = mtree, x = ppd$iPT, method = "lambda")
phylosig(tree = mtree, x = ppd$T, method = "lambda")
phylosig(tree = mtree, x = ppd$SM, method = "lambda")

###############################################################
# Compare muscle proportions in burrowers and subterrean only #
###############################################################

propbur <- prop_match[which(prop_match$lifestyle_simple == "burrower" | prop_match$lifestyle_simple == "subterranean"), c("DM","ePT", "iPT", "SM","T","ZM", "scaledDIG")]

burandsub <- rownames(prop_match)[which(prop_match$lifestyle_simple == "burrower" | prop_match$lifestyle_simple == "subterranean")]

bs <- gsub(x = burandsub, pattern = " ", replacement = "_")

bur_tre <- keep.tip(phy = sp_tre, tip=bs[-11])

#bur_tre <- read.tree(file = "burrowers_tree.txt")

prop_bur <- data.pruning(tree = bur_tre, dat = propbur, split.t = "_", split.d = " ", positions = c(1,2))

if (print.plots) {
pdf(paste("plots/", "Burrower_phylo_muscle_prop_onebar.pdf"), width = 6, height = 5)
par(mar = c(1,1,1,1))
plot(prop_bur$pruned.tree, x.lim = 250, adj = 0, cex = 0.85)
phydataplot(x = prop_bur$pruned.data[,-7], phy = prop_bur$pruned.tree, offset = 90, style = "bars", beside = F, col = cols)
dev.off()
}

if (print.plots) {
pdf(paste("plots/", "Burrower_phylo_muscle_prop.pdf"), width = 10, height = 5)
par(mar = c(0,0,1,1))
plot(prop_bur$pruned.tree, x.lim = 430)
phydataplot(x = prop_bur$pruned.data$SM, phy = prop_bur$pruned.tree, offset = 110, style = "bars", beside = F, col = cols[1], border ="white")
phydataplot(x = prop_bur$pruned.data$DM, phy = prop_bur$pruned.tree, offset =160, style = "bars", beside = F, col = cols[2], border ="white")
phydataplot(x = prop_bur$pruned.data$ZM, phy = prop_bur$pruned.tree, offset =230, style = "bars", beside = F, col = cols[3], border ="white")
phydataplot(x = prop_bur$pruned.data$T, phy = prop_bur$pruned.tree, offset =280, style = "bars", beside = F, col = cols[4], border ="white")
phydataplot(x = prop_bur$pruned.data$ePT, phy = prop_bur$pruned.tree, offset =340, style = "bars", beside = F, col = cols[5], border ="white")
phydataplot(x = prop_bur$pruned.data$iPT, phy = prop_bur$pruned.tree, offset =360, style = "bars", beside = F, col = cols[6], border ="white")
dev.off()

pdf(paste("plots/", "Burrower_phylo_muscle_prop_scaled.pdf"), width = 10, height = 5)
par(mar = c(0,0,0,1))
scaling <- 10
plot(prop_bur$pruned.tree, x.lim = 470)
phydataplot(x = scale(prop_bur$pruned.data$SM), phy = prop_bur$pruned.tree, offset = 120, style = "bars", beside = F, col = cols[1], border ="white", scaling = scaling)
phydataplot(x = scale(prop_bur$pruned.data$DM), phy = prop_bur$pruned.tree, offset = 155, style = "bars", beside = F, col = cols[2], border ="white", scaling = scaling)
phydataplot(x = scale(prop_bur$pruned.data$ZM), phy = prop_bur$pruned.tree, offset = 200, style = "bars", beside = F, col = cols[3], border ="white", scaling = scaling)
phydataplot(x = scale(prop_bur$pruned.data$T), phy = prop_bur$pruned.tree, offset = 270, style = "bars", beside = F, col = cols[4], border ="white", scaling = scaling)
phydataplot(x = scale(prop_bur$pruned.data$ePT), phy = prop_bur$pruned.tree, offset = 310, style = "bars", beside = F, col = cols[5], border ="white", scaling = scaling)
phydataplot(x = scale(prop_bur$pruned.data$iPT), phy = prop_bur$pruned.tree, offset = 355, style = "bars", beside = F, col = cols[6], border ="white", scaling = scaling)
phydataplot(x = prop_bur$pruned.data$scaledDIG, phy = prop_bur$pruned.tree, offset = 400, style = "bars", beside = F, col = cols[7], border ="white", scaling = scaling)
dev.off()
}

#############################################################
# Similar analyses and plot, but phylogenetically corrected #
#############################################################

prnd <- data.pruning(tree = sp_tre, dat = m, split.t = "_",split.d = " ", positions = c(1,2))

# PIC regressions

picsize <- pic(prnd$pruned.data$size, phy = prnd$pruned.tree)
picmass <- pic(prnd$pruned.data$mass, phy = prnd$pruned.tree)
picpcsa <- pic(prnd$pruned.data$pcsa, phy = prnd$pruned.tree)


modpic1 <- lm(picmass ~ picsize)
modpic2 <- lm(picpcsa ~ picsize)

print("Linear regression of phylogenetic contrasts of total muscle mass against contrasts of size with 95% confidence intervals")
summary(modpic1)
confint(modpic1)
print("Linear regression of phylogenetic contrasts of total muscle PCSA against contrasts of size, with confidence intervals")
summary(modpic2)
confint(modpic2)

# PGLS models

dat <- prnd$pruned.data[-which(is.na(prnd$pruned.data$lifestyle_simple)),]
dat$lifestyle_simple <- relevel(as.factor(dat$lifestyle_simple), ref = "subterranean")


pglsmass <- gls(mass ~ size + lifestyle_simple, data = dat, correlation = corBrownian(phy = prnd$pruned.tree))
post.hoc_mass <- glht(pglsmass, linfct=mcp(lifestyle_simple="Tukey"))
print("PGLS regression of total muscle mass against size + lifestyle, with post-hoc pairwise difference test")
summary(pglsmass)
summary(post.hoc_mass)

pglspcsa <- gls(pcsa ~ size + lifestyle_simple, data = dat, correlation = corBrownian(phy = prnd$pruned.tree))
post.hoc_pcsa <- glht(pglspcsa, linfct=mcp(lifestyle_simple="Tukey"))
print("PGLS regression of total muscle pcsa against size + lifestyle, with post-hoc pairwise difference test")
summary(pglspcsa)
summary(post.hoc_pcsa)

fchisel <- as.factor(prnd$pruned.data$chisel)

pglsmasschisel <- gls(mass ~ size + fchisel, data = prnd$pruned.data, correlation = corBrownian(phy = prnd$pruned.tree))
post.hoc_mass_chisel <- glht(pglsmasschisel, linfct=mcp(fchisel="Tukey"))
print("PGLS regression of total muscle mass against size + chisel tooth digging (true or false), with post-hoc pairwise difference test")
summary(pglsmasschisel)
summary(post.hoc_mass_chisel)

pglspcsachisel <- gls(pcsa ~ size + fchisel, data = prnd$pruned.data, correlation = corBrownian(phy = prnd$pruned.tree))
post.hoc_pcsa_chisel <- glht(pglspcsachisel, linfct=mcp(fchisel="Tukey"))
print("PGLS regression of total muscle pcsa against size + chisel tooth digging (true or false), with post-hoc pairwise difference test")
summary(pglspcsachisel)
summary(post.hoc_pcsa_chisel)


if (print.plots) {
pdf(paste("plots/","muscle_size_PGLS_regression.pdf"), height = 6, width = 10)
palette(viridis(length(unique(prnd$pruned.data$lifestyle_simple))))
# Plot size vs muscle mass and size vs pcsa 
layout(matrix(1:2, ncol=2))
plot(prnd$pruned.data$size, prnd$pruned.data$mass, bg = alpha(palette(), alpha=0.85)[as.factor(prnd$pruned.data$lifestyle_simple)], pch=c(21,24)[as.factor(prnd$pruned.data$chisel)], cex=2, xlab = "Log cranium centroid size", ylab = "Log total muscle mass")
abline(pglsmass$coefficients[1], pglsmass$coefficients[2])
abline(lm(prnd$pruned.data$mass~prnd$pruned.data$size), lty = 2)
plot.lm.CI(prnd$pruned.data$size, prnd$pruned.data$mass, line = F, alpha=0.5, border = F)

plot(prnd$pruned.data$size, prnd$pruned.data$pcsa, bg = alpha(palette(), alpha=0.85)[as.factor(prnd$pruned.data$lifestyle_simple)], pch=c(21,24)[as.factor(prnd$pruned.data$chisel)], cex=2, xlab = "Log cranium centroid size", ylab = "Log total muscle PCSA")
abline(pglspcsa$coefficients[1], pglspcsa$coefficients[2])
abline(lm(prnd$pruned.data$pcsa~prnd$pruned.data$size), lty = 2)
plot.lm.CI(prnd$pruned.data$size, prnd$pruned.data$pcsa, line = F, alpha=0.5, border = F)

legend("topleft", bty = "n", pt.bg = alpha(c(palette()[1:length(unique(na.omit(prnd$pruned.data$lifestyle_simple)))],NA), alpha=0.85), pch=c(rep(21,length(unique(na.omit(prnd$pruned.data$lifestyle_simple)))),24), pt.cex=2, legend = c(unique(na.omit(prnd$pruned.data$lifestyle_simple)), "Chisel-tooth digger"))
dev.off()
}

#######################################################################
# GMM Morphospaces, with lifestyles, and compute convergence measures #
#######################################################################

cras <- scan("CRANIUM_LM_FOR_TPS_Apodemus_sylvaticus_MB_847_2018_34_16_skull_1476_1494_2935_000948818.resampled.tps", what = "character")
mcra <- matrix(as.numeric(cras[-c(1,length(cras))]), ncol = 3, byrow=T)
scaled_mcra <- mcra/centsiz(mcra[1:20,])$centroid_size
centrd <- centcoord(scaled_mcra[1:20,])
trsc_mcra <- cbind(scaled_mcra[,1]-centrd[1], scaled_mcra[,2]-centrd[2], scaled_mcra[,3]-centrd[3])

mshp <- pAcc$mshape

sup <- pPsup(trsc_mcra[1:20,], mshp)
rotrest <- trsc_mcra %*% sup$rotation
line_dors <- c(10,70,23,21,22,71,11, NA, 85:97,85, NA, 72:84,72, NA, 37, 98:100,103:101,39)
line_lat <- c(58:69,1,20,26,28,29,30,32,36,42,44,46,48,50,52,54,56,58, NA,83:92,93)
#plot(rotrest[,c(1,3)], type ="n", asp=1)
#text(rotrest[,c(1,3)])
#lines(rotrest[line_dors,c(1,3)])

#plot(rotrest[,3], -rotrest[,2], type ="n", asp=1)
#text(rotrest[,3], -rotrest[,2])
#lines(rotrest[line_lat,3], -rotrest[line_lat,2])

# 1. Cranium, landmarks only

o <- cbind(Mc,d)

prnd.cranium <- array.pruning(sp_tre, o, A = pAcc$rotated, split.t = "_", split.d = "_", positions = c(1,2))

data <- prnd.cranium$data[,1:60]
metadata <- prnd.cranium$data[,-c(1:60)]
tree <- prnd.cranium$tree

chisel <- grep("(C)", metadata$LIFESTYLE_bis)

lifestyle_simple <- metadata$LIFESTYLE_bis
lifestyle_simple[which(metadata$LIFESTYLE_bis == "arboreal/gliding")] <- "arboreal"
lifestyle_simple[which(metadata$LIFESTYLE_bis == "semiaquatic/burrower")] <- "burrower"
lifestyle_simple[which(metadata$LIFESTYLE_bis == "subterranean (C)")] <- "subterranean"
lifestyle_simple[which(metadata$LIFESTYLE_bis == "burrower (C)")] <- "burrower"
unique(lifestyle_simple)

A <- prnd.cranium$array
mshp <- mshape(A)

dimnames(A)[[3]] <- tree$tip.label
pca.geomorph <- gm.prcomp(A= A, phy = tree)

#Measure phylogenetic signal in the first 2 PCs

lambdaPC1 <- phylosig(tree = tree, x = pca.geomorph$x[,1], method ="lambda", test = T)
lambdaPC2 <- phylosig(tree = tree, x = pca.geomorph$x[,2], method ="lambda", test = T)
KPC1 <- phylosig(tree = tree, x = pca.geomorph$x[,1], method ="K", test = T)
KPC2 <- phylosig(tree = tree, x = pca.geomorph$x[,2], method ="K", test = T)
print("Pagel's lambda and Blomberg's K along PC1 and PC2 for the cranium LM dataset")
print(lambdaPC1)
print(lambdaPC2)
print(KPC1)
print(KPC2)

if (print.plots) {
pdf(paste("plots/", if (plot.names) {"named_"}, "craniumLM_morphospace_lifestyle.pdf", sep = ""), width = 7, height = 7)
#layout(matrix(1:2, ncol=2))
#plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F,anc.states=F), cex=0.1)
#points(pca.geomorph$x[,1:2], bg = alpha(palette(), 0.85)[as.factor(metadata$LIFESTYLE)], pch=21, cex = 2)
#title(main = "Cranium (LM only) morphospace (lifestyles)")
if (plot.names) {layout(matrix(c(1,3,3,6,2,3,3,6,0,4,5,6), ncol=4, byrow=T))} else {layout(matrix(c(1,3,3,2,3,3,0,4,5), ncol=3, byrow=T))}

par(mar = c(0.5,0.5,0.5,0.5))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 2)
#maxdef <- tps2d(rotrest[,c(1,3)], matr = mshp[,c(1,3)], matt = shp$max[,c(1,3)])
#mindef <- tps2d(rotrest[,c(1,3)], matr = mshp[,c(1,3)], matt = shp$min[,c(1,3)])

plot(mshp[,c(3,1)], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,3], shp$min[,3], rotrest[,3])), ylim = minetmax(c(shp$max[,1], shp$min[,1], rotrest[,1])))
#lines(maxdef[line_dors,3], mindef[line_dors,1], col="red")
lines(mshp[c(7,9,NA,6,8,NA,11,14,NA,10,13),3], mshp[c(7,9,NA,6,8,NA,11,14,NA,10,13),1], lwd=3)
lines(rotrest[line_dors,3], rotrest[line_dors,1])
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = c(shp$max[i,1], mshp[i,1]), col="red", lwd=3)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = c(shp$min[i,1], mshp[i,1]), col="blue", lwd=3)}

plot(mshp[,3], -mshp[,2], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,3], shp$min[,3], rotrest[,3])), ylim = minetmax(-c(shp$max[,2], shp$min[,2], rotrest[,2])))
lines(rotrest[line_lat,3], -rotrest[line_lat,2])
lines(mshp[c(7,9,NA,6,8,NA,11,14,NA,10,13),3], -mshp[c(7,9,NA,6,8,NA,11,14,NA,10,13),2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = -c(shp$max[i,2], mshp[i,2]), col="red", lwd=3)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = -c(shp$min[i,2], mshp[i,2]), col="blue", lwd=3)}

par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F, anc.states=F, edge.color = alpha("black", 0.2)), col = alpha(c("grey40", "orangered2", "limegreen"), 0.5)[as.factor(metadata$Clade)], pch=20, cex = 1)
points(pca.geomorph$x[which(lifestyle_simple=="burrower"),1:2], pch = 21, cex=2, bg = alpha(c("grey40", "orangered2", "limegreen"), 0.85)[as.factor(metadata$Clade[which(lifestyle_simple=="burrower")])])
points(pca.geomorph$x[which(lifestyle_simple=="subterranean"),1:2], pch = 24, cex=2, bg = alpha(c("grey40", "orangered2", "limegreen"), 0.85)[as.factor(metadata$Clade[which(lifestyle_simple=="subterranean")])])
points(pca.geomorph$x[chisel,1:2], pch = c(21,24)[as.factor(lifestyle_simple[chisel])], cex=2, lwd=3, col = "black")

if (plot.names) {
labels <- abbreviate(rownames(pca.geomorph$x))
text(pca.geomorph$x[,1:2], labels = labels, pos = 1, cex = 0.6)}

par(mar = c(0.5,0.5,0.5,0.5))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 1)
plot(mshp[,c(3,1)], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,3], shp$min[,3], rotrest[,3])), ylim = minetmax(c(shp$max[,1], shp$min[,1], rotrest[,1])))
lines(rotrest[line_dors,3], rotrest[line_dors,1])
lines(mshp[c(7,9,NA,6,8,NA,11,14,NA,10,13),3], mshp[c(7,9,NA,6,8,NA,11,14,NA,10,13),1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = c(shp$max[i,1], mshp[i,1]), col="red", lwd=3)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = c(shp$min[i,1], mshp[i,1]), col="blue", lwd=3)}

plot(mshp[,3], -mshp[,2], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,3], shp$min[,3], rotrest[,3])), ylim = minetmax(-c(shp$max[,2], shp$min[,2], rotrest[,2])))
lines(rotrest[line_lat,3], -rotrest[line_lat,2])
lines(mshp[c(7,9,NA,6,8,NA,11,14,NA,10,13),3], -mshp[c(7,9,NA,6,8,NA,11,14,NA,10,13),2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = -c(shp$max[i,2], mshp[i,2]), col="red", lwd=3)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = -c(shp$min[i,2], mshp[i,2]), col="blue", lwd=3)}

if (plot.names) {
par(mar = c(0,0,0,0))
labels <- abbreviate(rownames(pca.geomorph$x))
plot(rep(1,length(labels)), 1:length(labels), axes = F, xlab = "", ylab = "", type = "n")
text(rep(1,length(labels)), 1:length(labels), labels = sort(paste(labels, names(labels), sep = "-")), pos = 1, cex = 0.5)}

#title(main = "Cranium (LM only) morphospace (clades)")
dev.off()
}

ConvCran <- convSig(phy = pca.geomorph$phy, traits = pca.geomorph$x[,1:2], focaltaxa=rownames(pca.geomorph$x)[which(lifestyle_simple=="burrower" | lifestyle_simple=="subterranean")], nsim = NSIM)
print("Stayton's C1-C4 indices for PC1-2 of the cranium LM dataset, for burrowers and subterranean")
print(ConvCran)

ConvCranSub <- convSig(phy = pca.geomorph$phy, traits = pca.geomorph$x[,1:2], focaltaxa=rownames(pca.geomorph$x)[which(lifestyle_simple=="subterranean")], nsim = NSIM)
print("Stayton's C1-C4 indices for PC1-2 of the cranium LM dataset, for subterranean only")
print(ConvCranSub)

# 2. Cranium, landmarks and curves

o <- cbind(MC,dd)

mshp <- pAC$mshape

prnd.C <- array.pruning(sp_tre, o, A = pAC$rotated, split.t = "_", split.d = "_", positions = c(1,2))

data <- prnd.C$data[,1:435]
metadata <- prnd.C$data[,-c(1:435)]
tree <- prnd.C$tree

chisel <- grep("(C)", metadata$LIFESTYLE_bis)

lifestyle_simple <- metadata$LIFESTYLE_bis
lifestyle_simple[which(metadata$LIFESTYLE_bis == "arboreal/gliding")] <- "arboreal"
lifestyle_simple[which(metadata$LIFESTYLE_bis == "semiaquatic/burrower")] <- "burrower"
lifestyle_simple[which(metadata$LIFESTYLE_bis == "subterranean (C)")] <- "subterranean"
lifestyle_simple[which(metadata$LIFESTYLE_bis == "burrower (C)")] <- "burrower"
unique(lifestyle_simple)

A <- prnd.C$array
mshp <- mshape(A)

dimnames(A)[[3]] <- tree$tip.label
pca.geomorph <- gm.prcomp(A= A, phy = tree)

#Measure phylogenetic signal in the first 2 PCs
cran_lambdaPC1 <- phylosig(tree = tree, x = pca.geomorph$x[,1], method ="lambda", test = T)
cran_lambdaPC2 <- phylosig(tree = tree, x = pca.geomorph$x[,2], method ="lambda", test = T)
cran_KPC1 <- phylosig(tree = tree, x = pca.geomorph$x[,1], method ="K", test = T)
cran_KPC2 <- phylosig(tree = tree, x = pca.geomorph$x[,2], method ="K", test = T)
print("Pagel's lambda and Blomberg's K along PC1 and PC2 for the cranium curve dataset")
print(cran_lambdaPC1)
print(cran_lambdaPC2)
print(cran_KPC1)
print(cran_KPC2)

if (print.plots) {
pdf(paste("plots/", if (plot.names) {"named_"}, "craniumCURVES_morphospace_lifestyle.pdf", sep = ""), width = 7, height = 7)
#layout(matrix(1:2, ncol=2))
#plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F,anc.states=F), cex=0.1)
#points(pca.geomorph$x[,1:2], bg = alpha(palette(), 0.85)[as.factor(metadata$LIFESTYLE)], pch=21, cex = 2)
#title(main = "Cranium morphospace (lifestyles)")
if (plot.names) {layout(matrix(c(1,3,3,6,2,3,3,6,0,4,5,6), ncol=4, byrow=T))} else {layout(matrix(c(1,3,3,2,3,3,0,4,5), ncol=3, byrow=T))}

par(mar = c(0.5,0.5,0.5,0.5))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 2)
plot(mshp[,c(3,1)], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,3], shp$min[,3])), ylim = minetmax(c(shp$max[,1], shp$min[,1])))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = c(shp$max[i,1], mshp[i,1]), col="red", lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = c(shp$min[i,1], mshp[i,1]), col="blue", lwd=2)}

plot(mshp[,3], -mshp[,2], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black", 0.5), xlim = minetmax(c(shp$max[,3], shp$min[,3])), ylim = minetmax(-c(shp$max[,2], shp$min[,2])))
lines(mshp[lincranlat,3], -mshp[lincranlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = -c(shp$max[i,2], mshp[i,2]), col="red", lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = -c(shp$min[i,2], mshp[i,2]), col="blue", lwd=2)}

par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F, anc.states=F, edge.color = alpha("black", 0.2)), col = alpha(c("grey40", "#FFC700", "#19975D"), 0.5)[as.factor(metadata$Clade)], pch=20, cex = 1)
points(pca.geomorph$x[which(lifestyle_simple=="burrower"),1:2], pch = 21, cex=2, bg = alpha(c("grey40", "#FFC700", "#19975D"), 0.85)[as.factor(metadata$Clade[which(lifestyle_simple=="burrower")])])
points(pca.geomorph$x[which(lifestyle_simple=="subterranean"),1:2], pch = 24, cex=2, bg = alpha(c("grey40", "#FFC700", "#19975D"), 0.85)[as.factor(metadata$Clade[which(lifestyle_simple=="subterranean")])])
points(pca.geomorph$x[chisel,1:2], pch = c(21,24)[as.factor(lifestyle_simple[chisel])], cex=2, lwd=3, col = "black")

if (plot.names) {
labels <- abbreviate(rownames(pca.geomorph$x))
text(pca.geomorph$x[,1:2], labels = labels, pos = 1, cex = 0.6)}

par(mar = c(0.5,0.5,0.5,0.5))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 1)
plot(mshp[,c(3,1)], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,3], shp$min[,3])), ylim = minetmax(c(shp$max[,1], shp$min[,1])))
lines(mshp[lincrandors,3], mshp[lincrandors,1], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = c(shp$max[i,1], mshp[i,1]), col="red", lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = c(shp$min[i,1], mshp[i,1]), col="blue", lwd=2)}

plot(mshp[,3], -mshp[,2], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,3], shp$min[,3])), ylim = minetmax(-c(shp$max[,2], shp$min[,2])))
lines(mshp[lincranlat,3], -mshp[lincranlat,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,3], mshp[i,3]), y = -c(shp$max[i,2], mshp[i,2]), col="red", lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,3], mshp[i,3]), y = -c(shp$min[i,2], mshp[i,2]), col="blue", lwd=2)}

if (plot.names) {
par(mar = c(0,0,0,0))
labels <- abbreviate(rownames(pca.geomorph$x))
plot(rep(1,length(labels)), 1:length(labels), axes = F, xlab = "", ylab = "", type = "n")
text(rep(1,length(labels)), 1:length(labels), labels = sort(paste(labels, names(labels), sep = "-")), pos = 1, cex = 0.5)}

dev.off()
}

ConvCurv <- convSig(phy = pca.geomorph$phy, traits = pca.geomorph$x[,1:2], focaltaxa=rownames(pca.geomorph$x)[which(lifestyle_simple=="burrower" | lifestyle_simple=="subterranean")], nsim = NSIM)
print("Stayton's indices for PC1-2 of the cranium curve dataset, for burrowers and subterranean")
print(ConvCurv)

ConvCurvSub <- convSig(phy = pca.geomorph$phy, traits = pca.geomorph$x[,1:2], focaltaxa=rownames(pca.geomorph$x)[which(lifestyle_simple=="subterranean")], nsim = NSIM)
print("Stayton's indices for PC1-2 of the cranium curve dataset, for subterranean only")
print(ConvCurvSub)

# 3. Mandible space

o <- cbind(Mmd,ddd)

prnd.md <- array.pruning(sp_tre, o, A = pAmc$rotated, split.t = "_", split.d = "_", positions = c(1,2))

data <- prnd.md$data[,1:36]
metadata <- prnd.md$data[,-c(1:36)]
tree <- prnd.md$tree

A <- prnd.md$array
mshp <- mshape(A)

chisel <- grep("(C)", metadata$LIFESTYLE_bis)

lifestyle_simple <- metadata$LIFESTYLE_bis
lifestyle_simple[which(metadata$LIFESTYLE_bis == "arboreal/gliding")] <- "arboreal"
lifestyle_simple[which(metadata$LIFESTYLE_bis == "semiaquatic/burrower")] <- "burrower"
lifestyle_simple[which(metadata$LIFESTYLE_bis == "subterranean (C)")] <- "subterranean"
lifestyle_simple[which(metadata$LIFESTYLE_bis == "burrower (C)")] <- "burrower"
unique(lifestyle_simple)

dimnames(A)[[3]] <- tree$tip.label
pca.geomorph <- gm.prcomp(A= A, phy = tree)

mand_lambdaPC1 <- phylosig(tree = tree, x = pca.geomorph$x[,1], method ="lambda", test = T)
mand_lambdaPC2 <- phylosig(tree = tree, x = pca.geomorph$x[,2], method ="lambda", test = T)
mand_KPC1 <- phylosig(tree = tree, x = pca.geomorph$x[,1], method ="K", test = T)
mand_KPC2 <- phylosig(tree = tree, x = pca.geomorph$x[,2], method ="K", test = T)
print("Pagel's lambda and Blomberg's K for PC1-2 of the mandible morphospace")
print(mand_lambdaPC1)
print(mand_lambdaPC2)
print(mand_KPC1)
print(mand_KPC2)

if (print.plots) {
pdf(paste("plots/", if (plot.names) {"named_"}, "mandible_morphospace_lifestyle.pdf", sep = ""), width = 7, height = 7)
if (plot.names) {layout(matrix(c(1,3,3,6,2,3,3,6,0,4,5,6), ncol=4, byrow=T))} else {layout(matrix(c(1,3,3,2,3,3,0,4,5), ncol=3, byrow=T))}

par(mar = c(0.5,0.5,0.5,0.5))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 2)
plot(mshp[,1], -mshp[,2], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,1], shp$min[,1])), ylim = minetmax(-c(shp$max[,2], shp$min[,2])))
lines(mshp[linmd,1], -mshp[linmd,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,1], mshp[i,1]), y = -c(shp$max[i,2], mshp[i,2]), col="red", lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,1], mshp[i,1]), y = -c(shp$min[i,2], mshp[i,2]), col="blue", lwd=2)}

plot(mshp[,1], mshp[,3], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,1], shp$min[,1])), ylim = minetmax(c(shp$max[,3], shp$min[,3])))
lines(mshp[linmd,1], mshp[linmd,3], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,1], mshp[i,1]), y = c(shp$max[i,3], mshp[i,3]), col="red", lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,1], mshp[i,1]), y = c(shp$min[i,3], mshp[i,3]), col="blue", lwd=2)}

par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F, anc.states=F, edge.color = alpha("black", 0.2)), col = alpha(c("grey40", "#FFC700", "#19975D"), 0.5)[as.factor(metadata$Clade)], pch=20, cex = 1)
points(pca.geomorph$x[which(lifestyle_simple=="burrower"),1:2], pch = 21, cex=2, bg = alpha(c("grey40", "#FFC700", "#19975D"), 0.85)[as.factor(metadata$Clade[which(lifestyle_simple=="burrower")])])
points(pca.geomorph$x[which(lifestyle_simple=="subterranean"),1:2], pch = 24, cex=2, bg = alpha(c("grey40", "#FFC700", "#19975D"), 0.85)[as.factor(metadata$Clade[which(lifestyle_simple=="subterranean")])])
points(pca.geomorph$x[chisel,1:2], pch = c(21,24)[as.factor(lifestyle_simple[chisel])], cex=2, lwd=3, col = "black")

if (plot.names) {
labels <- abbreviate(rownames(pca.geomorph$x))
text(pca.geomorph$x[,1:2], labels = labels, pos = 1, cex = 0.6)}

par(mar = c(0.5,0.5,0.5,0.5))
shp <- shp.deform(shp = mshp, PCA = pca.geomorph, axis = 1)
plot(mshp[,1], -mshp[,2], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,1], shp$min[,1])), ylim = minetmax(-c(shp$max[,2], shp$min[,2])))
lines(mshp[linmd,1], -mshp[linmd,2], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,1], mshp[i,1]), y = -c(shp$max[i,2], mshp[i,2]), col="red", lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,1], mshp[i,1]), y = -c(shp$min[i,2], mshp[i,2]), col="blue", lwd=2)}

plot(mshp[,1], mshp[,3], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = minetmax(c(shp$max[,1], shp$min[,1])), ylim = minetmax(c(shp$max[,3], shp$min[,3])))
lines(mshp[linmd,1], mshp[linmd,3], lwd=3)
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$max[i,1], mshp[i,1]), y = c(shp$max[i,3], mshp[i,3]), col="red", lwd=2)}
for (i in 1:dim(mshp)[1]) {lines(x= c(shp$min[i,1], mshp[i,1]), y = c(shp$min[i,3], mshp[i,3]), col="blue", lwd=2)}

if (plot.names) {
par(mar = c(0,0,0,0))
labels <- abbreviate(rownames(pca.geomorph$x))
plot(rep(1,length(labels)), 1:length(labels), axes = F, xlab = "", ylab = "", type = "n")
text(rep(1,length(labels)), 1:length(labels), labels = sort(paste(labels, names(labels), sep = "-")), pos = 1, cex = 0.5)}

dev.off()
}

ConvMd <- convSig(phy = pca.geomorph$phy, traits = pca.geomorph$x[,1:2], focaltaxa=rownames(pca.geomorph$x)[which(lifestyle_simple=="burrower" | lifestyle_simple=="subterranean")], nsim = NSIM)
print("Stayton's indices along PC1-2 of mandible morphospace, for burrowers and subterranean")
print(ConvMd)

ConvMdSub <- convSig(phy = pca.geomorph$phy, traits = pca.geomorph$x[,1:2], focaltaxa=rownames(pca.geomorph$x)[which(lifestyle_simple=="subterranean")], nsim = NSIM)
print("Stayton's indices along PC1-2 of mandible morphospace, for subterranean only")
print(ConvMdSub)

###########################################################################
# Muscles Morphospaces, with lifestyles, and compute convergence measures #
###########################################################################

prnd.muscles <- data.pruning(sp_tre, lsr_match, split.t = "_", split.d = " ", positions = c(1,2))
phy <- prnd.muscles$pruned.tree
traits <- prnd.muscles$pruned.data[, c("DM","ePT", "iPT", "SM","T","ZM")]
lfstl <- as.factor(prnd.muscles$pruned.data$lifestyle_simple)
clade <- as.factor(prnd.muscles$pruned.data$clade)
chisel <- grep("(C)", prnd.muscles$pruned.data$lifestyle)

pca.geomorph <- gm.prcomp(A= traits, phy = phy)

muscle_lambdaPC1 <- phylosig(tree = tree, x = pca.geomorph$x[,1], method ="lambda", test = T)
muscle_lambdaPC2 <- phylosig(tree = tree, x = pca.geomorph$x[,2], method ="lambda", test = T)
muscle_KPC1 <- phylosig(tree = tree, x = pca.geomorph$x[,1], method ="K", test = T)
muscle_KPC2 <- phylosig(tree = tree, x = pca.geomorph$x[,2], method ="K", test = T)
print("Pagel's lambda and Blomberg's K for PC1-2 of the muscle LSR morphospace")
print(muscle_lambdaPC1)
print(muscle_lambdaPC2)
print(muscle_KPC1)
print(muscle_KPC2)

if (print.plots) {
pdf(paste("plots/", if (plot.names) {"named_"}, "muscle_mass_morphospace_lifestyle.pdf", sep =""), width = 6, height= 6)
par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F,
anc.states=F, edge.color=alpha("black", 0.2)), col = alpha(c("grey40", "#FFC700", "#19975D"), 0.5)[clade], pch=20, cex = 1)
points(pca.geomorph$x[which(lfstl=="burrower"),1:2], pch = 21, cex=2, bg = alpha(c("grey40", "#FFC700", "#19975D"), 0.85)[clade[which(lfstl=="burrower")]])
points(pca.geomorph$x[which(lfstl=="subterranean"),1:2], pch = 24, cex=2, bg = alpha(c("grey40", "#FFC700", "#19975D"), 0.85)[clade[which(lfstl=="subterranean")]])
points(pca.geomorph$x[chisel,1:2], pch = c(21,24)[droplevels(lfstl[chisel])], cex=2, lwd=3, col = "black")

if (plot.names) {
labels <- abbreviate(rownames(pca.geomorph$x))
text(pca.geomorph$x[,1:2], labels = labels, pos = 1, cex = 0.6)}

par(fig=c(0,0.35,0.65,1), new=T)
plot.loads(x=pca.geomorph$rotation)
box()
dev.off()
}

ConvMusclesMass <- convSig(phy = pca.geomorph$phy, traits = pca.geomorph$x[,1:2], focaltaxa=rownames(pca.geomorph$x)[which(lifestyle_simple=="burrower" | lifestyle_simple=="subterranean")], nsim = NSIM)
print("Stayton's indices for PC1-2 of muscle LSR morphospace, for burrowers and subterranean")
print(ConvMusclesMass)

ConvMusclesMassSub <- convSig(phy = pca.geomorph$phy, traits = pca.geomorph$x[,1:2], focaltaxa=rownames(pca.geomorph$x)[which(lifestyle_simple=="subterranean")], nsim = NSIM)
print("Stayton's indices for PC1-2 of muscle LSR morphospace, for subterranean only")
print(ConvMusclesMassSub)

# PCA for muscle proportions NOT phylogenetically corrected
pca_prop <- prcomp(prop_match[, c("DM", "ePT", "iPT", "SM","T","ZM")])

if (print.plots) {
plot(pca_prop$x[,1:2], col = alpha(c("grey40", "#FFC700", "#19975D"), 0.5)[as.factor(prop_match$clade)], pch=20, cex = 1)
points(pca_prop$x[which(prop_match$lifestyle_simple == "burrower"),1:2], pch = 21, cex=2, lwd=2, bg = alpha(c("grey40", "#FFC700", "#19975D"), 0.85)[as.factor(prop_match$clade[which(prop_match$lifestyle_simple == "burrower")])])
points(pca_prop$x[which(prop_match$lifestyle_simple == "subterranean"),1:2], pch = 24, cex=2, lwd=2, bg = alpha(c("grey40", "#FFC700", "#19975D"), 0.85)[as.factor(prop_match$clade[which(prop_match$lifestyle_simple == "subterranean")])])

text(pca_prop$x[which(prop_match$lifestyle_simple == "burrower"),1:2], labels = rownames(prop_match)[which(prop_match$lifestyle_simple == "burrower")], cex = 0.6, pos = c(2,4))
text(pca_prop$x[which(prop_match$lifestyle_simple == "subterranean"),1:2], labels = rownames(prop_match)[which(prop_match$lifestyle_simple == "subterranean")], cex = 0.6, pos = c(2,4))
}

# Phylomorphospace for muscle PROPORTIONS
prnd.muscles <- data.pruning(sp_tre, prop_match, split.t = "_", split.d = " ", positions = c(1,2))
phy <- prnd.muscles$pruned.tree
traits <- prnd.muscles$pruned.data[, c("DM","ePT", "iPT", "SM","T","ZM")]
lfstl <- as.factor(prnd.muscles$pruned.data$lifestyle_simple)
clade <- as.factor(prnd.muscles$pruned.data$clade)
chisel <- grep("(C)", prnd.muscles$pruned.data$lifestyle)
pca.geomorph <- gm.prcomp(A= traits, phy = phy)

if (print.plots) {
pdf(paste("plots/", "muscle_proportions_morphospace_lifestyle.pdf"), width =6, height=6)
par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F,
anc.states=F, edge.color=alpha("black", 0.2)), col = alpha(c("grey40", "#FFC700", "#19975D"), 0.5)[clade], pch=20, cex = 1)
points(pca.geomorph$x[which(lfstl=="burrower"),1:2], pch = 21, cex=2, bg = alpha(c("grey40", "#FFC700", "#19975D"), 0.85)[clade[which(lfstl=="burrower")]])
points(pca.geomorph$x[which(lfstl=="subterranean"),1:2], pch = 24, cex=2, bg = alpha(c("grey40", "#FFC700", "#19975D"), 0.85)[clade[which(lfstl=="subterranean")]])
points(pca.geomorph$x[chisel,1:2], pch = c(21,24)[droplevels(lfstl[chisel])], cex=2, lwd=3, col = "black")

if (plot.names) {
labels <- abbreviate(rownames(pca.geomorph$x))
text(pca.geomorph$x[,1:2], labels = labels, pos = 1, cex = 0.6)}

par(fig=c(0,0.35,0.65,1), new=T)
plot.loads(x=pca.geomorph$rotation)
box()
dev.off()
}

ConvMusclesProp <- convSig(phy = pca.geomorph$phy, traits = pca.geomorph$x[,1:2], focaltaxa=rownames(pca.geomorph$x)[which(lifestyle_simple=="burrower" | lifestyle_simple=="subterranean")], nsim = NSIM)
print("Stayton's indices for PC1-2 of muscle proportions morphospace for burrowers and subterranean")
print(ConvMusclesProp)
ConvMusclesPropSub <- convSig(phy = pca.geomorph$phy, traits = pca.geomorph$x[,1:2], focaltaxa=rownames(pca.geomorph$x)[which(lifestyle_simple=="subterranean")], nsim = NSIM)
print("Stayton's indices for PC1-2 of muscle proportions morphospace for subterranean only")
print(ConvMusclesPropSub)

#########################################################
# Write results of convergence analyses in a table file #
#########################################################

Convdf <- data.frame(ConvCran, ConvCranSub, ConvCurv, ConvCurvSub, ConvMd, ConvMdSub, ConvMusclesMass, ConvMusclesMassSub, ConvMusclesProp, ConvMusclesPropSub)
colnames(Convdf) <- c("Cranium (LM)", "P-value", "Cranium (LM Sub only)", "P-value", "Cranium (LM and curves)", "P-value", "Cranium (LM and curves Sub only)", "P-value", "Mandible", "P-value", "Mandible (Sub only)", "P-value", "Muscle Mass", "P-value", "Muscle Mass (Sub only)", "P-value", "Muscle Prop", "P-value", "Muscle Prop (Sub only)", "P-value")

write.table(x = Convdf, file = "tables/convergence_indices_stayton.csv", row.names = T, col.names = T)


#####################################################
# Test that chisel tooth diggers have longer fibers #
#####################################################

fib_lgt <- read.csv("tables/adductor_fib_length.csv", h = T, sep = ",", dec = ".", row.names = 1)
shortnam <- shorten.labels(rownames(fib_lgt), positions = c(1,2), split = " ")
fib_lgt <- apply(X = fib_lgt, 2 , FUN = tapply, INDEX = shortnam, mean)

burrow_gen <- unique(lifestyles$Genre[c(grep("burrower", x = lifestyles$LIFESTYLE), grep("subterranean", x = lifestyles$LIFESTYLE))])
index_burrowers <- grep(paste(burrow_gen, collapse="|"), x = rownames(fib_lgt))

index_chisel <- grep(paste(chisels, collapse="|"), x = rownames(fib_lgt))

chiz <- rep(F, nrow(fib_lgt))
chiz[index_chisel] <- T
#boxplot(fib_lgt ~ chiz)
#boxplot(df_lsr_fib ~ chiz)
df_fib <-log(fib_lgt/gmean_mass[which(names(gmean_mass) %in% rownames(fib_lgt))])

sp_fib <- shorten.labels(rownames(fib_lgt), positions = c(1,2), split = " ")
sp_avg_lgt <- avg.matrix.fac(fib_lgt, INDEX = as.factor(sp_fib))
sisi <- c(sp_avg_csize[match(rownames(sp_avg_lgt), names(sp_avg_csize))])
sp_sizcorrect_lgt <- na.omit(sp_avg_lgt/sisi)

id_burr <- grep(paste(burrow_gen, collapse="|"), x = rownames(sp_sizcorrect_lgt))
id_chisel <- grep(paste(chisels, collapse="|"), x = rownames(sp_sizcorrect_lgt))
burrnonchiz <- id_burr[!id_burr %in% id_chisel]

chizz <- rep(F, nrow(sp_sizcorrect_lgt))
chizz[id_chisel] <- T
#boxplot(sp_sizcorrect_lgt ~ chizz)

if (print.plots & !plot.names) {
pdf(paste("plots/", "boxplot_fiber_length_chisel.pdf", sep = ""), width = 3.5, height = 6)
par(mar = c(15,5,0.5,0.5))
boxplot(
sp_sizcorrect_lgt[id_chisel,1],
sp_sizcorrect_lgt[-id_chisel,1],
sp_sizcorrect_lgt[id_chisel,2],
sp_sizcorrect_lgt[-id_chisel,2],
sp_sizcorrect_lgt[id_chisel,3],
sp_sizcorrect_lgt[-id_chisel,3],
sp_sizcorrect_lgt[id_chisel,4],
sp_sizcorrect_lgt[-id_chisel,4],
sp_sizcorrect_lgt[id_chisel,5],
sp_sizcorrect_lgt[-id_chisel,5],
sp_sizcorrect_lgt[id_chisel,6],
sp_sizcorrect_lgt[-id_chisel,6],
ylab = "Size-corrected fiber length",
col = rep(cols[c(2,5,6,1,4,3)],each=2),
names = c("Deep masseter (chisel)", "Deep masseter (non-chisel)", "External Pterygoid (chisel)", "External Pterygoid (non-chisel)","Internal Pterygoid (chisel)", "Internal Pterygoid (non-chisel)", "Superficial masseter (chisel)", "Superficial masseter (non-chisel)", "Temporal (chisel)", "Temporal (non-chisel)", "Zygomaticomandibularis (chisel)", "Zygomaticomandibularis (non-chisel)"),
las = 3)
dev.off()
}

print("T-tests of individual muscle fiber length differences between chisel-tooth diggers VS all other taxa")
for (i in 1:dim(sp_sizcorrect_lgt)[2]) {
print(t.test(sp_sizcorrect_lgt[id_chisel,i], sp_sizcorrect_lgt[-id_chisel,i]))
}

if (print.plots & !plot.names) {
pdf(paste("plots/", "boxplot_FL_scratch_vs_chisel.pdf", sep = ""), width = 3.5, height = 6)
par(mar = c(15,5,0.5,0.5))
boxplot(
sp_sizcorrect_lgt[id_chisel,1],
sp_sizcorrect_lgt[burrnonchiz,1],
sp_sizcorrect_lgt[id_chisel,2],
sp_sizcorrect_lgt[burrnonchiz,2],
sp_sizcorrect_lgt[id_chisel,3],
sp_sizcorrect_lgt[burrnonchiz,3],
sp_sizcorrect_lgt[id_chisel,4],
sp_sizcorrect_lgt[burrnonchiz,4],
sp_sizcorrect_lgt[id_chisel,5],
sp_sizcorrect_lgt[burrnonchiz,5],
sp_sizcorrect_lgt[id_chisel,6],
sp_sizcorrect_lgt[burrnonchiz,6],
ylab = "Size-corrected fiber length",
col = rep(cols[c(2,5,6,1,4,3)],each=2),
names = c("Deep masseter (chisel)", "Deep masseter (non-chisel)", "External Pterygoid (chisel)", "External Pterygoid (non-chisel)","Internal Pterygoid (chisel)", "Internal Pterygoid (non-chisel)", "Superficial masseter (chisel)", "Superficial masseter (non-chisel)", "Temporal (chisel)", "Temporal (non-chisel)", "Zygomaticomandibularis (chisel)", "Zygomaticomandibularis (non-chisel)"),
las = 3)
dev.off()
}

print("T-tests of individual muscle fiber length differences between chisel-tooth diggers VS non-chisel burrowing taxa")
for (i in 1:dim(sp_sizcorrect_lgt)[2]) {
print(t.test(sp_sizcorrect_lgt[id_chisel,i], sp_sizcorrect_lgt[burrnonchiz,i]))
}

######################################
# Combine morphospaces in one figure #
######################################

#prepare cranium objects
oC <- cbind(MC,dd)
prndC <- array.pruning(sp_tre, oC, A = pAC$rotated, split.t = "_", split.d = "_", positions = c(1,2))
dC <- prndC$data[,1:435]
mC <- prndC$data[,-c(1:435)]
treeC <- prndC$tree
chiselC <- grep("(C)", mC$LIFESTYLE_bis)
lC <- mC$LIFESTYLE_bis
lC[which(mC$LIFESTYLE_bis == "arboreal/gliding")] <- "arboreal"
lC[which(mC$LIFESTYLE_bis == "semiaquatic/burrower")] <- "burrower"
lC[which(mC$LIFESTYLE_bis == "subterranean (C)")] <- "subterranean"
lC[which(mC$LIFESTYLE_bis == "burrower (C)")] <- "burrower"
Acranium <- prndC$array
mshpC <- mshape(Acranium)
dimnames(Acranium)[[3]] <- treeC$tip.label
pca.geomorphC <- gm.prcomp(A=Acranium, phy = treeC)

#prepare mandible objects
oM <- cbind(Mmd,ddd)
prndM <- array.pruning(sp_tre, oM, A = pAmc$rotated, split.t = "_", split.d = "_", positions = c(1,2))
dM <- prndM$data[,1:36]
mM <- prndM$data[,-c(1:36)]
treeM <- prndM$tree
AM <- prndM$array
mshpM <- mshape(AM)
chiselM <- grep("(C)", mM$LIFESTYLE_bis)
lM <- mM$LIFESTYLE_bis
lM[which(mM$LIFESTYLE_bis == "arboreal/gliding")] <- "arboreal"
lM[which(mM$LIFESTYLE_bis == "semiaquatic/burrower")] <- "burrower"
lM[which(mM$LIFESTYLE_bis == "subterranean (C)")] <- "subterranean"
lM[which(mM$LIFESTYLE_bis == "burrower (C)")] <- "burrower"
dimnames(AM)[[3]] <- treeM$tip.label
pca.geomorphM <- gm.prcomp(A = AM, phy = treeM)

#prepare muscle objects
prndMU <- data.pruning(sp_tre, lsr_match, split.t = "_", split.d = " ", positions = c(1,2))
treeMU <- prndMU$pruned.tree
traitsMU <- prndMU$pruned.data[, c("DM","ePT", "iPT", "SM","T","ZM")]
lfstlMU <- as.factor(prndMU$pruned.data$lifestyle_simple)
cladeMU <- as.factor(prndMU$pruned.data$clade)
chiselMU <- grep("(C)", prndMU$pruned.data$lifestyle)
pca.geomorphMU <- gm.prcomp(A= traitsMU, phy = treeMU)

#Plot
if (print.plots & !plot.names) {
pdf(paste("plots/", "Combined_morphospaces_figure.pdf", sep=""), height=10, width=9)

#Set various parameters
layout(matrix(c(1,1,1,1,2,2,2,2,1,1,1,1,2,2,2,2,3,3,3,3,4,5,6,7,3,3,3,3,8,9,10,10), ncol=4))
colclade <- c("grey40", "orangered2", "limegreen")
alpha1 <- 0.6
alpha2 <- 0.9
pchbur <- 21
pchsub <- 24
cexall <- 1
cexbs <- 2
cexlab <- 1.8
cexmain <- 2

par(mar = c(4.5,4.5,4,0.5))
plot(pca.geomorphC, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F, anc.states=F, edge.color = alpha("black", 0.2)), col = alpha(colclade, alpha1)[as.factor(mC$Clade)], pch=20, cex = cexall, cex.lab = cexlab)
title(main = "A. Cranium morphospace", adj = 0, cex.main = cexmain, line = 1)
points(pca.geomorphC$x[which(lC=="burrower"),1:2], pch = pchbur, cex=cexbs, bg = alpha(colclade, alpha2)[as.factor(mC$Clade[which(lC=="burrower")])])
points(pca.geomorphC$x[which(lC=="subterranean"),1:2], pch = pchsub, cex=cexbs, bg = alpha(colclade, alpha2)[as.factor(mC$Clade[which(lC=="subterranean")])])
points(pca.geomorphC$x[chiselC,1:2], pch = c(pchbur,pchsub)[as.factor(lC[chiselC])], cex=cexbs, lwd=3, col = "black")
legend("bottomleft", bty = "n", pch = c(20,20,20,21,24,21), col = c(colclade, "black","black","black"), pt.cex = c(rep(cexall,3),rep(cexbs,3)), pt.bg = c(NA,NA,NA,"white","white","white"), pt.lwd = c(NA,NA,NA,1,1,3), legend = c("Hystricomorpha", "Supramyomorpha", "Sciuromorpha", "Burrower", "Subterranean", "Chisel-tooth digger"), cex = 1.2)

plot(pca.geomorphM, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F, anc.states=F, edge.color = alpha("black", 0.2)), col = alpha(colclade, alpha1)[as.factor(mC$Clade)], pch=20, cex = cexall, cex.lab = cexlab)
title(main = "B. Mandible morphospace", adj = 0, cex.main = cexmain, line = 1)
points(pca.geomorphM$x[which(lM=="burrower"),1:2], pch = pchbur, cex = cexbs, bg = alpha(colclade, alpha2)[as.factor(mM$Clade[which(lM=="burrower")])])
points(pca.geomorphM$x[which(lM=="subterranean"),1:2], pch = pchsub, cex = cexbs, bg = alpha(colclade, alpha2)[as.factor(mM$Clade[which(lM=="subterranean")])])
points(pca.geomorphM$x[chiselM,1:2], pch = c(pchbur, pchsub)[as.factor(lM[chiselM])], cex=cexbs, lwd=3, col = "black")

plot(pca.geomorphMU, phylo=T, phylo.par=list(tip.labels=F, node.labels=F, tip.labels=F, anc.states=F, edge.color = alpha("black", 0.2)), col = alpha(colclade, alpha1)[as.factor(cladeMU)], pch=20, cex = cexall, cex.lab = cexlab)
title(main = "C. Muscle mass morphospace", adj = 0, cex.main = cexmain, line = 1)
points(pca.geomorphMU$x[which(lfstlMU=="burrower"),1:2], pch = pchbur, cex = cexbs, bg = alpha(colclade, alpha2)[as.factor(cladeMU[which(lfstlMU=="burrower")])])
points(pca.geomorphMU$x[which(lfstlMU=="subterranean"),1:2], pch = pchsub, cex = cexbs, bg = alpha(colclade, alpha2)[as.factor(cladeMU[which(lfstlMU=="subterranean")])])
points(pca.geomorphMU$x[chiselMU,1:2], pch = c(pchbur, pchsub)[droplevels(as.factor(lfstlMU[chiselMU]))], cex = cexbs, lwd=3, col = "black")

par(mar = c(0.5,0.5,0.5,0.5))
shpC <- shp.deform(shp = mshpC, PCA = pca.geomorphC, axis = 1)
plot(mshpC[,c(3,1)], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1.5, col=alpha("black",0.5), xlim = maxetmin(c(shpC$max[,3], shpC$min[,3])), ylim = minetmax(c(shpC$max[,1], shpC$min[,1])), bty = "n")
#lines(mshpC[c(7,9,NA,6,8,NA,11,14,NA,10,13),3], mshpC[c(7,9,NA,6,8,NA,11,14,NA,10,13),1], lwd=2)
#lines(rotrest[line_dors,3], rotrest[line_dors,1])
for (i in 1:dim(mshpC)[1]) {lines(x=c(shpC$max[i,3], mshpC[i,3]), y = c(shpC$max[i,1], mshpC[i,1]), col="red", lwd=2)}
for (i in 1:dim(mshpC)[1]) {lines(x=c(shpC$min[i,3], mshpC[i,3]), y = c(shpC$min[i,1], mshpC[i,1]), col="blue", lwd=2)}
lines(mshpC[lincrandors,3], mshpC[lincrandors,1], lwd=2, col="grey60")
legend("bottomleft", legend = "PC1", bty="n", cex = 1.3)
legend("topleft", legend = "D.", bty="n", cex = 1.8, text.font = 2)
plot(mshpC[,3], -mshpC[,2], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = maxetmin(c(shpC$max[,3], shpC$min[,3])), ylim = minetmax(-c(shpC$max[,2], shpC$min[,2])), bty = "n")
#lines(rotrest[line_lat,3], -rotrest[line_lat,2])
#lines(mshpC[c(7,9,NA,6,8,NA,11,14,NA,10,13),3], -mshpC[c(7,9,NA,6,8,NA,11,14,NA,10,13),2], lwd=2)
for (i in 1:dim(mshpC)[1]) {lines(x= c(shpC$max[i,3], mshpC[i,3]), y = -c(shpC$max[i,2], mshpC[i,2]), col="red", lwd=3)}
for (i in 1:dim(mshpC)[1]) {lines(x= c(shpC$min[i,3], mshpC[i,3]), y = -c(shpC$min[i,2], mshpC[i,2]), col="blue", lwd=3)}
lines(mshpC[lincranlat,3], -mshpC[lincranlat,2], lwd=2, col="grey60")
legend("topleft", legend = "E.", bty="n", cex = 1.8, text.font = 2)

shpC <- shp.deform(shp = mshpC, PCA = pca.geomorphC, axis = 2)
plot(mshpC[,c(3,1)], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 1.5, col=alpha("black",0.5), xlim = maxetmin(c(shpC$max[,3], shpC$min[,3])), ylim = minetmax(c(shpC$max[,1], shpC$min[,1])), bty = "n")
#lines(mshpC[c(7,9,NA,6,8,NA,11,14,NA,10,13),3], mshpC[c(7,9,NA,6,8,NA,11,14,NA,10,13),1], lwd=2)
#lines(rotrest[line_dors,3], rotrest[line_dors,1])
for (i in 1:dim(mshpC)[1]) {lines(x=c(shpC$max[i,3], mshpC[i,3]), y = c(shpC$max[i,1], mshpC[i,1]), col="red", lwd=2)}
for (i in 1:dim(mshpC)[1]) {lines(x=c(shpC$min[i,3], mshpC[i,3]), y = c(shpC$min[i,1], mshpC[i,1]), col="blue", lwd=2)}
lines(mshpC[lincrandors,3], mshpC[lincrandors,1], lwd=2, col="grey60")
legend("bottomleft", legend = "PC2", bty="n", cex = 1.3)
legend("topleft", legend = "F.", bty="n", cex = 1.8, text.font = 2)
plot(mshpC[,3], -mshpC[,2], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = maxetmin(c(shpC$max[,3], shpC$min[,3])), ylim = minetmax(-c(shpC$max[,2], shpC$min[,2])), bty = "n")
#lines(rotrest[line_lat,3], -rotrest[line_lat,2])
#lines(mshpC[c(7,9,NA,6,8,NA,11,14,NA,10,13),3], -mshpC[c(7,9,NA,6,8,NA,11,14,NA,10,13),2], lwd=2)
for (i in 1:dim(mshpC)[1]) {lines(x= c(shpC$max[i,3], mshpC[i,3]), y = -c(shpC$max[i,2], mshpC[i,2]), col="red", lwd=3)}
for (i in 1:dim(mshpC)[1]) {lines(x= c(shpC$min[i,3], mshpC[i,3]), y = -c(shpC$min[i,2], mshpC[i,2]), col="blue", lwd=3)}
lines(mshpC[lincranlat,3], -mshpC[lincranlat,2], lwd=2, col="grey60")
legend("topleft", legend = "G.", bty="n", cex = 1.8, text.font = 2)

shpM <- shp.deform(shp = mshpM, PCA = pca.geomorphM, axis = 1)
plot(mshpM[,1], -mshpM[,2], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = maxetmin(c(shpM$max[,1], shpM$min[,1])), ylim = minetmax(-c(shpM$max[,2], shpM$min[,2])), bty = "n")
for (i in 1:dim(mshpM)[1]) {lines(x= c(shpM$max[i,1], mshpM[i,1]), y = -c(shpM$max[i,2], mshpM[i,2]), col="red", lwd=2)}
for (i in 1:dim(mshpM)[1]) {lines(x= c(shpM$min[i,1], mshpM[i,1]), y = -c(shpM$min[i,2], mshpM[i,2]), col="blue", lwd=2)}
lines(mshpM[linmd,1], -mshpM[linmd,2], lwd=2, col="grey60")
legend("bottomleft", legend = "PC1", bty="n", cex = 1.3)
legend("topleft", legend = "H.", bty="n", cex = 1.8, text.font = 2)

shpM <- shp.deform(shp = mshpM, PCA = pca.geomorphM, axis = 2)
plot(mshpM[,1], -mshpM[,2], asp = 1, xaxt ="n", yaxt = "n", xlab="",ylab="", pch = 20, cex = 2, col=alpha("black",0.5), xlim = maxetmin(c(shpM$max[,1], shpM$min[,1])), ylim = minetmax(-c(shpM$max[,2], shpM$min[,2])), bty = "n")

for (i in 1:dim(mshpM)[1]) {lines(x= c(shpM$max[i,1], mshpM[i,1]), y = -c(shpM$max[i,2], mshpM[i,2]), col="red", lwd=2)}
for (i in 1:dim(mshpM)[1]) {lines(x= c(shpM$min[i,1], mshpM[i,1]), y = -c(shpM$min[i,2], mshpM[i,2]), col="blue", lwd=2)}
lines(mshpM[linmd,1], -mshpM[linmd,2], lwd=2, col="grey60")
legend("bottomleft", legend = "PC2", bty="n", cex = 1.3)
legend("topleft", legend = "I.", bty="n", cex = 1.8, text.font = 2)

plot.loads(pca.geomorphMU$rotation, cex = 1.3)
legend("topleft", legend = "J. Muscle loadings", bty="n", cex = 1.8, text.font = 2)

dev.off()
}

#Create generic level phylogeny

gen_tre_lad <- ladderize(gen_tre, right=F)
lifestyle_gen <- disparat$LIFESTYLE_bis[match(gen_tre_lad$tip.label, disparat$Genre)]
burrowers <- rep(1, length(gen_tre_lad$tip.label))
burrowers[grep("burrower", lifestyle_gen)] <- 2
burrowers[grep("subterr", lifestyle_gen)] <- 2
edgwidth <- find.edgecols(phy = gen_tre_lad, factor = as.factor(burrowers), color = 1:2)

suborders <- as.factor(c(rep("Sciuromorpha", 28), 
		rep("Hystricomorpha", 42),
		rep("Supramyomorpha", 150)))
edg <- find.edgecols(phy = gen_tre_lad, factor = suborders, color = c("grey40", "orangered2", "limegreen"))
fnt <- rep(3, length(gen_tre_lad$tip.label))
fnt[grep("subterr", lifestyle_gen)] <- 4

wdth <- as.numeric(edgwidth[[2]])
wdth[which(wdth==2)] <- 2.5

if (print.plots) {
pdf(paste("plots/","generic_phylo_disparat.pdf", sep=""), height = 2, width = 8)
par(mar = c(0,0,0,0))
plot(gen_tre_lad, show.tip.label=T, edge.color=edg[[2]], edge.width=wdth, cex=0.3, font = fnt, direction = "upwards")
dev.off()
}

# Stop sinking if R was sinking
if (sinking) {sink()}

