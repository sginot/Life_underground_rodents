#Required package
library(ape)
library(phytools)
library(corrplot)
library(geomorph)
source("CUSTOM_R_FUNCTIONS.R")


# Opt to print plots or not:
print.plots <- F
if (print.plots & !dir.exists("plots")) {dir.create("plots")} #Create output folder if it does not exist

#Load phylogenetic trees
gen_tre <- read.tree("Disparat_generic_tree")

# Load adductor muscle data
datadduct <- read.csv2("tables/datadduct_taxo.csv", sep =",", h = T)

#Taxonomic data
taxo <- datadduct[,1:3]
lgen <- unlist(lapply(X = strsplit(x = datadduct$species, split = " "),
                      FUN = f <- function(x) {x[[1]]}))
taxo <- cbind(taxo, lgen)
utaxo <- data.frame(apply(X = taxo, 2 , FUN = tapply, INDEX = lgen, unique))

#Load data tables produced by previous script
mscl_lsr_mass <- read.csv("tables/adductor_LSR_mass.csv", h = T, sep = ",", dec = ".", row.names = 1)
fib_lgt <- read.csv("tables/adductor_fib_length.csv", h = T, sep = ",", dec = ".", row.names = 1)
lsr_fib_lgt <- read.csv("tables/adductor_LSR_fib_length.csv", h = T, sep = ",", dec = ".", row.names = 1)
mscl_lsr_pcsa <- read.csv("tables/adductor_LSR_pcsa.csv", h = T, sep = ",", dec = ".", row.names = 1)

#Compute generic average data (only changes species for which several specimens are present in the original data).
#1. Reduce row names to generic names
shortnam <- shorten.labels(rownames(mscl_lsr_pcsa), positions = 1, split = " ")
#2. Apply mean function across data frame based on generic name as a factor.
mscl_lsr_pcsa <- apply(X = mscl_lsr_pcsa, 2 , FUN = tapply, INDEX = shortnam, mean)

shortnam <- shorten.labels(rownames(mscl_lsr_mass), positions = 1, split = " ")
mscl_lsr_mass <- apply(X = mscl_lsr_mass, 2 , FUN = tapply, INDEX = shortnam, mean)

shortnam <- shorten.labels(rownames(lsr_fib_lgt), positions = 1, split = " ")
lsr_fib_lgt <- apply(X = lsr_fib_lgt, 2 , FUN = tapply, INDEX = shortnam, mean)

shortnam <- shorten.labels(rownames(fib_lgt), positions = 1, split = " ")
fib_lgt <- apply(X = fib_lgt, 2 , FUN = tapply, INDEX = shortnam, mean)

#Prune the generic level pcsa LSR data
prnd.pcsa <- data.pruning(gen_tre, mscl_lsr_pcsa, split.t = "_", split.d = " ", positions = 1)
#Just a check; not necessary
#par(mar = c(0,0,0,0))
#plot(prnd.pcsa$pruned.tree, type = "fan")
#Prune muscle mass data
prnd.mass <- data.pruning(gen_tre, mscl_lsr_mass, split.t = "_", split.d = " ", positions = 1)
#Prune fiber length data
prnd.lsr.fib <- data.pruning(gen_tre, lsr_fib_lgt, split.t = "_", split.d = " ", positions = 1)
prnd.fib <- data.pruning(gen_tre, fib_lgt, split.t = "_", split.d = " ", positions = 1)

if (print.plots) {
plot.mapd.traits(prnd = prnd.fib, 
                 file = "maps_fiber_lgt_gen.pdf", 
                 titles = c("DM", "PT", "SM", "T", "ZM"))
plot.mapd.traits(prnd = prnd.lsr.fib, 
                 file = "maps_LSR_fiber_lgt_gen.pdf", 
                 titles = c("DM", "PT", "SM", "T", "ZM"))
plot.mapd.traits(prnd = prnd.pcsa, 
                 file = "maps_LSR_PCSA_gen.pdf", 
                 titles = c("DM", "PT", "SM", "T", "ZM"))
plot.mapd.traits(prnd = prnd.mass, 
                 file = "maps_LSR_mass_gen.pdf",
                 titles = c("DM", "PT", "SM", "T", "ZM"))
}

#Phylogenetic PCA function
phyp <- function(prnd, 
                 tax = utaxo, 
                 PCs = c(1,2), 
                 map = T,
                 tax.level = "clade", 
                 make.plot = T,
                 add.labels = T,
                 add.loads = T) {
  
  dat <- prnd$pruned.data
  tre <- prnd$pruned.tree
  
  #Some duplicates may appear if several species were sampled in the same genus. This creates a bug du to matrix singularity.
  #Find duplicates and 1) remove duplicated tips of the tree, and 2) Average the traits values across duplicates.
  rn <- rownames(dat)
  
  if (any(table(rn) > 1)) {
    duplic <- names(which(table(rn) > 1))
    datnodup <- apply(dat, 2, tapply, INDEX = rn, mean)
  
    tipstorm <- rep(NA, length(duplic))
    for (i in 1:length(duplic)) {
      tipstorm[i] <- which(tre$tip.label == duplic[i])[1]
    }
  
    tre <- drop.tip(tip = tipstorm, tre)
    dat <- datnodup[match(tre$tip.label,rownames(datnodup)),]
  }

  #Get taxonomy for each observation
  nam <- rownames(dat)
  snam <- strsplit(nam, split = " ")
  g <- unlist(lapply(snam, FUN = f <- function(x) {x[1]}))
  ta <- tax[match(g, tax$lgen),]
  
  #Map taxonomic group as a discrete character on the tree
  taxfac <- as.factor(ta[,tax.level])
  taxfac <- setNames(taxfac, tre$tip.label)
  if (map) {tresimmap <- make.simmap(tre, x=taxfac)}
  
  #Compute PCA and extract scores
  phylpca <- phyl.pca(tree = tre, Y = dat)
  pcs <- scores(phylpca)
  lo <- phylpca$L[,PCs]
  varpct <- diag(phylpca$Eval)/sum(diag(phylpca$Eval))*100
  vv <- paste(names(varpct), round(varpct,2), "%")
  
  #Plotting
  if (make.plot) {
  #Define vector of colors
    cols <- setNames(palette()[1:length(levels(taxfac))], levels(taxfac))
    phylomorphospace(tree = if(map) {tresimmap} else {tre}, pcs[,PCs], ftype = "off", node.by.map = T, 
                   bty="n",
                   colors=cols,
                   node.size=c(0,1.2),
                   xlab = vv[PCs[1]],
                   ylab = vv[PCs[2]])
     
    if (!map) {points(pcs[,PCs], col=cols[taxfac], pch=20)}
                  
    if (add.labels) {text(pcs[,PCs], g, cex = 0.5, pos = 2)}
    
    if (add.loads) {
    par(fig = c(0.7,1,0.7,1), new = T, mar=c(0,0,4,2))
    plot.loads(lo)
    }
    }
  
if(map) {return(list(PCA = phylpca, Tree.taxo.mapped = tresimmap))}
else {return(list(PCA = phylpca))}}

#Next section creates morphospace, but takes a few minutes, optional.
#phypout <- phyp(prnd = prnd.mass, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = F)

#pdf(file = paste("plots", "LSR_fib_gen_phylospace12.pdf", sep = "/"))
#phyp(prnd = prnd.lsr.fib, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "LSR_fib_gen_phylospace34.pdf", sep = "/"))
#phyp(prnd = prnd.lsr.fib, tax = utaxo, PCs = c(3,4), tax.level = "clade", make.plot = T)
#dev.off()

#pdf(file = paste("plots", "fib_gen_phylospace12.pdf", sep = "/"))
#phyp(prnd = prnd.fib, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "fib_gen_phylospace34.pdf", sep = "/"))
#phyp(prnd = prnd.fib, tax = utaxo, PCs = c(3,4), tax.level = "clade", make.plot = T)
#dev.off()

#pdf(file = paste("plots", "pcsa_gen_phylospace12.pdf", sep = "/"))
#phyp(prnd = prnd.pcsa, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "pcsa_gen_phylospace34.pdf", sep = "/"))
#phyp(prnd = prnd.pcsa, tax = utaxo, PCs = c(3,4), tax.level = "clade", make.plot = T)
#dev.off()

#pdf(file = paste("plots", "mass_gen_phylospace12.pdf", sep = "/"))
#phyp(prnd = prnd.mass, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "mass_gen_phylospace34.pdf", sep = "/"))
#phyp(prnd = prnd.mass, tax = utaxo, PCs = c(3,4), tax.level = "clade", make.plot = T)
#dev.off()

#pdf(file = paste("plots", "muscle_spaces_PC12_nonames.pdf", sep = "/"), width = 15, height = 6)
#layout(matrix(1:3, ncol=3))
#par(mar=c(4,4,4,1))
#phyp(prnd = prnd.lsr.fib, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T, add.labels=F, add.loads=F)
#phyp(prnd = prnd.mass, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T, add.labels=F, add.loads=F)
#phyp(prnd = prnd.pcsa, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T, add.labels=F, add.loads=F)
#dev.off()


#Phylogenetic independent contrasts
pics.df <- function(prnd) {
  tre <- prnd$pruned.tree
  dat <- prnd$pruned.data
  pics <- apply(dat, 2, pic, phy = tre)
  return(pics)
}

fib_pics <- pics.df(prnd = prnd.fib)
lsr_fib_pics <- pics.df(prnd = prnd.lsr.fib)
pcsa_pics <- pics.df(prnd = prnd.pcsa)
mass_pics <- pics.df(prnd = prnd.mass)

write.csv(x = lsr_fib_pics, file = paste("tables", "gen_pics_LSR_fiber_length.csv", sep = "/"))
write.csv(x = fib_pics, file = paste("tables", "gen_pics_fiber_length.csv", sep = "/"))
write.csv(x = pcsa_pics, file = paste("tables", "gen_pics_pcsa.csv", sep = "/"))
write.csv(x = mass_pics, file = paste("tables", "gen_pics_muscle_mass.csv", sep = "/"))

############################
# Clade by clade phylo plots
############################

utaxo_mous <- utaxo[utaxo$clade == "MOUS",]
utaxo_squir <- utaxo[utaxo$clade == "SQUIR",]
utaxo_cteno <- utaxo[utaxo$clade == "CTENO",]

mass_mous <- mscl_lsr_mass[utaxo$clade == "MOUS",]
pcsa_mous <- mscl_lsr_pcsa[utaxo$clade == "MOUS",]
fib_mous <- lsr_fib_lgt[utaxo$clade == "MOUS",]

mass_squir <- mscl_lsr_mass[utaxo$clade == "SQUIR",]
pcsa_squir <- mscl_lsr_pcsa[utaxo$clade == "SQUIR",]
fib_squir <- lsr_fib_lgt[utaxo$clade == "SQUIR",]

mass_cteno <- mscl_lsr_mass[utaxo$clade == "CTENO",]
pcsa_cteno <- mscl_lsr_pcsa[utaxo$clade == "CTENO",]
fib_cteno <- lsr_fib_lgt[utaxo$clade == "CTENO",]

prnd.pcsa_mous <- data.pruning(gen_tre, pcsa_mous, split.t = "_", split.d = " ", positions = 1)
prnd.mass_mous <- data.pruning(gen_tre, mass_mous, split.t = "_", split.d = " ", positions = 1)
prnd.fib_mous <- data.pruning(gen_tre, fib_mous, split.t = "_", split.d = " ", positions = 1)

prnd.pcsa_squir <- data.pruning(gen_tre, pcsa_squir, split.t = "_", split.d = " ", positions = 1)
prnd.mass_squir <- data.pruning(gen_tre, mass_squir, split.t = "_", split.d = " ", positions = 1)
prnd.fib_squir <- data.pruning(gen_tre, fib_squir, split.t = "_", split.d = " ", positions = 1)

prnd.pcsa_cteno <- data.pruning(gen_tre, pcsa_cteno, split.t = "_", split.d = " ", positions = 1)
prnd.mass_cteno <- data.pruning(gen_tre, mass_cteno, split.t = "_", split.d = " ", positions = 1)
prnd.fib_cteno <- data.pruning(gen_tre, fib_cteno, split.t = "_", split.d = " ", positions = 1)

#Next section creates morphospace, but takes a few minutes, optional.
#pdf(file = paste("plots", "muscle_space_mouse.pdf", sep = "/"), width = 15, height = 6)
#palette(rainbow(12))
#layout(matrix(1:3, ncol=3))
#par(mar=c(4,4,4,1))
#phyp(prnd = prnd.fib_mous, tax = utaxo_mous, PCs = c(1,2), tax.level = "family", make.plot = T, add.labels=T, add.loads=F, map = F)
#title(main = "Fiber length")
#phyp(prnd = prnd.mass_mous, tax = utaxo_mous, PCs = c(1,2), tax.level = "family", make.plot = T, add.labels=T, add.loads=F, map = F)
#title(main = "Muscle mass")
#phyp(prnd = prnd.pcsa_mous, tax = utaxo_mous, PCs = c(1,2), tax.level = "family", make.plot = T, add.labels=T, add.loads=F, map = F)
#title(main = "PCSA")
#dev.off()

#pdf(file = paste("plots", "muscle_space_squir.pdf", sep = "/"), width = 15, height = 6)
#palette(rainbow(3))
#layout(matrix(1:3, ncol=3))
#par(mar=c(4,4,4,1))
#phyp(prnd = prnd.fib_squir, tax = utaxo_squir, PCs = c(1,2), tax.level = "family", make.plot = T, add.labels=T, add.loads=F, map = F)
#title(main = "Fiber length")
#phyp(prnd = prnd.mass_squir, tax = utaxo_squir, PCs = c(1,2), tax.level = "family", make.plot = T, add.labels=T, add.loads=F, map = F)
#title(main = "Muscle mass")
#phyp(prnd = prnd.pcsa_squir, tax = utaxo_squir, PCs = c(1,2), tax.level = "family", make.plot = T, add.labels=T, add.loads=F, map = F)
#title(main = "PCSA")
#dev.off()

#pdf(file = paste("plots", "muscle_space_cteno.pdf", sep = "/"), width = 15, height = 6)
#layout(matrix(1:3, ncol=3))
#par(mar=c(4,4,4,1))
#palette(rainbow(17))
#phyp(prnd = prnd.fib_cteno, tax = utaxo, PCs = c(1,2), tax.level = "family", make.plot = T, add.labels=T, add.loads=F, map = F)
#title(main = "Fiber length")
#phyp(prnd = prnd.mass_cteno, tax = utaxo_cteno, PCs = c(1,2), tax.level = "family", make.plot = T, add.labels=T, add.loads=F, map = F)
#title(main = "Muscle mass")
#phyp(prnd = prnd.pcsa_cteno, tax = utaxo_cteno, PCs = c(1,2), tax.level = "family", make.plot = T, add.labels=T, add.loads=F, map = F)
#title(main = "PCSA")
#dev.off()

##########################################################################################
#Remake PCAs with geomorph"s function which allows to project the phylogeny in the space.
##########################################################################################

#For PCSA

pca.geomorph <- gm.prcomp(A= prnd.pcsa$pruned.data, phy = prnd.pcsa$pruned.tree)

o <- match(rownames(prnd.pcsa$pruned.data),utaxo$lgen)

if (print.plots) {
pdf(paste("plots/", "phyloPCA_geomorph_PCSA_generic.pdf"), width = 5, height = 10)
layout(1:2)
par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F))
points(pca.geomorph$x[,1:2], bg = c(1:3)[as.factor(utaxo$clade[o])], pch=21, cex = 2)
#identify(pca.geomorph$x[,1:2], labels = utaxo$lgen[o])
u <- par("usr")
v <- c(grconvertX(u[1:2], "user", "ndc"), grconvertY(u[3:4], "user", "ndc"))
v <- c(v[1], v[1]+((v[2]-v[1])/3), v[3], v[3]+((v[4]-v[3])/3))
par(fig=v, new=TRUE, mar=c(0,0,0,0))
plot.loads(pca.geomorph$rotation[,1:2])

par(fig=c(0,1,0,0.5), new=TRUE, mar=c(4,4,0.5,0.5) )
plot(pca.geomorph, axis1 = 3, axis2= 4, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F))
points(pca.geomorph$x[,3:4], bg = c(1:3)[as.factor(utaxo$clade[o])], pch=21, cex = 2)
#identify(pca.geomorph$x[,3:4], labels = utaxo$lgen[o])
u <- par("usr")
v <- c(grconvertX(u[1:2], "user", "ndc"), grconvertY(u[3:4], "user", "ndc"))
v <- c(v[2]-((v[2]-v[1])/3), v[2], v[4]-((v[4]-v[3])/3), v[4])
par(fig=v, new=TRUE, mar=c(0,0,0,0))
plot.loads(pca.geomorph$rotation[,3:4])
dev.off()
}

#For fiber lengths

pca.geomorph <- gm.prcomp(A= prnd.lsr.fib$pruned.data, phy = prnd.lsr.fib$pruned.tree)

o <- match(rownames(prnd.lsr.fib$pruned.data),utaxo$lgen)

if (print.plots) {
pdf(paste("plots/", "phyloPCA_geomorph_fibers_generic.pdf"), width = 5, height = 10)
layout(1:2)
par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F))
points(pca.geomorph$x[,1:2], bg = c(1:3)[as.factor(utaxo$clade[o])], pch=21, cex = 2)
#identify(pca.geomorph$x[,1:2], labels = utaxo$lgen[o])
u <- par("usr")
v <- c(grconvertX(u[1:2], "user", "ndc"), grconvertY(u[3:4], "user", "ndc"))
v <- c(v[2]-((v[2]-v[1])/3), v[2], v[3], v[3]+((v[4]-v[3])/3))
par(fig=v, new=TRUE, mar=c(0,0,0,0))
plot.loads(pca.geomorph$rotation[,1:2])

par(fig=c(0,1,0,0.5), new=TRUE, mar=c(4,4,0.5,0.5) )
plot(pca.geomorph, axis1 = 3, axis2= 4, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F))
points(pca.geomorph$x[,3:4], bg = c(1:3)[as.factor(utaxo$clade[o])], pch=21, cex = 2)
#identify(pca.geomorph$x[,3:4], labels = utaxo$lgen[o])
u <- par("usr")
v <- c(grconvertX(u[1:2], "user", "ndc"), grconvertY(u[3:4], "user", "ndc"))
v <- c(v[2]-((v[2]-v[1])/3), v[2], v[3], v[3]+((v[4]-v[3])/3))
par(fig=v, new=TRUE, mar=c(0,0,0,0))
plot.loads(pca.geomorph$rotation[,3:4])
dev.off()
}

#For muscle mass

pca.geomorph <- gm.prcomp(A= prnd.mass$pruned.data, phy = prnd.mass$pruned.tree)

o <- match(rownames(prnd.mass$pruned.data),utaxo$lgen)

if (print.plots) {
pdf(paste("plots/", "phyloPCA_geomorph_musclemass_generic.pdf"), width = 5, height = 10)
layout(1:2)
par(mar = c(4,4,0.5,0.5))
plot(pca.geomorph, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F))
points(pca.geomorph$x[,1:2], bg = c(1:3)[as.factor(utaxo$clade[o])], pch=21, cex = 2)
#identify(pca.geomorph$x[,1:2], labels = utaxo$lgen[o])
u <- par("usr")
v <- c(grconvertX(u[1:2], "user", "ndc"), grconvertY(u[3:4], "user", "ndc"))
v <- c(v[1], v[1]+((v[2]-v[1])/3), v[4]-((v[4]-v[3])/3), v[4])
par(fig=v, new=TRUE, mar=c(0,0,0,0))
plot.loads(pca.geomorph$rotation[,1:2])

par(fig=c(0,1,0,0.5), new=TRUE, mar=c(4,4,0.5,0.5) )
plot(pca.geomorph, axis1 = 3, axis2= 4, phylo=T, phylo.par=list(tip.labels=T, node.labels=F, tip.labels=F,anc.states=F))
points(pca.geomorph$x[,3:4], bg = c(1:3)[as.factor(utaxo$clade[o])], pch=21, cex = 2)
#identify(pca.geomorph$x[,3:4], labels = utaxo$lgen[o])
u <- par("usr")
v <- c(grconvertX(u[1:2], "user", "ndc"), grconvertY(u[3:4], "user", "ndc"))
v <- c(v[1], v[1]+((v[2]-v[1])/3), v[3], v[3]+((v[4]-v[3])/3))
par(fig=v, new=TRUE, mar=c(0,0,0,0))
plot.loads(pca.geomorph$rotation[,3:4])
dev.off()
}



