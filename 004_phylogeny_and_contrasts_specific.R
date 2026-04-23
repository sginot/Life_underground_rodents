#Required package
library(ape)
library(phytools)
library(corrplot)
source("CUSTOM_R_FUNCTIONS.R")

# Opt to print plots or not:
print.plots <- F
if (print.plots & !dir.exists("plots")) {dir.create("plots")} #Create output folder if it does not exist

#Load phylogenetic trees
sp_tre <- read.tree("Disparat_specific_tree")

#Load data tables produced by previous script
df_sum_mass <- read.csv("tables/adductor_muscle_mass.csv", h = T, sep = ",", dec = ".", row.names = 1)
df_sum_pcsa <- read.csv("tables/adductor_pcsa.csv", h = T, sep = ",", dec = ".", row.names = 1)
mscl_lsr_mass <- read.csv("tables/adductor_LSR_mass.csv", h = T, sep = ",", dec = ".", row.names = 1)
fib_lgt <- read.csv("tables/adductor_fib_length.csv", h = T, sep = ",", dec = ".", row.names = 1)
lsr_fib_lgt <- read.csv("tables/adductor_LSR_fib_length.csv", h = T, sep = ",", dec = ".", row.names = 1)
mscl_lsr_pcsa <- read.csv("tables/adductor_LSR_pcsa.csv", h = T, sep = ",", dec = ".", row.names = 1)
dat_dig <- read.csv("tables/digastric_data.csv", h = T, sep = ",", dec = ".", row.names = 1)

# Compute geometric mean of muscle masses
geom_mean <- function(x) {(prod(x))^(1/length(x))}
gmean_mass <- apply(df_sum_mass, 1, geom_mean)

#Compute the total adductor muscle mass and pcsa per species
total_mass <- apply(df_sum_mass, 1, sum)
total_pcsa <- apply(df_sum_pcsa, 1, sum)

#Compute muscle proportions as percentages
muscle_prop <- df_sum_mass/total_mass*100

#Compute species average data (only changes species for which several specimens are present in the original data).
#1. Reduce row names to species binomial names
shortnam <- shorten.labels(rownames(mscl_lsr_pcsa), positions = c(1,2), split = " ")
#2. Apply mean function across data frame based on species name as a factor.
mscl_lsr_pcsa <- apply(X = mscl_lsr_pcsa, 2 , FUN = tapply, INDEX = shortnam, mean)

shortnam <- shorten.labels(rownames(mscl_lsr_mass), positions = c(1,2), split = " ")
mscl_lsr_mass <- apply(X = mscl_lsr_mass, 2 , FUN = tapply, INDEX = shortnam, mean)

shortnam <- shorten.labels(rownames(lsr_fib_lgt), positions = c(1,2), split = " ")
lsr_fib_lgt <- apply(X = lsr_fib_lgt, 2 , FUN = tapply, INDEX = shortnam, mean)

shortnam <- shorten.labels(rownames(fib_lgt), positions = c(1,2), split = " ")
fib_lgt <- apply(X = fib_lgt, 2 , FUN = tapply, INDEX = shortnam, mean)

shortnam <- shorten.labels(rownames(dat_dig), positions = c(1,2), split = " ")
dat_dig <- na.omit(apply(X = dat_dig, 2 , FUN = tapply, INDEX = shortnam, mean))

#Prune the species level pcsa LSR data
prnd.pcsa <- data.pruning(sp_tre, mscl_lsr_pcsa, split.t = "_", split.d = " ", positions = c(1,2))
#Just a check; not necessary
#par(mar = c(0,0,0,0))
#plot(prnd$pruned.tree, type = "fan")
#Prune muscle mass data
prnd.mass <- data.pruning(sp_tre, mscl_lsr_mass, split.t = "_", split.d = " ", positions = c(1,2))
prnd.dig <- data.pruning(sp_tre, dat_dig, split.t = "_", split.d = " ", positions = c(1,2))
#Prune fiber length data
prnd.lsr.fib <- data.pruning(sp_tre, lsr_fib_lgt, split.t = "_", split.d = " ", positions = c(1,2))
prnd.fib <- data.pruning(sp_tre, fib_lgt, split.t = "_", split.d = " ", positions = c(1,2))

if (print.plots) {
plot.mapd.traits(prnd = prnd.fib, 
                 file = "sp_maps_fiber_lgt.pdf", 
                 titles = colnames(prnd.fib$pruned.data))
plot.mapd.traits(prnd = prnd.lsr.fib, 
                 file = "sp_maps_LSR_fiber_lgt.pdf", 
                 titles = colnames(prnd.lsr.fib$pruned.data))
plot.mapd.traits(prnd = prnd.pcsa, 
                 file = "sp_maps_LSR_PCSA.pdf", 
                 titles = colnames(prnd.pcsa$pruned.data))
plot.mapd.traits(prnd = prnd.mass, 
                 file = "sp_maps_LSR_mass.pdf",
                 titles = colnames(prnd.mass$pruned.data))
plot.mapd.traits(prnd = prnd.dig, 
                 file = "sp_maps_dig.pdf", 
                 titles = colnames(prnd.dig$pruned.data))
}

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
dig_pics <- pics.df(prnd = prnd.dig)

if (print.plots) {
pdf(file = paste("plots","corrplots_muscle_variables.pdf", sep="/"), width = 20, height = 10)

layout(matrix(1:8, ncol = 4, byrow=T))

corrplot(cor(prnd.lsr.fib$pruned.data), 
         method = "ellipse", 
         type = "lower", 
         title = "LSR fiber length",
         mar = c(0,0,4,0))
corrplot(cor(prnd.fib$pruned.data), 
         method = "ellipse", 
         type = "lower", 
         title = "Raw fiber length",
         mar = c(0,0,4,0))
corrplot(cor(prnd.mass$pruned.data), 
         method = "ellipse", 
         type = "lower", 
         title = "LSR muscle mass",
         mar = c(0,0,4,0))
corrplot(cor(prnd.pcsa$pruned.data), 
         method = "ellipse", 
         type = "lower", 
         title = "LSR PCSA",
         mar = c(0,0,4,0))
         
corrplot(cor(lsr_fib_pics), 
         method = "ellipse", 
         type = "lower", 
         title = "PICs LSR fiber length",
         mar = c(0,0,4,0))
corrplot(cor(fib_pics), 
         method = "ellipse", 
         type = "lower", 
         title = "PICs Raw fiber length",
         mar = c(0,0,4,0))
corrplot(cor(mass_pics), 
         method = "ellipse", 
         type = "lower", 
         title = "PICs LSR muscle mass",
         mar = c(0,0,4,0))
corrplot(cor(pcsa_pics), 
         method = "ellipse", 
         type = "lower", 
         title = "PICs LSR PCSA",
         mar = c(0,0,4,0))
dev.off()
}

#Phylogenetic PCA function
phyp <- function(prnd, 
                 tax = utaxo, 
                 PCs = c(1,2), 
                 tax.level = "clade", 
                 make.plot = T) {
  
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
  s <- unlist(lapply(snam, FUN = f <- function(x) {paste(x[[1]], x[[2]])}))
  ta <- tax[match(g, tax$lgen),]
  
  #Map taxonomic group as a discrete character on the tree
  taxfac <- as.factor(ta[,tax.level])
  taxfac <- setNames(taxfac, tre$tip.label)
  tresimmap <- make.simmap(tre, x=taxfac)
  
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
    phylomorphospace(tresimmap, pcs[,PCs], ftype = "off", node.by.map = T, 
                   bty="n",
                   colors=cols,
                   node.size=c(0,1.2),
                   xlab = vv[PCs[1]],
                   ylab = vv[PCs[2]])
    text(pcs[,PCs], s, cex = 0.5, pos = 2)
    par(fig = c(0.7,1,0.7,1), new = T, mar=c(0,0,4,2))
    plot.loads(lo)
    }
  
return(list(PCA = phylpca, Tree.taxo.mapped = tresimmap))}

# Next section creates morphospace, but takes a few minutes, optional.

#phypout <- phyp(prnd = prnd.lsr.fib, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T)

#pdf(file = paste("plots", "LSR_fib_sp_phylospace12.pdf", sep = "/"))
#phyp(prnd = prnd.lsr.fib, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "LSR_fib_sp_phylospace34.pdf", sep = "/"))
#phyp(prnd = prnd.lsr.fib, tax = utaxo, PCs = c(3,4), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "fib_sp_phylospace12.pdf", sep = "/"))
#phyp(prnd = prnd.fib, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "fib_sp_phylospace34.pdf", sep = "/"))
#phyp(prnd = prnd.fib, tax = utaxo, PCs = c(3,4), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "pcsa_sp_phylospace12.pdf", sep = "/"))
#phyp(prnd = prnd.pcsa, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "pcsa_sp_phylospace34.pdf", sep = "/"))
#phyp(prnd = prnd.pcsa, tax = utaxo, PCs = c(3,4), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "mass_sp_phylospace12.pdf", sep = "/"))
#phyp(prnd = prnd.mass, tax = utaxo, PCs = c(1,2), tax.level = "clade", make.plot = T)
#dev.off()
#pdf(file = paste("plots", "mass_sp_phylospace34.pdf", sep = "/"))
#phyp(prnd = prnd.mass, tax = utaxo, PCs = c(3,4), tax.level = "clade", make.plot = T)
#dev.off()

write.csv(x = lsr_fib_pics, file = paste("tables", "pics_sp_LSR_fiber_length.csv", sep = "/"))
write.csv(x = fib_pics, file = paste("tables", "pics_sp_fiber_length.csv", sep = "/"))
write.csv(x = pcsa_pics, file = paste("tables", "pics_sp_pcsa.csv", sep = "/"))
write.csv(x = mass_pics, file = paste("tables", "pics_sp_muscle_mass.csv", sep = "/"))
write.csv(x = dig_pics, file = paste("tables", "pics_sp_dig.csv", sep = "/"))

