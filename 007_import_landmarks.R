source("Rfunctions1.txt")
source("CUSTOM_R_FUNCTIONS.R")
library(geomorph)

# Opt to print plots or not:
print.plots <- F
if (print.plots & !dir.exists("plots")) {dir.create("plots")} #Create output folder if it does not exist

# Opt to symmetrize cranium landmarks/curves or not:
symmetriz <- T

##############################
# Import landmarks from files#
##############################

input_dir <- "landmarks/"
fil <- list.files(input_dir, full.names = T)
fi <- list.files(input_dir, full.names = F)
sp_short <- shorten.labels(fi, c(1,2), split = "_")
specimens <- read.csv("specimens_landmarked.csv")
specimens <- specimens[order(specimens$file_name),]

imp_lm <- lapply(fil, scan, what = "charater")
lslm <- lapply(X = imp_lm, FUN = matrix, ncol = 4, byrow = T)

nlm <- unlist(lapply(X=lslm, FUN= myfun <- function(x) {length(grep(x[,1], pattern =  "Landmark"))}))

for (i in 1:length(lslm)) {
  coo <- apply(lslm[[i]][,2:4], 2, as.numeric)
  rownames(coo) <- lslm[[i]][,1]
  lslm[[i]] <- coo}

#########################################
#Extract cranium landmarks into an array#
#########################################

Ac <- array(dim = c(20, 3, length(lslm)), 
           dimnames = list(c(1:20), c("x","y","z")))

for (i in 1:length(lslm)) {
  Ac[,,i] <- lslm[[i]][1:20,]}

mssglm <- which(apply(is.na(Ac), FUN=sum, MARGIN=3)>0)

#Estimate missing cranium landmarks using thin plate splines
Acc <- estimate.missing(A = Ac, method="TPS")

#plot(Ac[,1, mssglm[1]], Ac[,3, mssglm[1]], asp=1)
#points(Acc[,c(1,3), mssglm[1]], col="red", pch=20)

#plot(Ac[,1, mssglm[2]], Ac[,3, mssglm[2]], asp=1)
#points(Acc[,c(1,3), mssglm[2]], col="red", pch=20)

#plot(Ac[,1, mssglm[3]], Ac[,3, mssglm[3]], asp=1)
#points(Acc[,c(1,3), mssglm[3]], col="red", pch=20)

#Alignment of cranium shapes
pAcc <- pgpa(Acc)
sizAcc <- pAcc$cent.size

sp_centroid_size <- data.frame(specimens$genus, specimens$sp, pAcc$cent.size)
gen_csize <- tapply(X=sp_centroid_size$pAcc.cent.size, INDEX=sp_centroid_size$specimens.genus, FUN = mean)

if (symmetriz) {
#Symmetrize cranium LM
ro <- c(1:2, 4, 3, 5,  7, 6, 9, 8, 11, 10, 12, 14, 13, 15:17, 19, 18,20)
symAcc <- apply(pAcc$rotated, MARGIN = 3, FUN = maksym, reorder = ro)
symAcc_array <- array(NA, dim=c(dim(symAcc)[1]/3, 3, dim(symAcc)[2]))
for (i in 1:dim(symAcc_array)[3]) {symAcc_array[,,i] <- matrix(symAcc[,i], ncol=3)}
pAcc <- pgpa(symAcc_array)
}

##################################################
#Extract curve landmarks (all located on cranium)#
##################################################

lsrn <- lapply(X = lslm, rownames)

#Find specimens with missing curves
mssgcrv <- rep(F, length(lsrn))

for (i in 1:length(lsrn)) {
  o <- lsrn[[i]]
  allcrv <- length(grep(o, pattern =  "Curve_segment:14")) > 0
  if (!allcrv) {mssgcrv[i] <- T}
}

Acrv <- array(dim = c(125, 3, length(lslm)-length(which(mssgcrv))), 
           dimnames = list(c(1:125), c("x","y","z")))

tl <- lslm[which(!mssgcrv)]

for (i in 1:dim(Acrv)[3]) {
  o <- tl[[i]]
  Acrv[,,i] <- o[grep(rownames(o), pattern="Curve"),]
}

##################################################
#Combine all cranium data (landmarks and curves) #
##################################################

AC <- array(dim = c(145, 3, length(lslm)-length(which(mssgcrv))), 
           dimnames = list(c(1:145), c("x","y","z")))
           
#Remove specimens with missing curves from the cranium LM array
tA <- Acc[,,-which(mssgcrv)]

#Fill array with cranium LM and curves
for (i in 1:dim(AC)[3]) {
  AC[1:20,,i] <- tA[,,i]
  AC[21:145,,i] <- Acrv[,,i]
}
  
#Remove specimens with missing curves from metadata table
spec <- specimens[-which(mssgcrv),]

#Alignment of shapes
pAC <- pgpa(AC)
sizAC <- pAC$cent.size

if (symmetriz) {
# Symmetrize cranium LM and curves
ro <- c(1:2, 4, 3, 5,  7, 6, 9, 8, 11, 10, 12, 14, 13, 15:17, 19, 18,20, 21:45, 96:145, 46:95)
symAC <- apply(pAC$rotated, MARGIN = 3, FUN = maksym, reorder = ro)
symAC_array <- array(NA, dim=c(dim(symAC)[1]/3, 3, dim(symAC)[2]))
for (i in 1:dim(symAC_array)[3]) {symAC_array[,,i] <- matrix(symAC[,i], ncol=3)}
pAC <- pgpa(symAC_array)
}

#######################################################################
# Alignment of entire configurations and curves, but based only on LM#
######################################################################

pAcrv <- pgpa.curves(A=tA,Acurv=AC)
sizAcrv <- pAcrv$cent.size

if (symmetriz) {
# Symmetrize cranium LM and curves
ro <- c(1:2, 4, 3, 5,  7, 6, 9, 8, 11, 10, 12, 14, 13, 15:17, 19, 18,20, 21:45, 96:145, 46:95)
symAcrv <- apply(pAcrv$curves_rot, MARGIN = 3, FUN = maksym, reorder = ro)
symAcrv_array <- array(NA, dim=c(dim(symAcrv)[1]/3, 3, dim(symAcrv)[2]))
for (i in 1:dim(symAcrv_array)[3]) {symAcrv_array[,,i] <- matrix(symAcrv[,i], ncol=3)}
pAcrv <- pgpa.curves(symAcrv_array[1:20,,], Acurv = symAcrv_array)
}

######################################################
#Extract mandible landmarks into an array and average#
#the two hemimandibles shapes.
######################################################

#Empty list
lsmd <- list()

# Fill list with matrices containing HEMImandible landmarks

for (i in 1:length(lslm)) {
  m <- lslm[[i]]
  wlm <- grep(pattern = "Landmark", rownames(m))
  mdm <- m[wlm,][-c(1:20),]
  
  if (dim(mdm)[1] == 12) {lsmd[[i]] <- list(mdm)} #If there is only one hemimandible
  else if (dim(mdm)[1] == 24) {lsmd[[i]] <- list(left = mdm[c(1:10,21:22),], right = mdm[c(11:20,23:24),])} #If both hemimandibles are present.
  #Note: LMs 21 to 24 were added later, and therefore do not follow the number sequence.
  else {lsmd[[i]] <- NA} #One specimen had both hemimandibles broken
}

onemd <- which(unlist(lapply(lapply(lsmd,names),length))==0) #Which specimens have only one hemimandible or none
remaining_md <- c("right", "right", "left", "right", "right", NA, "right", "right", "right", "left", "left") #Must be checked manually on 3D models

for (i in 1:length(onemd)) {
  md <- remaining_md[i]
  names(lsmd[[onemd[i]]]) <- md
} # Add the side identification to list elements which had only one hemimandible

# Align mandibles on a reference plane, mirror all right mandible into left mandibles, and average left and right shapes for each individual, in cases where both hemimandibles are present 

lsavmd <- list()

for (i in 1:length(lsmd)) {
  d <- lsmd[[i]]
  
  if (length(d) == 2) {
    l <- d$left
    r <- d$right
    la <- bookstein3d(l, 2, 7, 9) #Alignement based on LMs 2, 7 and 9
    ra <- bookstein3d(r, 2, 7, 9)
    
    dmir <- rep(NA,dim(l)[2])    
    for (j in 1:dim(l)[2]) {
      mra <- ra   
      mra[,j] <- ra[,j]*(-1)
      dmir[j] <- sqrt(sum((la-mra)^2, na.rm=T))
      }
    amir <- which.min(dmir)
    
    mirrd <- ra
    mirrd[,amir] <- ra[,amir]*(-1)
    
    avmd <- (mirrd+la)/2
    
    #replace NAs in one hemimandible by coordinates from the other (when present).
    sumNA <- is.na(la) + is.na(ra)
    avmd[sumNA == 1] <- na.omit(c(ra[sumNA == 1], la[sumNA == 1]))
    
    lsavmd[[i]] <- avmd
  }
  
  if (length(d) == 1) {
    if (is.na(d)) {lsavmd[[i]] <- NA}
    
    else if (names(d) == "right") {
      r <- d$right
      ra <- bookstein3d(r, 2, 7, 9)
      
      dmir <- rep(NA,dim(r)[2])    
      for (j in 1:dim(r)[2]) {
        mra <- ra   
        mra[,j] <- ra[,j]*(-1)
        dmir[j] <- sqrt(sum((ra-mra)^2, na.rm=T))
        }
      amir <- which.min(dmir)
      mirrd <- ra
      mirrd[,amir] <- ra[,amir]*(-1) 
      lsavmd[[i]] <- mra
    }
    
    else if (names(d) == "left") {
      l <- d$left
      la <- bookstein3d(l, 2, 7, 9)
      lsavmd[[i]] <- la
    }
  }
}

#Empty array
Am <- array(dim = c(12, 3, length(lslm)), 
           dimnames = list(c(1:12), c("x","y","z")))

for (i in 1:dim(Am)[3]) {Am[,,i] <- lsavmd[[i]]}

mdNA <- which(apply(Am, FUN = anyNA, MARGIN=3))

Amc <- estimate.missing(Am[,,-83]) #No mandibles digitized for id 83 (idurus zenkeri)

#layout(matrix(1:40, nrow=4))
#for (i in 1:40) {plot(Amc[,1,i], Amc[,2,i],asp=1);points(Am[,1,i], Am[,2,i],col=2, pch=20)}

pAmc <- pgpa(Amc)
sizAmc <- pAmc$cent.size

spem <- specimens[-83,]

mandi_csize <- data.frame(spem$genus, spem$sp, pAmc$cent.size)
gen_csize_mandi <- tapply(X=mandi_csize$pAmc.cent.size, INDEX=mandi_csize$spem.genus, FUN = mean)
 
#####################################################
#Convert rotated arrays into matrices of coordinates# 
#####################################################

# Only cranium landmarks
Mc <- matrix(NA, 
            ncol = dim(pAcc$rotated)[1] * dim(pAcc$rotated)[2], 
            nrow = dim(pAcc$rotated)[3])

for (i in 1:dim(pAcc$rotated)[3]) {
  o <- pAcc$rotated[,,i]
  Mc[i,] <- c(t(o))
}

rownames(Mc) <- paste(specimens$genus, specimens$sp, specimens$number, sep="_")

# Cranium landmarks and curves            
MC <- matrix(NA, 
            ncol = dim(pAC$rotated)[1] * dim(pAC$rotated)[2], 
            nrow = dim(pAC$rotated)[3])

for (i in 1:dim(pAC$rotated)[3]) {
  o <- pAC$rotated[,,i]
  MC[i,] <- c(t(o))
}

Mcrv <- matrix(NA, 
            ncol = dim(pAcrv$curves_rot)[1] * dim(pAcrv$curves_rot)[2], 
            nrow = dim(pAcrv$curves_rot)[3])

for (i in 1:dim(pAcrv$curves_rot)[3]) {
  o <- pAcrv$curves_rot[,,i]
  Mcrv[i,] <- c(t(o))
}


Mt <- matrix(NA, 
            ncol = dim(pAcrv$curves_rot)[1] * dim(pAcrv$curves_rot)[2], 
            nrow = dim(pAcrv$curves_rot)[3])

for (i in 1:dim(pAcrv$curves_rot)[3]) {
  o <- pAcrv$curves_rot[,,i]
  Mt[i,] <- c(o)
}


rownames(MC) <- rownames(Mcrv) <- paste(spec$genus, spec$sp, spec$number, sep="_")

# Mandible landmarks
Mmd <- matrix(NA, 
            ncol = dim(pAmc$rotated)[1] * dim(pAmc$rotated)[2], 
            nrow = dim(pAmc$rotated)[3])

for (i in 1:dim(pAmc$rotated)[3]) {
  o <- pAmc$rotated[,,i]
  Mmd[i,] <- c(t(o))
}

rownames(Mmd) <- paste(spem$genus, spem$sp, spem$number, sep="_")

#Combine all landmarks and curves
M <- cbind(Mcrv, Mmd[match(rownames(Mcrv), rownames(Mmd)),])
M2 <- cbind(MC, Mmd[match(rownames(MC), rownames(Mmd)),])

#################################
# Raw PCAs to check for outliers#
#################################

pcac <- prcomp(Mc)
pcaC <- prcomp(MC)
pcacrv <- prcomp(Mcrv)
pcamd <- prcomp(Mmd)
pcaM <- prcomp(M)
pcaM2 <- prcomp(M2)

if (print.plots) {
pdf(file = paste("plots/", "rawPCAs_for_outlier.pdf", sep=""), height = 12, width = 20)
layout(matrix(1:6, nrow = 2))
par(mar = c(4,4,3,1))
plot(pcac$x[,1:2], asp = 1, pch = 21, bg = "grey", main = "Cranium LMs only")
text(pcac$x[,1:2], labels = specimens$genus, cex = 0.7, pos = 2)

plot(pcaC$x[,1:2], asp = 1, pch = 21, bg = "grey", main = "Cranium LMs and curves, all used to align")
text(pcaC$x[,1:2], labels = spec$genus, cex = 0.7, pos = 2)

plot(pcacrv$x[,1:2], asp = 1, pch = 21, bg = "grey", main = "Cranium LMs and curves, only LMs used to align")
text(pcacrv$x[,1:2], labels = spec$genus, cex = 0.7, pos = 2)

plot(pcamd$x[,1:2], asp = 1, pch = 21, bg = "grey", main = "Mandible LMs only")
text(pcamd$x[,1:2], labels = spem$genus, cex = 0.7, pos = 2)

plot(pcaM$x[,1:2], asp = 1, pch = 21, bg = "grey", main = "All LMs and curves, only cranium LMs used to align")
text(pcaM$x[,1:2], labels = spec$genus, cex = 0.7, pos = 2)

plot(pcaM2$x[,1:2], asp = 1, pch = 21, bg = "grey", main = "All LMs and curves, cranium LMs and curves used to align")
text(pcaM2$x[,1:2], labels = spec$genus, cex = 0.7, pos = 2)
dev.off()
}

# Split global dataset into "modules":

cran <- 1:60
vault <- 61:135
darcd <- c(136:210,286:360)
varcd <- c(211:285,361:435)
mand <- 436:471

lsMmodu <- list(M[,cran], M[,vault], M[,darcd], M[,varcd], M[,mand])

# Correlate centroid size with PCs to check for allometric effects

size_correl <- cor(cbind(sizAcrv, pcaM$x))[,1]

if (print.plots) {
pdf(file = paste("plots/","Correlation_size_PCs.pdf", sep=""), width = 12, height = 7)
layout(matrix(1:6, ncol=3, byrow=T))
for (i in 1:6) {
plot(sizAcrv, pcaM$x[,i], xlab = "Centroid Size", ylab = paste("PC",i), pch=20)
legend("topright", legend = round(size_correl[i+1],2), bty="n", cex=2)}
dev.off()
}

