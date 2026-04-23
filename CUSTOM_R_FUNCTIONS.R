#----------------------------------------------
# Remove the n last characters from a vector of character strings
#-----------------------------------------------

rm.last.chr <- function(x, n) {
  nc <- nchar(x)
  brk <- nc - n
  substr(x, 1, brk)
}

#----------------------------------------------
# Keep the n last characters from a vector of character strings
#-----------------------------------------------

keep.last.chr <- function(x, n) {
  nc <- nchar(x)
  brk <- nc - n
  substr(x, brk + 1, nc)
}


#-------------------------------------------------------
# Add CI polygon to a biplot
#-------------------------------------------------------

plot.lm.CI <- function(x, y, 
                       border = NULL,
                       color = "gray", 
                       alpha = NA, 
                       level = 0.95,
                       line = T) {
  
  mod <- lm(y ~ x)
  
  minx <- min(x, na.rm = T)
  maxx <- max(x, na.rm = T)
  rg <- abs(minx - maxx)
  
  newx <- seq(minx, 
              maxx, 
              by = rg/1000) 
  
  prd <- predict(object = mod,
                 newdata = data.frame(x = newx),
                 interval = "confidence",
                 level = level)
  
  if (line) {
    clip(minx, maxx, min(y)/2, max(y)*2)
    abline(mod, col = color)}
  
  if (!is.na(alpha)) {
    require(scales)
    polygon(x = c(newx, rev(newx)),
          y = c(prd[,2], rev(prd[,3])),
          col = alpha(color, alpha = alpha), 
          border = border)
  } else {
    polygon(x = c(newx, rev(newx)),
            y = c(prd[,2], rev(prd[,3])),
            col = color)
  }
  
} 

#--------------------------------------------------------
# Check if a number is even (return TRUE or FALSE)
#--------------------------------------------------------

is.even <- function(x, verbose = T) {
  if (verbose) {
    print(x/2 == round(x/2))
    x/2 == round(x/2)
  } else {
    x/2 == round(x/2)
  }
}

#--------------------------------------------------------
# Match two vectors of different length based on the names 
# of their elements, output data frame with common elements
# for each vector
#--------------------------------------------------------

vec.mtch <- function(vec1, vec2) {
  if (is.null(names(vec1)) | is.null(names(vec2))) {
    stop("Both vectors must have names")
  }
  
  s1 <- vec1[names(vec1) %in% names(vec2)]
  s1 <- s1[order(names(s1))]
  s2 <- vec2[names(vec2) %in% names(vec1)]
  s2 <- s2[order(names(s2))]
  
  d <- data.frame(s1, s2)
  rownames(d) <- names(s1)
  colnames(d)[1] <- deparse(substitute(vec1))
  colnames(d)[2] <- deparse(substitute(vec2))
  d
}


#--------------------------------------------------------
# Match two matrices of different row numbers based on row
# names, output data frame with common rows for each matrix
#--------------------------------------------------------

mat.mtch <- function(mat1, mat2) {
  if (is.null(rownames(mat1)) | is.null(rownames(mat2))) {
    stop("Both matrices must have row names")
  }
  
  s1 <- as.matrix(mat1[rownames(mat1) %in% rownames(mat2),])
  s1 <- as.matrix(s1[order(rownames(s1)),])
  s2 <- as.matrix(mat2[rownames(mat2) %in% rownames(mat1),])
  s2 <- as.matrix(s2[order(rownames(s2)),])
  
  d <- data.frame(s1, s2)
  rownames(d) <- rownames(s1)
  d
}


#--------------------------------------------------------
# Function to shorten labels (e.g. tree tips or rownames), 
# based on strplit of the vectors into the different components. 
# Positions argument defines which element(s) of the 
# strsplit output to keep (e.g. only generic / specific name)
#--------------------------------------------------------

shorten.labels <- function(x, positions, split) {
  splt <- strsplit(x = x, split = split)
  newlab <- list()  
  for (i in 1:length(splt)) {
    newlab[[i]] <- paste(splt[[i]][positions], collapse = " ")
  }
  return(unlist(newlab))
}

#--------------------------------------------------------
# Function to average matrix based on factor, 
# i.e. apply / tapply combo
#--------------------------------------------------------

avg.matrix.fac <- function(mat, MARGIN = 2, INDEX) {
  mav <- apply(mat, 2, tapply, INDEX = INDEX, mean, na.rm=T)
return(mav)}

#--------------------------------------------------------
#Function to concomitantly prune data frames and phylogenies based on matching species (row names and tip labels, respectively). Output should be a tree and a data frame with perfectly matching tips / rows.
#--------------------------------------------------------

data.pruning <- function(tree, dat, split.t, split.d, positions) {
  treetip <- tree$tip.label
  shorttips <- shorten.labels(x = treetip, positions, split.t)
  shortnam <- shorten.labels(rownames(dat), positions, split.d)
    
  tip_match <- which(shorttips %in% shortnam)
  sp_match <- shorttips[tip_match]
  
  if (any(grepl(pattern="sp\\.", x = shortnam))) {
  msg_sp <- shortnam[grep(pattern="sp\\.", x = shortnam)]
  gen_msgsp <- gsub(pattern=" sp.", replacement="", x= msg_sp)
  possibles <- shorttips[grep(pattern = paste(gen_msgsp, collapse = "|"), x=shorttips)]
  real_poss <- possibles[-which(possibles %in% sp_match)]
  gg <- unique(unlist(lapply(strsplit(real_poss, split = " "), FUN = fun <- function(x){x[[1]]})))
    for (i in gg) {
    to_find <- real_poss[grep(i, x = real_poss)[1]]
    shorttips[grep(pattern = to_find, x = shorttips)] <- paste(i, "sp.", sep = " ")
    }
  tip_match <- which(shorttips %in% shortnam)
  sp_match <- shorttips[tip_match]
  }
  
  pdat <- dat[match(sp_match, shortnam),]
  rownames(pdat) <- sp_match
  ptree <- keep.tip(phy = tree, tip = tip_match)
  ptree$tip.label <- sp_match
  
  return(list(pruned.tree = ptree, pruned.data = pdat))
} 


array.pruning <- function(tree, dat, A, split.t, split.d, positions) {
  treetip <- tree$tip.label
  shorttips <- shorten.labels(x = treetip, positions, split.t)
  shortnam <- shorten.labels(rownames(dat), positions, split.d)
  
  tip_match <- which(shorttips %in% shortnam)
  sp_match <- shorttips[tip_match]
  
  pdat <- dat[match(sp_match, shortnam),]
  pA <- A[,,match(sp_match, shortnam)]
  rownames(pdat) <- sp_match
  ptree <- keep.tip(phy = tree, tip = tip_match)
  ptree$tip.label <- sp_match
  
  return(list(tree = ptree, data = pdat, array = pA))
} 


#---------------------------------------------
#Function to compute and plot character maps
#---------------------------------------------

plot.mapd.traits <- function(prnd, file, output = "plots", titles = NULL, width = 15) {
  tre <- prnd$pruned.tree
  dat <- prnd$pruned.data
  l_maps <- list()
  #Compute character maps along tree
  for (i in 1:dim(dat)[2]) {l_maps[[i]] <- contMap(tre, dat[,i], plot = F)}
 
  if (dim(dat)[2] %% 2 == 0) {mat <- matrix(1:dim(dat)[2], ncol = 2)}
  else {mat <- matrix(c(1:dim(dat)[2], 0), ncol = 2)}
  #Plot all contMap trees
  pdf(file = paste(output, file, sep = "/"), width = width, height = 4 * length(l_maps))
  layout(mat)
  par(mar = c(0,1,0,1))
  for (i in 1:dim(dat)[2]) {
    plot(l_maps[[i]], type = "fan", cex = 0.5)
    mtext(side = 3, line = -1, text = titles[i])
  }
  
  dev.off()
}

#----------------------------------------
# Function to plot pca loadings (eg to add to a pca biplot)
#----------------------------------------

plot.loads <- function(x, cex = 1) {

plot(0, 0, pch = 20, cex = 3, xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
abline(v=0,h=0, col="grey")

for (i in 1:dim(x)[1]) {
lines(x = c(0, x[i,1]), y = c(0, x[i,2]))
text(x = x[i,1], y = x[i,2], labels = rownames(x)[i], cex = cex)
}
}


#----------------------------------------
# Function to define colors for edges of a phylogeny
#----------------------------------------

find.edgecols <- function(phy, factor, color = rainbow(length(unique(factor)))) {

tips <- phy$tip.label
edges <- phy$edge

edge.fac <- rep(NA, dim(edges)[1])
cols <- rep("black", dim(edges)[1])

lvl <- levels(factor)

for (i in 1:length(lvl)){

index <- which(factor == lvl[i])
wedge <- which.edge(phy = phy, group = index)

edge.fac[wedge] <- lvl[i]
cols[wedge] <- color[i]
}

return(list(edge.fac = edge.fac, cols = cols))
}

#-------------------------------------------------------------------
# Function to compute congruence matrix for (3D) data matrix and PICs
#-------------------------------------------------------------------

congru <- function(M, phy, plot = T, split.t = "_", split.d= "_", positions = 1) {

hack <- array(NA, dim = c(dim(M)[2]/3, 3, dim(M)[1]))

for (i in 1:dim(M)[1]) {hack[,,i] <- matrix(M[i,], ncol=3, byrow=T)}

congru <- dotcorr(hack)

m <- abs(congru)
m[!lower.tri(m, diag = T)] <- NA

prnd <- data.pruning(tree = phy, dat = M, split.t, split.d, positions)

prnd.M <- prnd$pruned.data

PICs <- pics.df(prnd)

hackPIC <- array(NA, dim = c(dim(PICs)[2]/3, 3, dim(PICs)[1]))

for (i in 1:dim(PICs)[1]) {hackPIC[,,i] <- matrix(PICs[i,], ncol=3, byrow=T)}

congruPIC <- dotcorr(hackPIC)

mpic <- abs(congruPIC)
mpic[lower.tri(mpic, diag = T)] <- NA

if (plot) {
  par(mar = c(5,5,0.5,5))
  image(z = m, col=hcl.colors(n=12, "viridis"), xaxt="n", yaxt="n", ylab = "", xlab = "")
  image(z = mpic, col=hcl.colors(n=12, "viridis"), xaxt="n", yaxt="n", ylab = "", xlab = "", add =T)
}

return(list(m, mpic))}

#-------------------------------------------------------------------
# Function to compute 2B PLS for (3D) data matrix and PICs
#-------------------------------------------------------------------

plspic <- function(M, phy, parti, split.t = "_", split.d= "_", positions = 1, PIC = T) {

hack <- array(NA, dim = c(dim(M)[2]/3, 3, dim(M)[1]))

for (i in 1:dim(M)[1]) {hack[,,i] <- matrix(M[i,], ncol=3, byrow=T)}

integ <- integration.test(A = hack, partition.gp = parti)

if (PIC) { # Account for phylogeny with PIC approach
prnd <- data.pruning(tree = phy, dat = M, split.t, split.d, positions)

prnd.M <- prnd$pruned.data

PICs <- pics.df(prnd)

hackPIC <- array(NA, dim = c(dim(PICs)[2]/3, 3, dim(PICs)[1]))

for (i in 1:dim(PICs)[1]) {hackPIC[,,i] <- matrix(PICs[i,], ncol=3, byrow=T)}

integphylo <- integration.test(A = hackPIC, partition.gp = parti)

return(list(integ = integ, integphylo = integphylo))}

if (!PIC) { # Account for phylogeny with geomorph approach

prnd <- data.pruning(tree = phy, dat = M, split.t, split.d, positions)

prnd.M <- prnd$pruned.data

hackprnd <- array(NA, dim = c(dim(prnd.M)[2]/3, 3, dim(prnd.M)[1]))

dimnames(hackprnd) <- list(NULL, NULL, rownames(prnd.M))

for (i in 1:dim(prnd.M)[1]) {hackprnd[,,i] <- matrix(prnd.M[i,], ncol=3, byrow=T)}

integphylo <- phylo.integration(A = hackprnd, partition.gp = parti, phy = prnd$pruned.tree)

return(list(integ = integ, integphylo = integphylo))}
}

#-------------------------------------------------------------------
# Function to compute and plot 2B PLS results
#-------------------------------------------------------------------

plotplspic <- function(
  M, phy, fact = NULL,
  first.block, second.block, 
  first.GMM = F, second.GMM = T, 
  split.t = "_", split.d= "_", positions = 1, 
  print.progress = F, 
  PLS = 1, 
  PIC = T, 
  main = NULL, 
  col.line = "red", col.pch = "black",
  which.curves1 = NULL, which.curves2 = NULL, which.LM = NULL,
  view1 = c(1,2), view2 = c(2,3)) {


A1 <- M[,first.block]
A2 <- M[,second.block]

if (first.GMM) {
  hack <- array(NA, dim = c(dim(A1)[2]/3, 3, dim(A1)[1]))
  for (i in 1:dim(A1)[1]) {hack[,,i] <- matrix(A1[i,], ncol=3, byrow=T)}
  A1 <- hack}

if (second.GMM) {
  hack <- array(NA, dim = c(dim(A2)[2]/3, 3, dim(A2)[1]))
  for (i in 1:dim(A2)[1]) {hack[,,i] <- matrix(A2[i,], ncol=3, byrow=T)}
  A2 <- hack}

pls <- two.b.pls(A1 = A1, A2 = A2, print.progress = print.progress)

preds <- shape.predictor(A2, A1, Intercept=F, method = "PLS", pred1 = max(pls$XScores[,PLS]), pred2 = min(pls$XScores[,PLS]))

if (PIC) {
  prnd <- data.pruning(tree = phy, dat = M, split.t, split.d, positions)
  prnd.M <- prnd$pruned.data
  PICs <- pics.df(prnd)
  A1PIC <- PICs[,first.block]
  A2PIC <- PICs[,second.block]

  plspic <- two.b.pls(A1 = A1PIC, A2 = A2PIC, print.progress = print.progress)}

if (!PIC) {
  prnd <- data.pruning(tree = phy, dat = M, split.t, split.d, positions)
  prnd.M <- prnd$pruned.data
  prnd.names <- prnd$pruned.tree$tip.label
  prnd.A1 <- prnd.M[,first.block]
  prnd.A2 <- prnd.M[,second.block]

  if (first.GMM) {
    hack <- array(NA, dim = c(dim(prnd.A1)[2]/3, 3, dim(prnd.A1)[1]))
    for (i in 1:dim(prnd.A1)[1]) {hack[,,i] <- matrix(prnd.A1[i,], ncol=3, byrow=T)}
    prnd.A1 <- hack
    dimnames(prnd.A1) <- list(NULL, NULL, prnd.names)}

  if (second.GMM) {
    hack <- array(NA, dim = c(dim(prnd.A2)[2]/3, 3, dim(prnd.A2)[1]))
    for (i in 1:dim(prnd.A2)[1]) {hack[,,i] <- matrix(prnd.A2[i,], ncol=3, byrow=T)}
    prnd.A2 <- hack
    dimnames(prnd.A2) <- list(NULL, NULL, prnd.names)}

  plspic <- phylo.integration(A = prnd.A1, A2 = prnd.A2, print.progress = print.progress, phy = prnd$pruned.tree)

  predsphylo <- shape.predictor(prnd.A2, prnd.A1, Intercept=F, method = "PLS", pred1 = max(plspic$XScores[,PLS]), pred2 = min(plspic$XScores[,PLS]))

}

x <- pls$XScores[,PLS]
y <- pls$YScores[,PLS]

xpic <- plspic$XScores[,PLS]
ypic <- plspic$YScores[,PLS]

pc <- prcomp(cbind(x, y))$x[, 1]
px <- predict(lm(x ~ pc))
py <- predict(lm(y ~ pc))

pcpic <- prcomp(cbind(xpic, ypic))$x[, 1]
pxpic <- predict(lm(xpic ~ pcpic))
pypic <- predict(lm(ypic ~ pcpic))

varp <- round(pls$svd$d/sum(pls$svd$d)*100, 2)
varpic <- round(plspic$svd$d/sum(plspic$svd$d)*100, 2)

if (!is.null(fact)) {col.pch <- c(1:nlevels(fact))[fact]}

layout(matrix(c(1,1,2,2,3,3,4,4,5,6,7,8), ncol=3))

par(mar = c(4,4,4,0))
plot(x, y, pch = 20, xlab = paste("X PLS",PLS, varp[PLS], "%"), ylab = paste("Y PLS",PLS), col = col.pch)
abline(lm(py ~ px), col = col.line)
minmax <- c(which.min(y), which.max(y))
text(x[minmax], y[minmax], labels = names(x)[minmax], pos = c(4,2))
legend("bottomright", legend = c(paste("r-PLS=", round(pls$r.pls,2)), paste("P=", pls$P.value)), cex=0.7)

if (!is.null(fact)) {legend("topleft", legend = levels(fact), pch = 20, col = 1:nlevels(fact), bty = "n", cex = 0.7)}

title(main=main, outer = T, line = -2)

col.pch <- "black"
if (!PIC & !is.null(fact)) {
factpic <- fact[match(names(xpic), names(fact))]
col.pch <- c(1:nlevels(factpic))[factpic]}

par(mar = c(4,4,4,0), cex.axis=0.5)
plot(xpic, ypic, pch = 20, xlab = paste("PICs X PLS",PLS, varpic[PLS], "%"), ylab = paste("PICs Y PLS",PLS), col = col.pch)
abline(lm(pypic ~ pxpic), col = col.line)
minmax <- c(which.min(ypic), which.max(ypic))
text(xpic[minmax], ypic[minmax], labels = names(xpic)[minmax], pos = c(4,2))
legend("bottomright", legend = c(paste("r-PLS=", round(plspic$r.pls,2)), paste("P=", plspic$P.value)), cex=0.7)

par(mar = c(4,4,4,4))
barplot(pls$left.pls.vectors[,PLS], names.arg = rownames(pls$left.pls.vectors))

par(mar = c(4,4,4,4))
barplot(plspic$left.pls.vectors[,PLS], names.arg = rownames(plspic$left.pls.vectors))

  if (!is.null(which.curves1) & !is.null(which.curves2) & !is.null(which.LM)) {
  
par(mar=c(1,1,1,1))
    plot(preds$pred1[, view1[1]], preds$pred1[, view1[2]], asp=1, type = "n", axes = F, xlab = "", ylab = "")
    points(preds$pred1[which.LM, view1[1]], preds$pred1[which.LM, view1[2]], pch = 20, col = "black")
    lines(preds$pred1[which.curves1, view1[1]], preds$pred1[which.curves1, view1[2]], col = "black")
    points(preds$pred2[which.LM, view1[1]], preds$pred2[which.LM, view1[2]], pch = 20, col = "red")
    lines(preds$pred2[which.curves1, view1[1]], preds$pred2[which.curves1, view1[2]], col = "red")

par(mar=c(1,1,1,1))
    plot(preds$pred1[, view2[1]], -preds$pred1[, view2[2]], asp=1, type = "n", axes = F, xlab = "", ylab = "")
    points(preds$pred1[which.LM, view2[1]], -preds$pred1[which.LM, view2[2]], pch = 20, col = "black")
    lines(preds$pred1[which.curves2, view2[1]], -preds$pred1[which.curves2, view2[2]], col = "black")
    points(preds$pred2[which.LM, view2[1]], -preds$pred2[which.LM, view2[2]], pch = 20, col = "red")
    lines(preds$pred2[which.curves2, view2[1]], -preds$pred2[which.curves2, view2[2]], col = "red")
        
    if (!PIC) {
par(mar=c(1,1,1,1))
    plot(predsphylo$pred1[, view1[1]], predsphylo$pred1[, view1[2]], asp=1, type = "n", axes = F, xlab = "", ylab = "")
    points(predsphylo$pred1[which.LM, view1[1]], predsphylo$pred1[which.LM, view1[2]], pch = 20, col = "black")
    lines(predsphylo$pred1[which.curves1, view1[1]], predsphylo$pred1[which.curves1, view1[2]], col = "black")
    points(predsphylo$pred2[which.LM, view1[1]], predsphylo$pred2[which.LM, view1[2]], pch = 20, col = "red")
    lines(predsphylo$pred2[which.curves1, view1[1]], predsphylo$pred2[which.curves1, view1[2]], col = "red")} else {plot(1,1, axes = F, type = "n", xlab = "", ylab = "")}
    
    if (!PIC) {    
par(mar=c(1,1,1,1))
    plot(predsphylo$pred1[, view2[1]], -predsphylo$pred1[, view2[2]], asp=1, type = "n", axes = F, xlab = "", ylab = "")
    points(predsphylo$pred1[which.LM, view2[1]], -predsphylo$pred1[which.LM, view2[2]], pch = 20, col = "black")
    lines(predsphylo$pred1[which.curves2, view2[1]], -predsphylo$pred1[which.curves2, view2[2]], col = "black")
    points(predsphylo$pred2[which.LM, view2[1]], -predsphylo$pred2[which.LM, view2[2]], pch = 20, col = "red")
    lines(predsphylo$pred2[which.curves2, view2[1]], -predsphylo$pred2[which.curves2, view2[2]], col = "red")} else {plot(1,1, axes = F, type = "n", xlab = "", ylab = "")}

  }
  
  else {
par(mar=c(1,1,1,1))
    plot(preds$pred1[, view1[1]], preds$pred1[, view1[2]], asp=1, type = "n", axes = F, xlab = "", ylab = "")
    points(preds$pred1[, view1[1]], preds$pred1[, view1[2]], pch = 20, col = "black")
    points(preds$pred2[, view1[1]], preds$pred2[, view1[2]], pch = 20, col = "red")

par(mar=c(1,1,1,1))
    plot(preds$pred1[, view2[1]], -preds$pred1[, view2[2]], asp=1, type = "n", axes = F, xlab = "", ylab = "")
    points(preds$pred1[, view2[1]], -preds$pred1[, view2[2]], pch = 20, col = "black")
    points(preds$pred2[, view2[1]], -preds$pred2[, view2[2]], pch = 20, col = "red")
    
    if (!PIC) { 
par(mar=c(1,1,1,1))
    plot(predsphylo$pred1[, view1[1]], predsphylo$pred1[, view1[2]], asp=1, type = "n", axes = F, xlab = "", ylab = "")
    points(predsphylo$pred1[, view1[1]], predsphylo$pred1[, view1[2]], pch = 20, col = "black")
    points(predsphylo$pred2[, view1[1]], predsphylo$pred2[, view1[2]], pch = 20, col = "red")
    } else {plot(1,1, axes = F, type = "n", xlab = "", ylab = "")}  
        
    if (!PIC) { 
par(mar=c(1,1,1,1))
    plot(predsphylo$pred1[, view2[1]], -predsphylo$pred1[, view2[2]], asp=1, type = "n", axes = F, xlab = "", ylab = "")
    points(predsphylo$pred1[, view2[1]], -predsphylo$pred1[, view2[2]], pch = 20, col = "black")
    points(predsphylo$pred2[, view2[1]], -predsphylo$pred2[, view2[2]], pch = 20, col = "red")
    } else {plot(1,1, axes = F, type = "n", xlab = "", ylab = "")}    
  }
}

#-------------------------------
# Function to symmetrize a shape
#-------------------------------

maksym <- function(shp, axis = 1, reorder = NULL) {

ishp <- shp
ishp[,axis] <- -ishp[,axis]
symet <- (shp+ishp[reorder,])/2

symet}

#---------------------------------------------------------------------
# Function for alignment of entire configurations and curves, but based only on landmarks (curve points not used in the alignement itself)
#---------------------------------------------------------------------

pgpa.curves <- function(A, Acurv) { #A = Array with reference configurations, landmarks only; Acurv = Array with identical landmarks as in A, plus associated curves.
  p<-dim(A)[1]
  k<-dim(A)[2]
  n<-dim(A)[3]
  
  if (dim(Acurv)[3] != n) {stop("Both arrays must have the same number of specimens")}
  
  temp2<-temp1<-array(NA, dim=c(p,k,n))
  temp3<-temp4<-array(NA, dim=dim(Acurv))
  Siz<-numeric(n)
  ctroid<-list()
  
  for (i in 1:n)
       {Acs<-centsiz(A[,,i])
        Siz[i]<-Acs[[1]]
        trs <- trans1(Acs[[2]])
        temp1[,,i]<-trs
        ctr <- attr(trs, "scaled:center")
        ctroid[[i]] <- ctr
        scl <- Acurv[,,i]/Acs[[1]]
        trscurv <- t(t(scl) - ctr)
        temp3[,,i] <- trscurv}
        
  Qm1<-dist(t(matrix(temp1,k*p,n)))
  Q<-sum(Qm1)
  iter<-0
  
  while (abs(Q)>0.00001) {
  for(i in 1:n){
       M<-mshape(temp1[,,-i])
       sup <- pPsup(temp1[,,i],M)
       temp2[,,i]<- sup[[1]]
       rot <- sup$rotation
       temp4[,,i] <- temp3[,,i] %*% rot
       }
    Qm2<-dist(t(matrix(temp2,k*p,n)))
    Q<-sum(Qm1)-sum(Qm2)
    Qm1<-Qm2
    iter=iter+1
    temp1<-temp2
    temp3<-temp4}
 list("rotated"=temp2, "curves_rot"=temp4, "it.number"=iter, "Q"=Q, "intereucl.dist"=Qm2, "mshape"=centsiz(mshape(temp2)), "cent.size"=Siz)}
 
#------------------------------------
#Function to compute Sum of variances
#------------------------------------

SOV.fac <- function(x, fac, bootstrap = F, it = NULL) {

nlev <- length(levels(fac))
SOV <- rep(NA, nlev)
  # Make empty vector to store results

if (bootstrap) {mbSOV <- matrix(NA, ncol = nlev, nrow = it)}

for (i in 1:nlev) {
 
    submorphospace <- x[which(fac == levels(fac)[i]),]
      # Select data to compute variance: all observations of category i x all variables.   
    vars <- apply(submorphospace, 2, var)
      # Compute variances for each individual PC
    SOV[i] <- sum(vars)
    
    if (bootstrap) {
    bootSOV <- rep(NA, it)
      for (j in 1:it) {
      bootspace <- submorphospace[sample(nrow(submorphospace), replace = T),]
      bvars <- apply(bootspace, 2, var)
      bootSOV[j] <- sum(bvars)
      mbSOV[,i] <- bootSOV
      }
    }  
  }
names(SOV) <- levels(fac)
if (!bootstrap) {return(SOV)} else {return(list(SOV = SOV, bsSOV = mbSOV))}
}

#---------------------------------
#Function to compute Sum of ranges
#---------------------------------

SOR.fac <- function(x, fac, bootstrap = F, it = NULL) {

nlev <- length(levels(fac))
SOR <- rep(NA, nlev)
  # Make empty vector to store results

if (bootstrap) {mbSOR <- matrix(NA, ncol = nlev, nrow = it)}

for (i in 1:nlev) {
 
    submorphospace <- x[which(fac == levels(fac)[i]),]
      # Select data to compute variance: all observations of category i x all variables.   
    ranges <- apply(apply(submorphospace, 2, range), 2, diff)
    # Find minimum and maximum values for each PC, then substract them from 
    # each other to obtain the range value
    SOR[i] <- sum(ranges)
    
    if (bootstrap) {
    bootSOR <- rep(NA, it)
      for (j in 1:it) {
      bootspace <- submorphospace[sample(nrow(submorphospace), replace = T),]
      branges <- apply(apply(bootspace, 2, range), 2, diff)
      bootSOR[j] <- sum(branges)
      mbSOR[,i] <- bootSOR
      }
    }  
  }

names(SOR) <- levels(fac)
if (!bootstrap) {return(SOR)} else {return(list(SOR = SOR, bsSOR = mbSOR))}
}

#------------------------------------------------
#Function to return both min and max of a vector
#------------------------------------------------

minetmax <- function(x) {c(min(x), max(x))}
maxetmin <- function(x) {c(max(x), min(x))}
#------------------------------------------------
#Function to create shapes deformed along PC axes
#------------------------------------------------

shp.deform <- function(shp, PCA, axis = 1, magnification = 1) {
dims <- dim(shp)
v <- as.vector(t(shp))
vmax <- v + magnification * max(PCA$x[,axis]) * PCA$rotation[,axis]
vmin <- v + magnification * min(PCA$x[,axis]) * PCA$rotation[,axis]

shpmax <- matrix(vmax, nrow = dims[1], ncol = dims[2], byrow = T)
shpmin <- matrix(vmin, nrow = dims[1], ncol = dims[2], byrow = T)

list(max = shpmax, min = shpmin)}

#---------------------------------
# Function to compute and plot lda
#---------------------------------

ldawithplot <- function(mat, grouping, axes = c(1,2), pos.inset = c(0.7,1,0.7,1), pos.legend = "topleft") {

library(MASS)

x <- mat
grouping <- grouping

if (anyNA(grouping)) {
x <- x[-which(is.na(grouping)),]
grouping <- grouping[-which(is.na(grouping))]
}

if (anyNA(x)) {
grouping <- grouping[-which(is.na(x))]
x <- x[-which(is.na(x)),]
}

grouping <- droplevels(grouping)
ngrp <- nlevels(grouping)

discrim <- lda(x = x, grouping = grouping)
pdm <-predict(discrim)

clr <- c(1:ngrp)[grouping]

plot(pdm$x[,axes], pch = 21, bg = clr, asp = 1)
legend(pos.legend, legend = levels(grouping), bty ="n", pt.bg = c(1:ngrp), pch = 21)

scal <- discrim$scaling

par(fig = pos.inset, new = T, mar=c(0,0,4,2))
plot(scal[,axes], asp = 1, type = "n", 
xlab = "", ylab = "", 
bty = "n", xaxt = "n", yaxt = "n")
abline(v=0,h=0, col="grey", lty = 2)
for (i in 1:dim(scal)[1]) {lines(x = c(0, scal[i,axes[1]]), y = c(0, scal[i,axes[2]]), lwd = 2, col = "gray")}
text(scal[,axes], labels = rownames(scal))

list(discrim = discrim, predictions = pdm)}

#--------------------------------------------------------------------
# Conservative matching, trying to keep specimens which are "sp.", or a different species of a genus that's included in the phylogeny
#---------------------------------------------------------------------

sp.match.conservative <- function(x, y) {
m <- na.omit(match(x,y))
idx <- attr(m, which = "na.action")
msng <- x[idx]
n <- numeric()
  for (i in 1:length(msng)) {
  g <- strsplit(msng[i], split = "_")[[1]][1]
  gmtch <- grep(pattern = g, x = y)
    if (!any(gmtch %in% m)) {n <- c(n, gmtch[1])}
    else {n <- c(n,NA)}
  }

list(tokeep = c(m,na.omit(n)), namestomodify = na.omit(n))}



#---------------------------------------------------------------------

