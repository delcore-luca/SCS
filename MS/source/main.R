#####################
## SCS - main code ##
#####################
##
## Author: Luca Del Core

cat("\nInstall/load packages")
inst.pkgs <- installed.packages()

## required packages:
l.pkgs <- c("Matrix", 
            "splines", 
            "splines2", 
            "nloptr", 
            "scatterplot3d", 
            "scales",
            "parallel",
            "vegan",
            "beeswarm",
            "xtable")
## check if packages are installed
lapply(l.pkgs, function(pkg){
  if(!(pkg %in% rownames(inst.pkgs))){
    install.packages(pkg)
  }
})

lapply(l.pkgs, function(pkg){library(pkg,character.only=TRUE)})


##############
## currTime ##
##############

currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- substr(currTime, start = 1, stop = 13)

#############
## folders ##
#############

sourcePath <- paste("./source/")
isMatrices <- "./input/MS/matrices/"

StatsFolder <- "./input/MS/stats/"
allResPath <- "./results/"
resPath <- "./results/MS/"

currAllResPath <- paste(allResPath, "/", currTime, sep = "")
currResPath <- paste(resPath, "/", currTime, sep = "")
currOutPath <- paste(currResPath, "/output/", sep = "")
currFigPath <- paste(currResPath, "/figures/", sep = "")

ifelse(!dir.exists(file.path(allResPath)), dir.create(file.path(allResPath)), FALSE)
ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
ifelse(!dir.exists(file.path(currResPath)), dir.create(file.path(currResPath)), FALSE)
ifelse(!dir.exists(file.path(currOutPath)), dir.create(file.path(currOutPath)), FALSE)
ifelse(!dir.exists(file.path(currFigPath)), dir.create(file.path(currFigPath)), FALSE)


cat("\nload functions")
## load functions:
source(paste(sourcePath,"fun.R" , sep = ""))

cat("\nload data")
metadata <- read.csv(paste(StatsFolder, "/metadata_MS.csv", sep = ""), sep = ",", dec = ".", header = T, row.names = 1, stringsAsFactors = FALSE)
ISmat <- read.csv(paste(isMatrices, "/ISmat_MS.tsv", sep = ""), sep = "\t", header = T, row.names = NULL, check.names = FALSE) 
ISmat <- Matrix(as.matrix(ISmat), sparse = TRUE)

## compute Shannon entropy:
metadata$h.idx <- NA
h.idx <- apply(ISmat, 2, function(s){
  s <- s[s>0];
  p <- s/sum(s);
  h <- -sum(p*log(p));
  return(h)
})
## compute sequencing depth:
metadata$seqDepth <- NA
seq.depth <- apply(ISmat, 2, function(s){
  s <- s[s>0];
  sd <- sum(s)
  return(sd)
})
## compute nIS:
metadata$nIS<- NA
nIS <- apply(ISmat, 2, function(s){
  nIS <- sum(s>0);
  return(nIS)
})

metadata[names(h.idx), "h.idx"] <- h.idx
metadata[names(seq.depth), "seqDepth"] <- seq.depth
metadata[names(nIS), "nIS"] <- nIS
yX <- metadata

## Table 1:
xtable(x = summary(yX[which(yX$VectorID == "PGK"), c("DNAngUsed", "VCN" ,"N.mice.pool", "seqDepth", "nIS")]))
xtable(x = summary(yX[which(yX$VectorID == "LV.SF.LTR"), c("DNAngUsed", "VCN" ,"N.mice.pool", "seqDepth", "nIS")]))

cat("\nSCS rescaling...")
## confounding factors:
x.L <- list(DNA = as.numeric(yX[,"DNAngUsed"]),
            VCN = as.numeric(yX[,"VCN"]),
            PS = as.numeric(yX[,"N.mice.pool"]),
            SD = as.numeric(yX[,"seqDepth"]))
nKTS.L <- list(kts.1 = 2, kts.2 = 2, kts.3 = 2, kts.4 = 2)
sp.ord.L <- list(sp.ord.1 = 2, sp.ord.2 = 2, sp.ord.3 = 2, sp.ord.4 = 2)
confounders.L <- names(x.L)

## Generating matrix of additional factors of interest:
vecMrkComb <- as.matrix(unique.matrix(expand.grid(yX$VectorID, as.character(yX$CellMarker)), MARGIN = 1))
rownames(vecMrkComb) <- 1:nrow(vecMrkComb)
AF.mats <- get.AF(yX)


## Frequentist Model Averaging (FMA):
res.FMA <- fma(yX, x.L, nKTS.L, sp.ord.L, AF.mats, RELTOL = 1e-4, graphic = F)
yX$h.idx.res <- NA
yX$h.idx.res <- as.numeric(res.FMA$residuals)
b_MLE <- res.FMA$par
b_MLE_cfd <- res.FMA$par.cfd
s2_MLE <- res.FMA$s2_MLE

###########
## Plots ##
###########

pdf(file = paste(currFigPath, "h_confounders.pdf", sep = ""), width = 12, height = 8)
par(mar = c(5,5,2,2), mfrow = c(2,2))
plot(yX$DNAngUsed, yX$h.idx, cex.axis = 2, cex.lab = 2, xlab = "DNA (ng)", ylab = "entropy", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("a", region="figure", pos="topleft", cex = 3, font = 2)
plot(yX$VCN, yX$h.idx, cex.axis = 2, cex.lab = 2, xlab = "VCN", ylab = "entropy", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("b", region="figure", pos="topleft", cex = 3, font = 2)
plot(yX$N.mice.pool, yX$h.idx, cex.axis = 2, cex.lab = 2, xlab = "pool size (n. mice)", ylab = "entropy", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("c", region="figure", pos="topleft", cex = 3, font = 2)
plot(yX$seqDepth, yX$h.idx, cex.axis = 2, cex.lab = 2, xlab = "seq. depth (cells)", ylab = "entropy", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("d", region="figure", pos="topleft", cex = 3, font = 2)
dev.off()

pdf(file = paste(currFigPath, "frequentistPosterior.pdf", sep = ""), width = 7, height = 5)
par(mar = c(5,5,2,2), mfrow = c(1,1))
barplot(res.FMA$probs, ylab = "posterior probability", xlab = "model", 
        cex.axis = 2, cex.lab = 2, col = "black", font = 2, cex.names = 2,
        names.arg = 1:15)
fig_label("a", region="figure", pos="topleft", cex = 3, font = 2)
dev.off()


pdf(file = paste(currFigPath, "frequentistPosterior_VarProbs.pdf", sep = ""), width = 7, height = 5)
par(mar = c(5,5,2,2), mfrow = c(1,1))
barplot(res.FMA$varProbs, ylab = "posterior probability", xlab = "variable", 
        cex.axis = 2, cex.lab = 2, col = "black", font = 2, cex.names = 2,
        # names.arg = names(res.FMA$varProbs)
        names.arg = c("DNA", "VCN", "PS", "SD")
)
fig_label("b", region="figure", pos="topleft", cex = 3, font = 2)
dev.off()

#################
## Fitting SCS ##
#################
##

C <- get.C(H = res.FMA$H, 
           AT = res.FMA$AT, 
           CM = res.FMA$CM, 
           C.con = res.FMA$C.con,
           b = res.FMA$par, 
           s2 = res.FMA$s2_MLE, 
           nSim = 1000)
z.alpha <- qnorm(p = 1 - .05/2, mean = 0, sd = 1)

mrk_lab <- c("a", "b", "c", "d")
names(mrk_lab) <- as.vector(unique(yX$CellMarker))
pdf(file = paste(currFigPath, "rescaledFitted_FMA_smoothing.pdf", sep = ""), width = 12, height = 8)
par(mar = c(5,5,2,2), mfrow = c(2,2))
for (mrk in as.vector(unique(yX$CellMarker))) {
  
  xLim <- range(yX[, "TimePoint"])
  yLim <- range(yX[, "h.idx.res"])
  
  vecMrk <- which(vecMrkComb[,2] == mrk)[1]
  
  yX_vecMrk <- yX[which(yX$VectorID == vecMrkComb[vecMrk,1] & yX$CellMarker == vecMrkComb[vecMrk,2]),]
  x <- yX_vecMrk$TimePoint
  y.pred <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept","af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)])
  
  y.pred.sd <- z.alpha*sqrt(s2_MLE * diag((AF.mats$H.pred%*%AF.mats$AT) %*% C[colnames(AF.mats$AT),colnames(AF.mats$AT)] %*% t(AF.mats$H.pred%*%AF.mats$AT)))[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)]
  
  y.pred.sd.minus <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept","af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)] - y.pred.sd)
  y.pred.sd.plus <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept","af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)] + y.pred.sd)
  
  plot(yX_vecMrk$TimePoint, yX_vecMrk$h.idx.res, xlim = xLim, ylim = yLim,
       pch = 19, xlab = "time (days)", ylab = "SCS-entropy", cex.axis = 2, cex.lab = 2, # pch = pchPCR[yX_vecMrk$PCRMethod]
       main = mrk, cex.main = 2, col = alpha("red", .5), cex = 2)
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred, lwd = 3, col = "red")
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred.sd.minus, lwd = 3, col = "red", lty = 2)
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred.sd.plus, lwd = 3, col = "red", lty = 2)
  
  vecMrk <- which(vecMrkComb[,2] == mrk)[2]
  
  yX_vecMrk <- yX[which(yX$VectorID == vecMrkComb[vecMrk,1] & yX$CellMarker == vecMrkComb[vecMrk,2]),]
  x <- yX_vecMrk$TimePoint
  y.pred <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept", "af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)])
  
  y.pred.sd <- z.alpha*sqrt(s2_MLE * diag((AF.mats$H.pred%*%AF.mats$AT) %*% C[colnames(AF.mats$AT),colnames(AF.mats$AT)] %*% t(AF.mats$H.pred%*%AF.mats$AT)))[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)]
  
  y.pred.sd.minus <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept","af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)] - y.pred.sd)
  y.pred.sd.plus <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept","af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)] + y.pred.sd)
  
  points(yX_vecMrk$TimePoint, yX_vecMrk$h.idx.res, 
         pch = 17, col = alpha("blue", .5), cex = 2)
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred, lwd = 3, col = "blue")
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred.sd.minus, lwd = 3, col = "blue", lty = 2)
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred.sd.plus, lwd = 3, col = "blue", lty = 2)
  
  fig_label(mrk_lab[mrk], region="figure", pos="topleft", cex = 3, font = 2)
}
dev.off()



## Fit standard LM:
fit <- lm(formula = "log(h.idx) ~ bs(x = TimePoint, knots = quantile(TimePoint, probs = c(.33,.66)), degree = 2)*VectorID*CellMarker", data = yX)

mrks <- as.vector(unique(yX$CellMarker))
vcts <- as.vector(unique(yX$VectorID))
newdat_all <- do.call(rbind,lapply(mrks, function(mrk){
  do.call(rbind, lapply(vcts, function(vct){
    return(data.frame(TimePoint = seq(min(yX[which(yX$CellMarker == mrk & yX$VectorID == vct),]$TimePoint), 
                                      max(yX[which(yX$CellMarker == mrk & yX$VectorID == vct),]$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = vct))
  }))
}))

conf_interval_all <- exp(predict(fit, 
                                 newdata=newdat_all, 
                                 interval="confidence",
                                 level = 0.95))

## Plot results of standard LM:
pdf(file = paste(currFigPath, "standardLM.pdf", sep = ""), width = 12, height = 8)
par(mar = c(5,5,2,2), mfrow = c(2,2))
for (mrk in unique(yX$CellMarker)) {
  yX_mrk <- yX[which(yX$CellMarker == mrk),]
  
  
  yX_mrk_vID <- yX[which(yX$CellMarker == mrk & yX$VectorID == "LV.SF.LTR"),]
  newdat = data.frame(TimePoint = seq(min(yX_mrk_vID$TimePoint), max(yX_mrk_vID$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = "LV.SF.LTR")
  newdat$pred = exp(predict(fit, newdata = newdat))
  conf_interval <- exp(predict(fit, 
                               newdata=newdat, 
                               interval="confidence",
                               level = 0.95))
  
  plot(h.idx ~ TimePoint, data = yX_mrk_vID, cex = 2,
       main = mrk, cex.main = 2,
       pch = 19, cex.axis = 2, cex.lab = 2, col = alpha("red", .5),
       xlim = range(yX$TimePoint),
       ylim = range(yX_mrk$h.idx, conf_interval_all[,2:3]), xlab = "time (days)", ylab = "OBS-entropy")
  with(newdat, lines(x = TimePoint, y = pred, col = "red", lwd = 3))
  
  lines(newdat$TimePoint, conf_interval[,2], col="red", lty=2, lwd = 3)
  lines(newdat$TimePoint, conf_interval[,3], col="red", lty=2, lwd = 3)
  
  ####
  yX_mrk_vID <- yX[which(yX$CellMarker == mrk & yX$VectorID == "PGK"),]
  newdat = data.frame(TimePoint = seq(min(yX_mrk_vID$TimePoint), max(yX_mrk_vID$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = "PGK")
  newdat$pred = exp(predict(fit, newdata = newdat))
  conf_interval <- exp(predict(fit, 
                               newdata=newdat, 
                               interval="confidence",
                               level = 0.95))
  points(h.idx ~ TimePoint, data = yX_mrk_vID, cex = 2,
         pch = 17, col = alpha("blue", .5))
  with(newdat, lines(x = TimePoint, y = pred, col = "blue", lwd = 3))
  
  lines(newdat$TimePoint, conf_interval[,2], col="blue", lty=2, lwd = 3)
  lines(newdat$TimePoint, conf_interval[,3], col="blue", lty=2, lwd = 3)
  fig_label(mrk_lab[mrk], region="figure", pos="topleft", cex = 3, font = 2)
  
}
dev.off()

## Plot legends:
pdf(file = paste(currFigPath, "vectorLegends.pdf", sep = ""), width = 16, height = 4)
par(mar = c(5,5,2,2), mfrow = c(2,2))
plot.new()
legend(x = "center", legend = c("PGK", "LTR"), col = c("blue", "red"), pch = c(17, 20), cex = 4, lwd = 7, lty = 1, horiz = T)
dev.off()

########## additional plots #########

#################################################################################################
############################# COMPARISON WITH OTHER METHODS #####################################
#################################################################################################

cMin <- min(colSums(ISmat)) ## select the rarefaction level
set.seed(1) ## set seed of replicability
nRep <- 1000 ## set the number of random sub-samplings
cl <- makeCluster(getOption("cl.cores", round(detectCores()/2))) ## select the number of cores
clusterExport(cl = cl, varlist = list("ISmat", "cMin", "rrarefy2", "SRS_2")) # export variables to the cluster
# export packages to the cluster:
clusterEvalQ(cl, library("Matrix")) 
clusterEvalQ(cl, library("vegan"))

######## RAREFACTION ########
h.idx.RAR_allRep <- parApply(cl = cl, X = matrix(1:nRep), MARGIN = 1, FUN = function(r){
  message(r)
  ISmat_JY_RAR <- t(Matrix(rrarefy2(x = t(ISmat), cMin)));
  h.idx.RAR <- apply(ISmat_JY_RAR, 2, function(s){
    s <- s[s>0];
    p <- s/sum(s);
    h <- -sum(p*log(p));
    return(h)
  });
  return(h.idx.RAR)
})
h.idx.RAR_avg <- apply(h.idx.RAR_allRep, 1, mean)
yX$h.idx.res.RAR <- NA
yX[names(h.idx.RAR_avg), "h.idx.res.RAR"] <- h.idx.RAR_avg
h.idx.RAR_sd <- apply(h.idx.RAR_allRep, 1, sd)
yX$h.idx.res.RAR.sd <- NA
yX[names(h.idx.RAR_sd), "h.idx.res.RAR.sd"] <- h.idx.RAR_sd


######### SRS #########
h.idx.SRS_allRep <- parApply(cl = cl, X = matrix(1:nRep), MARGIN = 1, FUN = function(r){
  message(r)
  ISmat_JY_SRS <- SRS_2(ISmat, cMin)
  
  h.idx.SRS <- apply(ISmat_JY_SRS, 2, function(s){
    s <- s[s>0];
    p <- s/sum(s);
    h <- -sum(p*log(p));
    return(h)
  });
  
  return(h.idx.SRS)
})

h.idx.SRS_avg <- apply(h.idx.SRS_allRep, 1, mean)
yX$h.idx.res.SRS <- NA
yX[names(h.idx.SRS_avg), "h.idx.res.SRS"] <- h.idx.SRS_avg
h.idx.SRS_sd <- apply(h.idx.SRS_allRep, 1, sd)
yX$h.idx.res.SRS.sd <- NA
yX[names(h.idx.SRS_sd), "h.idx.res.SRS.sd"] <- h.idx.SRS_sd

stopCluster(cl)


## Fit standard LM on RAR-rescaled entropy:
fit <- lm(formula = "log(h.idx.res.RAR) ~ bs(x = TimePoint, knots = quantile(TimePoint, probs = c(.33,.66)), degree = 2)*VectorID*CellMarker", data = yX)

mrks <- as.vector(unique(yX$CellMarker))
vcts <- as.vector(unique(yX$VectorID))
newdat_all <- do.call(rbind,lapply(mrks, function(mrk){
  do.call(rbind, lapply(vcts, function(vct){
    return(data.frame(TimePoint = seq(min(yX[which(yX$CellMarker == mrk & yX$VectorID == vct),]$TimePoint), 
                                      max(yX[which(yX$CellMarker == mrk & yX$VectorID == vct),]$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = vct))
  }))
}))

conf_interval_all <- exp(predict(fit, 
                                 newdata=newdat_all, 
                                 interval="confidence",
                                 level = 0.95))

## Plot results of standard LM:
pdf(file = paste(currFigPath, "standardLM-RAR.pdf", sep = ""), width = 12, height = 8)
par(mar = c(5,5,2,2), mfrow = c(2,2))
for (mrk in unique(yX$CellMarker)) {
  yX_mrk <- yX[which(yX$CellMarker == mrk),]
  
  
  yX_mrk_vID <- yX[which(yX$CellMarker == mrk & yX$VectorID == "LV.SF.LTR"),]
  newdat = data.frame(TimePoint = seq(min(yX_mrk_vID$TimePoint), max(yX_mrk_vID$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = "LV.SF.LTR")
  newdat$pred = exp(predict(fit, newdata = newdat))
  conf_interval <- exp(predict(fit, 
                               newdata=newdat, 
                               interval="confidence",
                               level = 0.95))
  
  plot(h.idx.res.RAR ~ TimePoint, data = yX_mrk_vID, cex = 2,
       main = mrk, cex.main = 2,
       pch = 19, cex.axis = 2, cex.lab = 2, col = alpha("red", .5),
       xlim = range(yX$TimePoint),
       ylim = range(yX_mrk$h.idx.res.RAR, conf_interval_all[,2:3]), xlab = "time (days)", ylab = "RAR-entropy")
  with(newdat, lines(x = TimePoint, y = pred, col = "red", lwd = 3))
  
  lines(newdat$TimePoint, conf_interval[,2], col="red", lty=2, lwd = 3)
  lines(newdat$TimePoint, conf_interval[,3], col="red", lty=2, lwd = 3)
  
  ####
  yX_mrk_vID <- yX[which(yX$CellMarker == mrk & yX$VectorID == "PGK"),]
  newdat = data.frame(TimePoint = seq(min(yX_mrk_vID$TimePoint), max(yX_mrk_vID$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = "PGK")
  newdat$pred = exp(predict(fit, newdata = newdat))
  conf_interval <- exp(predict(fit, 
                               newdata=newdat, 
                               interval="confidence",
                               level = 0.95))
  points(h.idx.res.RAR ~ TimePoint, data = yX_mrk_vID, cex = 2,
         pch = 17, col = alpha("blue", .5))
  with(newdat, lines(x = TimePoint, y = pred, col = "blue", lwd = 3))
  
  lines(newdat$TimePoint, conf_interval[,2], col="blue", lty=2, lwd = 3)
  lines(newdat$TimePoint, conf_interval[,3], col="blue", lty=2, lwd = 3)
  fig_label(mrk_lab[mrk], region="figure", pos="topleft", cex = 3, font = 2)
  
}
dev.off()

## Fit standard LM on SRS-rescaled entropy:
fit <- lm(formula = "log(h.idx.res.SRS) ~ bs(x = TimePoint, knots = quantile(TimePoint, probs = c(.33,.66)), degree = 2)*VectorID*CellMarker", data = yX)

mrks <- as.vector(unique(yX$CellMarker))
vcts <- as.vector(unique(yX$VectorID))
newdat_all <- do.call(rbind,lapply(mrks, function(mrk){
  do.call(rbind, lapply(vcts, function(vct){
    return(data.frame(TimePoint = seq(min(yX[which(yX$CellMarker == mrk & yX$VectorID == vct),]$TimePoint), 
                                      max(yX[which(yX$CellMarker == mrk & yX$VectorID == vct),]$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = vct))
  }))
}))

conf_interval_all <- exp(predict(fit, 
                                 newdata=newdat_all, 
                                 interval="confidence",
                                 level = 0.95))

## Plot results of standard LM:
pdf(file = paste(currFigPath, "standardLM-SRS.pdf", sep = ""), width = 12, height = 8)
par(mar = c(5,5,2,2), mfrow = c(2,2))
for (mrk in unique(yX$CellMarker)) {
  yX_mrk <- yX[which(yX$CellMarker == mrk),]
  
  
  yX_mrk_vID <- yX[which(yX$CellMarker == mrk & yX$VectorID == "LV.SF.LTR"),]
  newdat = data.frame(TimePoint = seq(min(yX_mrk_vID$TimePoint), max(yX_mrk_vID$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = "LV.SF.LTR")
  newdat$pred = exp(predict(fit, newdata = newdat))
  conf_interval <- exp(predict(fit, 
                               newdata=newdat, 
                               interval="confidence",
                               level = 0.95))
  
  plot(h.idx.res.SRS ~ TimePoint, data = yX_mrk_vID, cex = 2,
       main = mrk, cex.main = 2,
       pch = 19, cex.axis = 2, cex.lab = 2, col = alpha("red", .5),
       xlim = range(yX$TimePoint),
       ylim = range(yX_mrk$h.idx.res.SRS, conf_interval_all[,2:3]), xlab = "time (days)", ylab = "SRS-entropy")
  with(newdat, lines(x = TimePoint, y = pred, col = "red", lwd = 3))
  
  lines(newdat$TimePoint, conf_interval[,2], col="red", lty=2, lwd = 3)
  lines(newdat$TimePoint, conf_interval[,3], col="red", lty=2, lwd = 3)
  
  ####
  yX_mrk_vID <- yX[which(yX$CellMarker == mrk & yX$VectorID == "PGK"),]
  newdat = data.frame(TimePoint = seq(min(yX_mrk_vID$TimePoint), max(yX_mrk_vID$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = "PGK")
  newdat$pred = exp(predict(fit, newdata = newdat))
  conf_interval <- exp(predict(fit, 
                               newdata=newdat, 
                               interval="confidence",
                               level = 0.95))
  points(h.idx.res.SRS ~ TimePoint, data = yX_mrk_vID, cex = 2,
         pch = 17, col = alpha("blue", .5))
  with(newdat, lines(x = TimePoint, y = pred, col = "blue", lwd = 3))
  
  lines(newdat$TimePoint, conf_interval[,2], col="blue", lty=2, lwd = 3)
  lines(newdat$TimePoint, conf_interval[,3], col="blue", lty=2, lwd = 3)
  fig_label(mrk_lab[mrk], region="figure", pos="topleft", cex = 3, font = 2)
  
}
dev.off()

#################################
## Comparisons of correlations ##
#################################

cors <- cor(yX[,c("DNAngUsed", "VCN", "N.mice.pool", "seqDepth", "h.idx", "h.idx.res","h.idx.res.RAR","h.idx.res.SRS")], method = "spearman")[c("h.idx", "h.idx.res","h.idx.res.RAR","h.idx.res.SRS"), c("DNAngUsed", "VCN", "N.mice.pool", "seqDepth")]
corPvals <- t(apply(matrix(rownames(cors)), 1, function(r){apply(matrix(colnames(cors)), 1, function(c){as.numeric(cor.test(x = yX[,c], y = yX[,r], method = "spearman")$p.value)})}))
dimnames(corPvals) <- dimnames(cors)

rownames(cors) <- c("OBS","SCS","RAR","SRS")
rownames(corPvals) <- rownames(cors)
cors <- cors[c("OBS", "RAR", "SRS", "SCS"),]
corPvals <- corPvals[c("OBS", "RAR", "SRS", "SCS"),]

pdf(file = paste(currFigPath, "Comparisons_correlations.pdf", sep = ""), width = 12, height = 4)
par(mar = c(5,6,2,2), mfrow = c(1,4))
bp <- barplot(abs(cors[,1]), cex.axis = 2, cex.names = 2, names.arg = rownames(cors[,]), col = "black", ylab = expression(rho(h,DNA)), cex.lab = 2)
lab <- attr(regexpr("(?<=\\.)0+", corPvals[,1], perl = TRUE), "match.length")
lab[lab < 0] <- 0
lab <- apply(matrix(lab), 1, FUN = function(l){if(l > 0){
  return(paste(rep("*", l), collapse = ""))
}else{
  return("")
}})
text(bp, abs(cors[,1])*.9, 
     labels = lab, col = "red", font = 2, srt = 90, cex = 4 ,pos = 1)
fig_label("a", region="figure", pos="topleft", cex = 3, font = 2)
bp <- barplot(abs(cors[,2]), cex.axis = 2, cex.names = 2, names.arg = rownames(cors[,]), col = "black", ylab = expression(rho(h,VCN)), cex.lab = 2)
lab <- attr(regexpr("(?<=\\.)0+", formatC(corPvals[,2], format = "f", digits = 10), perl = TRUE), "match.length")
lab[lab < 0] <- 0
lab <- apply(matrix(lab), 1, FUN = function(l){if(l > 0){
  return(paste(rep("*", l), collapse = ""))
}else{
  return("")
}})
text(bp, abs(cors[,2])*.8, 
     labels = lab, col = "red", font = 2, srt = 90, cex = 4, pos = 1)
fig_label("b", region="figure", pos="topleft", cex = 3, font = 2)
bp <- barplot(abs(cors[,3]), cex.axis = 2, cex.names = 2, names.arg = rownames(cors[,]), col = "black", ylab = expression(rho(h,SD)), cex.lab = 2)
lab <- attr(regexpr("(?<=\\.)0+", formatC(corPvals[,3], format = "f", digits = 10), perl = TRUE), "match.length")
lab[lab < 0] <- 0
lab <- apply(matrix(lab), 1, FUN = function(l){if(l > 0){
  return(paste(rep("*", l), collapse = ""))
}else{
  return("")
}})
text(bp, abs(cors[,3])*.7, 
     labels = lab, col = "red", font = 2, srt = 90, cex = 4, pos = 1)
fig_label("c", region="figure", pos="topleft", cex = 3, font = 2)

bp <- barplot(abs(cors[,4]), cex.axis = 2, cex.names = 2, names.arg = rownames(cors[,]), col = "black", ylab = expression(rho(h,SD)), cex.lab = 2)
lab <- attr(regexpr("(?<=\\.)0+", formatC(corPvals[,4], format = "f", digits = 20), perl = TRUE), "match.length")
lab[lab < 0] <- 0
lab <- apply(matrix(lab), 1, FUN = function(l){if(l > 0){
  return(paste(rep("*", l), collapse = ""))
}else{
  return("")
}})
text(bp, abs(cors[,3])*.7, 
     labels = lab, col = "red", font = 2, srt = 90, cex = 4, pos = 1)
fig_label("c", region="figure", pos="topleft", cex = 3, font = 2)
dev.off()


pdf(file = paste(currFigPath, "Comparisons_boxplots.pdf", sep = ""), width = 7, height = 5)
boxplot(yX[, c("h.idx", "h.idx.res")], lwd = 3, cex.axis = 2, cex.lab = 2, cex = 2, names = c("OBS", "SCS"), outline = F)
beeswarm(yX[, c("h.idx", "h.idx.res")], add=T, pch = 20, cex = 1.5, col = "red")
dev.off()

pdf(file = paste(currFigPath, "Comparisons_boxplots_allMethods.pdf", sep = ""), width = 15, height = 5)
par(mar = c(5,5,2,2))
boxplot(yX[, c("h.idx", "h.idx.res", "h.idx.res.RAR", "h.idx.res.SRS")], lwd = 3, cex.axis = 2, cex.lab = 2, cex = 2, names = c("OBS", "SCS" ,"RAR", "SRS"), outline = F, ylab = "entropy")
beeswarm(yX[, c("h.idx", "h.idx.res", "h.idx.res.RAR", "h.idx.res.SRS")], add=T, pch = 20, cex = 1.5, col = "red")
dev.off()
