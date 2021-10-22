
#############
## get.C
#############
##
##  Description: Get C matrix for computing the covariance matrix s^2*DCD'.
##
get.C <- function(H, AT, CM, C.con, b, s2, nSim){
  
  colnames(H)[1] <- "intercept"
  intercept.idxs <- which(colnames(H) %in% "intercept")
  constrained.idxs <- which(!(colnames(H) %in% "intercept"))
  
  intercept.idxs2 <- which(colnames(AT) %in% "intercept")
  constrained.idxs2 <- which(!(colnames(AT) %in% "intercept"))
  
  b <- b[c(intercept.idxs2, constrained.idxs2)]
  AT <- AT[c(intercept.idxs, constrained.idxs), c(intercept.idxs2, constrained.idxs2)]
  CM <- CM[,c(intercept.idxs, constrained.idxs)]
  H <- H[,c(intercept.idxs, constrained.idxs)]
  C.con <- C.con[,c(intercept.idxs2, constrained.idxs2)]
  
  h.eta <- as.numeric(H %*% AT %*% b)
  
  C.as <- lapply(1:nSim, function(ss){
    eta.s <- rnorm(n = length(h.eta), mean = h.eta, sd = sqrt(s2))
    J.s <- qprog(q = t(H%*%AT)%*%H%*%AT,
                 c = as.numeric(t(t(H%*%AT)%*%eta.s)),
                 amat = C.con,
                 b = rep(0, nrow(C.con)))$face
    
    if(length(J.s) > 0){
      C.s <- get.Cs(J.s, H, AT)
      return(list(C.s = C.s,
                  J.s = J.s))
    }else{
      return(NULL)
    }
    
  })
  
  C.as <-  C.as[!as.vector(unlist(lapply(C.as, function(l){is.null(l)})))]
  list_name <- lapply(C.as, function(l){l$J.s})
  pjs <- table(unlist(lapply(list_name, paste, collapse = " ")))/length(list_name)
  
  C <- Reduce("+", lapply(C.as, function(l){(l$C.s)*(pjs[paste(l$J.s, collapse = " ")])}))
  
  rownames(C) <- colnames(C) <- c(intercept.idxs2, constrained.idxs2)
  
  C <- C[as.character(sort(as.numeric(rownames(C)), decreasing = F)), as.character(sort(as.numeric(rownames(C)), decreasing = F))]
  rownames(C) <- colnames(C) <- names(b)
  
  return(C)
}

#############
## get.Cs
#############
##
##  Description: Get C matrix corresponding to face J.
##
get.Cs <- function(J, H, AT){
  X <- H %*% AT
  m <- ncol(X)
  Xj <- X[,c(1:4, J+4)]
  
  C <- diag(0, m, m)
  C[J+4, J+4] <- solve(t(Xj) %*% Xj)[-c(1:4),-c(1:4)]
  
  C[J+4, 1:4] <- solve(t(Xj) %*% Xj)[-c(1:4), c(1:4)]
  C[1:4, J+4] <- solve(t(Xj) %*% Xj)[c(1:4),-c(1:4)]
  
  C[1:4, 1:4] <- solve(t(Xj) %*% Xj)[1:4,1:4]
  
  return(C)
}

#############
## fig_label
#############
##
##  Description: Add a label to a plot into the external area.
##
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

###########
## get.AF
###########
##
## Description: Get all the matrices (design and constraints)
##              for the additional factors of interest
##
get.AF <- function(yX, sp.ord = 2, nKTS = 2){
  addFactors <- lapply(as.vector(unique(yX$CellMarker)), function(mrk){
    
    res.vID <- lapply(as.vector(unique(yX$VectorID)), function(vID){
      p <- nKTS + sp.ord
      
      yX_vecMrk <- yX[which(yX$VectorID == vID & yX$CellMarker == mrk),]
      x <- yX_vecMrk$TimePoint
      
      kts <- tail(head(seq(from = 0, to = max(x), length.out = nKTS + 2), -1), -1)
      bdKts <- c(0, max(x))
      
      ## design:
      H <- bs(x, intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)  
      rownames(H) <- rownames(yX_vecMrk)
      colnames(H) <- rep(paste(vID, mrk, collapse = "-"), ncol(H))
      H <- as.matrix(as.data.frame(H))
      ## design for prediction:
      H.pred <- bs(seq(min(x), max(x), length.out = 100), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)  
      colnames(H.pred) <- rep(paste(vID, mrk, collapse = "-"), ncol(H.pred))
      H.pred <- as.matrix(as.data.frame(H.pred))
      
      ## derivative of pline basis:
      dM <- tail(dbs(sort(x, decreasing = F), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts, derivs = 1), 1)
      AT <- matrix(data = 0, nrow = ncol(H), ncol = ncol(H) - 1)
      AT[tail(which(dM != 0), 1), head(which(dM != 0), 1)] <- -dM[head(which(dM != 0), 1)]/dM[tail(which(dM != 0), 1)]
      AT[-tail(which(dM != 0), 1),] <- diag(1, p - 1)
      
      ## linear inequality constraints for monotonicity:
      CC.mnt <- diag(1, p)
      CC.mnt <- diff(CC.mnt)
      CC.mnt <- rbind(c(1, rep(0, p - 1)),
                      CC.mnt)
      
      return(list(H = H,
                  H.pred = H.pred,
                  AT = AT,
                  CC.mnt = CC.mnt,
                  rn = rownames(H),
                  cn = colnames(H)
      ))
    })
    
    if(mrk != "MNC"){
      H <- cbind(1, bdiag(lapply(res.vID, function(v){v$H})))
      rownames(H) <- as.vector(unlist(lapply(res.vID, function(v){rownames(v$H)})))
      colnames(H) <- c("intercept", as.vector(unlist(lapply(res.vID, function(v){colnames(v$H)}))))
      
      H.pred <- cbind(1, bdiag(lapply(res.vID, function(v){v$H.pred})))
      
      colnames(H.pred) <- c("intercept", as.vector(unlist(lapply(res.vID, function(v){colnames(v$H.pred)}))))
      
      AT <- bdiag(1, bdiag(lapply(res.vID, function(v){v$AT})))
      colnames(AT) <- rep("af", ncol(AT))
      
      CC.mnt <- bdiag(0, bdiag(lapply(res.vID, function(v){v$CC.mnt})))
    }else{
      H <- bdiag(lapply(res.vID, function(v){v$H}))
      rownames(H) <- as.vector(unlist(lapply(res.vID, function(v){rownames(v$H)})))
      colnames(H) <- as.vector(unlist(lapply(res.vID, function(v){colnames(v$H)})))
      
      H.pred <- bdiag(lapply(res.vID, function(v){v$H.pred}))
      colnames(H.pred) <- as.vector(unlist(lapply(res.vID, function(v){colnames(v$H.pred)})))
      
      AT <- bdiag(lapply(res.vID, function(v){v$AT}))
      colnames(AT) <- rep("af", ncol(AT))
      
      CC.mnt <- bdiag(lapply(res.vID, function(v){v$CC.mnt}))
    }
    
    
    return(list(H = H,
                H.pred = H.pred,
                AT = AT,
                CC.mnt = CC.mnt,
                rn = rownames(H),
                cn = colnames(H)
    ))
  })
  
  AT.af <- bdiag(lapply(addFactors, function(v){v$AT}))
  colnames(AT.af) <- c(rep("af", (ncol(AT.af) + 1)/length(addFactors) - 1), rep(c("intercept", rep("af", (ncol(AT.af) + 1)/length(addFactors) - 1)), times = length(addFactors) - 1))
  
  
  H.af <- bdiag(lapply(addFactors, function(v){v$H}))
  rownames(H.af) <- as.vector(unlist(lapply(addFactors, function(v){rownames(v$H)})))
  colnames(H.af) <- as.vector(unlist(lapply(addFactors, function(v){colnames(v$H)})))
  H.af <- H.af[rownames(yX),]
  C.af <- bdiag(lapply(addFactors, function(v){v$CC.mnt}))
  
  H.pred.af <- bdiag(lapply(addFactors, function(v){v$H.pred}))
  colnames(H.pred.af) <- as.vector(unlist(lapply(addFactors, function(v){colnames(v$H.pred)})))
  
  return(list(H = H.af,
              H.pred = H.pred.af,
              AT = AT.af,
              CC.mnt = C.af,
              rn = rownames(H.af),
              cn = colnames(H.af)
  ))
}

###########
## get.CFDs
###########
##
## Description: Get all the matrices (design and constraints)
##              for the confounders
##
get.CFDs <- function(x.L, nKTS.L, sp.ord.L){
  H_CM.L <- lapply(1:length(x.L), function(v){
    
    x <- x.L[[v]]
    nKTS <- nKTS.L[[v]]
    sp.ord <- sp.ord.L[[v]]
    
    kts <- tail(head(seq(from = min(x), to = max(x), length.out = nKTS + 2), -1), -1)
    bdKts <- range(x)
    
    p <- nKTS + sp.ord
    H <- bs(x, intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)
    colnames(H) <- rep(names(x.L)[v], ncol(H))
    
    H.sorted <- bs(sort(x, decreasing = F), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)
    colnames(H.sorted) <- rep(names(x.L)[v], ncol(H))
    
    if(v %% 2 == 1){
      H.pred <- bs(rep(seq(min(x), max(x), length.out = 100), each = 100), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)
    }else{
      H.pred <- bs(rep(seq(min(x), max(x), length.out = 100), times = 100), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)
    }
    
    dM <- tail(dbs(sort(x, decreasing = F), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts, derivs = 1), 1)
    
    CC.mnt <- diag(1, p)
    CC.mnt <- diff(CC.mnt)
    CC.mnt <- tail(CC.mnt, 1)
    CC.mnt <- rbind(c(1, rep(0, p - 1)),
                    CC.mnt)
    
    CC.cncv <- diag(1, p)
    CC.cncv <- diff(diff(CC.cncv))
    
    CM <- as.matrix(rbind(-1*CC.mnt, CC.cncv))
    
    AT <- rbind(diag(1, p - 1),
                c(rep(0,p-2), -dM[,p-1]/dM[,p]))
    colnames(AT) <- rep(names(x.L)[v], ncol(AT))
    
    return(list(H = H,
                H.sorted = H.sorted,
                H.pred = H.pred,
                CM = CM,
                AT = AT))
  })
  
  AT.cfd <- as.matrix(bdiag(1, bdiag(lapply(H_CM.L, function(v){v$AT}))))
  colnames(AT.cfd) <- c("intercept", as.vector(unlist(lapply(H_CM.L, function(v){colnames(v$AT)}))))
  H.cfd <- cbind(1, Reduce(cbind, lapply(H_CM.L, function(v){v$H})))
  CM.cfd <- as.matrix(cbind(0, bdiag(lapply(H_CM.L, function(v){v$CM}))))
  
  return(list(H = H.cfd,
              H.sorted = cbind(1, Reduce(cbind, lapply(H_CM.L, function(v){v$H.sorted}))),
              H.pred = cbind(1, Reduce(cbind, lapply(H_CM.L, function(v){v$H.pred}))),
              CM = CM.cfd,
              AT = AT.cfd))
}

#######
## scs
#######
##
## Description: Fit the SCS for a particular set of confounders and factors of interest.
##
scs <- function(x.L, y, nKTS.L, sp.ord.L, confounders.L, AF.mats, RELTOL = 1e-4, graphic = F){
  
  H.af <- AF.mats$H
  AT.af <- AF.mats$AT
  C.af <- AF.mats$CC.mnt
  
  CFD.mats <- get.CFDs(x.L, nKTS.L, sp.ord.L)
  
  AT.cfd <- CFD.mats$AT
  H.cfd <- CFD.mats$H
  CM.cfd <- CFD.mats$CM
  
  
  H <- as.matrix(cbind(H.cfd, H.af))
  CM <- as.matrix(bdiag(CM.cfd, C.af))
  AT <- as.matrix(bdiag(AT.cfd, AT.af))
  colnames(AT) <- c(colnames(AT.cfd),
                    colnames(AT.af)) 
  
  # CM <- rbind(CM, cbind(0, diag(-1, ncol(H)-1)))
  C.con <- -CM %*% AT
  C.con <- C.con[rowSums(abs(C.con)) != 0,]
  C.con <- unique.matrix(C.con)
  
  D.mat <- t(H%*%AT)%*%H%*%AT
  
  res.opt <- solve.QP(D = D.mat,
                      dvec = as.numeric(t(t(H%*%AT)%*%log(y))),
                      Amat = t(C.con), # t(-CM %*% AT),
                      bvec = rep(0, nrow(C.con))) # rep(0, nrow(CM)))
  b_MLE <- res.opt$solution
  names(b_MLE) <- colnames(H%*%AT)
  
  nDF <-  qprog(q = as.matrix(D.mat),
                c = as.numeric(t(t(H%*%AT)%*%log(y))),
                amat = -CM %*% AT,
                b = rep(0, nrow(CM)))$df
  
  cpar <- 1.5
  s2_MLE <- as.numeric(t(log(y) - H%*%AT%*%b_MLE)%*%(log(y) - H%*%AT%*%b_MLE)/(nrow(H) - cpar*nDF))
  AIC_MLE <- 2*nDF - 2*sum(dnorm(x = log(y), mean = as.numeric(H%*%AT%*%b_MLE), sd = sqrt(s2_MLE), log = T)) + 2*nDF*(nDF+1)/(nrow(H) - nDF)
  BIC_MLE <- log(nrow(H))*nDF - 2*sum(dnorm(x = log(y), mean = as.numeric(H%*%AT%*%b_MLE), sd = sqrt(s2_MLE), log = T))
  
  b_MLE_cfd <- b_MLE[1:ncol(CFD.mats$AT)]
  
  if(graphic & length(x.L) >= 2){
    H.pred <- CFD.mats$H.pred
    
    par(mar = c(5,2,2,2))
    sc3 <- scatterplot3d(x.L[[1]], 
                         x.L[[2]], 
                         y, pch = 20, 
                         cex = 3, 
                         color = alpha("black", .5),
                         main = paste("AICc = ", round(AIC_MLE,2), sep = ""), cex.main = 2,
                         xlab = confounders.L[1], 
                         ylab = confounders.L[2], 
                         zlab = "entropy", 
                         cex.axis = 1.5, 
                         cex.lab = 1.5)
    
    # b_MLE_cfd.avgInt <- mean(b_MLE[names(b_MLE) %in% "intercept"])
    # sc3$points3d(rep(seq(min(x.L[[1]]), max(x.L[[1]]), length.out = 100), each = 100),
    #              rep(seq(min(x.L[[2]]), max(x.L[[2]]), length.out = 100), times = 100),
    #              exp(H.pred %*% AT.cfd %*% b_MLE_cfd.avgInt), col = alpha("green", .2), pch = 20, cex = .3)
    
    sc3$points3d(rep(seq(min(x.L[[1]]), max(x.L[[1]]), length.out = 100), each = 100),
                 rep(seq(min(x.L[[2]]), max(x.L[[2]]), length.out = 100), times = 100),
                 exp(H.pred %*% AT.cfd %*% b_MLE_cfd), col = alpha("green", .5), pch = 20, cex = .3)
  }
  
  ####
  res <- list()
  res$par <- b_MLE
  res$residuals <- exp(log(y) - H.cfd[,which(colnames(H.cfd) %in% confounders.L)]%*%(AT.cfd%*%b_MLE_cfd)[which(colnames(H.cfd) %in% confounders.L)])
  H.sorted <- CFD.mats$H.sorted
  H.diff <- matrix(data = tail(H.sorted,1), nrow = nrow(H.cfd), ncol = ncol(H.cfd), byrow = T) - H.cfd
  res$residualsFromMax <- exp(log(y) + H.diff[,which(colnames(H.cfd) %in% confounders.L)]%*%(AT.cfd%*%b_MLE_cfd)[which(colnames(H.cfd) %in% confounders.L)])
  res$AIC <- AIC_MLE
  res$BIC <- BIC_MLE
  res$nDF <- nDF
  res$s2 <- s2_MLE
  
  return(res)
}

#######
## fma
#######
##
## Description: Compute the Frequentist Model Averaging (FMA) estimator 
##              and all the corresponding statistics.
##
fma <- function(yX, x.L, nKTS.L, sp.ord.L, AF.mats, RELTOL, graphic = F){
  
  allAICRes <- apply(head(expand.grid(c(T,F), c(T,F), c(T,F), c(T,F)), -1), 1, 
                     FUN = function(model){
                       scs(x.L[t(model)], 
                           y = as.numeric(yX$h.idx), 
                           nKTS.L = nKTS.L[t(model)], 
                           sp.ord.L = sp.ord.L[t(model)], 
                           confounders.L = names(x.L[t(model)]), 
                           AF.mats = AF.mats, 
                           RELTOL = RELTOL, 
                           graphic = graphic)
                     })
  
  all_pars <- matrix(data = 0, nrow = length(allAICRes), ncol = length(allAICRes[[1]]$par))
  colnames(all_pars) <- names(allAICRes[[1]]$par)
  apply(matrix(1:length(allAICRes)), 1, function(model){all_pars[model, which(colnames(all_pars) %in% names(allAICRes[[model]]$par))] <<- allAICRes[[model]]$par})
  
  probs <- exp(-as.vector(unlist(lapply(allAICRes, function(l){l$BIC})))/2)/sum(exp(-as.vector(unlist(lapply(allAICRes, function(l){l$BIC})))/2))
  varProbs <- apply(head(expand.grid(c(T,F), c(T,F), c(T,F), c(T,F)), -1), 2, function(v){sum(probs[v])})
  names(varProbs) <- c("DNA", "VCN", "PS", "SD")
  b_FMA <- colSums(all_pars * t(t(probs))[,rep(1, ncol(all_pars))])
  
  H.af <- AF.mats$H
  AT.af <- AF.mats$AT
  C.af <- AF.mats$CC.mnt
  
  CFD.mats <- get.CFDs(x.L, nKTS.L, sp.ord.L)
  
  AT.cfd <- CFD.mats$AT
  H.cfd <- CFD.mats$H
  CM.cfd <- CFD.mats$CM
  
  H <- as.matrix(cbind(H.cfd, H.af))
  CM <- as.matrix(bdiag(CM.cfd, C.af))
  AT <- as.matrix(bdiag(AT.cfd, AT.af))
  colnames(AT) <- c("intercept",
                    tail(colnames(AT.cfd), -1),
                    colnames(AT.af)) 

  C.con <- -CM %*% AT
  C.con <- C.con[rowSums(abs(C.con)) != 0,]
  C.con <- unique.matrix(C.con)
  
  nDF <- mean(as.vector(unlist(lapply(allAICRes, function(l){l$nDF}))))
  s2_MLE <- mean(as.vector(unlist(lapply(allAICRes, function(l){l$s2}))))
  
  res <- list()
  res$par <- b_FMA
  res$nDF <- nDF
  res$s2_MLE <- s2_MLE
  
  H.cfd <- CFD.mats$H
  AT.cfd <- CFD.mats$AT
  b_MLE_cfd <- b_FMA[1:ncol(AT.cfd)]
  res$par.cfd <- b_MLE_cfd
  

  res$residuals <- exp(log(yX$h.idx) - H.cfd[,which(colnames(H.cfd) %in% confounders.L)]%*%(AT.cfd%*%b_MLE_cfd)[which(colnames(H.cfd) %in% confounders.L)])
  
  H.sorted <- get.CFDs(x.L, nKTS.L, sp.ord.L)$H.sorted
  H.diff <- matrix(data = tail(H.sorted,1), nrow = nrow(H.cfd), ncol = ncol(H.cfd), byrow = T) - H.cfd
  res$residualsFromMax <- exp(log(yX$h.idx) + H.diff[,which(colnames(H.cfd) %in% confounders.L)]%*%(AT.cfd%*%b_MLE_cfd)[which(colnames(H.cfd) %in% confounders.L)])
  res$probs <- probs
  res$varProbs <- varProbs
  res$H <- H
  res$CM <- CM
  res$AT <- AT
  res$C.con <- C.con
  
  return(res)
}

#########
## SRS_2
#########
##
## Description: reimplementation of SRS's SRS() function which allows the use of sparse matrices.
##
SRS_2 <- function (data, Cmin) 
{
  if (Cmin > min(colSums(data))) {
    print("ERROR: Cmin > minimum library size. Please select a Cmin that is <= the minimum library size of the dataset.")
  }
  else {
    if (Cmin < 0) {
      print("ERROR: Cmin < 0. Please select a Cmin >= 0.")
    }
    else {
      if (Cmin%%1 > 0) {
        print("ERROR: Please select a Cmin without decimal places.")
      }
      else {
        counter = 1
        for (i in seq(1, ncol(data), 1)) {
          if (i == 1) {
            fixed_factor <- (data[, i]/(sum(data[, i])/Cmin))
            assign(paste(colnames(data)[i], sep = ""), fixed_factor)
            fixed_factor_1 <- Matrix(get(colnames(data)[i]))
            colnames(fixed_factor_1)[i] <- colnames(data)[i]
          }
          else {
            fixed_factor <- (data[, i]/(sum(data[, i])/Cmin))
            assign(paste(colnames(data)[i], sep = ""), fixed_factor)
            fixed_factor_1 <- cbind(fixed_factor_1, fixed_factor)
            colnames(fixed_factor_1)[{
              counter = counter + 1
            }] <- colnames(data)[i]
          }
        }
        
        revtrunc_fixed_factor_1 <- floor(fixed_factor_1)
        
        diff_counts <- Cmin - colSums(revtrunc_fixed_factor_1)
        
        revtrunc <- function(x) {
          sign(x) * (x - floor(x))
        }
        revtrunc_fixed_factor <- (round(revtrunc(fixed_factor_1), 1e+07))
        
        x <- revtrunc_fixed_factor
        counter = 1
        for (i in seq(1, ncol(x), 1)) {
          if (i == 1) {
            if (diff_counts[i] == 0) {
              fixed_factor <- revtrunc_fixed_factor_1[, i]
              assign(paste(colnames(data)[i], sep = ""), fixed_factor)
              fixed_factor_1 <- Matrix(get(colnames(data)[i]))
              colnames(fixed_factor_1)[i] <- colnames(data)[i]
            }
            else {
              maxN <- function(x, N = diff_counts[i]) {
                len <- length(x)
                if (N > len) {
                  warning("N greater than length(x).  Setting N=length(x)")
                  N <- length(x)
                }
                sort(x, partial = len - N + 1)[len - N + 1]
              }
              max <- which(x[, i] == maxN(x[, i]), arr.ind = TRUE)
              
              normalization_value <- diff_counts[i] - sum(x[, i] > unique(x[, i][max]))
              
              lowest_level_choise <- matrix(which(x[, i] == unique(maxN(x[, i]))))
              
              if (sum(revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]]) == 0) {
                lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise[, 1]), normalization_value, replace = F)))
                y <- t(Matrix(rep(0, length(x[, 1]))))
                y[lowest_level] = 1
                
              }
              else {
                sub_int <- matrix(subset(lowest_level_choise, (revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]] >= 1) == TRUE))
                sub_int_bind <- matrix(cbind(sub_int,  revtrunc_fixed_factor_1[, i][sub_int[,  1]]), ncol = 2)
                
                colnames(sub_int_bind) <- c("V1", "V2")
                
                sub_int_bind_ordered <- matrix(sub_int_bind[order(sub_int_bind[, "V2"], decreasing = TRUE), ], ncol = 2)
                
                colnames(sub_int_bind_ordered) <- c("V1", "V2")
                
                sub_int_bind_ordered_V1 <- sub_int_bind_ordered[, "V1"]
                sub_int_bind_ordered_V2 <- sub_int_bind_ordered[, "V2"]
                if ((length(unique(sub_int_bind_ordered_V2)) == 
                     1 & length(sub_int_bind_ordered_V2) > 
                     as.vector(normalization_value))) {
                  lowest_level <- as.numeric(as.vector(sample(as.factor(sub_int_bind_ordered_V1), normalization_value, replace = F)))
                  y <- t(Matrix(rep(0, length(x[, 1]))))
                  y[lowest_level] = 1
                  
                }
                else {
                  if (length(sub_int_bind_ordered_V1) > 
                      normalization_value) {
                    maxN_1 <- function(x, N = normalization_value) {
                      len <- length(x)
                      if (N > len) {
                        warning("N greater than length(x).  Setting N=length(x)")
                        N <- length(x)
                      }
                      sort(x, partial = len - N + 1)[len - N + 1]
                    }
                    max_1 <- which(Matrix(sub_int_bind_ordered_V2)[, 1] == maxN_1(Matrix(sub_int_bind_ordered_V2)[, 1]), arr.ind = TRUE)
                    
                    normalization_value_1 <- normalization_value - sum(Matrix(sub_int_bind_ordered_V2)[,  1] > unique(Matrix(sub_int_bind_ordered_V2)[, 1][max_1]))
                    
                    lowest_level_choise_1 <- Matrix(which(Matrix(sub_int_bind_ordered_V2)[, 1] == unique(maxN_1(Matrix(sub_int_bind_ordered_V2)[, 1]))))
                    
                    lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise_1[, 1]), normalization_value_1, replace = F)))
                    lowest_level <- sub_int_bind_ordered_V1[lowest_level]
                    
                    lowest_level_1 <- sub_int_bind_ordered_V1[(Matrix(sub_int_bind_ordered_V2)[, 1] > unique(Matrix(sub_int_bind_ordered_V2)[, 1][max_1]))]
                    lowest_level <- c(lowest_level_1, lowest_level)
                    y <- t(Matrix(rep(0, length(x[, 1]))))
                    y[lowest_level] = 1
                  }
                  else {
                    if (length(sub_int_bind_ordered_V1) < 
                        normalization_value) {
                      sub_int_zeros <- Matrix(subset(lowest_level_choise, (revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]] < 1) == TRUE))
                      length(t(sub_int_zeros))
                      lowest_level_2 <- as.numeric(as.vector(sample(as.factor(sub_int_zeros[, 1]), (normalization_value - length(sub_int_bind_ordered_V1)), replace = F)))
                      
                      lowest_level_3 <- c(sub_int_bind_ordered_V1, lowest_level_2)
                      y <- t(Matrix(rep(0, length(x[, 1]))))
                      y[lowest_level_3] = 1
                      
                    }
                    else {
                      y <- t(Matrix(rep(0, length(x[, 1]))))
                      y[sub_int_bind_ordered_V1] = 1
                      
                    }
                  }
                }
              }
              SRS <- t(revtrunc_fixed_factor_1[, i] + ceiling(x[, i] > unique(x[, i][max])) + y)
              
              assign(paste(colnames(data)[i], sep = ""), SRS)
              fixed_factor_1 <- get(colnames(data)[i])
              colnames(fixed_factor_1)[i] <- colnames(data)[i]
            }
          }
          else {
            if (diff_counts[i] == 0) {
              fixed_factor <- revtrunc_fixed_factor_1[, i]
              assign(paste(colnames(data)[i], sep = ""), fixed_factor)
              fixed_factor_1 <- cbind(fixed_factor_1, fixed_factor)
              colnames(fixed_factor_1)[{
                counter = counter + 1
              }] <- colnames(data)[i]
            }
            else {
              maxN <- function(x, N = diff_counts[i]) {
                len <- length(x)
                if (N > len) {
                  warning("N greater than length(x).  Setting N=length(x)")
                  N <- length(x)
                }
                sort(x, partial = len - N + 1)[len - N + 1]
              }
              max <- which(x[, i] == maxN(x[, i]), arr.ind = TRUE)
              
              normalization_value <- diff_counts[i] - sum(x[, i] > unique(x[, i][max]))
              
              lowest_level_choise <- matrix(which(x[, i] == unique(maxN(x[, i]))))
              
              length(t(lowest_level_choise))
              if (sum(revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]]) == 0) {
                lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise[, 1]), normalization_value, replace = F)))
                y <- t(Matrix(rep(0, length(x[, 1]))))
                y[lowest_level] = 1
                
              }
              else {
                sub_int <- matrix(subset(lowest_level_choise, (revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]] >= 1) == TRUE))
                sub_int_bind <- Matrix(cbind(sub_int, revtrunc_fixed_factor_1[, i][sub_int[, 1]]))
                colnames(sub_int_bind) <- c("V1", "V2")
                sub_int_bind_ordered <- matrix(sub_int_bind[order(sub_int_bind[, "V2"], decreasing = TRUE), ], ncol = 2)
                colnames(sub_int_bind_ordered) <- c("V1", "V2")
                
                sub_int_bind_ordered_V1 <- sub_int_bind_ordered[, "V1"]
                sub_int_bind_ordered_V2 <- sub_int_bind_ordered[, "V2"]
                if ((length(unique(sub_int_bind_ordered_V2)) == 1 & length(sub_int_bind_ordered_V2) > as.vector(normalization_value))) {
                  lowest_level <- as.numeric(as.vector(sample(as.factor(sub_int_bind_ordered_V1), normalization_value, replace = F)))
                  y <- t(Matrix(rep(0, length(x[, 1]))))
                  y[lowest_level] = 1
                  
                }
                else {
                  if (length(sub_int_bind_ordered_V1) > normalization_value) {
                    maxN_1 <- function(x, N = normalization_value) {
                      len <- length(x)
                      if (N > len) {
                        warning("N greater than length(x).  Setting N=length(x)")
                        N <- length(x)
                      }
                      sort(x, partial = len - N + 1)[len - N + 1]
                    }
                    max_1 <- which(Matrix(sub_int_bind_ordered_V2)[, 1] == maxN_1(Matrix(sub_int_bind_ordered_V2)[, 1]), arr.ind = TRUE)
                    
                    normalization_value_1 <- normalization_value - sum(Matrix(sub_int_bind_ordered_V2)[, 1] > unique(Matrix(sub_int_bind_ordered_V2)[, 1][max_1]))
                    
                    lowest_level_choise_1 <- Matrix(which(Matrix(sub_int_bind_ordered_V2)[, 1] == unique(maxN_1(Matrix(sub_int_bind_ordered_V2)[, 1]))))
                    
                    lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise_1[, 1]), normalization_value_1, replace = F)))
                    lowest_level <- sub_int_bind_ordered_V1[lowest_level]
                    
                    lowest_level_1 <- sub_int_bind_ordered_V1[(Matrix(sub_int_bind_ordered_V2)[, 1] > unique(Matrix(sub_int_bind_ordered_V2)[, 1][max_1]))]
                    lowest_level <- c(lowest_level_1, lowest_level)
                    y <- t(Matrix(rep(0, length(x[, 1]))))
                    y[lowest_level] = 1
                  }
                  else {
                    if (length(sub_int_bind_ordered_V1) < normalization_value) {
                      sub_int_zeros <- Matrix(subset(lowest_level_choise, (revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]] < 1) == TRUE))
                      length(t(sub_int_zeros))
                      lowest_level_2 <- as.numeric(as.vector(sample(as.factor(sub_int_zeros[, 1]), (normalization_value - length(sub_int_bind_ordered_V1)), replace = F)))
                      lowest_level_3 <- c(sub_int_bind_ordered_V1, lowest_level_2)
                      y <- t(Matrix(rep(0, length(x[, 1]))))
                      y[lowest_level_3] = 1
                    }
                    else {
                      y <- t(Matrix(rep(0, length(x[, 1]))))
                      y[sub_int_bind_ordered_V1] = 1
                    }
                  }
                }
              }
              SRS <- t(revtrunc_fixed_factor_1[, i] + ceiling(x[, i] > unique(x[, i][max])) + y)
              assign(paste(colnames(data)[i], sep = ""), SRS)
              fixed_factor_1 <- cbind(fixed_factor_1, SRS)
              colnames(fixed_factor_1)[{
                counter = counter + 1
              }] <- colnames(data)[i]
            }
          }
        }
        SRS_output <- fixed_factor_1
        gc()
        return(SRS_output)
      }
    }
  }
}

############
## rrarefy2
############
##
## Description: reimplementation of VEGAN's rrarefy() function which allows the use of sparse matrices.
##
rrarefy2 <- function (x, sample) 
{
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function is meaningful only for integers (counts)")
  if (!is.integer(x)) 
    x <- round(x)
  if (ncol(x) == 1) 
    x <- t(x)
  if (length(sample) > 1 && length(sample) != nrow(x)) 
    stop(gettextf("length of 'sample' and number of rows of 'x' do not match"))
  sample <- rep(sample, length = nrow(x))
  for (i in 1:nrow(x)) {
    x[i, ] <- .Call(do_rrarefy, x[i, ], sample[i])
  }
  return(Matrix(x,sparse = T))
}
environment(rrarefy2) <- environment(rrarefy)
