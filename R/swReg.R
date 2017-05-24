#' Perform forward stagewise regression
#' 
#' \code{swReg} performs forward stagewise regression.
#' 
#' @param X Matrix of predictor variables.
#' @param Y Outcome variable (vector).
#' @param stepsize Value with which the standardized coefficients are updated in each stage (iteration).
#' @param threshold Convergence criterion: When none of the predictors has a correlation greater than or equal to the threshold, the algorithm converges.
#' @details The function performs forward stagewise regression, as described in Efron, B., Hastie, T., Johnstone, I., & Tibshirani, R. (2004). Least angle regression. The Annals of Statistics, 32(2), 407-499.
#' @return The function returns a list with the following elements:
#' unstandardized.coef = bUnstand, 
#' standardized.coef = b,
#' iteration = iteration, 
#' stepsize = stepsize, 
#' threshold = threshold, 
#' coef.path = dataframe with standardized coefficient values for each predictor variable, at each stage or iteration of the algorithm
#' data = list with original data (X andY)
#' @export
#' @examples ## Example using Boston housing data:
#'  library(MASS)
#'  X <- as.matrix(Boston[,-14])
#'  Y <- Boston$medv
#'    
#'  ## Run forward stagewise regression:
#'  swReg.2 <- swReg(X, Y, st = .2, thr = .2)
#'  swReg.1 <- swReg(X, Y, st = .1, thr = .1)
#'  swReg.01 <- swReg(X, Y, st = .01, thr = .01)
#'  swReg.001 <- swReg(X, Y, st = .001, thr = .001)
#'    
#'  ## Plot coefficient paths:
#'  par(mfrow=c(2,2))
#'  plot.swReg(swReg.2)
#'  plot.swReg(swReg.1)
#'  plot.swReg(swReg.01)
#'  plot.swReg(swReg.001)
#'  
#'  
swReg <- function (X, Y, stepsize = 0.1, threshold = 0.1) 
{
  sX <- scale(X)
  sY <- scale(Y)
  b <- vector(length = ncol(X))
  rY <- sY
  sgn <- 0
  lastsgn <- 0
  index <- 0
  lastindex <- 0
  iteration <- 0
  path <- list()
  abs(max(t(rY) %*% sX)) > threshold & (index != lastindex | sgn == lastsgn)
  # if there are >1 variables with the maximum correlation, this gives warnings:  
  while(abs(max(t(rY) %*% sX)) > threshold & (index != lastindex | sgn == lastsgn)) {
    iteration <- iteration + 1
    lastindex <- index
    lastsgn <- sgn
    index <- which(abs(t(rY) %*% sX) == max(abs(t(rY) %*% sX)))
    sgn <- sign((t(rY) %*% sX)[index])
    b[index] <- b[index] + sgn*stepsize
    rY <- sY - sX %*% b
    path[[iteration]] <- b
  }
  bUnstand <- (b/attr(sX, "scaled:scale"))*attr(sY, "scaled:scale") 
  intercept <- attr(sY, "scaled:center") - bUnstand %*% attr(sX, "scaled:center")
  bUnstand <- c(intercept, bUnstand)
  names(bUnstand) <- c("(Intercept)", colnames(X))
  
  return(list(unstandardized.coef = bUnstand, standardized.coef = b,
              iteration = iteration, stepsize = stepsize, threshold = threshold, 
              coef.path = data.frame(matrix(
                unlist(path), ncol = ncol(X), nrow = iteration, byrow = TRUE, 
                dimnames = list(1:iteration, colnames(X)))),
              data = list(X = X, Y = Y))
  )  
}

#' Plot forward stagewise regression path
#' 
#' \code{plot.swReg} plots non-zero coefficients for each iterations of forward stagewise regression.
#' 
#' @param x an R object resulting from application of \code{swReg}
#' @param legend logical. Should a legend be printed? Defaults to TRUE
#' @param ... not used
#' @return A plot with iteration numbers on the x-axis and coefficient values on the y-axis,
#' @export
#' @method plot swReg
plot.swReg <- function(x, legend = TRUE, ...) {
  plot(x$coef.path[,1], type = "l", 
       ylim = c(min(x$coef.path), max(x$coef.path)), 
       ylab = "coefficient", xlab = "iteration", 
       main = "Standardized coefficient paths", 
       sub = paste("stepsize =", x$stepsize))
  for (i in 2:ncol(x$data$X)) {lines(x$coef.path[,i], col = i)}
  if(legend) {
    legend("topleft", legend = colnames(x$data$X), cex = .5, 
           lty = 1, col = 1:ncol(x$coef.path), y.intersp = .25, 
           x.intersp = .25, bty = "n", seg.len = .5)
  }
}


#' Get predictions for a given swReg model
#' 
#' \code{predict.swReg} returns predictions for a given swReg model
#' 
#' @param object an R object resulting from application of \code{swReg}
#' @param newdata a dataframe or matrix with predictor variable values.
#' @param ... currently not used.
#' @return A vector of predictions for \code{newdata}.
#' @export
predict.swReg <- function(object, newdata = object$data$X, ...) {
  cbind(rep(1, times = nrow(newdata)), as.matrix(newdata)) %*% object$unstandardized.coef
}
                            
                            
#' Get k-fold cross validated predictions for a given swReg model
#' 
#' \code{swReg.xval} returns cross validated predictions for the swReg model specified.
#' 
#' @param object an R object resulting from application of \code{swReg}
#' @param k integer. Number of folds to be used.
#' @return A list with iteration numbers on the x-axis and coefficient values on the y-axis,
#' @export
swReg.xval <- function(object, k = 10) {
  # get observation ids for in folds:
  ids <- peperr::resample.indices(n = nrow(object$data$X), sample.n = 10, method = "cv")
  # create a list for gathering results:
  models <- list()
  xvaldatasets <- list()
  # for each fold k:
  for (i in 1:k){
    # make datasets of training and test observations:
    xvaldatasets[[i]] <- list(train = list(X = object$data$X[ids$sample.index[[i]],],
                                           Y = object$data$Y[ids$sample.index[[i]]]),
                              test = list(X = object$data$X[ids$not.in.sample[[i]],],
                                          Y = object$data$Y[ids$not.in.sample[[i]]]))
    # fit models on training observations
    models[[i]] <- swReg(xvaldatasets[[i]]$train$X, xvaldatasets[[i]]$train$Y, 
                         stepsize = object$stepsize, threshold = object$threshold)
    # make predictions for test observations:
    xvaldatasets[[i]]$test$xvalpredY <- predict.swReg(models[[i]], newdata = xvaldatasets[[i]]$test$X)
  }
  # combine X, Y, cv predictions and foldnumber
  dataset <- data.frame(object$data$X, fold = NA, Y = object$data$Y, Ycvpred = NA)
  coefficients <- list()
  for(i in 1:k){
    dataset[,"fold"][ids$not.in.sample[[i]]] <- i
    dataset[,"Ycvpred"][ids$not.in.sample[[i]]] <- xvaldatasets[[i]]$test$xvalpredY 
  }
  dataset <- cbind(dataset, error =  dataset[,"Ycvpred"] - dataset[,"Y"])
  return(dataset = dataset)
}