#' Perform forward stagewise regression
#' 
#' \code{swReg} performs forward stagewise regression.
#' 
#' @param X Matrix of predictor variables.
#' @param Y Outcome variable (vector).
#' @param stepsize numeric. Value with which the (un)standardized coefficients 
#' are updated in each iteration (a.k.a. epsilon).
#' @param standardizeY logical. Should the response be standardized prior to
#' application of the algorithm? If \code{TRUE}, stepsize can be interpreted
#' as a correlation. If \code{sd(Y) > 1}, setting \code{standardizeY = FALSE} 
#' will increase the number of steps needed for convergence.  
#' @details The function performs incremental forward stagewise regression, as 
#' described in Hastie, Tisbhirani & Friedman (2009). Hastie, T., Tibshirani, 
#' R., & Friedman, J. (2009). The Elements of Statistical Learning.
#' @return The function returns a list with the following elements:
#' unstandardized.coef = bUnstand, 
#' standardized.coef = b,
#' iteration = iteration, 
#' stepsize = stepsize (epsilon), 
#' coef.path = dataframe with standardized coefficient values for each predictor variable, at each stage or iteration of the algorithm
#' data = list with original data (X andY)
#' @export
#' @examples ## Example using Boston housing data:
#'  library(MASS)
#'  X <- as.matrix(Boston[,-14])
#'  Y <- Boston$medv
#'    
#'  ## Run forward stagewise regression:
#'  swReg1 <- swReg(X, Y, stepsize = .01)
#'  swReg2 <- swReg(X, Y, stepsize = .01, standardizeY = FALSE)
#'    
swReg <- function (X, Y, stepsize = 0.1, standardizeY = TRUE) 
{
  # ESL, print 10, page 86:
  #
  # Algorithm 3.4 Incremental Forward Stagewise Regression-FS_epsilon
  # Step 1. Start with the residual r equal to y and beta_1, beta_2, ..., beta_p = 0. All the
  # predictors are standardized to have mean zero and unit norm.
  n <- nrow(X)
  p <- ncol(X)
  sX <- scale(X)
  b <- rep(0, times = p)
  if(standardizeY) {r <- scale(Y)} else {r <- Y}
  sgn <- index <- iteration <- 0
  path <- list()
  cur_max_cor <- stepsize
  # Step 4. Repeat Steps 2 and 3 until the residuals are uncorrelated with all the predictors:
  while(cur_max_cor >= stepsize) {
    iteration <- iteration + 1
    lastindex <- index
    lastsgn <- sgn
    # Step 2. Find the predictor x_j most correlated with r:
    cur_cors <- (t(r - mean(r)) %*% sX)/n
    cur_max_cor <- max(abs(cur_cors))
    j <- which(abs(cur_cors) == cur_max_cor)[1]
    # Step 3. Update beta_j  <-  beta_j + delta_j , where delta_j = epsilon x sign[<x_j, r>] and epsilon > 0 is a small
    # step size, and set r <- r ??? delta_j x_j.
    sgn <- sign(cur_cors[j])
    b[j] <- b[j] + sgn*stepsize
    r <- r - sgn*stepsize*sX[,j]
    path[[iteration]] <- b
  }
  if(standardizeY) {
    b <- b*sd(Y)
  }
  bUnstand <- b/attr(sX, "scaled:scale") 
  intercept <- mean(Y) - attr(sX, "scaled:center") %*% bUnstand
  bUnstand <- c(intercept, bUnstand)
  names(bUnstand) <- c("(Intercept)", colnames(X))
  
  result<- list(unstandardized.coef = bUnstand, standardized.coef = b,
              iteration = iteration, stepsize = stepsize, standardizeY = standardizeY, 
              coef.path = data.frame(matrix(
                unlist(path), ncol = p, nrow = iteration, byrow = TRUE, 
                dimnames = list(1:iteration, colnames(X)))),
              data = list(X = X, Y = Y))
  class(result) <- "swReg"
  return(result)
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
#' @examples ## Example using Boston Housing data:
#'  library(MASS)
#'  X <- as.matrix(Boston[,-14])
#'  Y <- Boston$medv
#'  swReg1 <- swReg(X, Y, stepsize = .01)
#'  swReg2 <- swReg(X, Y, stepsize = .01, standardizeY = FALSE)
#'    
#'  ## Plot coefficient paths:
#'  par(mfrow=c(1,2))
#'  plot.swReg(swReg1)
#'  plot.swReg(swReg2)
#'  
plot.swReg <- function(x, legend = TRUE, ...) {
  if(x$standardizeY) {
    main = "Completely standardized coeffcient paths"
  } else {
    main = "Coefficient paths (standardized X)"
  }
  plot(x$coef.path[,1], type = "l", 
       ylim = c(min(x$coef.path), max(x$coef.path)), 
       ylab = "coefficient", xlab = "iteration", 
       main = main, 
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
#' @method predict swReg
#' @examples ## Example using Boston Housing data:
#'  library(MASS)
#'  X <- as.matrix(Boston[1:400,-14])
#'  Y <- Boston$medv[1:400]
#'  swReg1 <- swReg(X, Y, stepsize = .01)
#'  swReg2 <- swReg(X, Y, stepsize = .01, standardizeY = FALSE)
#'  plot.swReg(swReg1, newdata = Boston[401:502,-14])
#'  plot.swReg(swReg2, newdata = Boston[401:502,-14])
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
#' @examples ## Example using Boston Housing data:
#'  library(MASS)
#'  X <- as.matrix(Boston[,-14])
#'  Y <- Boston$medv
#'  swReg2 <- swReg(X, Y, stepsize = .01)
#'  xval2 <- swReg.xval(swReg2)
#'  mean(xval2$error^2) # MSE
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
                         stepsize = object$stepsize, standardizeY = object$standardizeY)
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