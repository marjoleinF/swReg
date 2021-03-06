% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swReg.R
\name{swReg}
\alias{swReg}
\title{Perform forward stagewise regression}
\usage{
swReg(X, Y, stepsize = 0.1, standardizeY = TRUE)
}
\arguments{
\item{X}{Matrix of predictor variables.}

\item{Y}{Outcome variable (vector).}

\item{stepsize}{numeric. Value with which the (un)standardized coefficients 
are updated in each iteration (a.k.a. epsilon).}

\item{standardizeY}{logical. Should the response be standardized prior to
application of the algorithm? If \code{TRUE}, stepsize can be interpreted
as a correlation. If \code{sd(Y) > 1}, setting \code{standardizeY = FALSE} 
will increase the number of steps needed for convergence.}
}
\value{
The function returns a list with the following elements:
unstandardized.coef = bUnstand, 
standardized.coef = b,
iteration = iteration, 
stepsize = stepsize (epsilon), 
coef.path = dataframe with standardized coefficient values for each predictor variable, at each stage or iteration of the algorithm
data = list with original data (X andY)
}
\description{
\code{swReg} performs forward stagewise regression.
}
\details{
The function performs incremental forward stagewise regression, as 
described in Hastie, Tisbhirani & Friedman (2009). Hastie, T., Tibshirani, 
R., & Friedman, J. (2009). The Elements of Statistical Learning.
}
\examples{
## Example using Boston housing data:
 library(MASS)
 X <- as.matrix(Boston[,-14])
 Y <- Boston$medv
   
 ## Run forward stagewise regression:
 swReg1 <- swReg(X, Y, stepsize = .01)
 swReg2 <- swReg(X, Y, stepsize = .01, standardizeY = FALSE)
   
}
