\name{Whitening_Lasso}
\alias{Whitening_Lasso}
\title{
Whitening Lasso
}
\description{
The function implements the approach described in the paper Zhu et al. (2020) given in the references.
}
\usage{
Whitening_Lasso(X, Y, Sigma, gamma = 0.95, maxsteps = 2000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Design matrix of the linear model.
}
  \item{Y}{
Response variable of the linear model.
}
  \item{Sigma}{
Correlation matrix of the rows of the design matrix. If not specified, the function \code{\link{Sigma_Estimation}} will be used to estimate this matrix.}
  \item{gamma}{
Parameter \eqn{\gamma} defined in the paper Zhu et al. (2020) given in the references. Its default value is 0.95.}
  \item{maxsteps}{
Integer specifying the maximum number of steps for the generalized Lasso algorithm. Its default value is 2000.}
}

\value{
Returns a list with the following components
\item{lambda }{different values of the parameter \eqn{\lambda} considered.}
\item{beta }{matrix of the estimations of \eqn{\beta} for all the \eqn{\lambda} considered.}
\item{trans_mat }{tranformation matrix which is the inverse of the square root of the estimation of the correlation matrix of the rows of the design matrix X.}
\item{beta.min }{estimation of \eqn{\beta} which minimize the MSE.}
\item{mse }{MSE for all the \eqn{\lambda} considered.}
}

\references{
W. Zhu, C. Levy-Leduc, N. Ternes. "A variable selection approach for highly correlated predictors in 
high-dimensional genomic data". arXiv:2007.10768.
}

\author{
Wencan Zhu [aut, cre], Celine Levy-Leduc [ctb], Nils Ternes [ctb]
}

\examples{
data(X)
data(Y)
Sigma_est <- Sigma_Estimation(X)$mat
res_wlasso <- Whitening_Lasso(X=X, Y=Y, Sigma=Sigma_est, maxsteps=100, gamma=0.9)
beta_est <- res_wlasso$beta
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Whitening Lasso }% use one of  RShowDoc("KEYWORDS")

