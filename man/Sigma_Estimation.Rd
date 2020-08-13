\name{Sigma_Estimation}
\alias{Sigma_Estimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the correlation matrix
}
\description{
This function estimates the correlation matrix of the rows of a design matrix in linear models
in the case where the correlation matrix has a block-wise structure.
}
\usage{
Sigma_Estimation(X)
}
\arguments{
  \item{X}{
Design matrix of the linear model considered.
}
}

\value{
\item{mat }{Estimation of the correlation matrix of the rows of X}
\item{alpha }{Estimation of the coefficients of the blocks of the correlation matrix}
\item{group_act }{Indices of the active or non active variables}
%% ...
}

\author{
Wencan Zhu [aut, cre], Celine Levy-Leduc [ctb], Nils Ternes [ctb]
}

\references{
W. Zhu, C. Levy-Leduc, N. Ternes. "A variable selection approach for highly correlated predictors in 
high-dimensional genomic data". arXiv:2007.10768.
}

\examples{
data(X)
Sigma_Estimation(X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~estimation of the correlation matrix }% use one of  RShowDoc("KEYWORDS")



