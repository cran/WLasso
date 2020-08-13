\name{top_thresh}
\alias{top_thresh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Thresholding to a given threshold of the smallest values}
\description{
This function keeps only the K largest values of the vector \code{sorted_vect} and sets the others to the smallest value among the K largest.
}
\usage{
top_thresh(x, thresh, sorted_vect)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector to threshold
}
  \item{thresh}{
threshold}
  \item{sorted_vect}{
vector x sorted in descending order}
}


\value{
This function returns the thresholded vector.
}

\author{
Wencan Zhu [aut, cre], Celine Levy-Leduc [ctb], Nils Ternes [ctb]
}

\references{
W. Zhu, C. Levy-Leduc, N. Ternes. "A variable selection approach for highly correlated predictors in 
high-dimensional genomic data". arXiv:2007.10768.
}

\examples{
x=sample(1:10,10)
sorted_vect=sort(x,decreasing=TRUE)
thresh=3
top_thresh(x,thresh,sorted_vect)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~thresholding }% use one of  RShowDoc("KEYWORDS")

