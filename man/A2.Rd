\name{A2}
\alias{A2}
\title{A2
}
\description{
Calculates the standard error and constructs a confidence interval for the A statistic for two correlated samples using bootstrap methods.
}
\usage{
A2(y1, y2, weights = FALSE, w = 0, n.bootstrap = 1999,
conf.level = .95, method = "BCA", seed = 1)
}
\arguments{
  \item{y1}{Scores for group 1 (vector).}
  \item{y2}{Scores for group 2 (vector).}
  \item{weights}{Whether to weight cases (default = FALSE).}
  \item{w}{Weights for cases in group 1 (optional) (vector, default is 0).}
  \item{n.bootstrap}{Number of bootstrap samples (scalar, default = 1999).}
  \item{conf.level}{Confidence level (scalar, default = .95).}
  \item{method}{Method used to construct confidence interval (character, default = "BCA", user can also call "percentile").}
  \item{seed}{Random number seed (scalar, default = 1).}
}
\value{
Nothing; displays the A statistic, its estimated standard error, and the confidence interval.
}
\references{Ruscio (2008) & Ruscio and Mullen (2012) & Ruscio and Gera (2013)
}
\author{John Ruscio
}
\examples{
y1 <- c(6, 7, 8, 7, 9, 6, 5, 4, 7, 8, 7, 6, 9, 5, 4)
y2 <- c(7, 5, 6, 7, 6, 4, 3, 5, 4, 5, 4, 5, 7, 4, 5)
A2(y1, y2)
}

