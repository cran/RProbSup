\name{Ord2}
\alias{Ord2}
\title{Ord2
}
\description{Calculates the confidence interval for the ordinal comparison of the A statistic for two or more correlated samples.
}
\usage{
Ord2(y, weights = FALSE, increase = FALSE, n.bootstrap = 1999,
conf.level = .95, method = "BCA", seed = 1)
}
\arguments{
  \item{y}{Matrix of cases (rows) by scores (column 1) and group codes (column 2) (matrix).}
  \item{weights}{Weight of each case. Set to TRUE to weight cases; if so, column 3 contains case weights (default = FALSE).}
  \item{increase}{Set to TRUE if scores are predicted to increase with group codes (default = FALSE).}
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
x1 <- rnorm(25)
x2 <- x1 - rnorm(25, mean = 1)
x3 <- x2 - rnorm(25, mean = 1)
y <- cbind(x1, x2, x3)
Ord2(y)
}
