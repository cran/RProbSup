\name{AAD1}
\alias{AAD1}
\title{AAD1
}
\description{Calculates the confidence interval for the A statistic for the average absolute deviation for two or more groups.
}
\usage{
AAD1(y, r = 0, weights = FALSE, n.bootstrap = 1999,
conf.level = .95, method = "BCA", seed = 1)
}
\arguments{
  \item{y}{Matrix of cases (rows) by scores (column 1) and group codes (column 2) (matrix).}
  \item{r}{Vector of proportions (default = 0, represents equal proportions) (vector).}
  \item{weights}{Weight of each case. Set to TRUE to weight cases; if so, column 3 contains case weights (default = FALSE).}
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
y <- cbind(c(x1, x2, x3), c(rep(1, 25), rep(2, 25), rep(3, 25)))
AAD1(y)
}
