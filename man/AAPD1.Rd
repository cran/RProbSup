\name{AAPD1}
\alias{AAPD1}
\title{AAPD1
}
\description{Calculates the confidence interval for the A statistic for the average absolute paired deviation for two or more groups.
}
\usage{
AAPD1(y, weights = FALSE, n.bootstrap = 1999,
conf.level = .95, ci.method = 1, seed = 1)
}
\arguments{
  \item{y}{Matrix of cases (rows) by scores (column 1) and group codes (column 2) (matrix).}
  \item{weights}{Weight of each case. Set to TRUE to weight cases; if so, column 3 contains case weights (default = FALSE).}
  \item{n.bootstrap}{Number of bootstrap samples (scalar, default = 1999).}
  \item{conf.level}{Confidence level (scalar, default = .95).}
  \item{ci.method}{Method used to construct confidence interval (scalar, default = 1 (for BCA), user can also call 2 (for percentile).}
  \item{seed}{Random number seed (scalar, default = 1).}
}
\value{
A vector containing the A statistic, its estimated standard error, and the upper and lower bounds of the confidence interval.
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
AAPD1(y)
}
