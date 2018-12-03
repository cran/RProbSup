\name{CalcIK2}
\alias{CalcIK2}
\title{CalcIK2
}
\description{Calculates the A statistic while singling out one group for two or more correlated samples. Note: This function is not meant to be called by the user, but it is called by IK2.
}
\usage{
CalcIK2(y, ref = 1, weights = FALSE)
}
\arguments{
  \item{y}{Matrix of cases (rows) by scores (column 1) and group codes (column 2) (matrix).}
  \item{ref}{Reference group (to compare to all others) (scalar, default = 1).}
  \item{weights}{Weight of each case. Set to TRUE to weight cases; if so, column 3 contains case weights (default = FALSE).}
}
\value{
 \item{a}{The A statistic.}
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
CalcIK2(y)
}
