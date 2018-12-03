\name{CalcAAD1}
\alias{CalcAAD1}
\title{CalcAAD1
}
\description{Calculates the A statistic for the average absolute deviation for two or more groups. Note: This function is not meant to be called by the user, but it is called by AAD1.
}
\usage{
CalcAAD1(y, r = 0, weights = FALSE)
}
\arguments{
  \item{y}{Matrix of cases (rows) by scores (column 1) and group codes (column 2) (matrix).}
  \item{r}{Vector of proportions (default = 0, represents equal proportions) (vector).}
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
y <- cbind(c(x1, x2, x3), c(rep(1, 25), rep(2, 25), rep(3, 25)))
CalcAAD1(y)
}
