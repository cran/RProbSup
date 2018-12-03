\name{CalcAAPD2}
\alias{CalcAAPD2}
\title{CalcAAPD2
}
\description{Calculates the A statistic for the average absolute paired deviation for two or more correlated samples. Note: This function is not meant to be called by the user, but it is called by AAPD2.
}
\usage{
CalcAAPD2(y, weights = FALSE)
}
\arguments{
  \item{y}{Matrix of cases (rows) by scores (column 1) and group codes (column 2) (matrix).}
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
AAPD2(y)
}
