\name{RemoveMissing}
\alias{RemoveMissing}
\title{RemoveMissing
}
\description{
Checks for missing data and performs listwise deletion if any is detected.
}
\usage{
RemoveMissing(data)
}
\arguments{
  \item{data}{For a between subjects design, a matrix of cases (rows) by scores (column 1) and group codes (column 2). For a within subjects design, a matrix of scores with each sample in its own column (matrix).}
}
\value{
Data matrix with any missing data removed using listwise deletion of cases.
}
\references{Ruscio (2008) & Ruscio and Mullen (2012) & Ruscio and Gera (2013)
}
\author{John Ruscio
}
\examples{
x1 <- c(rnorm(25), NA)
x2 <- x1 - rnorm(26, mean = 1)
x3 <- x2 - rnorm(26, mean = 1)
data <- cbind(c(x1, x2, x3), c(rep(1, 26), rep(2, 26), rep(3, 26)))
A(data, 1, 2)
}

