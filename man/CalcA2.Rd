\name{CalcA2}
\alias{CalcA2}
\title{CalcA2
}
\description{Calculates the A statistic for 2 correlated samples.
}
\usage{
CalcA2(y1, y2, weights = FALSE, w = 0)
}
\arguments{
  \item{y1}{Scores for variable 1 (vector).}
  \item{y2}{Scores for variable 2 (vector).}
  \item{weights}{Whether to weight cases (default = FALSE).}
  \item{w}{Weights (optional) (vector, default is 0).}
}
\value{
 \item{a}{The A statistic.}
}
\references{Ruscio (2008) & Ruscio and Mullen (2012) & Ruscio and Gera (2013)
}
\author{John Ruscio
}
\examples{
y1 <- c(6, 7, 8, 7, 9, 6, 5, 4, 7, 8, 7, 6, 9, 5, 4)
y2 <- c(7, 5, 6, 7, 6, 4, 3, 5, 4, 5, 4, 5, 7, 4, 5)
CalcA2(y1, y2)
}

