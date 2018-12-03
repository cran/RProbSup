\name{CalcA1}
\alias{CalcA1}
\title{CalcA1
}
\description{Calculates the A statistic for 2 groups.
}
\usage{
CalcA1(y1, y2, weights = FALSE, w1 = 0, w2 = 0)
}
\arguments{
  \item{y1}{Scores for group 1 (vector).}
  \item{y2}{Scores for group 2 (vector).}
  \item{weights}{Whether to weight cases (default = FALSE).}
  \item{w1}{Weights for cases in group 1 (optional) (vector, default is 0).}
  \item{w2}{Weights for cases in group 2 (optional) (vector, default is 0).}
}
\value{
  \item{a}{The A statistic.}
}
\references{Ruscio (2008) & Ruscio and Mullen (2012) & Ruscio and Gera (2013)
}
\author{John Ruscio
}
\examples{
#Example used in Ruscio and Mullen (2012)
y1 <- c(6, 7, 8, 7, 9, 6, 5, 4, 7, 8, 7, 6, 9, 5, 4)
y2 <- c(4, 3, 5, 3, 6, 2, 2, 1, 6, 7, 4, 3, 2, 4, 3)
CalcA1(y1, y2)
}
