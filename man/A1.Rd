\name{A1}
\alias{A1}
\title{A1
}
\description{Calculates the standard error and constructs a confidence interval for the A statistic using bootstrap methods.
}
\usage{
A1(y1, y2, weights = FALSE, w1 = 0, w2 = 0, n.bootstrap = 1999,
conf.level = .95, method = "BCA", seed = 1)
}
\arguments{
  \item{y1}{Scores for group 1 (vector).}
  \item{y2}{Scores for group 2 (vector).}
  \item{weights}{Whether to weight cases (default = FALSE).}
  \item{w1}{Weights for cases in group 1 (optional) (vector, default is 0).}
  \item{w2}{Weights for cases in group 2 (optional) (vector, default is 0).}
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
#Example used in Ruscio and Mullen (2012)
y1 <- c(6, 7, 8, 7, 9, 6, 5, 4, 7, 8, 7, 6, 9, 5, 4)
y2 <- c(4, 3, 5, 3, 6, 2, 2, 1, 6, 7, 4, 3, 2, 4, 3)
A1(y1, y2)
}

