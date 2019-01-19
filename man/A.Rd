\name{A}
\alias{A}
\title{A
}
\description{
Calculates probability of superiority (A), its standard error, and a confidence interval.
}
\usage{
A(data, design = 1, statistic = 1, weights = FALSE,
w = 0, w1 = 0, w2 = 0, increase = FALSE, ref = 1, r = 0,
n.bootstrap = 1999, conf.level = .95, ci.method = 1, seed = 1)
}
\arguments{
  \item{data}{For a between subjects design, a matrix of cases (rows) by scores (column 1) and group codes (column 2). For a within subjects design, a matrix of scores with each sample in its own column (matrix).}
  \item{design}{Design of experiment (scalar, default = 1 (for between subjects design), user can also call 2 (for within subjects design)).}
  \item{statistic}{Statistic to be calculated (scalar, default = 1 (A), user can also call 2 (A.AAD), 3 (A.AAPD), 4 (A.IK), or 5 (A.Ord)).}
  \item{weights}{Whether to assign weights to cases (default = FALSE); if set to TRUE, data contains case weights in final column.}
  \item{w}{Weights for cases (vector; default = 0).}
  \item{w1}{Weights for cases in group 1 (vector; default = 0).}
  \item{w2}{Weights for cases in group 2 (vector; default = 0).}
  \item{increase}{Set to TRUE if scores are predicted to increase with group codes (default = FALSE).}
  \item{ref}{Reference group (to compare to all others) (scalar, default = 1).}
  \item{r}{Vector of proportions (vector, default = 0, represents equal proportions).}
  \item{n.bootstrap}{Number of bootstrap samples (scalar, default = 1999).}
  \item{conf.level}{Confidence level (default = .95).}
  \item{ci.method}{Method used to construct confidence interval (scalar, default = 1 (for BCA), user can also call 2 (for percentile)).}
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
data <- cbind(c(x1, x2, x3), c(rep(1, 25), rep(2, 25), rep(3, 25)))
A(data, 1, 2)
}

