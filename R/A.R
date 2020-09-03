################################################################################
A <- function(data, design = 1, statistic = 1, weights = FALSE, w = 0,
              w1 = 0, w2 = 0, increase = FALSE, ref = 1, r = 0,
              n.bootstrap = 1999, conf.level = .95, ci.method = 1, seed = 1) {
# Calculates probability of superiority (A), its standard error, and a
# confidence interval.
#
# Args:
#   data        : For a between subjects design, a matrix of cases (rows) by
#                 scores (column 1) and group codes (column 2). For a within
#                 subjects design, a matrix of scores with each sample in its
#                 own column (matrix).
#   design      : Design of experiment (scalar, default = 1 (for between
#                 subjects design), user can also call 2 (for within subjects
#                 design)).
#   statistic   : Statistic to be calculated (scalar, default = 1 (A),
#                 user can also call 2 (A.AAD), 3 (A.AAPD), 4 (A.IK), or
#                 5 (A.Ord)).
#   weights     : Whether to assign weights to cases (default = FALSE); if
#                 set to TRUE, data contains case weights in final column.
#   increase    : Set to TRUE if scores are predicted to increase with group
#                 codes (default = FALSE).
#   ref         : Reference group (to compare to all others) (scalar,
#                 default = 1).
#   r           : Vector of proportions (vector, default = 0, represents equal
#                 proportions).
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for
#                 percentile)).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   Nothing; displays the A statistic, its estimated standard error, and the
#   confidence interval.
  set.seed(seed)
  data <- RemoveMissing(data)
  if ((design == 1) & (statistic == 1)) {
    x <- A1(data[(data[, 2] == 1), 1], data[(data[, 2] == 2), 1], weights,
            w1, w2, n.bootstrap, conf.level, ci.method, seed)
  }
  if ((design == 1) & (statistic == 2)) {
    x <- AAD1(data, r, weights, n.bootstrap, conf.level, ci.method, seed)
  }
  if ((design == 1) & (statistic == 3)) {
    x <- AAPD1(data, weights, n.bootstrap, conf.level, ci.method, seed)
  }
  if ((design == 1) & (statistic == 4)) {
    x <- IK1(data, ref, weights, n.bootstrap, conf.level, ci.method, seed)
  }
  if ((design == 1) & (statistic == 5)) {
    x <- Ord1(data, weights, increase, n.bootstrap, conf.level, ci.method, seed)
  }
  if ((design == 2) & (statistic == 1)) {
    x <- A2(data[, 1], data[, 2], weights, w, n.bootstrap, conf.level, ci.method, seed)
  }
  if ((design == 2) & (statistic == 2)) {
    x <- AAD2(data, r, weights, n.bootstrap, conf.level, ci.method, seed)
  }
  if ((design == 2) & (statistic == 3)) {
    x <- AAPD2(data, weights, n.bootstrap, conf.level, ci.method, seed)
  }
  if ((design == 2) & (statistic == 4)) {
    x <- IK2(data, ref, weights, n.bootstrap, conf.level, ci.method, seed)
  }
  if ((design == 2) & (statistic == 5)) {
    x <- Ord2(data, weights, increase, n.bootstrap, conf.level, ci.method, seed)
  }
  cat("     A: ", round(x[1], 3), "\n")
  cat("    SE: ", round(x[2], 3), "\n")
  cat(100 * conf.level, "% CI:  ", round(x[3], 3), " to ", round(x[4], 3),
    "\n", sep = "")
  if (ci.method == 1) {
    cat("        Constructed using BCA method with B = ", n.bootstrap,
        " bootstrap samples\n", sep = "")
  }
  if (ci.method == 2) {
    cat("        Constructed using percentile method with B = ", n.bootstrap,
        " bootstrap samples\n", sep = "")
  }
}

################################################################################
RemoveMissing <- function(data) {
#
#
# Args :
#   data: between- or within-subjects data set (matrix).
#
# Returns :
#   Data matrix with any missing data removed using listwise deletion of cases.
#
  n <- dim(data)[1]
  k <- dim(data)[2]
  complete <- rep(TRUE, n)
  for (i in 1:n)
    for (j in 1:k)
      if (is.na(data[i, j])) complete[i] <- FALSE
   data2 <- data[complete, ]
   n2 <- dim(data2)[1]
   cat("\n", n2 ," out of ", n ," cases (", round(100 * n2 / n, 2),
       "%) retained for analysis; ", n - n2," cases (",
       round(100 * (n - n2) / n, 2), "%) contained missing data.\n", sep = "")
   return(data2)
}

################################################################################
A1 <- function(y1, y2, weights = FALSE, w1 = 0, w2 = 0, n.bootstrap = 1999,
               conf.level = .95, ci.method = 1, seed = 1) {
# Calculates the standard error and constructs a confidence interval for the A
# statistic for two groups using bootstrap methods.
#
# Args :
#   y1          : Scores for group 1 (vector).
#   y2          : Scores for group 2 (vector).
#   weights     : Whether to assign weights to cases (default = FALSE).
#   w1          : Weights for cases in group 1 (vector; default = 0).
#   w2          : Weights for cases in group 2 (vector; default = 0).
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (scalar, default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for percentile).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   A vector containing the A statistic, its estimated standard error, and the
#   upper and lower bounds of the confidence interval.
#
  n1 <- length(y1)
  n2 <- length(y2)
  if (!weights) {
  	w1 <- rep(0, n1)
  	w2 <- rep(0, n2)
  }
  set.seed(seed)
  a.obs <- CalcA1(y1, y2, weights, w1, w2)
  alpha <- 1 - conf.level
  ci.lower <- ci.upper <- pi
  a.boot <- rep(0, n.bootstrap)
  for (i in 1:n.bootstrap) {
    cases1 <- sample(1:n1, n1, replace = TRUE)
    cases2 <- sample(1:n2, n2, replace = TRUE)
    a.boot[i] <- CalcA1(y1[cases1], y2[cases2], weights, w1[cases1], w2[cases2])
  }
  a.boot <- sort(a.boot)
  if (min(a.boot) == max(a.boot))
    ci.lower <- ci.upper <- a.boot[1]
  if ((a.obs < min(a.boot)) | (a.obs > max(a.boot))) {
    ci.lower <- a.boot[round((alpha / 2) * n.bootstrap)]
    ci.upper <- a.boot[round((1 - alpha / 2) * n.bootstrap)]
  }
  if ((ci.lower == pi) & (ci.upper == pi)) {
    z0 <- qnorm(mean(a.boot < a.obs))
    jk <- rep(0, (n1 + n2))
    for (i in 1:n1)
      jk[i] <- CalcA1(y1[-i], y2, weights, w1[-i], w2)
    for (i in 1:n2)
      jk[n1 + i] <- CalcA1(y1, y2[-i], weights, w1, w2[-i])
    diff <- mean(jk) - jk
    a <- sum(diff ^ 3) / (6 * (sum(diff ^ 2)) ^ 1.5)
    alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2)) / (1 - a * (z0 +
                    qnorm(alpha/2))))
    alpha2 <- pnorm(z0 + (z0 - qnorm(alpha/2)) / (1 - a * (z0 -
                    qnorm(alpha/2))))
    if (is.na(alpha1))
      alpha1 <- alpha / 2
    if (is.na(alpha2))
      alpha2 <- 1 - alpha / 2
    if (round(alpha1 * n.bootstrap) < 1) {
      ci.lower <- a.boot[1]
    } else {
      ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    }
  }
  if (ci.method == 2) {
    alpha1 <- (1 - conf.level) / 2
    alpha2 <- 1 - alpha1
    ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
    ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
  }
  se.a <- sd(a.boot)
  return (c(a.obs, se.a, ci.lower, ci.upper))
}

################################################################################
AAD1 <- function(y, r = 0, weights = FALSE, n.bootstrap = 1999,
                 conf.level = .95, ci.method = 1, seed = 1) {
# Calculates the confidence interval for the A statistic for the average
# absolute deviation for two or more groups.
#
# Args:
#   y           : Matrix of cases (rows) by scores (column 1) and group codes
#                 (column 2) (matrix).
#   r           : Vector of proportions (default = 0, represents equal
#                 proportions) (vector).
#   weights     : Whether to assign weights to cases (default = FALSE); if
#                 set to TRUE, y contains case weights in column 3.
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (scalar, default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for
#                 percentile)).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   A vector containing the A statistic, its estimated standard error, and the
#   upper and lower bounds of the confidence interval.
#
  set.seed(seed)
  gr <- sort(unique(y[,2]))
  k <- length(gr)
  ns <- rep(0, k)
  for (i in 1:k)
    ns[i] <- sum(y[,2] == gr[i])
  y.bs <- y
  a.boot <- rep(0, n.bootstrap)
  alpha <- 1 - conf.level
  ci.lower <- ci.upper <- pi
  a.obs <- CalcAAD1(y, r, weights)
  for (i in 1:n.bootstrap) {
    for (j in 1:k)
      y.bs[(y[,2] == gr[j]), 1] <- sample(y[(y[,2] == gr[j]), 1],
                                          replace = TRUE)
    a.boot[i] <- CalcAAD1(y.bs, r, weights)
  }
  a.boot <- sort(a.boot)
  if (min(a.boot) == max(a.boot))
    ci.lower <- ci.upper <- a.boot[1]
  if ((a.obs < min(a.boot)) | (a.obs > max(a.boot))) {
    ci.lower <- a.boot[round((alpha / 2) * n.bootstrap)]
    ci.upper <- a.boot[round((1 - alpha / 2) * n.bootstrap)]
  }
  if ((ci.lower == pi) & (ci.upper == pi)) {
    z0 <- qnorm(mean(a.boot < a.obs))
    n <- dim(y)[1]
    jk <- rep(0, n)
    for (i in 1:n)
      jk[i] <- CalcAAD1(y[-i, ], r, weights)
    diff <- mean(jk) - jk
    a <- sum(diff ^ 3) / (6 * (sum(diff ^ 2)) ^ 1.5)
    alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2)) / (1 -
                    a * (z0 + qnorm(alpha/2))))
    alpha2 <- pnorm(z0 + (z0 - qnorm(alpha/2)) / (1 -
                    a * (z0 - qnorm(alpha/2))))
    if (is.na(alpha1))
      alpha1 <- alpha / 2
    if (is.na(alpha2))
      alpha2 <- 1 - alpha / 2
    if (round(alpha1 * n.bootstrap) < 1) {
      ci.lower <- a.boot[1]
    } else {
      ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    }
  }
  if (ci.method == 2) {
    alpha1 <- (1 - conf.level) / 2
    alpha2 <- 1 - alpha1
    ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
    ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
  }
  se.a <- sd(a.boot)
  return (c(a.obs, se.a, ci.lower, ci.upper))
}

################################################################################
AAPD1 <- function(y, weights = FALSE, n.bootstrap = 1999, conf.level = .95,
                  ci.method = 1, seed = 1) {
# Calculates the confidence interval for the A statistic for the average
# absolute paired deviation for two or more groups.
#
# Args:
#   y           : Matrix of cases (rows) by scores (column 1) and group codes
#                 (column 2) (matrix).
#   weights     : Whether to assign weights to cases (default = FALSE); if
#                 set to TRUE, y contains case weights in column 3.
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for
#                 percentile)).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   A vector containing the A statistic, its estimated standard error, and the
#   upper and lower bounds of the confidence interval.
#
  set.seed(seed)
  gr <- sort(unique(y[,2]))
  k <- length(gr)
  ns <- rep(0, k)
  for (i in 1:k)
    ns[i] <- sum(y[,2] == gr[i])
  y.bs <- y
  a.boot <- rep(0, n.bootstrap)
  alpha <- 1 - conf.level
  ci.lower <- ci.upper <- pi
  a.obs <- CalcAAPD1(y, weights)
  for (i in 1:n.bootstrap) {
    for (j in 1:k)
      y.bs[(y[,2] == gr[j]), 1] <- sample(y[(y[,2] == gr[j]), 1],
                                          replace = TRUE)
    a.boot[i] <- CalcAAPD1(y.bs, weights)
  }
  a.boot <- sort(a.boot)
  if (min(a.boot) == max(a.boot))
    ci.lower <- ci.upper <- a.boot[1]
  if ((a.obs < min(a.boot)) | (a.obs > max(a.boot))) {
    ci.lower <- a.boot[round((alpha / 2) * n.bootstrap)]
    ci.upper <- a.boot[round((1 - alpha / 2) * n.bootstrap)]
  }
  if ((ci.lower == pi) & (ci.upper == pi)) {
    z0 <- qnorm(mean(a.boot < a.obs))
    n <- dim(y)[1]
    jk <- rep(0, n)
    for (i in 1:n)
      jk[i] <- CalcAAPD1(y[-i, ], weights)
    diff <- mean(jk) - jk
    a <- sum(diff ^ 3) / (6 * (sum(diff ^ 2)) ^ 1.5)
    alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2)) / (1 -
                    a * (z0 + qnorm(alpha/2))))
    alpha2 <- pnorm(z0 + (z0 - qnorm(alpha/2)) / (1 -
                    a * (z0 - qnorm(alpha/2))))
    if (is.na(alpha1))
      alpha1 <- alpha / 2
    if (is.na(alpha2))
      alpha2 <- 1 - alpha / 2
    if (round(alpha1 * n.bootstrap) < 1) {
      ci.lower <- a.boot[1]
    } else {
      ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    }
  }
  if (ci.method == 2) {
    alpha1 <- (1 - conf.level) / 2
    alpha2 <- 1 - alpha1
    ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
    ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
  }
  se.a <- sd(a.boot)
  return (c(a.obs, se.a, ci.lower, ci.upper))
}

################################################################################
IK1 <- function(y, ref = 1, weights = FALSE, n.bootstrap = 1999,
                conf.level = .95, ci.method = 1, seed = 1) {
# Calculates the confidence interval for the A statistic while singling out one
# group for two or more groups.
#
# Args:
#   y           : Matrix of cases (rows) by scores (column 1) and group codes
#                 (column 2) (matrix).
#   ref         : Reference group (to compare to all others) (scalar,
#                 default = 1).
#   weights     : Whether to assign weights to cases (default = FALSE); if
#                 set to TRUE, y contains case weights in column 3.
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for
#                 percentile)).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   A vector containing the A statistic, its estimated standard error, and the
#   upper and lower bounds of the confidence interval.
#
  set.seed(seed)
  gr <- sort(unique(y[,2]))
  k <- length(gr)
  ns <- rep(0, k)
  for (i in 1:k)
    ns[i] <- sum(y[,2] == gr[i])
  y.bs <- y
  a.boot <- rep(0, n.bootstrap)
  alpha <- 1 - conf.level
  ci.lower <- ci.upper <- pi
  a.obs <- CalcIK1(y, ref, weights)
  for (i in 1:n.bootstrap) {
    for (j in 1:k)
      y.bs[(y[,2] == gr[j]), 1] <- sample(y[(y[,2] == gr[j]), 1],
                                          replace = TRUE)
    a.boot[i] <- CalcIK1(y.bs, ref, weights)
  }
  a.boot <- sort(a.boot)
  if (min(a.boot) == max(a.boot))
    ci.lower <- ci.upper <- a.boot[1]
  if ((a.obs < min(a.boot)) | (a.obs > max(a.boot))) {
    ci.lower <- a.boot[round((alpha / 2) * n.bootstrap)]
    ci.upper <- a.boot[round((1 - alpha / 2) * n.bootstrap)]
  }
  if ((ci.lower == pi) & (ci.upper == pi)) {
    z0 <- qnorm(mean(a.boot < a.obs))
    n <- dim(y)[1]
    jk <- rep(0, n)
    for (i in 1:n)
      jk[i] <- CalcIK1(y[-i, ], ref, weights)
    diff <- mean(jk) - jk
    a <- sum(diff ^ 3) / (6 * (sum(diff ^ 2)) ^ 1.5)
    alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2)) / (1 - a * (z0 +
                    qnorm(alpha/2))))
    alpha2 <- pnorm(z0 + (z0 - qnorm(alpha/2)) / (1 - a * (z0 -
                    qnorm(alpha/2))))
    if (is.na(alpha1))
      alpha1 <- alpha / 2
    if (is.na(alpha2))
      alpha2 <- 1 - alpha / 2
    if (round(alpha1 * n.bootstrap) < 1) {
    	  ci.lower <- a.boot[1]
    } else {
    	  ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    }
  }
  if (ci.method == 2) {
    alpha1 <- (1 - conf.level) / 2
    alpha2 <- 1 - alpha1
    ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
    ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
  }
  se.a <- sd(a.boot)
  return (c(a.obs, se.a, ci.lower, ci.upper))
}

################################################################################
Ord1 <- function(y, weights = FALSE, increase = FALSE, n.bootstrap = 1999,
                 conf.level = .95, ci.method = 1, seed = 1) {
# Calculates the confidence interval for the ordinal comparison of the A
# statistic for two or more groups.
#
# Args:
#   y           : Matrix of cases (rows) by scores (column 1) and group codes
#                 (column 2) (matrix).
#   weights     : Whether to assign weights to cases (default = FALSE); if
#                 set to TRUE, y contains case weights in column 3.
#   increase    : Set to TRUE if scores are predicted to increase with group
#                 codes (default = FALSE).
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for
#                 percentile)).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   A vector containing the A statistic, its estimated standard error, and the
#   upper and lower bounds of the confidence interval.
#
  set.seed(seed)
  gr <- sort(unique(y[,2]))
  k <- length(gr)
  ns <- rep(0, k)
  for (i in 1:k)
    ns[i] <- sum(y[,2] == gr[i])
  y.bs <- y
  a.boot <- rep(0, n.bootstrap)
  alpha <- 1 - conf.level
  ci.lower <- ci.upper <- pi
  a.obs <- CalcOrd1(y, weights, increase)
  for (i in 1:n.bootstrap) {
    for (j in 1:k)
      y.bs[(y[,2] == gr[j]), 1] <- sample(y[(y[,2] == gr[j]), 1],
                                          replace = TRUE)
    a.boot[i] <- CalcOrd1(y.bs, weights, increase)
  }
  a.boot <- sort(a.boot)
  if (min(a.boot) == max(a.boot))
    ci.lower <- ci.upper <- a.boot[1]
  if ((a.obs < min(a.boot)) | (a.obs > max(a.boot))) {
      ci.lower <- a.boot[round((alpha / 2) * n.bootstrap)]
      ci.upper <- a.boot[round((1 - alpha / 2) * n.bootstrap)]
  }
  if ((ci.lower == pi) & (ci.upper == pi)) {
    z0 <- qnorm(mean(a.boot < a.obs))
    n <- dim(y)[1]
    jk <- rep(0, n)
    for (i in 1:n)
      jk[i] <- CalcOrd1(y[-i, ], weights, increase)
    diff <- mean(jk) - jk
    a <- sum(diff ^ 3) / (6 * (sum(diff ^ 2)) ^ 1.5)
    alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2)) / (1 - a * (z0 +
                    qnorm(alpha/2))))
    alpha2 <- pnorm(z0 + (z0 - qnorm(alpha/2)) / (1 - a * (z0 -
                    qnorm(alpha/2))))
    if (is.na(alpha1))
      alpha1 <- alpha / 2
    if (is.na(alpha2))
      alpha2 <- 1 - alpha / 2
    if (round(alpha1 * n.bootstrap) < 1) {
    	  ci.lower <- a.boot[1]
     } else {
       ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
       ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
     }
  }
  if (ci.method == 2) {
    alpha1 <- (1 - conf.level) / 2
    alpha2 <- 1 - alpha1
    ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
    ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
  }
  se.a <- sd(a.boot)
  return (c(a.obs, se.a, ci.lower, ci.upper))
}

################################################################################
A2 <- function(y1, y2, weights = FALSE, w = 0, n.bootstrap = 1999,
               conf.level = .95, ci.method = 1, seed = 1) {
# Calculates the standard error and constructs a confidence interval for the A
# statistic for two correlated samples using bootstrap methods.
#
# Args :
#   y1          : Scores for sample 1 (vector).
#   y2          : Scores for sample 2 (vector).
#   weights     : Whether to assign weights to cases (default = FALSE).
#   w           : Weights for cases (vector; default = 0).
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (scalar, default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for
#                 percentile)).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   A vector containing the A statistic, its estimated standard error, and the
#   upper and lower bounds of the confidence interval.
#
  n <- length(y1)
  if (!weights)
  	w <- rep(0, n)
  set.seed(seed)
  a.obs <- CalcA2(y1, y2, weights, w)
  alpha <- 1 - conf.level
  ci.lower <- ci.upper <- pi
  a.boot <- rep(0, n.bootstrap)
  for (i in 1:n.bootstrap) {
    cases <- sample(1:n, replace = TRUE)
    a.boot[i] <- CalcA2(y1[cases], y2[cases], weights, w[cases])
  }
  a.boot <- sort(a.boot)
  if (min(a.boot) == max(a.boot))
    ci.lower <- ci.upper <- a.boot[1]
  if ((a.obs < min(a.boot)) | (a.obs > max(a.boot))) {
    ci.lower <- a.boot[round((alpha / 2) * n.bootstrap)]
    ci.upper <- a.boot[round((1 - alpha / 2) * n.bootstrap)]
  }
  if ((ci.lower == pi) & (ci.upper == pi)) {
    z0 <- qnorm(mean(a.boot < a.obs))
    jk <- rep(0, (n * 2))
    for (i in 1:n)
      jk[i] <- CalcA2(y1[-i], y2[-i], weights, w[-i])
    for (i in 1:n)
      jk[n + i] <- CalcA2(y1[-i], y2[-i], weights, w[-i])
    diff <- mean(jk) - jk
    a <- sum(diff ^ 3) / (6 * (sum(diff ^ 2)) ^ 1.5)
    alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2)) / (1 - a * (z0 +
                    qnorm(alpha/2))))
    alpha2 <- pnorm(z0 + (z0 - qnorm(alpha/2)) / (1 - a * (z0 -
                    qnorm(alpha/2))))
    if (is.na(alpha1))
      alpha1 <- alpha / 2
    if (is.na(alpha2))
      alpha2 <- 1 - alpha / 2
    if (round(alpha1 * n.bootstrap) < 1) {
      ci.lower <- a.boot[1]
    } else {
      ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    }
  }
  if (ci.method == 2) {
    alpha1 <- (1 - conf.level) / 2
    alpha2 <- 1 - alpha1
    ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
    ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
  }
  se.a <- sd(a.boot)
  return (c(a.obs, se.a, ci.lower, ci.upper))
}

################################################################################
AAD2 <- function(y, r = 0, weights = FALSE, n.bootstrap = 1999,
                 conf.level = .95, ci.method = 1, seed = 1) {
# Calculates the confidence interval for the A statistic for the average
# absolute deviation for two or more correlated samples.
#
# Args:
#   y           : Matrix of scores with each sample in its own column (matrix).
#   r           : Vector of proportions (vector, default = 0, represents equal
#                 proportions).
#   weights     : Whether to assign weights to cases (default = FALSE); if
#                 set to TRUE, y contains case weights in final column.
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (scalar, default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for
#                 percentile)).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   A vector containing the A statistic, its estimated standard error, and the
#   upper and lower bounds of the confidence interval.
#
  set.seed(seed)
  n <- dim(y)[1]
  k <- dim(y)[2]
  if (weights)
    k <- k - 1
  y.bs <- y
  a.boot <- rep(0, n.bootstrap)
  alpha <- 1 - conf.level
  ci.lower <- ci.upper <- pi
  a.obs <- CalcAAD2(y, r, weights)
  for (i in 1:n.bootstrap) {
    y.bs[, 1:k] <- y[sample(1:n, replace = TRUE), 1:k]
    a.boot[i] <- CalcAAD2(y.bs, r, weights)
  }
  a.boot <- sort(a.boot)
  if (min(a.boot) == max(a.boot))
    ci.lower <- ci.upper <- a.boot[1]
  if ((a.obs < min(a.boot)) | (a.obs > max(a.boot))) {
    ci.lower <- a.boot[round((alpha / 2) * n.bootstrap)]
    ci.upper <- a.boot[round((1 - alpha / 2) * n.bootstrap)]
  }
  if ((ci.lower == pi) & (ci.upper == pi)) {
    z0 <- qnorm(mean(a.boot < a.obs))
    jk <- rep(0, n)
    for (i in 1:n)
      jk[i] <- CalcAAD2(y[-i, ], r, weights)
    diff <- mean(jk) - jk
    a <- sum(diff ^ 3) / (6 * (sum(diff ^ 2)) ^ 1.5)
    alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2)) / (1 - a * (z0 +
                    qnorm(alpha/2))))
    alpha2 <- pnorm(z0 + (z0 - qnorm(alpha/2)) / (1 - a * (z0 -
                    qnorm(alpha/2))))
    if (is.na(alpha1))
      alpha1 <- alpha / 2
    if (is.na(alpha2))
      alpha2 <- 1 - alpha / 2
    if (round(alpha1 * n.bootstrap) < 1) {
      ci.lower <- a.boot[1]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    } else {
      ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    }
  }
  if (ci.method == 2) {
    alpha1 <- (1 - conf.level) / 2
    alpha2 <- 1 - alpha1
    ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
    ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
  }
  se.a <- sd(a.boot)
  return (c(a.obs, se.a, ci.lower, ci.upper))
}

################################################################################
AAPD2 <- function(y, weights = FALSE, n.bootstrap = 1999, conf.level = .95,
                  ci.method = 1, seed = 1) {
# Calculates the confidence interval for the A statistic for the average
# absolute paired deviation for two or more correlated samples.
#
# Args:
#   y           : Matrix of cases (rows) by scores (column 1) and group codes
#                 (column 2) (matrix).
#   weights     : Whether to assign weights to cases (default = FALSE); if
#                 set to TRUE, y contains case weights in final column.
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for
#                 percentile)).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   A vector containing the A statistic, its estimated standard error, and the
#   upper and lower bounds of the confidence interval.
#
  set.seed(seed)
  n <- dim(y)[1]
  k <- dim(y)[2]
  if (weights)
    k <- k - 1
  y.bs <- y
  a.boot <- rep(0, n.bootstrap)
  alpha <- 1 - conf.level
  ci.lower <- ci.upper <- pi
  a.obs <- CalcAAPD2(y, weights)
  for (i in 1:n.bootstrap) {
    y.bs[, 1:k] <- y[sample(1:n, replace = TRUE), 1:k]
    a.boot[i] <- CalcAAPD2(y.bs, weights)
  }
  a.boot <- sort(a.boot)
  if (min(a.boot) == max(a.boot))
    ci.lower <- ci.upper <- a.boot[1]
  if ((a.obs < min(a.boot)) | (a.obs > max(a.boot))) {
    ci.lower <- a.boot[round((alpha / 2) * n.bootstrap)]
    ci.upper <- a.boot[round((1 - alpha / 2) * n.bootstrap)]
  }
  if ((ci.lower == pi) & (ci.upper == pi)) {
    z0 <- qnorm(mean(a.boot < a.obs))
    jk <- rep(0, n)
    for (i in 1:n)
      jk[i] <- CalcAAPD2(y[-i, ], weights)
    diff <- mean(jk) - jk
    a <- sum(diff ^ 3) / (6 * (sum(diff ^ 2)) ^ 1.5)
    alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2)) / (1 - a * (z0 +
                    qnorm(alpha/2))))
    alpha2 <- pnorm(z0 + (z0 - qnorm(alpha/2)) / (1 - a * (z0 -
                    qnorm(alpha/2))))
    if (is.na(alpha1))
      alpha1 <- alpha / 2
    if (is.na(alpha2))
      alpha2 <- 1 - alpha / 2
    if (round(alpha1 * n.bootstrap) < 1) {
      ci.lower <- a.boot[1]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    } else {
      ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    }
  }
  if (ci.method == 2) {
    alpha1 <- (1 - conf.level) / 2
    alpha2 <- 1 - alpha1
    ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
    ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
  }
  se.a <- sd(a.boot)
  return (c(a.obs, se.a, ci.lower, ci.upper))
}

################################################################################
IK2 <- function(y, ref = 1, weights = FALSE, n.bootstrap = 1999,
                conf.level = .95, ci.method = 1, seed = 1) {
# Calculates the confidence interval for the A statistic while singling out one
# group for two or more correlated samples.
#
# Args:
#   y           : Matrix of cases (rows) by scores (column 1) and group codes
#                 (column 2) (matrix).
#   ref         : reference group (to compare to all others) (scalar,
#                 default = 1).
#   weights     : Whether to assign weights to cases (default = FALSE); if
#                 set to TRUE, y contains case weights in final column.
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for
#                 percentile)).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   A vector containing the A statistic, its estimated standard error, and the
#   upper and lower bounds of the confidence interval.
#
  set.seed(seed)
  n <- dim(y)[1]
  k <- dim(y)[2]
  if (weights)
    k <- k - 1
  y.bs <- y
  a.boot <- rep(0, n.bootstrap)
  alpha <- 1 - conf.level
  ci.lower <- ci.upper <- pi
  a.obs <- CalcIK2(y, ref, weights)
  for (i in 1:n.bootstrap) {
    y.bs[, 1:k] <- y[sample(1:n, replace = TRUE), 1:k]
    a.boot[i] <- CalcIK2(y.bs, ref, weights)
  }
  a.boot <- sort(a.boot)
  if (min(a.boot) == max(a.boot))
    ci.lower <- ci.upper <- a.boot[1]
  if ((a.obs < min(a.boot)) | (a.obs > max(a.boot))) {
    ci.lower <- a.boot[round((alpha / 2) * n.bootstrap)]
    ci.upper <- a.boot[round((1 - alpha / 2) * n.bootstrap)]
  }
  if ((ci.lower == pi) & (ci.upper == pi)) {
    z0 <- qnorm(mean(a.boot < a.obs))
    n <- dim(y)[1]
    jk <- rep(0, n)
    for (i in 1:n)
      jk[i] <- CalcIK2(y[-i, ], ref, weights)
    diff <- mean(jk) - jk
    a <- sum(diff ^ 3) / (6 * (sum(diff ^ 2)) ^ 1.5)
    alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2)) / (1 - a * (z0 +
                    qnorm(alpha/2))))
    alpha2 <- pnorm(z0 + (z0 - qnorm(alpha/2)) / (1 - a * (z0 -
                    qnorm(alpha/2))))
    if (is.na(alpha1))
      alpha1 <- alpha / 2
    if (is.na(alpha2))
      alpha2 <- 1 - alpha / 2
    if (round(alpha1 * n.bootstrap) < 1) {
    	  ci.lower <- a.boot[1]
    } else {
    	  ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    }
  }
  if (ci.method == 2) {
    alpha1 <- (1 - conf.level) / 2
    alpha2 <- 1 - alpha1
    ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
    ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
  }
  se.a <- sd(a.boot)
  return (c(a.obs, se.a, ci.lower, ci.upper))
}

################################################################################
Ord2 <- function(y, weights = FALSE, increase = FALSE, n.bootstrap = 1999,
                 conf.level = .95, ci.method = 1, seed = 1) {
# Calculates the confidence interval for the ordinal comparison of the A
# statistic for two or more correlated samples.
#
# Args:
#   y           : Matrix of cases (rows) by scores (column 1) and group codes
#                 (column 2) (matrix).
#   weights     : Whether to assign weights to cases (default = FALSE); if
#                 set to TRUE, y contains case weights in final column.
#   increase    : Set to TRUE if scores are predicted to increase with group
#                 codes (default = FALSE).
#   n.bootstrap : Number of bootstrap samples (scalar, default = 1999).
#   conf.level  : Confidence level (default = .95).
#   ci.method   : Method used to construct confidence interval (scalar,
#                 default = 1 (for BCA), user can also call 2 (for
#                 percentile)).
#   seed        : Random number seed (scalar, default = 1).
#
# Returns :
#   A vector containing the A statistic, its estimated standard error, and the
#   upper and lower bounds of the confidence interval.
#
  set.seed(seed)
  n <- dim(y)[1]
  k <- dim(y)[2]
  if (weights)
    k <- k - 1
  y.bs <- y
  a.boot <- rep(0, n.bootstrap)
  alpha <- 1 - conf.level
  ci.lower <- ci.upper <- pi
  a.obs <- CalcOrd2(y, weights, increase)
  for (i in 1:n.bootstrap) {
    y.bs[, 1:k] <- y[sample(1:n, replace = TRUE), 1:k]
    a.boot[i] <- CalcOrd2(y.bs, weights, increase)
  }
  a.boot <- sort(a.boot)
  if (min(a.boot) == max(a.boot))
    ci.lower <- ci.upper <- a.boot[1]
  if ((a.obs < min(a.boot)) | (a.obs > max(a.boot))) {
    ci.lower <- a.boot[round((alpha / 2) * n.bootstrap)]
    ci.upper <- a.boot[round((1 - alpha / 2) * n.bootstrap)]
  }
  if ((ci.lower == pi) & (ci.upper == pi)) {
    z0 <- qnorm(mean(a.boot < a.obs))
    jk <- rep(0, n)
    for (i in 1:n)
      jk[i] <- CalcOrd2(y[-i, ], weights, increase)
    diff <- mean(jk) - jk
    a <- sum(diff ^ 3) / (6 * (sum(diff ^ 2)) ^ 1.5)
    alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2)) / (1 - a * (z0 +
                    qnorm(alpha/2))))
    alpha2 <- pnorm(z0 + (z0 - qnorm(alpha/2)) / (1 - a * (z0 -
                    qnorm(alpha/2))))
    if (is.na(alpha1))
      alpha1 <- alpha / 2
    if (is.na(alpha2))
      alpha2 <- 1 - alpha / 2
    if (round(alpha1 * n.bootstrap) < 1) {
      ci.lower <- a.boot[1]
    } else {
    	  ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
      ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
    }
  }
  if (ci.method == 2) {
    alpha1 <- (1 - conf.level) / 2
    alpha2 <- 1 - alpha1
    ci.lower <- a.boot[round(alpha1 * n.bootstrap)]
    ci.upper <- a.boot[round(alpha2 * n.bootstrap)]
  }
  se.a <- sd(a.boot)
  return (c(a.obs, se.a, ci.lower, ci.upper))
}

################################################################################
CalcA1 <- function(y1, y2, weights = FALSE, w1 = 0, w2 = 0) {
# Calculates the A statistic for 2 groups.
#
# Args:
#   y1     : Scores for group 1 (vector).
#   y2     : Scores for group 2 (vector).
#   weights: Whether to weight cases (default = FALSE).
#   w1     : Weights for group 1 (optional) (vector, default is 0).
#   w2     : Weights for group 2 (optional) (vector, default is 0).
#
# Returns:
# a : The A statistic
#
  n1 <- length(y1)
  n2 <- length(y2)
  if (sum(w1) == 0) {
    r1 <- sum(rank(c(y1, y2))[1:n1])
    a <- (r1 / n1 - (n1 + 1) / 2) / n2
  } else {
    num <- den <- 0
    if (n2 < n1) {
      yt <- y1
      y1 <- y2
      y2 <- yt
      wt <- w1
      w1 <- w2
      w2 <- wt
    }
    for (i in 1:min(n1,n2)) {
      num <- num + sum(w1[i] * w2 * ((y1[i] > y2) + .5 * (y1[i] == y2)))
      den <- den + sum(w1[i] * w2)
    }
    a <- num / den
    if (n2 < n1)
      a <- 1 - a
  }
  return(a)
}

################################################################################
CalcAAD1 <- function(y, r = 0, weights = FALSE) {
# Calculates the A statistic for the average absolute deviation for two or more
# groups. Note: This function is not meant to be called by the user, but is
# called by AAD1.
#
# Args:
#   y       : Matrix of cases (rows) by scores (column 1) and group codes
#             (column 2) (matrix).
#   r       : Vector of proportions (vector, default = 0, represents equal
#             proportions).
#   weights : Whether to assign weights to cases (default = FALSE); if
#             set to TRUE, y contains case weights in column 3.
#
# Returns:
#   a : The A statistic
#
  gr <- sort(unique(y[,2]))
  k <- length(gr)
  if ((sum(r) == 0) | (length(r) != k))
    r <- rep(1/k, k)
  if (sum(r) != 1)
    r <- r / sum(r)
  a <- rep(0, k)
  for (i in 1:k) {
    g1 <- y[(y[,2] == gr[i]), 1]
    for (j in 1:k) {
      if (i != j) {
        g2 <- y[(y[,2] == gr[j]), 1]
          if (!weights) {
            a[i] <- a[i] + (r[j] / (1 - r[i])) * CalcA1(g1, g2)
          } else {
            w1 <- y[(y[,2] == gr[i]), 3]
            w2 <- y[(y[,2] == gr[j]), 3]
            a[i] <- a[i] + (r[j] / (1 - r[i])) * CalcA1(g1, g2, weights = TRUE,
                    w1, w2)
          }
      }
    }
  }
  return(mean(abs(a - .5)) + .50)
}


################################################################################
CalcAAPD1 <- function(y, weights = FALSE) {
# Calculates the A statistic for the average absolute paired deviation for two
# or more groups. Note: This function is not meant to be called by the user, but
# it is called by AAPD1.
#
# Args:
#   y       : Matrix of cases (rows) by scores (column 1) and group codes
#             (column 2) (matrix).
#   weights : Whether to assign weights to cases (default = FALSE); if
#             set to TRUE, y contains case weights in column 3.
#
# Returns:
#   a : The A statistic
#
  gr <- sort(unique(y[,2]))
  k <- length(gr)
  a <- rep(0, k * (k - 1) / 2)
  m <- 0
  for (i in 1:(k - 1))
    for (j in (i + 1):k) {
      m <- m + 1
      g1 <- y[(y[,2] == gr[i]), 1]
      g2 <- y[(y[,2] == gr[j]), 1]
      if (!weights) {
        a[m] <- CalcA1(g1, g2)
      } else {
        w1 <- y[(y[,2] == gr[i]), 3]
        w2 <- y[(y[,2] == gr[j]), 3]
        a[m] <- CalcA1(g1, g2, weights = TRUE, w1, w2)
      }
    }
  return(mean(abs(a - .5)) + .50)
}

################################################################################
CalcIK1 <- function(y, ref = 1, weights = FALSE) {
# Calculates the A statistic while singling out one group for two or more groups.
# Note: This function is not meant to be called by the user, but it is called by
# IK1.
#
# Args:
#   y       : Matrix of cases (rows) by scores (column 1) and group codes
#             (column 2) (matrix).
#   ref     : Reference group (to compare to all others) (scalar, default = 1).
#   weights : Whether to assign weights to cases (default = FALSE); if
#             set to TRUE, y contains case weights in column 3.
#
# Returns:
#   a : The A statistic
#
  g1 <- y[(y[,2] == ref), 1]
  g2 <- y[(y[,2] != ref), 1]
  if (!weights) {
    a <- CalcA1(g1, g2)
  } else {
    w1 <- y[(y[,2] == ref), 3]
    w2 <- y[(y[,2] != ref), 3]
    a <- CalcA1(g1, g2, weights = TRUE, w1, w2)
  }
  return(a)
}

################################################################################
CalcOrd1 <- function(y, weights = FALSE, increase = FALSE) {
# Calculates the ordinal comparison of the A statistic for two or more groups.
# Note: This function is not meant to be called by the user, but it is called by
# AOrd1.
#
# Args:
#   y        : Matrix of cases (rows) by scores (column 1) and group codes
#              (column 2) (matrix).
#   weights  : Whether to assign weights to cases (default = FALSE); if
#              set to TRUE, y contains case weights in column 3.
#   increase : Set to TRUE if scores are predicted to increase with group codes
#              (default = FALSE).
#
# Returns:
#   a : The A statistic
#
  gr <- sort(unique(y[,2]))
  k <- length(gr)
  a <- rep(0, k - 1)
  for (i in 1:(k - 1)) {
    g1 <- y[(y[,2] == gr[i]), 1]
    g2 <- y[(y[,2] == gr[i + 1]), 1]
    if (!weights) {
    	  a[i] <- CalcA1(g1, g2)
    } else {
      w1 <- y[(y[,2] == gr[i]), 3]
      w2 <- y[(y[,2] == gr[i + 1]), 3]
      a[i] <- CalcA1(g1, g2, weights = TRUE, w1, w2)
    }
  }
  if (increase) {
    return (1 - mean(a))
  } else {
    return(mean(a))
  }
}

################################################################################
CalcA2 <- function(y1, y2, weights = FALSE, w = 0) {
# Calculates the A statistic for 2 correlated samples.
#
# Args:
#   y1     : Scores for sample 1 (vector).
#   y2     : Scores for sample 2 (vector).
#   weights: Whether to weight cases (default = FALSE).
#   w      : Weights for cases (optional) (vector, default is 0).
#
# Returns:
# a : The A statistic
#
  n <- length(y1)
  if (!weights)
    w <- rep(1, n)
  num <- sum(w * ((y1 > y2) + .5 * (y1 == y2)))
  den <- sum(w)
  a <- num / den
  return(a)
}

################################################################################
CalcAAD2 <- function(y, r = 0, weights = FALSE) {
# Calculates the A statistic for the average absolute deviation for two or more
# correlated samples. Note: This function is not meant to be called by the user,
# but it is called by AAD2.
#
# Args:
#   y       : Matrix of scores with each sample in its own column (matrix).
#   r       : Vector of proportions (vector, default = 0, represents equal
#             proportions).
#   weights : Whether to assign weights to cases (default = FALSE); if
#             set to TRUE, y contains case weights in final column.
#
# Returns:
#   a : The A statistic
#
  k <- dim(y)[2]
  if (weights)
    k <- k - 1
  if ((sum(r) == 0) | (length(r) != k))
    r <- rep(1/k, k)
  if (sum(r) != 1)
    r <- r / sum(r)
  a <- rep(0, k)
  for (i in 1:k) {
    v1 <- y[, i]
    for (j in 1:k) {
      if (i != j) {
        v2 <- y[, j]
        if (!weights) {
          a[i] <- a[i] + (r[j] / (1 - r[i])) * CalcA2(v1, v2)
        } else {
          w <- y[, k + 1]
          a[i] <- a[i] + (r[j] / (1 - r[i])) * CalcA2(v1, v2,
                  weights = TRUE, w)
        }
      }
    }
  }
  return(mean(abs(a - .5)) + .50)
}

################################################################################
CalcAAPD2 <- function(y, weights = FALSE) {
# Calculates the A statistic for the average absolute paired deviation for two
# or more correlated samples. Note: This function is not meant to be called by
# the user, but it is called by AAPD2.
#
# Args:
#   y       : Matrix of scores with each sample in its own column (matrix).
#   weights : Whether to assign weights to cases (default = FALSE); if
#             set to TRUE, y contains case weights in final column.
#
# Returns:
#   a : The A statistic
#
  k <- dim(y)[2]
  if (weights)
    k <- k - 1
  a <- rep(0, k * (k - 1) / 2)
  m <- 0
  for (i in 1:(k - 1))
  for (j in (i + 1):k) {
    m <- m + 1
    v1 <- y[, i]
    v2 <- y[, j]
    if (!weights) {
      a[m] <- CalcA2(v1, v2)
    } else {
       w <- y[, k + 1]
       a[m] <- CalcA2(v1, v2, weights = TRUE, w)
    }
  }
  return(mean(abs(a - .5)) + .50)
}

################################################################################
CalcIK2 <- function(y, ref = 1, weights = FALSE) {
# Calculates the A statistic while singling out one group for two or more
# correlated samples. Note: This function is not meant to be called by the user,
# but it is called by IK2.
#
# Args:
#   y       : Matrix of scores with each sample in its own column (matrix).
#   ref     : reference group (to compare to all others) (scalar, default = 1).
#   weights : Whether to assign weights to cases (default = FALSE); if
#             set to TRUE, y contains case weights in final column.
#
# Returns:
#   a : The A statistic
#
  k <- dim(y)[2]
  if (weights)
    k <- k - 1
  v1 <- rep(as.vector(y[, (1:k)[ref]]), k - 1)
  v2 <- as.vector(y[, (1:k)[-ref]])
  if (!weights) {
    a <- CalcA2(v1, v2)
  } else {
    w <- rep(y[, k + 1], k - 1)
    a <- CalcA2(v1, v2, weights = TRUE, w)
  }
  return(a)
}

################################################################################
CalcOrd2 <- function(y, weights = FALSE, increase = FALSE) {
# Calculates the ordinal comparison of the A statistic for two or more
# correlated samples. Note: This function is not meant to be called by the user,
# but it is called by AOrd2.
#
# Args:
#   y        : Matrix of scores with each sample in its own column (matrix).
#   weights  : Whether to assign weights to cases (default = FALSE); if
#              set to TRUE, y contains case weights in final column.
#   increase : Set to TRUE if scores are predicted to increase with group codes
#              (default = FALSE).
#
# Returns:
#   a : The A statistic
#
  n <- dim(y)[1]
  k <- dim(y)[2]
  if (weights)
    k <- k - 1
  a <- rep(0, k - 1)
  for (i in 1:(k - 1)) {
    v1 <- y[, i]
    v2 <- y[, i + 1]
    if (!weights) {
    	  a[i] <- CalcA2(v1, v2)
    } else {
      w <- y[, k + 1]
      a[i] <- CalcA2(v1, v2, weights = TRUE, w)
    }
  }
  if (increase) {
    return (1 - mean(a))
  } else {
    return(mean(a))
  }
}
utils::globalVariables(c("y"))
