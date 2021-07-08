# ciBinomial roxy [sinew] ----
#' @seealso
#'  \code{\link[stats]{Normal}},\code{\link[stats]{uniroot}}
#' @rdname varBinomial
#' @export
#' @importFrom stats qnorm uniroot
# ciBinomial function [sinew] ----
ciBinomial <- function(x1, x2, n1, n2, alpha = .05, adj = 0, scale = "Difference") {
  # check input arguments
  checkScalar(n1, "integer", c(1, Inf))
  checkScalar(n2, "integer", c(1, Inf))
  checkScalar(x1, "integer", c(0, n1))
  checkScalar(x2, "integer", c(0, n2))
  checkScalar(alpha, "numeric", c(0, 1), c(FALSE, FALSE))
  checkScalar(adj, "integer", c(0, 1))
  checkScalar(scale, "character")
  scale <- match.arg(tolower(scale), c("difference", "rr", "or"))
  checkLengths(x1, x2)

  if (scale == "difference") {
    delta <- x1 / n1 - x2 / n2

    if (delta == -1) {
      lower <- -1
    }
    else if (testBinomial(delta0 = -.9999, x1 = x1, x2 = x2, n1 = n1, n2 = n2, adj = adj)
    < stats::qnorm(alpha / 2)) {
      lower <- -1
    }
    else {
      lower <- stats::uniroot(bpdiff,
        interval = c(-.9999, delta), x1 = x1, x2 = x2,
        n1 = n1, n2 = n2, adj = adj, alpha = alpha
      )$root
    }

    if (delta == 1) {
      upper <- 1
    }
    else if (testBinomial(delta0 = .9999, x1 = x1, x2 = x2, n1 = n1, n2 = n2, adj = adj)
    > -stats::qnorm(alpha / 2)) {
      upper <- 1
    }
    else {
      upper <- stats::uniroot(bpdiff,
        interval = c(delta, .9999), lower.tail = TRUE,
        x1 = x1, x2 = x2, n1 = n1, n2 = n2, adj = adj,
        alpha = alpha
      )$root
    }
  }
  else if (scale == "rr") {
    if (x1 == 0) {
      lower <- 0
    } else {
      lower <- stats::uniroot(bpdiff,
        interval = c(-20, 20), x1 = x1, x2 = x2,
        n1 = n1, n2 = n2, adj = adj, scale = "RR", alpha = alpha
      )$root
      lower <- exp(lower)
    }
    if (x2 == 0) {
      upper <- Inf
    } else {
      upper <- stats::uniroot(bpdiff,
        interval = c(-20, 20), lower.tail = TRUE,
        x1 = x1, x2 = x2, n1 = n1, n2 = n2, adj = adj, scale = "RR",
        alpha = alpha
      )$root
      upper <- exp(upper)
    }
  }
  else {
    if (x1 == 0 || x2 == n2) {
      lower <- -Inf
    } else {
      lower <- stats::uniroot(bpdiff,
        interval = c(-10, 10), x1 = x1, x2 = x2,
        n1 = n1, n2 = n2, adj = adj, scale = scale, alpha = alpha
      )$root
      lower <- exp(lower)
    }
    if (x2 == 0 || x1 == n1) {
      upper <- Inf
    } else {
      upper <- stats::uniroot(bpdiff,
        interval = c(-10, 10), lower.tail = TRUE,
        x1 = x1, x2 = x2, n1 = n1, n2 = n2, adj = adj, scale = scale,
        alpha = alpha
      )$root
      upper <- exp(upper)
    }
  }
  data.frame(lower = lower, upper = upper)
}

# nBinomial roxy [sinew] ----
#' @rdname varBinomial
#' @export
#' @importFrom stats qnorm pnorm
# nBinomial function [sinew] ----
nBinomial <- function(p1, p2, alpha = 0.025, beta = 0.1, delta0 = 0, ratio = 1, sided = 1, outtype = 1, scale = "Difference", n = NULL) {
  checkVector(p1, "numeric", c(0, 1), c(FALSE, FALSE))
  checkVector(p2, "numeric", c(0, 1), c(FALSE, FALSE))
  checkScalar(sided, "integer", c(1, 2))
  checkScalar(alpha, "numeric", c(0, 1 / sided), c(FALSE, FALSE))
  checkVector(beta, "numeric", c(0, 1 - alpha / sided), c(
    FALSE,
    FALSE
  ))
  checkVector(delta0, "numeric")
  checkVector(ratio, "numeric", c(0, Inf), c(FALSE, FALSE))
  checkScalar(outtype, "integer", c(1, 3))
  checkScalar(scale, "character")
  if (!is.null(n)) {
    checkVector(n, "numeric")
  }
  scale <- match.arg(tolower(scale), c(
    "difference", "rr",
    "or", "lnor"
  ))
  if (is.null(n)) {
    checkLengths(p1, p2, beta, delta0, ratio, allowSingle = TRUE)
  } else {
    checkLengths(n, p1, p2, delta0, ratio, allowSingle = TRUE)
  }
  len <- max(sapply(list(p1, p2, beta, delta0, ratio), length))
  if (len > 1) {
    if (length(p1) == 1) {
      p1 <- rep(p1, len)
    }
    if (length(p2) == 1) {
      p2 <- rep(p2, len)
    }
    if (length(alpha) == 1) {
      alpha <- rep(alpha, len)
    }
    if (length(beta) == 1) {
      beta <- rep(beta, len)
    }
    if (length(delta0) == 1) {
      delta0 <- rep(delta0, len)
    }
    if (length(ratio) == 1) {
      ratio <- rep(ratio, len)
    }
  }
  if (max(delta0 == 0) > 0 && max(p1[delta0 == 0] == p2[delta0 ==
    0]) > 0) {
    stop("p1 may not equal p2 when delta0 is zero")
  }
  z.beta <- stats::qnorm(1 - beta)
  sided[sided != 2] <- 1
  z.alpha <- stats::qnorm(1 - alpha / sided)
  d0 <- (delta0 == 0)
  if (scale == "difference") {
    if (min(abs(p1 - p2 - delta0)) < 1e-11) {
      stop("p1 - p2 may not equal delta0 when scale is \"Difference\"")
    }
    a <- 1 + ratio
    b <- -(a + p1 + ratio * p2 + delta0 * (ratio + 2))
    c <- delta0^2 + delta0 * (2 * p1 + a) + p1 + ratio *
      p2
    d <- -p1 * delta0 * (1 + delta0)
    v <- (b / (3 * a))^3 - b * c / 6 / a^2 + d / 2 / a
    u <- (sign(v) + (v == 0)) * sqrt((b / 3 / a)^2 - c / 3 / a)
    w <- (pi + acos(v / u^3)) / 3
    p10 <- 2 * u * cos(w) - b / 3 / a
    p20 <- p10 - delta0
    p10[d0] <- (p1[d0] + ratio[d0] * p2[d0]) / (1 + ratio[d0])
    p20[d0] <- p10[d0]
    sigma0 <- sqrt((p10 * (1 - p10) + p20 * (1 - p20) / ratio) *
      (ratio + 1))
    sigma1 <- sqrt((p1 * (1 - p1) + p2 * (1 - p2) / ratio) *
      (ratio + 1))
    if (is.null(n)) {
      n <- ((z.alpha * sigma0 + z.beta * sigma1) / (p1 - p2 - delta0))^2
      if (outtype == 2) {
        return(data.frame(cbind(n1 = n / (ratio + 1), n2 = ratio *
          n / (ratio + 1))))
      }
      else if (outtype == 3) {
        return(data.frame(cbind(
          n = n, n1 = n / (ratio + 1),
          n2 = ratio * n / (ratio + 1), alpha = alpha,
          sided = sided, beta = beta, Power = 1 - beta,
          sigma0 = sigma0, sigma1 = sigma1, p1 = p1,
          p2 = p2, delta0 = delta0, p10 = p10, p20 = p20
        )))
      }
      else {
        return(n = n)
      }
    }
    else {
      pwr <- stats::pnorm(-(stats::qnorm(1 - alpha / sided) - sqrt(n) *
        ((p1 - p2 - delta0) / sigma0)) * sigma0 / sigma1)
      if (outtype == 2) {
        return(data.frame(cbind(n1 = n / (ratio + 1), n2 = ratio *
          n / (ratio + 1), Power = pwr)))
      }
      else if (outtype == 3) {
        return(data.frame(cbind(
          n = n, n1 = n / (ratio + 1),
          n2 = ratio * n / (ratio + 1), alpha = alpha,
          sided = sided, beta = 1 - pwr, Power = pwr,
          sigma0 = sigma0, sigma1 = sigma1, p1 = p1,
          p2 = p2, delta0 = delta0, p10 = p10, p20 = p20
        )))
      }
      else {
        return(Power = pwr)
      }
    }
  }
  else if (scale == "rr") {
    RR <- exp(delta0)
    if (min(abs(p1 / p2 - RR)) < 1e-07) {
      stop("p1/p2 may not equal exp(delta0) when scale=\"RR\"")
    }
    a <- (1 + ratio)
    b <- -(RR * (1 + ratio * p2) + ratio + p1)
    c <- RR * (p1 + ratio * p2)
    p10 <- (-b - sqrt(b^2 - 4 * a * c)) / 2 / a
    p20 <- p10 / RR
    p10[d0] <- (p1[d0] + ratio[d0] * p2[d0]) / (1 + ratio[d0])
    p20[d0] <- p10[d0]
    sigma0 <- sqrt((ratio + 1) * (p10 * (1 - p10) + RR^2 *
      p20 * (1 - p20) / ratio))
    sigma1 <- sqrt((ratio + 1) * (p1 * (1 - p1) + RR^2 *
      p2 * (1 - p2) / ratio))
    if (is.null(n)) {
      n <- ((z.alpha * sigma0 + z.beta * sigma1) / (p1 - p2 * RR))^2
      if (outtype == 2) {
        return(data.frame(cbind(
          n1 = n / (ratio + 1),
          n2 = ratio * n / (ratio + 1)
        )))
      }
      else if (outtype == 3) {
        return(data.frame(cbind(
          n = n, n1 = n / (ratio + 1),
          n2 = ratio * n / (ratio + 1), alpha = alpha,
          sided = sided, beta = beta, Power = 1 - beta,
          sigma0 = sigma0, sigma1 = sigma1, p1 = p1,
          p2 = p2, delta0 = delta0, p10 = p10, p20 = p20
        )))
      }
      else {
        return(n = n)
      }
    }
    else {
      pwr <- stats::pnorm(-(stats::qnorm(1 - alpha / sided) - sqrt(n) *
        ((p1 - p2 * RR) / sigma0)) * sigma0 / sigma1)
      if (outtype == 2) {
        return(data.frame(cbind(n1 = n / (ratio + 1), n2 = ratio *
          n / (ratio + 1), Power = pwr)))
      }
      else if (outtype == 3) {
        return(data.frame(cbind(
          n = n, n1 = n / (ratio + 1),
          n2 = ratio * n / (ratio + 1), alpha = alpha,
          sided = sided, beta = 1 - pwr, Power = pwr,
          sigma0 = sigma0, sigma1 = sigma1, p1 = p1,
          p2 = p2, delta0 = delta0, p10 = p10, p20 = p20
        )))
      }
      else {
        return(Power = pwr)
      }
    }
  }
  else {
    OR <- exp(-delta0)
    if (min(abs(p1 / (1 - p1) / p2 * (1 - p2) * OR - 1)) < 1e-07) {
      stop("p1/(1-p1)/p2*(1-p2) may not equal exp(delta0) when scale=\"OR\"")
    }
    a <- OR - 1
    b <- 1 + ratio * OR + (1 - OR) * (ratio * p2 + p1)
    c <- -(ratio * p2 + p1)
    p10 <- (-b + sqrt(b^2 - 4 * a * c)) / 2 / a
    p20 <- OR * p10 / (1 + p10 * (OR - 1))
    p10[d0] <- (p1[d0] + ratio[d0] * p2[d0]) / (1 + ratio[d0])
    p20[d0] <- p10[d0]
    sigma0 <- sqrt((ratio + 1) * (1 / p10 / (1 - p10) + 1 / p20 / (1 - p20) / ratio))
    sigma1 <- sqrt((ratio + 1) * (1 / p1 / (1 - p1) + 1 / p2 / (1 - p2) / ratio))
    if (is.null(n)) {
      n <- ((z.alpha * sigma0 + z.beta * sigma1) / log(OR / p2 * (1 - p2) * p1 / (1 - p1)))^2
      if (outtype == 2) {
        return(data.frame(cbind(n1 = n / (ratio + 1), n2 = ratio *
          n / (ratio + 1))))
      }
      else if (outtype == 3) {
        return(data.frame(cbind(
          n = n, n1 = n / (ratio + 1),
          n2 = ratio * n / (ratio + 1), alpha = alpha,
          sided = sided, beta = beta, Power = 1 - beta,
          sigma0 = sigma0, sigma1 = sigma1, p1 = p1,
          p2 = p2, delta0 = delta0, p10 = p10, p20 = p20
        )))
      }
      else {
        return(n = n)
      }
    }
    else {
      pwr <- stats::pnorm(-(stats::qnorm(1 - alpha / sided) - sqrt(n) *
        (log(OR / p2 * (1 - p2) * p1 / (1 - p1)) / sigma0)) *
        sigma0 / sigma1)
      if (outtype == 2) {
        return(data.frame(cbind(n1 = n / (ratio + 1), n2 = ratio *
          n / (ratio + 1), Power = pwr)))
      }
      else if (outtype == 3) {
        return(data.frame(cbind(
          n = n, n1 = n / (ratio + 1),
          n2 = ratio * n / (ratio + 1), alpha = alpha,
          sided = sided, beta = 1 - pwr, Power = pwr,
          sigma0 = sigma0, sigma1 = sigma1, p1 = p1,
          p2 = p2, delta0 = delta0, p10 = p10, p20 = p20
        )))
      }
      else {
        return(Power = pwr)
      }
    }
  }
}

# simBinomial roxy [sinew] ----
#' @rdname varBinomial
#' @export
#' @importFrom stats rbinom
# simBinomial function [sinew] ----
simBinomial <- function(p1, p2, n1, n2, delta0 = 0, nsim = 10000, chisq = 0, adj = 0, scale = "Difference") {
  # check input arguments
  # Note: delta0, chisq, adj, and scale checked by testBinomial() function
  checkVector(p1, "numeric", c(0, 1), c(FALSE, FALSE))
  checkVector(p2, "numeric", c(0, 1), c(FALSE, FALSE))
  checkScalar(n1, "integer", c(1, Inf))
  checkScalar(n2, "integer", c(1, Inf))
  checkScalar(nsim, "integer", c(1, Inf))
  checkLengths(p1, p2)

  x1 <- stats::rbinom(prob = p1, size = n1, n = nsim)
  x2 <- stats::rbinom(prob = p2, size = n2, n = nsim)
  scale <- match.arg(tolower(scale), c("difference", "rr", "or", "lnor"))

  testBinomial(
    x1 = x1, x2 = x2, n1 = n1, n2 = n2, delta0 = delta0, adj = adj,
    chisq = chisq, scale = scale
  )
}

# testBinomial roxy [sinew] ----
#' @rdname varBinomial
#' @export
#' @importFrom stats rbinom
# testBinomial function [sinew] ----
testBinomial <- function(x1, x2, n1, n2, delta0 = 0, chisq = 0, adj = 0, scale = "Difference", tol = .1e-10) {
  # check input arguments
  checkVector(n1, "integer", c(1, Inf))
  checkVector(n2, "integer", c(1, Inf))
  checkVector(x1, "integer", c(0, Inf))
  checkVector(x2, "integer", c(0, Inf))
  checkVector(chisq, "integer", c(0, 1))
  checkVector(adj, "integer", c(0, 1))
  checkScalar(scale, "character")
  scale <- match.arg(tolower(scale), c("difference", "rr", "or", "lnor"))
  checkScalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  checkLengths(n1, n2, x1, x2, delta0, chisq, adj, allowSingle = TRUE)
  checkVector(n1 - x1, "integer", c(0, Inf))
  checkVector(n2 - x2, "integer", c(0, Inf))

  # make all vector arguments the same length (don't extend n1, n2)
  len <- max(sapply(list(n1, n2, x1, x2, chisq, adj, delta0), length))

  if (len > 1) {
    if (length(x1) == 1) p1 <- rep(x1, len)
    if (length(x2) == 1) p2 <- rep(x2, len)
    if (length(chisq) == 1) alpha <- rep(chisq, len)
    if (length(adj) == 1) beta <- rep(adj, len)
    if (length(delta0) == 1) delta0 <- rep(delta0, len)
  }

  ntot <- n1 + n2
  xtot <- x1 + x2

  # risk difference test - from Miettinen and Nurminen eqn (9)
  if (scale == "difference") {
    checkVector(delta0, "numeric", c(-1, 1), c(FALSE, FALSE))
    L2 <- (n1 + 2 * n2) * delta0 - ntot - xtot
    L1 <- (n2 * delta0 - ntot - 2 * x2) * delta0 + xtot
    L0 <- x2 * delta0 * (1 - delta0)
    q <- (L2 / (3 * ntot))^3 - L1 * L2 / 6 / ntot^2 + L0 / 2 / ntot
    p <- (sign(q) + (q == 0)) * sqrt((L2 / (3 * ntot))^2 - L1 / (3 * ntot))
    a <- q / p^3
    a[a > 1] <- 1
    a <- (pi + acos(a)) / 3
    R0 <- 2 * p * cos(a) - L2 / 3 / ntot
    R1 <- R0 + delta0
    V <- (R1 * (1 - R1) / n1 + R0 * (1 - R0) / n2)
    # V=0 only if no or all events and delta0=0
    # in which case test statistic will be 0
    V[V <= tol] <- 1
    # compute z-statistic
    z <- x1 / n1 - x2 / n2 - delta0
  }
  # relative risk test - from Miettinen and Nurminen eqn (10)
  else if (scale == "rr") {
    checkVector(delta0, "numeric", c(-Inf, Inf), c(FALSE, FALSE))
    delta0 <- exp(delta0) # value of 0 input represents equal rates
    A <- delta0 * ntot
    B <- -(n1 * delta0 + x1 + n2 + x2 * delta0)
    C <- xtot
    R0 <- (-B - sqrt(B^2 - 4 * A * C)) / 2 / A
    R1 <- R0 * delta0
    V <- R1 * (1 - R1) / n1 + delta0^2 * R0 * (1 - R0) / n2
    # V=0 only if no events or (all events and delta0=1)
    # in which case test statistic will be 0
    V[V <= 0 | is.na(V)] <- 1
    z <- x1 / n1 - x2 / n2 * delta0
  }
  # odds-ratio and log-odds-ratio
  else {
    checkVector(delta0, "numeric", c(-Inf, Inf), c(FALSE, FALSE))
    delta0 <- exp(delta0) # change from log odds-ratio to odds-ratio
    A <- n2 * (delta0 - 1)
    B <- n1 * delta0 + n2 - xtot * (delta0 - 1)
    C <- -xtot
    R0 <- (-B + sqrt(B^2 - 4 * A * C)) / 2 / A
    R1 <- R0 * delta0 / (1 + R0 * (delta0 - 1))
    # next 2 lines deal with length(delta0)>length(xtot) or length(ntot)
    xtem <- (delta0 == delta0) * xtot
    ntem <- (delta0 == delta0) * ntot
    R0[delta0 == 1] <- xtem[delta0 == 1] / ntem[delta0 == 1]
    R1[delta0 == 1] <- R0[delta0 == 1]

    # odds-ratio test - from Miettinen and Nurminen eqn (13)
    if (scale == "or") {
      V <- 1 / (1 / n1 / R1 / (1 - R1) + 1 / n2 / R0 / (1 - R0))
      V[xtot == 0 | xtot == ntot] <- 1
      z <- n1 * (x1 / n1 - R1)
    }
    # log-odds ratio - based on asymptotic distribution of log-odds
    # see vignette
    else if (scale == "lnor") {
      V <- 1 / n1 / R1 / (1 - R1) + 1 / n2 / R0 / (1 - R0)
      V[xtot == 0 | xtot == ntot] <- 1
      z <- log(x1 / (n1 - x1) / x2 * (n2 - x2) / delta0)
      z[xtot == 0 | xtot == ntot] <- 0
      z[xtot > 0 & xtot < ntot] <- z[xtot > 0 & xtot < ntot]
    }
  }

  # do continuity correction, where required
  one <- rep(TRUE, max(length(V), length(adj)))
  adj <- (adj == 1) & one
  V[adj] <- V[adj] * ntot / (ntot - 1)

  # square z where required to get chi-square
  one <- rep(TRUE, max(length(z), length(chisq)))
  chisq <- chisq & one
  if (length(z) == 1) z <- z * one
  z[chisq] <- z[chisq]^2 / V[chisq]
  z[ !chisq] <- z[ !chisq] / sqrt(V[ !chisq])

  z
}


# bpdiff function [sinew] ----
bpdiff <- function(delta, x1, x2, n1, n2, alpha = .05, adj = 0, scale = "Difference", lower.tail = FALSE) {
  stats::pnorm(testBinomial(x1, x2, n1, n2, delta0 = delta, adj = adj, scale = scale),
    lower.tail = lower.tail
  ) - alpha / 2
}
