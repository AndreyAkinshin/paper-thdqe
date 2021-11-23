library(ggplot2)
library(evd)
library(latex2exp)

source("thdqe.R")

draw.beta <- function(trimmed = F) {
  n <- 10
  p <- 0.5
  a <- (n + 1) * p
  b <- (n + 1) * (1 - p)

  if (trimmed) {
    hdi <- getBetaHdi(a, b, 1 / sqrt(n))
    L <- hdi[1]
    R <- hdi[2]
    filename <- "tbeta"
    title <- paste0("Truncated Beta(", a, ", ", b, ") PDF")
  } else {
    L <- 0
    R <- 1
    filename <- "beta"
    title <- paste0("Beta(", a, ", ", b, ") PDF")
  }
  
  scale <- 1 / (pbeta(R, a, b) - pbeta(L, a, b))
  step <- 0.001
  x <- c(L, seq(L, R, by = step), R)
  y <- c(0, dbeta(seq(L, R, by = step), a, b), 0) * scale
  df <- data.frame(x, y)
  
  x.segm <- 1:(n-1)/n
  x.segm <- x.segm[x.segm >= L & x.segm <= R]
  df.segm <- data.frame(
    x1 = x.segm,
    y1 = 0,
    x2 = x.segm,
    y2 = dbeta(x.segm, a, b) * scale
  )
  
  wx <- (c(x.segm, R) + c(L, x.segm)) / 2
  wy <- dbeta(wx, a, b) * scale / 2
  wy[wy < 0.1] <- wy[wy < 0.1] + 0.2
  wi <- floor(wx * 10) + 1
  wt <- paste0("W", wi, "")
  df.w <- data.frame(x = wx, y = wy, text = wt)
  
  p <- ggplot(df, aes(x, y)) +
    geom_line() +
    geom_segment(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      data = df.segm,
      linetype = "dashed") +
    geom_text(data = df.w, mapping = aes(x, y, label = text)) +
    scale_x_continuous(breaks = 0:10/10, limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = "density", title = title) +
    theme_bw()
  
  show(p)
}

draw.hdi <- function(sum = 10, width = 0.3, ps = c(0.05, 0.5, 0.95)) {
  build.df <- function(sum, width, p, comment) {
    a <- sum * p
    b <- sum * (1 - p)
    hdi <- getBetaHdi(a, b, width)
    
    x <- seq(0, 1, by = 0.01)
    y <- dbeta(x, a, b)
    if (y[1] > 6 || y[length(y)] > 6) {
      x <- seq(0.01, 0.99, by = 0.01)
      y <- dbeta(x, a, b)
    }
    
    max.value <- ifelse(y[1] < y[2] & y[length(y) - 1] > y[length(y)], max(y), 6)
    visible <- y <= max.value
    inside <- x >= hdi[1] & x <= hdi[2]
    y <- pmin(y, max.value)
    data.frame(x, y, p, visible, inside, comment)
  }

  df <- rbind(
    build.df(sum, width, 0.05, "L"),
    build.df(sum, width, 0.50, "M"),
    build.df(sum, width, 0.95, "R")
  )
  df$comment <- factor(df$comment, levels = c("L", "M", "R"))
  levels(df$comment) <- c(
    L = TeX("\\textit{(Left border)} $\\alpha = 0.5,\\; \\beta = 9.5$"),
    M = TeX("\\textit{(Middle)} $\\alpha = 5,\\; \\beta = 5$"),
    R = TeX("\\textit{(Right border)} $\\alpha = 9.5,\\; \\beta = 0.5$")
  )
  p <- ggplot(df, aes(x, y)) +
    geom_area(data = df[df$inside,], aes(x, y), fill = "#999999", alpha = 0.4) +
    geom_line(data = df[df$visible,], aes(x, y)) +
    facet_wrap(vars(comment), labeller=label_parsed) +
    labs(
      title = TeX(paste0("HDI of beta distribution ($\\alpha + \\beta = ", sum, "$, interval width = $", width, "$)")),
      y = "density"
    ) +
    theme_bw()
  show(p)
}

# Contaminated normal distribution
rcnorm <- function(n, mean = 0, sd = 1, eps = 0.01, c = 1000000)
  ifelse(runif(n) > eps, rnorm(n, mean, sd), rnorm(n, mean, sd * sqrt(c)))

# Quantile estimators
hf7qe <- function(x, probs) as.numeric(quantile(x, probs))
hdqe <- function(x, probs) as.numeric(hdquantile(x, probs))
thdqe <- function(x, probs) as.numeric(thdquantile(x, probs))

simulation1 <- function() {
  gen <- function() {
    x <- rcnorm(7)
    c(
      hf7 = hf7qe(x, 0.5),
      hd = hdqe(x, 0.5),
      thd = thdqe(x, 0.5))
  }
  set.seed(1729)
  df.raw <- data.frame(t(replicate(10000, gen())))
  probs <- c(seq(0, 0.05, by = 0.01), seq(0.95, 1, by = 0.01))
  df <- data.frame(
    quantile = probs,
    HF7 = quantile(df.raw$hf7, probs),
    HD = quantile(df.raw$hd, probs),
    "THD-SQRT" = quantile(df.raw$thd, probs)
  )
  rownames(df) <- c()
  print(kable(df, col.names = c("quantile", "HF7", "HD", "THD-SQRT")))
}

simulation2 <- function() {
  calc.efficiency <- function(distribution, n, p) {
    true.value <- distribution$q(p)
    calc.mse <- function() {
      df <- data.frame(t(replicate(200, {
        x <- distribution$r(n)
        c(hf7 = (hf7qe(x, p) - true.value)^2,
          hd = (hdqe(x, p) - true.value)^2,
          thd = (thdqe(x, p) - true.value)^2)
      })))
      c(
        hf7 = mean(df$hf7),
        hd = mean(df$hd),
        thd = mean(df$thd)
      )
    }
    df.mse <- data.frame(t(replicate(101, calc.mse())))
    list(
      distribution = distribution$title,
      p = p,
      n = n,
      hd = median(df.mse$hf7) / median(df.mse$hd),
      thd = median(df.mse$hf7) / median(df.mse$thd)
    )
  }
  
  d.unif <- list(title = "Uniform(a=0, b=1)", r = runif, q = qunif)
  d.tri_0_2_1 <- list(title = "Triangular(a=0, b=2, c=1)", r = function(n) rtri(n, 0, 2, 1), q = function(p) qtri(p, 0, 2, 1))
  d.tri_0_2_02 <- list(title = "Triangular(a=0, b=2, c=0.2)", r = function(n) rtri(n, 0, 2, 0.2), q = function(p) qtri(p, 0, 2, 0.2))
  d.beta2_4 <- list(title = "Beta(a=2, b=4)", r = function(n) rbeta(n, 2, 4), q = function(p) qbeta(p, 2, 4))
  d.beta2_10 <- list(title = "Beta(a=2, b=10)", r = function(n) rbeta(n, 2, 10), q = function(p) qbeta(p, 2, 10))
  
  d.norm <- list(title = "Normal(m=0, sd=1)", r = rnorm, q = qnorm)
  d.weibull1_2 <- list(title = "Weibull(scale=1, shape=2)", r = function(n) rweibull(n, 2), q = function(p) qweibull(p, 2))
  d.student3 <- list(title = "Student(df=3)", r = function(n) rt(n, 3), q = function(p) qt(p, 3))
  d.gumbel <- list(title = "Gumbel(loc=0, scale=1)", r = rgumbel, q = qgumbel)
  d.exp <- list(title = "Exp(rate=1)", r = rexp, q = qexp)
  
  d.cauchy <- list(title = "Cauchy(x0=0, gamma=1)", r = rcauchy, q = qcauchy)
  d.pareto1_05 <- list(title = "Pareto(loc=1, shape=0.5)", r = function(n) rpareto(n, 1, 0.5), q = function(p) qpareto(p, 1, 0.5))
  d.pareto1_2 <- list(title = "Pareto(loc=1, shape=2)", r = function(n) rpareto(n, 1, 2), q = function(p) qpareto(p, 1, 2))
  d.lnorm0_1 <- list(title = "LogNormal(mlog=0, sdlog=1)", r = function(n) rlnorm(n, 0, 1), q = function(p) qlnorm(p, 0, 1))
  d.lnorm0_2 <- list(title = "LogNormal(mlog=0, sdlog=2)", r = function(n) rlnorm(n, 0, 2), q = function(p) qlnorm(p, 0, 2))
  
  d.lnorm0_3 <- list(title = "LogNormal(mlog=0, sdlog=3)", r = function(n) rlnorm(n, 0, 3), q = function(p) qlnorm(p, 0, 3))
  d.weibull1_03 <- list(title = "Weibull(shape=0.3)", r = function(n) rweibull(n, 0.3), q = function(p) qweibull(p, 0.3))
  d.weibull1_05 <- list(title = "Weibull(shape=0.5)", r = function(n) rweibull(n, 0.5), q = function(p) qweibull(p, 0.5))
  d.frechet1 <- list(title = "Frechet(shape=1)", r = function(n) rfrechet(n, shape = 1), q = function(p) qfrechet(p, shape = 1))
  d.frechet3 <- list(title = "Frechet(shape=3)", r = function(n) rfrechet(n, shape = 3), q = function(p) qfrechet(p, shape = 3))
  
  ds <- list(
    d.unif, d.tri_0_2_1, d.tri_0_2_02, d.beta2_4, d.beta2_10,
    d.norm, d.weibull1_2, d.student3, d.gumbel, d.exp,
    d.cauchy, d.pareto1_05, d.pareto1_2, d.lnorm0_1, d.lnorm0_2,
    d.lnorm0_3, d.weibull1_03, d.weibull1_05, d.frechet1, d.frechet3
  )
  ns <- c(5) # c(3, 5, 10, 20, 40)
  ps <- seq(0.01, 0.99, by = 0.01)
  input <- expand.grid(d = ds, n = ns, p = ps)
  efficienyThreshold <- 3
  start.time <- Sys.time()
  if (!file.exists("efficiency.csv")) {
    df <- do.call("rbind",
                  lapply(1:nrow(input),
                         function(i) calc.efficiency(input$d[i][[1]], input$n[i], input$p[i])))
    df <- data.frame(df)
    df <- df %>% gather("estimator", "efficiency", -names(df)[1:3])
    df$distribution <- unlist(df$distribution)
    df$distribution <- factor(df$distribution, levels = unique(df$distribution))
    df$n <- unlist(df$n)
    df$p <- unlist(df$p)
    df$estimator <- factor(df$estimator, levels = c("hd", "thd"))
    df$efficiency <- unlist(df$efficiency)
    df$efficiency <- pmin(df$efficiency, efficienyThreshold)
    
    write.csv(df, "efficiency.csv", row.names = F)
  } else {
    df <- read.csv("efficiency.csv")
    df$distribution <- factor(df$distribution, levels = unique(df$distribution))
    df$estimator <- factor(df$estimator, levels = c("hd", "thd"))
  }
  end.time <- Sys.time()
  
  draw <- function(n) {
    p <- ggplot(df[df$n == n,], aes(x = p, y = efficiency, col = estimator)) +
      facet_wrap(vars(distribution), ncol = 5) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      xlim(0, 1) +
      ylim(0, efficienyThreshold) +
      scale_color_manual(values = c("#D55E00", "#56B4E9"), labels = c("HD", "THD-SQRT")) +
      geom_line(size = 0.5) +
      labs(
        title = paste0("Relative efficiency of quantile estimators (n = ", n, ")"),
        x = "Quantile",
        y = "Statistical efficiency",
        col = "Estimator"
      ) +
      theme_bw() +
      theme(legend.position = "bottom", text = element_text(size = 8))
    show(p)
  }
  draw
  # print(end.time - start.time)
}