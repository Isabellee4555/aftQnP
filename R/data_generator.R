#' data generator
#' @import dplyr
#' @noRd
data_generator <- function(n=1000, beta=c(1, -0.3, 0.5), gamma=c( 2, -2, 1.5), evt_p = 0.5, seed_ = NULL)
{
  delta <- deltaL <- deltaR <- NULL
  y <- y_L <- y_R <- deltaI <- NULL
  if(!is.null(seed_)){set.seed(seed_)}
  X1 <- rbinom(n, 1, 0.5)
  X2 <- runif(n, 0, 3)
  X3 <- runif(n, 0, 5)
  X <- cbind(X1, X2, X3)

  Z1 <- rep(1,n)
  Z2 <- X2
  Z3 <- X3
  Z <- cbind(Z1, Z2, Z3)

  kappa <- rweibull(n, shape=3, scale=1)
  y_true <- as.vector(kappa * exp(X %*% beta))

  alpha <- as.vector(exp(Z %*% gamma) / (1 + exp(Z %*% gamma)))
  cured <-  ifelse(runif(n) > alpha, 1, 0)

  width <-  0.8
  kappa_lower <- kappa - sapply(kappa, function(k) runif(1, 0, min(k, width)))
  kappa_upper <- kappa + runif(length(kappa), 0, width)
  kappa_lower[cured == 1] <- kappa_upper[cured == 1] <- 2.5 + runif(sum(cured), -1, 0)#qweibull(0.9999, 3, 1)
  y_lower <- as.vector(kappa_lower * exp(X %*% beta))
  y_upper <- as.vector(kappa_upper * exp(X %*% beta))

  df <- data.frame(cbind(X, Z, y_true, alpha, cured, kappa, kappa_lower, kappa_upper, y_lower, y_upper))
  df <- df[order(cured), ]

  n_event <- n - sum(cured)
  props <- c(evt_p, rep((1 - evt_p)/3,3))
  breaks <- round(cumsum(props) * n_event)

  df$delta <- 0
  if(breaks[1] != 0){
    df$delta[1:breaks[1]] <- 1
  }

  df$deltaL <- 0
  df$deltaL[(breaks[1] + 1):breaks[2]] <- 1

  df$deltaR <- 0
  df$deltaR[(breaks[2] + 1):breaks[3]] <- 1
  df$deltaR[(breaks[4] + 1):n] <- 1

  df$deltaI <- 0
  df$deltaI[(breaks[3] + 1):breaks[4]] <- 1


  df$y <- 0
  df$y[df$delta == 1] <- df$y_true[df$delta == 1]
  df$y[df$deltaL == 1] <- df$y_upper[df$deltaL == 1]
  df$y[df$deltaR == 1] <- df$y_lower[df$deltaR == 1]

  df$y_L <- 0
  df$y_L[df$deltaI == 1] <- df$y_lower[df$deltaI == 1]
  df$y_R <- 0
  df$y_R[df$deltaI == 1] <- df$y_upper[df$deltaI == 1]

  data <-  tibble(df) %>% mutate(y_L = ifelse(delta == 1 | deltaR==1, y, y_L),
                               y_R = ifelse(delta == 1 | deltaL==1, y, y_R)) %>% select(X1, X2, X3, Z1,Z2,Z3,
                                                                                        delta, deltaL, deltaR, deltaI,
                                                                                        y_L, y_R)

  return(data)
}
