#' Plot Predicted Hazard Function
#' @description a function that plots predicted hazard rates from the semi-parametric AFT mixture cure model.
#' @param fit an object of class "aftsur". A model fitted by the \code{aftsur} function.
#' @param x a vector of latency covariates that is used to estimate hazard rates. If \code{x = NULL}, the baseline hazard will be predicted.
#' @return a \code{ggplot2} object
#' @examples
#' require(survival)
#' # load data
#' data("simu_data")
#' # create a formula using a Surv
#' formula_aft <- Surv(y_L, y_R, type = "interval2") ~ X1 + X2 - 1
#' # fit a model
#' aft_fit <- aftsur(formula = formula_aft, cure_var = ~ Z1 + Z2 + Z3, offset = TRUE, data = simu_data)
#' plot_hz(aft_fit)
#' @export
#' @import ggplot2
plot_hz <- function(fit, x = NULL){
  X_base <-  if(is.null(x)){
    rep(0, length(fit$beta))
  }else{x}
  const <- exp(-X_base%*%fit$beta)
  max_event_t <- max(fit$k[fit$delta==1]* c(exp(X_base%*%fit$beta)))
  t_grid <- seq(0, max_event_t, 0.01)
  h_vals <- h0(t_grid%*%const, fit$theta, fit$basis_params)*c(const)
  asy_h_sd <- sqrt(diag(compute_hazard_covariance(fit, compute_covariance_matrix(fit,return_F = T), t_grid%*%const)))
  lty_ <- colors <- c("Hazard" = 1, "95% ASY.CI"=2)
  colors <- c("Hazard" = "red", "95% ASY.CI"="grey30")
  fills <- c("Hazard" = "red", "95% ASY.CI"="grey30")
  p <- ggplot(data = data.frame(x = t_grid))+
    geom_line(aes(x, y  = h_vals, col = "Hazard", lty = "Hazard")) +
    geom_line(aes(x, y  = h_vals + 1.959964*asy_h_sd, col = "95% ASY.CI", lty = "95% ASY.CI")) +
    geom_line(aes(x, y  = h_vals - 1.959964*asy_h_sd, col = "95% ASY.CI", lty = "95% ASY.CI")) +
    geom_ribbon(aes(x = x, ymin = pmax(h_vals - 1.959964*asy_h_sd, 0),ymax = pmin(h_vals + 1.959964*asy_h_sd, max(h_vals))),
                fill= "Orange", alpha = 0.1)+
    theme_minimal() + xlab("Survival time (t)") + ylab(expression("Hazard rate"~ " "~h[i](t))) +
    scale_x_continuous(breaks = seq(round(max_event_t/5, ifelse(max_event_t>=3, 0, 1)),
                                    max_event_t,round(max_event_t/5, ifelse(max_event_t>=3, 0, 1))),
                       limits = c(0, max_event_t), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, max(h_vals), round(max(h_vals)/5, ifelse(max(h_vals)>=3, 0, 1))),
                       limits = c(0, max(h_vals)), expand = c(0, 0)) +
    scale_color_manual(name = "same", values = colors)  +
    scale_linetype_manual(name = "same",values = lty_) +
    theme(legend.position = c(0.2, 0.75),
          legend.title = element_blank(),
          legend.background = element_rect(fill="white",color = "white"),
          legend.text = element_text(face="bold", size = 14,family = "serif"),
          axis.title.x = element_text( size = 14,family = "serif"),
          axis.title.y = element_text( size = 14,family = "serif"),
          plot.title = element_text(face = "bold", size = 16, hjust = 0, family = "serif"))
  suppressWarnings(print(p))
}

