#' Predicted Survival Function Plot
#' @description plot for predicted survival probabilities from the semi-parametric AFT mixture cure model.
#' @param fit a model fitted by the \code{aftsur} function.
#' @param x a vector of latency covariates that is used to estimate survival probabilities. If \code{x = NULL}, the baseline survival will be estimated.
#' @export
#' @import ggplot2
plot_sur <- function(fit, x = NULL){
  X_base <-  if(is.null(x)){
    rep(0, length(fit$beta))
  }else{x}
  const <- exp(-X_base%*%fit$beta)
  max_event_t <- max(fit$k[fit$delta==1]* c(exp(X_base%*%fit$beta)))
  t_grid <- seq(0, max_event_t, 0.01)
  H_vals <-  H0(t_grid%*%const, fit$theta, fit$basis_params)
  S_vals <- exp(-H_vals)
  variance_H <- diag(compute_cumu_hazard_covariance(fit, compute_covariance_matrix(fit,return_F = T), t_grid%*%const))
  se_S <- S_vals^2 * variance_H
  lty_ <- colors <- c("Survival" = 1, "95% ASY.CI"=2)
  colors <- c("Survival" = "red", "95% ASY.CI"="grey30")
  fills <- c("Survival" = "red", "95% ASY.CI"="grey30")
  p <- ggplot(data = data.frame(x = t_grid))+
  geom_line(aes(x, y  = S_vals, col = "Survival", lty = "Survival")) +
  geom_line(aes(x, y  = S_vals + 1.96*sqrt(se_S), col = "95% ASY.CI", lty = "95% ASY.CI")) +
  geom_line(aes(x, y  = S_vals - 1.96*sqrt(se_S), col = "95% ASY.CI", lty = "95% ASY.CI")) +
  geom_ribbon(aes(x = x, ymin = pmax(S_vals - 1.96*sqrt(se_S), 0),ymax = pmin(S_vals + 1.96*sqrt(se_S), 1)),
              fill= "Orange", alpha = 0.1)+
  theme_minimal() + xlab("Survival time (t)") + ylab(expression("Survival prob."~ " "~S[i](t))) +
  scale_x_continuous(breaks = seq(round(max_event_t/5, ifelse(max_event_t>=3, 0, 1)),
                                  max_event_t,round(max_event_t/5, ifelse(max_event_t>=3, 0, 1))),
                     limits = c(0, max_event_t), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, max(S_vals), round(max(S_vals)/5, 1)),
                     limits = c(0, max(S_vals)), expand = c(0, 0)) +
  scale_color_manual(name = "same", values = colors)  +
  scale_linetype_manual(name = "same",values = lty_) +
  theme(legend.position = c(0.8, 0.75),
        legend.title = element_blank(),
        legend.background = element_rect(fill="white",color = "white"),
        legend.text = element_text(face="bold", size = 14,family = "serif"),
        axis.title.x = element_text( size = 14,family = "serif"),
        axis.title.y = element_text( size = 14,family = "serif"),
        plot.title = element_text(face = "bold", size = 16, hjust = 0, family = "serif"))
  suppressWarnings(print(p))
}

