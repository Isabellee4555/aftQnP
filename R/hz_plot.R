#' Estimated Hazard Plot
#'
#' @param fit a aftsur object
#' @export


hz_plot <- function(fit){
  plot(x = fit$predicth$x, y = fit$predicth$h, type = "l",
       xlim = c(0, max(fit$predicth$x)), ylim = c(0, max(fit$predicth$h)),
       main = "Baseline Hazard Estimates",ylab = expression(h[0](t)), xlab = expression(kappa))
}





