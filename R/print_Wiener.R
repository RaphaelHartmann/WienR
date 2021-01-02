#' @export
print.Diffusion_deriv <- function(x, ...) {
  der <- c("dt", "da", "dv", "dw", "grad")
  cll <- as.character(x$call[1])
  fnc <- ifelse(grepl("PDF", cll), "PDF", "CDF")
  der <- ifelse(grepl("dt", cll), "t", 
                ifelse(grepl("da", cll), "a", 
                       ifelse(grepl("dv", cll), "v", 
                              ifelse(grepl("dw", cll), "w", "all params")))
  )
  
  grepl("grad", cll)
  cat(paste0("\nPartial derivative of ", fnc, " with respect to ", der, "\n\n"))
  print(x$deriv)
  cat("\n---------------------------\n")
}

#' @export
print.Diffusion_pdf <- function(x, ...) {
  cat(paste0("\nFirst-passage time PDF\n\n"))
  print(x$value)
  cat("\n---------------------------\n")
  cat(paste0("\nLog of first-passage time PDF\n\n"))
  print(x$logvalue)
  cat("\n---------------------------\n\n")
}

#' @export
print.Diffusion_cdf <- function(x, ...) {
  cat(paste0("\nFirst-passage time CDF\n\n"))
  print(x$value)
  cat("\n---------------------------\n")
  cat(paste0("\nLog of first-passage time CDF\n\n"))
  print(x$logvalue)
  cat("\n---------------------------\n\n")
}