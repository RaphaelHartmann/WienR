#' @export
print.Wiener_deriv <- function(x, ...) {
  der <- c("dt", "da", "dv", "dw", "grad")
  cll <- as.character(x$call[1])
  fnc <- ifelse(grepl("PDF", cll), "PDF", "CDF")
  der <- ifelse(grepl("dt", cll), "t", 
                ifelse(grepl("da", cll), "a", 
                       ifelse(grepl("dv", cll), "v", 
                              ifelse(grepl("dw", cll), "w", "all params")))
  )
  
  grepl("grad", cll)
  cat(paste0("\nDerivation of ", fnc, " with respect to ", der, "\n\n"))
  print(x$deriv)
  cat("\n---------------------------\n")
  cat(paste0("\nDerivation of log(", fnc, ") with respect to ", der, "\n\n"))
  print(x$derivln)
  cat("\n---------------------------\n\n")
}

#' @export
print.Wiener_pdf <- function(x, ...) {
  cat(paste0("\nWiener PDF\n\n"))
  print(x$value)
  cat("\n---------------------------\n")
  cat(paste0("\nlog of Wiener PDF\n\n"))
  print(x$logvalue)
  cat("\n---------------------------\n\n")
}

#' @export
print.Wiener_cdf <- function(x, ...) {
  cat(paste0("\nWiener CDF\n\n"))
  print(x$value)
  cat("\n---------------------------\n")
  cat(paste0("\nlog of Wiener CDF\n\n"))
  print(x$logvalue)
  cat("\n---------------------------\n\n")
}