
#' Partial derivative of the Wiener diffusion cumulative distribution function with respect to the relative starting point
#'
#' Calculate the partial derivative of the Wiener diffusion cumulative distribution function with respect to the relative starting point w.
#' @param t First-passage time. Numeric vector.
#' @param response Response boundary. Character vector with \code{"upper"} and \code{"lower"} as possible values. Alternatively a numeric vector with
#'   \code{1}=lower and \code{2}=upper.
#' @param a Upper barrier. Numeric vector.
#' @param v Drift rate. Numeric vector.
#' @param w Relative starting point. Numeric vector.
#' @param precision Optional numeric value. Precision of the partial derivative. Numeric value. Default is \code{NULL}, which takes default value 1e-12.
#' @param K Optional. Number of iterations to calculate the infinite sums. Numeric value (integer). Default is \code{NULL}.
#'   \itemize{
#'     \item \code{precision = NULL} and \code{K = NULL}: Default \code{precision = 1e-12} used to calculate internal K.
#'     \item \code{precision != NULL} and \code{K = NULL}: \code{precision} is used to calculate internal K,
#'     \item \code{precision = NULL} and \code{K != NULL}: \code{K} is used as internal K,
#'     \item \code{precision != NULL} and \code{K != NULL}: if internal K calculated through \code{precision} is smaller than \code{K}, \code{K} is used.
#'   }
#'   We recommend using either default (\code{precision = K = NULL}) or only \code{precision}.
#' @param n.threads Optional numerical or logical value. Number of threads to use. If not provided (or 1 or \code{FALSE}) parallelization is not used. If set to \code{TRUE} then all available threads are used.
#' @return A list of the class \code{Wiener_deriv} containing
#'   \itemize{
#'     \item \code{deriv}: the derivatives of the CDF with respect to w,
#'     \item \code{derivln}: the derivatives of the log-transformed CDF with respect to w,
#'     \item \code{call}: the function call.
#'   }
#' @examples
#' dwWienerCDF(t = 1.2, response = "upper", a = 1.1, v = 13, w = .6, precision = NULL, K = NULL)
#' @author Raphael Hartmann
#' @useDynLib "WienR", .registration=TRUE
#' @export
dwWienerCDF <- function(t,
                        response,
                        a,
                        v,
                        w,
                        precision = NULL,
                        K = NULL,
                        n.threads = FALSE) {




  # ---- VALUE CHECKS ---- #

  # general checks
  lengths <- c(length(t), length(response), length(a), length(v), length(w))
  max_len <- max(lengths)
  if(any(lengths != max_len & lengths != 1)) stop("t, response, a, v, and w must have same length (except length one)")
  if(length(t) != max_len) t <- rep(t, max_len)
  if(length(response) != max_len) response <- rep(response, max_len)
  if(length(a) != max_len) a <- rep(a, max_len)
  if(length(v) != max_len) v <- rep(v, max_len)
  if(length(w) != max_len) w <- rep(w, max_len)

  # t a v w checks
  if(!is.numeric(t) | !is.numeric(a) | !is.numeric(v) | !is.numeric(w)) stop("t, a, v, and w must be numeric")
  if(any(t <= 0) | any(a <= 0) | any(w <= 0)) stop("t, a, and w must be positive")
  if(any(w >= 1)) stop("w must be lower than one")

  # response checks
  if(!is.character(response) & !is.numeric(response)) stop("response must be a character with the values \"upper\" and/or \"lower\" OR numerics with the values 1=\"lower\" or 2=\"upper\"")
  if(!all(response %in% c("upper", "lower")) & !all(response %in% c(1,2)) ) stop("response must cannot include values other than \"upper\" and/or \"lower\" OR 1=\"lower\" or 2=\"upper\"")
  resps <- ifelse(response == "lower" | response == 1, 0, 1)

  # K checks
  if(!is.numeric(K) & !is.null(K)) stop("K must either be NULL or some numeric value")
  if(!is.null(K)) {
    if(length(K)!=1) stop("K must be of length one")
    if(K %% 1 != 0) stop("K must be an integer") else K <- as.integer(round(K))
  }

  # precision checks
  if(!is.numeric(precision) & !is.null(precision)) stop("precision must either be NULL or some numeric value")
  if(length(precision)!=1 & !is.null(precision)) stop("precision must be of length one")

  PRECISION_FLAG <- TRUE
  if(is.null(precision)) PRECISION_FLAG <- FALSE

  if(is.null(K)) K <- 0
  if(is.null(precision)) precision <- 0

  # thread checks
  if(!is.numeric(n.threads) & !is.logical(n.threads)) stop("n.threads must either be numerical or logical")
  if(is.numeric(n.threads)) if(n.threads %% 1 != 0) stop("n.threads must be an integer") else n.threads <- as.integer(n.threads)
  if(is.logical(n.threads)) n.threads <- ifelse(n.threads == TRUE, 99999, 0)
  if(n.threads < 2) n.threads <- 0



  # --- C++ FUNCTION CALL ---- #

  out <- .Call("dwpWiener",
               as.numeric(t),
               as.numeric(a),
               as.numeric(v),
               as.numeric(w),
               as.numeric(precision),
               as.integer(resps),
               as.integer(K),
               as.integer(max_len),
               as.integer(n.threads),
               as.logical(PRECISION_FLAG)
  )

  #print(out)

  derivative <- list(deriv = out$deriv, derivln = out$deriv_ln, call = match.call())

  # output
  class(derivative) <- "Wiener_deriv"
  return(derivative)

}
