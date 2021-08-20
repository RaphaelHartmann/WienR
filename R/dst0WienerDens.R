
#' Partial derivative of the first-passage time probability density function of the diffusion model with respect to the inter-trial variability of the non-decision time
#'
#' Calculates the partial derivative of the first-passage time probability density function of the diffusion model with respect to the inter-trial variability of the non-decision time st0.
#' @param t First-passage time. Numeric vector.
#' @param response Response boundary. Character vector with \code{"upper"} and \code{"lower"} as possible values. Alternatively a numeric vector with
#'   \code{1}=lower and \code{2}=upper.
#' @param a Upper barrier. Numeric vector.
#' @param v Drift rate. Numeric vector.
#' @param w Relative starting point. Numeric vector.
#' @param t0 Non-decision time. Numeric vector
#' @param sv Inter-trial variability of drift rate. Numeric vector. Standard deviation of a normal distribution \code{N(v, sv)}.
#' @param sw Inter-trial variability of relative starting point. Numeric vector. Range of uniform distribution \code{U(w-0.5*sw, w+0.5*sw)}.
#' @param st0 Inter-trial variability of non-decision time. Numeric vector. Range of uniform distribution \code{U(t0, t0+st0)}.
#' @param precision Optional numeric value. Precision of the partial derivative. Numeric value. Default is \code{NULL}, which takes default value 1e-12.
#' @param K Optional numeric value. Number of iterations to calculate the infinite sums. Numeric value (integer). Default is \code{NULL}.
#'   \itemize{
#'     \item \code{precision = NULL} and \code{K = NULL}: Default \code{precision = 1e-12} used to calculate internal K.
#'     \item \code{precision != NULL} and \code{K = NULL}: \code{precision} is used to calculate internal K,
#'     \item \code{precision = NULL} and \code{K != NULL}: \code{K} is used as internal K,
#'     \item \code{precision != NULL} and \code{K != NULL}: if internal K calculated through \code{precision} is smaller than \code{K}, \code{K} is used.
#'   }
#'   We recommend using either default (\code{precision = K = NULL}) or only \code{precision}.
#' @param n.threads Optional numerical or logical value. Number of threads to use. If not provided (or 1 or \code{FALSE}) parallelization is not used. If set to \code{TRUE} then all available threads are used.
#' @return A list of the class \code{Diffusion_deriv} containing
#'   \itemize{
#'     \item \code{deriv}: the derivatives of the PDF with respect to a,
#'     \item \code{call}: the function call,
#'     \item \code{err}: the absolute error.
#'   }
#' @examples
#' dst0WienerPDF(t = 1.2, response = "upper", a = 1.1, v = 13, w = .6, st0 = .2)
#' @author Raphael Hartmann
#' @useDynLib "WienR", .registration=TRUE
#' @export
dst0WienerPDF <- function(t,
                         response,
                         a,
                         v,
                         w,
                         t0 = 0,
                         sv = 0,
                         sw = 0,
                         st0,
                         precision = NULL,
                         K = NULL,
                         n.threads = FALSE) {




  # ---- VALUE CHECKS ---- #

  # general checks
  lengths <- c(length(t), length(response), length(a), length(v), length(w), length(t0), length(sv), length(sw), length(st0))
  max_len <- max(lengths)
  if(any(lengths != max_len & lengths != 1)) stop("t, response, a, v, w, t0, sv, sw, and st0 must have same length (except length one)")
  if(length(t) != max_len) t <- rep(t, max_len)
  if(length(response) != max_len) response <- rep(response, max_len)
  if(length(a) != max_len) a <- rep(a, max_len)
  if(length(v) != max_len) v <- rep(v, max_len)
  if(length(w) != max_len) w <- rep(w, max_len)
  if(length(t0) != max_len) t0 <- rep(t0, max_len)
  if(length(sv) != max_len) sv <- rep(sv, max_len)
  if(length(sw) != max_len) sw <- rep(sw, max_len)
  if(length(st0) != max_len) st0 <- rep(st0, max_len)

  # t a v w t0 sw sv st0 checks
  if(!is.numeric(t) | !is.numeric(a) | !is.numeric(v) | !is.numeric(w) | !is.numeric(t0) | !is.numeric(sv) | !is.numeric(sw) | !is.numeric(st0)) stop("t, a, v, w, t0, sv, sw, and st0 must be numeric")
  if(any(t <= 0) | any(a <= 0) | any(w <= 0)) stop("t, a, and w must be strictly positive")
  if(any(t0 < 0) | any(sw < 0) | any(sv < 0) | any(st0 < 0)) stop("t0, sw, sv, and st0 must be positive or zero")
  if(any(w >= 1)) stop("w must be lower than one")
  if(any(w-0.5*sw <= 0) | any(w+0.5*sw >= 1)) stop("w-0.5*sw must be greater than zero and w+0.5*sw must be lower than one")
  if(any(st0<=0)) stop("st0 must be larger than zero for this partial derivative (d/dst0)")

  # response checks
  if(!is.character(response) & !is.numeric(response)) stop("response must be a character with the values \"upper\" and/or \"lower\" OR numerics with the values 1=\"lower\" or 2=\"upper\"")
  if(!all(response %in% c("upper", "lower")) & !all(response %in% c(1,2)) ) stop("response cannot include values other than \"upper\" and/or \"lower\" OR 1=\"lower\" or 2=\"upper\"")
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

  out <- .Call("dDiffusion7",
               as.numeric(t),
               as.numeric(a),
               as.numeric(v),
               as.numeric(t0),
               as.numeric(w),
               as.numeric(sw),
               as.numeric(sv),
               as.numeric(st0),
               as.numeric(precision),
               as.integer(resps),
               as.integer(K),
               as.integer(max_len),
               as.integer(n.threads),
               as.integer(7),
               as.logical(PRECISION_FLAG)
  )



  #print(out)

  derivative <- list(deriv = out$deriv, call = match.call(), err = out$err)

  # output
  class(derivative) <- "Diffusion_deriv"
  return(derivative)

}
