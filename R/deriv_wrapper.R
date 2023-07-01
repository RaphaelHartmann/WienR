
#' Wrapper function for the partial derivative of the first-passage time probability density function of the diffusion model
#' 
#' Calculates the partial derivative of the first-passage time probability density function of the diffusion model with respect to one of t, a, v, w, t0, sv, sw, or st0, or calculate the gradient.
#' @param wrt partial derivative w.r.t. one of the following:
#'   \itemize{
#'     \item \code{"t"} the first-passage time,
#'     \item \code{"a"} the upper barrier,
#'     \item \code{"v"} the drift rate,
#'     \item \code{"w"} the relative starting point,
#'     \item \code{"t0"} the non-decision time,
#'     \item \code{"sv"} the inter-trial variability of drift rate,
#'     \item \code{"sw"} the inter-trial variability of relative starting point,
#'     \item \code{"st0"} the inter-trial variability of non-decision time, or
#'     \item \code{"grad"} all the above but t.
#'   }
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
#' @param K Optional. Number of iterations to calculate the infinite sums. Numeric value (integer). Default is \code{NULL}.
#'   \itemize{
#'     \item \code{precision = NULL} and \code{K = NULL}: Default \code{precision = 1e-12} used to calculate internal K.
#'     \item \code{precision != NULL} and \code{K = NULL}: \code{precision} is used to calculate internal K,
#'     \item \code{precision = NULL} and \code{K != NULL}: \code{K} is used as internal K,
#'     \item \code{precision != NULL} and \code{K != NULL}: if internal K calculated through \code{precision} is smaller than \code{K}, \code{K} is used.
#'   }
#'   We recommend using either default (\code{precision = K = NULL}) or only \code{precision}.
#' @param n.threads Optional numerical or logical value. Number of threads to use. If not provided (or 1 or \code{FALSE}) parallelization is not used. If set to \code{TRUE} then all available threads are used.
#' @param n.evals Optional. Number of maximal function evaluations in the numeric integral if sv, sw, and/or st0 are not zero. Default is \code{6000} and \code{0} implies no limit and the 
#'   numeric integration goes on until the specified \code{precision} is guaranteed.
#' @return A list of the class \code{Diffusion_deriv} containing
#'   \itemize{
#'     \item \code{deriv}: the derivatives of the PDF with respect to the chosen \code{wrt},
#'     \item \code{call}: the function call,
#'     \item \code{err}: the absolute error. Only provided if sv, sw, or st0 is non-zero. If numerical integration is used, the precision cannot always be guaranteed.
#'   }
#' @references
#' Hartmann, R., & Klauer, K. C. (2021). Partial derivatives for the first-passage time distribution in Wiener diffusion models. \emph{Journal of Mathematical Psychology, 103}, 102550. \doi{10.1016/j.jmp.2021.102550}
#' @examples
#' ddWDM(wrt = "a", t = 1.2, response = "upper", a = 1.1, v = 2, w = .6, precision = NULL, K = NULL)
#' @author Raphael Hartmann
#' @useDynLib "WienR", .registration=TRUE
#' @export
ddWDM <- function(wrt,
                  t,
                  response,
                  a,
                  v,
                  w,
                  t0 = 0,
                  sv = 0,
                  sw = 0,
                  st0 = 0,
                  precision = NULL,
                  K = NULL,
                  n.threads = FALSE,
                  n.evals = 6000) {
  
  if(length(wrt) > 1) stop("Argument \"wrt\" must be a single string, not a vector")
  if(!wrt %in% c("t", "a", "v", "w", "t0", "sv", "sw", "st0", "grad")) stop("Argument \"wrt\" not valid.")
  
  if(identical(wrt, "t")) return(dtWienerPDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "a")) return(daWienerPDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "v")) return(dvWienerPDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "w")) return(dwWienerPDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "t0")) return(dt0WienerPDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "sw")) return(dswWienerPDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "sv")) return(dsvWienerPDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "st0")) return(dst0WienerPDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "grad")) return(gradWienerPDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  
}




#' Wrapper function for the partial derivative of the first-passage time cumulative distribution function of the diffusion model
#' 
#' Calculates the partial derivative of the first-passage time cumulative distribution function of the diffusion model with respect to one of a, v, w, t0, sv, sw, or st0, or calculate the gradient.
#' @param wrt partial derivative w.r.t. one of the following:
#'   \itemize{
#'     \item \code{"a"} the upper barrier,
#'     \item \code{"v"} the drift rate,
#'     \item \code{"w"} the relative starting point,
#'     \item \code{"t0"} the non-decision time,
#'     \item \code{"sv"} the inter-trial variability of drift rate,
#'     \item \code{"sw"} the inter-trial variability of relative starting point,
#'     \item \code{"st0"} the inter-trial variability of non-decision time, or
#'     \item \code{"grad"} all the above.
#'   }
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
#' @param K Optional. Number of iterations to calculate the infinite sums. Numeric value (integer). Default is \code{NULL}.
#'   \itemize{
#'     \item \code{precision = NULL} and \code{K = NULL}: Default \code{precision = 1e-12} used to calculate internal K.
#'     \item \code{precision != NULL} and \code{K = NULL}: \code{precision} is used to calculate internal K,
#'     \item \code{precision = NULL} and \code{K != NULL}: \code{K} is used as internal K,
#'     \item \code{precision != NULL} and \code{K != NULL}: if internal K calculated through \code{precision} is smaller than \code{K}, \code{K} is used.
#'   }
#'   We recommend using either default (\code{precision = K = NULL}) or only \code{precision}.
#' @param n.threads Optional numerical or logical value. Number of threads to use. If not provided (or 1 or \code{FALSE}) parallelization is not used. If set to \code{TRUE} then all available threads are used.
#' @param n.evals Optional. Number of maximal function evaluations in the numeric integral if sv, sw, and/or st0 are not zero. Default is \code{6000} and \code{0} implies no limit and the 
#'   numeric integration goes on until the specified \code{precision} is guaranteed.
#' @return A list of the class \code{Diffusion_deriv} containing
#'   \itemize{
#'     \item \code{deriv}: the derivatives of the CDF with respect to the chosen \code{wrt},
#'     \item \code{call}: the function call,
#'     \item \code{err}: the absolute error. Only provided if sv, sw, or st0 is non-zero. If numerical integration is used, the precision cannot always be guaranteed.
#'   }
#' @references
#' Hartmann, R., & Klauer, K. C. (2021). Partial derivatives for the first-passage time distribution in Wiener diffusion models. \emph{Journal of Mathematical Psychology, 103}, 102550. \doi{10.1016/j.jmp.2021.102550}
#' @examples
#' dpWDM(wrt = "a", t = 1.2, response = "upper", a = 1.1, v = 2, w = .6, precision = NULL, K = NULL)
#' @author Raphael Hartmann
#' @useDynLib "WienR", .registration=TRUE
#' @export
dpWDM <- function(wrt,
                  t,
                  response,
                  a,
                  v,
                  w,
                  t0 = 0,
                  sv = 0,
                  sw = 0,
                  st0 = 0,
                  precision = NULL,
                  K = NULL,
                  n.threads = FALSE,
                  n.evals = 6000) {
  
  if(length(wrt) > 1) stop("Argument \"wrt\" must be a single string, not a vector")
  if(!wrt %in% c("a", "v", "w", "t0", "sv", "sw", "st0", "grad")) stop("Argument \"wrt\" not valid.")
  
  if(identical(wrt, "a")) return(daWienerCDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "v")) return(dvWienerCDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "w")) return(dwWienerCDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "t0")) return(dt0WienerCDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "sw")) return(dswWienerCDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "sv")) return(dsvWienerCDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "st0")) return(dst0WienerCDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  if(identical(wrt, "grad")) return(gradWienerCDF(t = t, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, precision = precision, K = K, n.threads = n.threads, n.evals = n.evals))
  
}