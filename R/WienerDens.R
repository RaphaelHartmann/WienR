
#' First-passage time probability density function of the diffusion model
#'
#' Calculate the first-passage time probability density function of the diffusion model.
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
#' @param precision Optional numeric value. Precision of the PDF. Numeric value. Default is \code{NULL}, which takes default value 1e-12.
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
#' @return A list of the class \code{Diffusion_pdf} containing
#'   \itemize{
#'     \item \code{pdf}: the PDF,
#'     \item \code{logpdf}: the log-transformed PDF,
#'     \item \code{call}: the function call,
#'     \item \code{err}: the absolute error. Only provided if sv, sw, or st0 is non-zero. If numerical integration is used, the precision cannot always be guaranteed.
#'   }
#' @references
#' Blurton, S. P., Kesselmeier, M., & Gondan, M. (2017). The first-passage time distribution for the diffusion model with variable drift. \emph{Journal of Mathematical Psychology, 76}, 7–12. \doi{10.1016/j.jmp.2016.11.003}
#' 
#' Foster, K., & Singmann, H. (2021). Another Approximation of the First-Passage Time Densities for the Ratcliff Diffusion Decision Model. \emph{arXiv preprint arXiv:2104.01902}.
#' 
#' Gondan, M., Blurton, S. P., & Kesselmeier, M. (2014). Even faster and even more accurate first-passage time densities and distributions for the Wiener diffusion model. \emph{Journal of Mathematical Psychology, 60}, 20–22. \doi{10.1016/j.jmp.2014.05.002}
#'
#' Hartmann, R., & Klauer, K. C. (2021). Partial derivatives for the first-passage time distribution in Wiener diffusion models. \emph{Journal of Mathematical Psychology, 103}, 102550. \doi{10.1016/j.jmp.2021.102550}
#' 
#' Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate calculations for first-passage times in Wiener diffusion models. \emph{Journal of Mathematical Psychology, 53}(4), 222–230. \doi{10.1016/j.jmp.2009.02.003}
#'
#' Wabersich, D., & Vandekerckhove, J. (2014). The RWiener Package: an R Package Providing Distribution Functions for the Wiener Diffusion Model. \emph{The R Journal, 6(1)}, 49. \doi{10.32614/rj-2014-005}
#' @examples
#' WienerPDF(t = 1.2, response = "upper", a = 1.1, v = 2, w = .6, precision = NULL, K = NULL)
#' @author Raphael Hartmann
#' @useDynLib "WienR", .registration=TRUE
#' @export
WienerPDF <- function(t,
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

  # num. integral evaluation checks
  if(any(sw!=0) | any(st0!=0)) {
    if(!is.numeric(n.evals)) stop("n.evals must numeric")
    if(n.evals %% 1 != 0 | n.evals < 0) stop("n.evals must be an integer and larger or equal to 0")
  }


  # --- C++ FUNCTION CALL ---- #

  indW <- which(sw==0 & st0==0)
  if(length(indW)==0) indD <- 1:max_len else indD <- (1:max_len)[-indW]

  out <- list(pdf = rep(NA, max_len), logpdf = rep(NA, max_len), err = rep(precision, max_len))

  if (length(indW) > 0) {
    tt <- t[indW]-t0[indW]
    temp <- .Call("dWiener",
                  as.numeric(ifelse(tt<0, 0, tt)),
                  as.numeric(a[indW]),
                  as.numeric(v[indW]),
                  as.numeric(w[indW]),
                  as.numeric(sv[indW]),
                  as.numeric(precision),
                  as.integer(resps[indW]),
                  as.integer(K),
                  as.integer(length(indW)),
                  as.integer(n.threads),
                  as.logical(PRECISION_FLAG)
    )
    out$pdf[indW] <- temp$pdf
    out$logpdf[indW] <- temp$logpdf
  }
  if (length(indD) > 0){
    temp <- .Call("dDiffusion7",
                  as.numeric(t[indD]),
                  as.numeric(a[indD]),
                  as.numeric(v[indD]),
                  as.numeric(t0[indD]),
                  as.numeric(w[indD]),
                  as.numeric(sw[indD]),
                  as.numeric(sv[indD]),
                  as.numeric(st0[indD]),
                  as.numeric(precision),
                  as.integer(resps[indD]),
                  as.integer(K),
                  as.integer(length(indD)),
                  as.integer(n.threads),
                  as.integer(0),
                  as.integer(n.evals),
                  as.logical(PRECISION_FLAG)
    )
    out$pdf[indD] <- temp$pdf
    out$logpdf[indD] <- temp$logpdf
    out$err[indD] <- temp$err
  }


  #print(out)

  outcome <- list(value = out$pdf, logvalue = out$logpdf, call = match.call())
  if (length(indD) > 0) outcome$err = out$err
  
  # output
  class(outcome) <- "Diffusion_pdf"
  return(outcome)

}


#' @rdname WienerPDF
#' @examples
#' dWDM(t = 1.2, response = "upper", a = 1.1, v = 2, w = .6, precision = NULL, K = NULL)
#' @export
dWDM <- WienerPDF