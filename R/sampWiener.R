
#' Random sampling from the Wiener diffusion model
#'
#' Draws random samples from the (truncated) first-passage time distribution of the Wiener diffusion model with up to 7 parameters.
#' @param N Number of samples. Numeric value (integer).
#' @param a Upper barrier. Numeric value.
#' @param v Drift rate. Numeric value.
#' @param w Relative starting point. Numeric value.
#' @param t0 Non-decision time. Numeric value.
#' @param sv Inter-trial variability of drift rate. Numeric value. Standard deviation of a normal distribution \code{N(v, sv)}.
#' @param sw Inter-trial variability of relative starting point. Numeric value. Range of uniform distribution \code{U(w-0.5*sw, w+0.5*sw)}.
#' @param st0 Inter-trial variability of non-decision time. Numeric value. Range of uniform distribution \code{U(t0, t0+st0)}.
#' @param response Response boundary. Character string, either \code{"upper"}, \code{"lower"}, or \code{"both"}. Alternatively a numeric value,
#'   \code{2}=upper, \code{1}=lower, or \code{0}=both. Default is "both".
#' @param bound Boundary for the first-passage time. Numeric value. Default is Inf.
#' @param method Sampling method. Either "ars", "rs", "its", or "p-ars". The method "ars" stands for adaptive rejection sampling, "rs" stands for rejection sampling,
#'   "its" stands for inverse transform sampling, and "p-ars" stands for pseudo adaptive rejection sampling. Default is "ars".
#' @param precision Optional numeric value. Precision of the infinite series approximation. Numeric value. Default is \code{NULL}, which takes default value 1e-12.
#' @param n.threads Optional numeric or logic value. Number of threads to use. If not provided (\code{FALSE} or 1) parallelization is not used. If set to 
#'   \code{TRUE} then all available threads are used.
#' @param ars_list Optional list for method "ars". For \code{response} "lower" or "upper" a list with upper hull, lower hull etc. is needed. 
#'   For \code{response} "both" a list with two lists must be provided. The corresponding list is produced when using the "ars" method 
#'   and the argument \code{ARS_STORE = TRUE}. Do not make the list yourself and do not mix the lists for the corresponding boundaries.
#' @param ARS_STORE Optional flag for method "ars". If \code{TRUE} saves upper hull, lower hull and some more values, which are updated at each rejection step,
#'   as a list. The list can then be used with the "ars" method in the argument \code{ars_list} to make the new sampling faster. If the first-passage times
#'   were sampled only from one boundary then the list will contain upper hull, etc. and if they were sampled from both boundaries then the list consists of two
#'   lists, each containing upper hull, etc. for the respective boundary.
#' @details The following \code{methods} can be used:
#'   \itemize{
#'     \item \code{"ars"}: adaptive rejection sampling method. This method builds on Gilks and Wild (1992) as well as Hartmann and Klauer (in press). The former provides
#'       a method for an adaptive rejection sampling method which assumes that the density is log-concave. This method is fastest for cases where \code{sv = 0}. This is
#'       the only method where the integral needs to be calculated. The advantage, though, is that after each rejection the upper and lower hull functions will be adapted
#'       (or updated), which leads to fewer and fewer rejections in the proceeding sampling steps.
#'     \item \code{"rs"}: rejection sampling method by Drugowitsch (2016). This method uses different proposal distributions in different conditions.
#'     \item \code{"its"}: inverse transform (a.k.a. probability integral transform) sampling method. A random sample u is sampled from a uniform
#'       distribution and the corresponding first-passage time, for which CDF(t) = u, is approximated.
#'     \item \code{"p-ars"}: pseudo-adaptive rejection sampling. A variation of "ars". In this method the hull functions will be adapted until the current sample is drawn,
#'       but the information from this adaptation will be discarded for the next sample.
#'   }
#'   Note: The speed of the methods do not depend on \code{t0} or \code{st0}.
#'   
#'   \code{ars_store}, one of the returned list objects if method \code{"ars"} and \code{ARS_STORE = TRUE}, consists of twelve vectors and three scalars:
#'   \itemize{
#'     \item \code{hstore_x}: vector of alpha values -- change of variable \code{alpha = (log(t)-start)/scale}, where t is the first-passage time -- relevant for the upper and
#'       lower hull functions.
#'     \item \code{hstore_h}: vector of log-density of change of variable \code{A = (log(T)-start)/scale} at the alpha points \code{hstore_x}
#'     \item \code{hstore_dh}: vector of partial derivative of log-density of A with respect to alpha.
#'     \item \code{upperstore_z}: vector of alpha values at which the piece-wise linear upper hull transitions from one linear segment to the next.
#'     \item \code{upperstore_slope}: same as \code{hstore_dh}. Gives the slope of the piece-wise linear functions for the upper hull.
#'     \item \code{upperstore_absc}: same as \code{hstore_h}. Gives the evaluation of the function h() at \code{hstore_x}, where the piece-wise linear function touches h().
#'     \item \code{upperstore_center}: same as \code{hstore_x}. Gives the alpha values, where the piece-wise linear function touches h().
#'     \item \code{lowerstore_z}: same as \code{hstore_x} but with an additional leading element (=-Inf) in the vector.
#'     \item \code{lowerstore_slope}: vector of zeros since not needed.
#'     \item \code{lowerstore_absc}: vector of zeros since not needed.
#'     \item \code{lowerstore_center:}: vector of zeros since not needed.
#'     \item \code{startstore}: scalar representing the "start" value for the change of variable \code{A = (log(T)-start)/scale}.
#'     \item \code{scalestore}: scalar representing the "scale" value for the change of variable \code{A = (log(T)-start)/scale}.
#'     \item \code{normstore}: scalar. Gives the value of h() at alpha = 0.
#'     \item \code{sstore}: vector of values at \code{log(s_k(hstore_x))}, where s_k() is the function defined in equation 3 in Gilks and Wild (1992).
#'   }
#' @return A list of the class \code{Diffusion_samp} containing
#'   \itemize{
#'     \item \code{q}: first-passage time sample(s),
#'     \item \code{response}: response(s) "lower" and/or "upper",
#'     \item \code{call}: the function call,
#'     \item \code{ars_store}: if \code{ARS_STORE = TRUE} is used with the method "ars" then either a list with upper hull, etc. is stored (either from the upper
#'       or lower boundary) or a list of two lists with corresponding upper hull, etc. is stored (from both boundaries) and can be used as function argument to
#'       (\code{ars_list}) for further sampling with the same parameters.
#'   }
#' @references
#' Drugowitsch, J. (2016). Fast and accurate Monte Carlo sampling of first-passage times from Wiener diffusion models. \emph{Scientific Reports, 6}(1). \doi{10.1038/srep20490}
#' 
#' Gilks, W. R., & Wild, P. (1992). Adaptive Rejection Sampling for Gibbs Sampling. \emph{Applied Statistics, 41}(2), 337. \doi{10.2307/2347565}
#' 
#' Hartmann, R., & Klauer, K. C. (in press). Partial Derivatives for the First-Passage Time Distribution in Wiener Diffusion Models. \emph{Journal of Mathematical Psychology}.
#' 
#' @examples
#' sample_list1 <- sampWiener(N = 100000, a = 1, v = .3, w = .5)
#' hist(sample_list1$q, 200)
#' 
#' sample_list2 <- sampWiener(N = 100000, a = 1, v = .3, w = .5, ARS_STORE = TRUE)
#' hist(sample_list2$q, 200)
#' sample_list2$ars_store
#' 
#' sample_list3 <- sampWiener(N = 100000, a = 1, v = .3, w = .5, ars_list = sample_list2$ars_store)
#' hist(sample_list3$q, 200)
#' 
#' @author Raphael Hartmann
#' @useDynLib "WienR", .registration=TRUE
#' @export
#' @importFrom stats runif
sampWiener <- function(N,
                       a,
                       v,
                       w, 
                       t0 = 0,
                       sv = 0,
                       sw = 0,
                       st0 = 0,
                       response = "both",
                       bound = Inf,
                       method = "ars",
                       precision = NULL,
                       n.threads = 1,
                       ars_list = NULL,
                       ARS_STORE = FALSE) {




  # ---- VALUE CHECKS ---- #

  # general checks
  lengths <- c(length(response), length(a), length(v), length(w), length(t0), length(sv), length(sw), length(st0), length(bound), length(method), length(n.threads))
  if(any(lengths > 1)) stop("N, a, v, w, t0, sv, sw, st0, response, bound, method, n.threads must be of length one.")

  # t a v w t0 sw sv st0 checks
  if(!is.numeric(a) | !is.numeric(v) | !is.numeric(w) | !is.numeric(t0) | !is.numeric(sv) | !is.numeric(sw) | !is.numeric(st0)) stop("a, v, w, t0, sv, sw, and st0 must be numeric")
  if(a <= 0 | w <= 0) stop("a and w must be strictly positive")
  if(t0 < 0 | sw < 0 | sv < 0 | st0 < 0) stop("t0, sw, sv, and st0 must be positive or zero")
  if(w >= 1) stop("w must be lower than one")
  if(w-0.5*sw <= 0 | w+0.5*sw >= 1) stop("w-0.5*sw must be greater than zero and w+0.5*sw must be lower than one")

  # response checks
  if(!is.character(response) & !is.numeric(response)) stop("response must be a character with the values \"upper\", \"lower\", or \"both\" OR numerics with the values 2=\"upper\", 1=\"lower\" or 0=\"both\"")
  if(length(response)!=1) stop("response must be of length one")
  if(!all(response %in% c("both", "lower", "upper")) & !all(response %in% c(1,2,0)) ) stop("response cannot include values other than \"upper\", \"lower\", and \"both\" OR 2=\"upper\", 1=\"lower\", and 0=\"both\".")
  if(is.character(response)) {
    resp <- which(c("both", "lower", "upper") %in% response) - 1
  } else {
    resp <- response
  } 

  # precision checks
  if(!is.numeric(precision) & !is.null(precision)) stop("precision must either be NULL or some numeric value")
  if(length(precision)!=1 & !is.null(precision)) stop("precision must be of length one")

  PRECISION_FLAG <- TRUE
  if(is.null(precision)) PRECISION_FLAG <- FALSE

  if(is.null(precision)) precision <- 1e-12

  # thread checks
  if(!is.numeric(n.threads) & !is.logical(n.threads)) stop("n.threads must either be numerical or logical")
  if(is.numeric(n.threads)) if(n.threads %% 1 != 0) stop("n.threads must be an integer") else n.threads <- as.integer(n.threads)
  if(is.logical(n.threads)) n.threads <- ifelse(n.threads == TRUE, 99999, 0)
  if(n.threads < 2) n.threads <- 0
  
  # ars_list check
  if(!is.null(ars_list)) {
    ars_vector <- NULL
    ars_nrs <- NULL
    if(!is.list(ars_list)) stop("ars_list must be a list")
    list_names <- c("hstore_x", "hstore_h", "hstore_dh", "lowerstore_z", "lowerstore_slope", "lowerstore_absc", "lowerstore_center", "upperstore_z", 
                    "upperstore_slope", "upperstore_absc", "upperstore_center", "startstore", "scalestore", "normstore", "sstore")
    if(length(ars_list) != 2) {
      
      if(!all(names(ars_list) %in% list_names)) {
        warning("The provided ars_list seems to be incorrect. ars_list will not be used.")
        ars_vector <- 0
      } else {
        ars_vector <- as.numeric(unlist(ars_list))
        ars_nrs <- c(1, length(ars_list$hstore_x), length(ars_list$lowerstore_z), length(ars_list$upperstore_z), length(ars_list$sstore))
      }

    } else {
      
      if(!all(names(ars_list[[1]]) %in% list_names) | !all(names(ars_list[[2]]) %in% list_names)) {
        warning("The provided ars_list seems to be incorrect. ars_list will not be used.")
        ars_vector <- 0
      } else {
        ars_vector <- as.numeric(c(unlist(ars_list[[1]]), unlist(ars_list[[2]])))
        ars_nrs <- c(2, 
                     length(ars_list[[1]]$hstore_x), length(ars_list[[1]]$lowerstore_z), length(ars_list[[1]]$upperstore_z), length(ars_list[[1]]$sstore),
                     length(ars_list[[2]]$hstore_x), length(ars_list[[2]]$lowerstore_z), length(ars_list[[2]]$upperstore_z), length(ars_list[[2]]$sstore))
      }
      
    }
    
  } else {ars_vector = 0; ars_nrs <- 0}
  
  # methods check
  method_names <- c("ars", "rs", "its", "p-ars")
  if(! method %in% method_names) stop("method not valid.")
  choice <- which(method_names == method)
  


  # --- C++ FUNCTION CALL ---- #

  samp_list <- .Call("randWiener", 
                     a, 
                     v, 
                     w, 
                     sv, 
                     sw, 
                     precision, 
                     bound, 
                     ars_vector, 
                     as.integer(resp),
                     as.integer(0), 
                     as.integer(N), 
                     as.integer(n.threads), 
                     as.integer(choice), 
                     as.integer(ars_nrs), 
                     as.integer(1), 
                     as.integer(ARS_STORE)
  )
  
  
  # add t0
  if(t0 != 0) {
    temp <- t0 + ifelse(st0 != 0, st0 * runif(N), 0)
    samp_list$q <- samp_list$q + temp
  }
  
  
  # responses
  responses <- ifelse(samp_list$resp == 1, "lower", "upper")


  # output
  output <- list(q = samp_list$q, response = responses, call = match.call())
  if(ARS_STORE == TRUE & method == "ars") {
    if(response=="both") {
      ars_list_save <- list(ars.store.upp = samp_list$ars.store.upp, ars.store.low = samp_list$ars.store.low)
    } else {
      ars_list_save <- samp_list$ars.store
    }
    output$ars_store <- ars_list_save
  }
  class(output) <- "Diffusion_samp"
  return(output)

}
