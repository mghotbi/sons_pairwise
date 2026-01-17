#' Compute SONS (Shared OTUs and Similarity) pairwise metrics for two communities
#'
#' This function implements the OTU-based similarity metrics described in:
#' Schloss PD & Handelsman J (2006) "Introducing SONS..."
#' Specifically:
#' - Classic incidence-based Jaccard and Sørenson: Eq. (3)–(4)
#' - Estimated fraction of sequences belonging to shared OTUs: U_est, V_est: Eq. (5)–(6)
#' - Abundance-based Jaccard and Sørenson derived from U_est and V_est: Eq. (7)–(8)
#' - Yue & Clayton structure similarity (theta): Eq. (9)
#'
#' Important notes about inputs:
#' - x and y must be aligned vectors of OTU counts (same OTUs in same order).
#' - Counts must be nonnegative; integer counts are expected but numeric is allowed.
#' - OTUs absent from a community have count 0.
#'
#' Handling f2 = 0:
#' Eq. (5)–(6) include ratios of the form f1/(2*f2). When f2 = 0 (no shared doubletons),
#' the original paper does not specify a correction. This implementation provides:
#' - f2_correction = "none": sets the correction term to 0 (conservative).
#' - f2_correction = "add1": continuity correction using 2*(f2+1) in the denominator.
#'
#' Bootstrapping:
#' - The original paper suggests bootstrap CIs for Eq. (5)–(8).
#' - Here we implement a straightforward multinomial bootstrap of sequence allocation
#'   within each community (conditioned on observed relative abundances).
#'
#' @param x Numeric vector of OTU counts for community A.
#' @param y Numeric vector of OTU counts for community B.
#' @param f2_correction Character. How to handle cases where f_.2 or f_2. equals 0 in Eq. (5)–(6).
#'   One of c("none","add1"). "none" sets the correction term to 0; "add1" uses a continuity correction
#'   by replacing (2*f2) with (2*(f2+1)) when f2 == 0.
#' @param bootstrap Logical; if TRUE, compute bootstrap intervals for selected metrics.
#' @param B Integer; number of bootstrap replicates (used if bootstrap = TRUE).
#' @param conf Numeric in (0,1); confidence level for bootstrap intervals (e.g., 0.95).
#' @param seed Optional integer seed for reproducible bootstrapping.
#'
#' @return A list with:
#' - inputs summary (n_total, m_total, etc.)
#' - incidence-based metrics (S1_obs, S2_obs, S12_obs, Jclas, Lclas)
#' - abundance-based overlap estimators (U_est, V_est, Jabund, Labund)
#' - structure similarity (theta_yc) and distance (1 - theta_yc)
#' - optional bootstrap intervals
#'
#' @examples
#' x <- c(10, 5, 1, 0, 2)
#' y <- c( 0, 3, 1, 7, 1)
#' sons_pairwise(x, y, f2_correction = "none")
sons_pairwise <- function(x,
                          y,
                          f2_correction = c("none", "add1"),
                          bootstrap = FALSE,
                          B = 999,
                          conf = 0.95,
                          seed = NULL) {
  f2_correction <- match.arg(f2_correction)
  
  # ----------------------------
  # 0) Defensive checks
  # ----------------------------
  if (!is.numeric(x) || !is.numeric(y)) stop("x and y must be numeric vectors.")
  if (length(x) != length(y)) stop("x and y must have the same length (aligned OTUs).")
  if (any(is.na(x)) || any(is.na(y))) stop("x and y must not contain NA.")
  if (any(x < 0) || any(y < 0)) stop("Counts must be nonnegative.")
  if (conf <= 0 || conf >= 1) stop("conf must be in (0, 1).")
  if (!is.logical(bootstrap) || length(bootstrap) != 1) stop("bootstrap must be TRUE/FALSE.")
  if (!is.numeric(B) || length(B) != 1 || B < 1) stop("B must be a positive integer.")
  
  # Convert to plain numeric (handles integer, etc.)
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  # Total sequences in each community (ntotal, mtotal in paper)
  ntotal <- sum(x)
  mtotal <- sum(y)
  
  # Presence/absence vectors
  pres_x <- x > 0
  pres_y <- y > 0
  
  # Shared OTUs are those with Xi>0 AND Yi>0
  shared <- pres_x & pres_y
  
  # ----------------------------
  # 1) Incidence-based richness counts and similarities: Eq. (3)–(4)
  # ----------------------------
  S1_obs  <- sum(pres_x)   # S1: number of OTUs observed in A
  S2_obs  <- sum(pres_y)   # S2: number of OTUs observed in B
  S12_obs <- sum(shared)   # S12: number of OTUs shared (observed)
  
  # Classic Jaccard: Jclas = S12 / (S1 + S2 - S12)
  denom_j <- (S1_obs + S2_obs - S12_obs)
  Jclas <- if (denom_j == 0) NA_real_ else S12_obs / denom_j
  
  # Classic Sørenson: Lclas = 2*S12 / (S1 + S2)
  denom_l <- (S1_obs + S2_obs)
  Lclas <- if (denom_l == 0) NA_real_ else (2 * S12_obs) / denom_l
  
  # ----------------------------
  # 2) Abundance-based overlap estimators U_est and V_est: Eq. (5)–(6)
  # ----------------------------
  x_shared <- x[shared]
  y_shared <- y[shared]
  
  if (ntotal == 0 || mtotal == 0) {
    U_est <- NA_real_
    V_est <- NA_real_
  } else if (length(x_shared) == 0) {
    U_est <- 0
    V_est <- 0
  } else {
    base_U <- sum(x_shared / ntotal)
    base_V <- sum(y_shared / mtotal)
    
    f_dot1 <- sum(y_shared == 1)  # f_.1
    f_dot2 <- sum(y_shared == 2)  # f_.2
    f_1dot <- sum(x_shared == 1)  # f_1.
    f_2dot <- sum(x_shared == 2)  # f_2.
    
    sumU_ind <- sum((x_shared / ntotal) * (y_shared == 1))
    sumV_ind <- sum((y_shared / mtotal) * (x_shared == 1))
    
    corr_U <- 0
    if (mtotal > 0) {
      if (f_dot2 > 0) {
        mult_U <- f_dot1 / (2 * f_dot2)
      } else {
        mult_U <- if (f2_correction == "none") 0 else f_dot1 / (2 * (f_dot2 + 1))
      }
      corr_U <- ((mtotal - 1) / mtotal) * mult_U * sumU_ind
    }
    
    corr_V <- 0
    if (ntotal > 0) {
      if (f_2dot > 0) {
        mult_V <- f_1dot / (2 * f_2dot)
      } else {
        mult_V <- if (f2_correction == "none") 0 else f_1dot / (2 * (f_2dot + 1))
      }
      corr_V <- ((ntotal - 1) / ntotal) * mult_V * sumV_ind
    }
    
    U_est <- base_U + corr_U
    V_est <- base_V + corr_V
    
    U_est <- min(max(U_est, 0), 1)
    V_est <- min(max(V_est, 0), 1)
  }
  
  # Abundance-based Jaccard and Sørenson from Eq. (7)–(8)
  if (is.na(U_est) || is.na(V_est)) {
    Jabund <- NA_real_
    Labund <- NA_real_
  } else {
    denom_jab <- (U_est + V_est - U_est * V_est)
    Jabund <- if (denom_jab == 0) NA_real_ else (U_est * V_est) / denom_jab
    
    denom_lab <- (U_est + V_est)
    Labund <- if (denom_lab == 0) NA_real_ else (2 * U_est * V_est) / denom_lab
  }
  
  # ----------------------------
  # 3) Yue & Clayton community structure similarity: Eq. (9)
  # ----------------------------
  if (ntotal == 0 || mtotal == 0) {
    theta_yc <- NA_real_
  } else {
    pA <- x / ntotal
    pB <- y / mtotal
    
    cross <- sum(pA * pB)
    denom_theta <- (sum(pA^2) + sum(pB^2) - cross)
    
    theta_yc <- if (denom_theta == 0) NA_real_ else cross / denom_theta
    theta_yc <- min(max(theta_yc, 0), 1)
  }
  
  dist_theta_yc <- if (is.na(theta_yc)) NA_real_ else (1 - theta_yc)
  
  # ----------------------------
  # 4) Optional bootstrap confidence intervals (basic, reproducible)
  # ----------------------------
  boot <- NULL
  if (bootstrap) {
    if (!is.null(seed)) set.seed(seed)
    
    if (ntotal == 0 || mtotal == 0) {
      boot <- list(note = "Bootstrap skipped because ntotal==0 or mtotal==0.")
    } else {
      alpha <- (1 - conf) / 2
      qs <- c(alpha, 1 - alpha)
      
      compute_core <- function(xb, yb) {
        tmp <- sons_pairwise(xb, yb, f2_correction = f2_correction, bootstrap = FALSE)
        c(
          U_est = tmp$U_est,
          V_est = tmp$V_est,
          Jabund = tmp$Jabund,
          Labund = tmp$Labund,
          theta_yc = tmp$theta_yc
        )
      }
      
      pA <- x / ntotal
      pB <- y / mtotal
      
      mat <- replicate(B, {
        xb <- as.numeric(rmultinom(1, size = ntotal, prob = pA))
        yb <- as.numeric(rmultinom(1, size = mtotal, prob = pB))
        compute_core(xb, yb)
      })
      
      boot_est <- t(mat)
      
      boot <- list(
        B = B,
        conf = conf,
        f2_correction = f2_correction,
        intervals = apply(boot_est, 2, quantile, probs = qs, na.rm = TRUE),
        boot_estimates = boot_est
      )
    }
  }
  
  # ----------------------------
  # 5) Return results
  # ----------------------------
  out <- list(
    n_total = ntotal,
    m_total = mtotal,
    S1_obs = S1_obs,
    S2_obs = S2_obs,
    S12_obs = S12_obs,
    Jclas = Jclas,
    Lclas = Lclas,
    U_est = U_est,
    V_est = V_est,
    Jabund = Jabund,
    Labund = Labund,
    theta_yc = theta_yc,
    dist_theta_yc = dist_theta_yc,
    f2_correction = f2_correction
  )
  
  if (bootstrap) out$bootstrap <- boot
  
  out
}


#' Convenience wrapper: compute SONS metrics from an OTU table for two named samples
#'
#' @param otu_mat Matrix/data.frame with OTUs in rows and samples in columns (counts).
#' @param sample_a Column name for community A.
#' @param sample_b Column name for community B.
#' @param ... Passed to sons_pairwise (e.g., f2_correction="add1", bootstrap=TRUE, B=999, seed=1).
#'
#' @return Output of sons_pairwise for the two selected samples.
sons_from_otu_table <- function(otu_mat, sample_a, sample_b, ...) {
  if (is.data.frame(otu_mat)) otu_mat <- as.matrix(otu_mat)
  if (!is.matrix(otu_mat)) stop("otu_mat must be a matrix or data.frame.")
  if (is.null(colnames(otu_mat))) stop("otu_mat must have column names (sample IDs).")
  if (!(sample_a %in% colnames(otu_mat))) stop("sample_a not found in otu_mat colnames.")
  if (!(sample_b %in% colnames(otu_mat))) stop("sample_b not found in otu_mat colnames.")
  
  x <- otu_mat[, sample_a]
  y <- otu_mat[, sample_b]
  
  sons_pairwise(x, y, ...)
}
