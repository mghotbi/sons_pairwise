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
#' Bootstrapping:
#' - The original paper suggests bootstrap CIs for Eq. (5)–(8). :contentReference[oaicite:1]{index=1}
#' - Here we implement a straightforward multinomial bootstrap of sequence allocation
#'   within each community (conditioned on observed relative abundances).
#'
#' @param x Numeric vector of OTU counts for community A.
#' @param y Numeric vector of OTU counts for community B.
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
#' sons_pairwise(x, y)
sons_pairwise <- function(x,
                          y,
                          bootstrap = FALSE,
                          B = 999,
                          conf = 0.95,
                          seed = NULL) {
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
  
  # Shared OTUs are those with Xi>0 AND Yi>0 (belongs to D12 in Eq. 5–6) :contentReference[oaicite:2]{index=2}
  shared <- pres_x & pres_y
  
  # ----------------------------
  # 1) Incidence-based richness counts and similarities: Eq. (3)–(4) :contentReference[oaicite:3]{index=3}
  # ----------------------------
  S1_obs  <- sum(pres_x)        # S1: number of OTUs observed in A
  S2_obs  <- sum(pres_y)        # S2: number of OTUs observed in B
  S12_obs <- sum(shared)        # S12: number of OTUs shared (observed)
  
  # Classic Jaccard: Jclas = S12 / (S1 + S2 - S12)
  denom_j <- (S1_obs + S2_obs - S12_obs)
  Jclas <- if (denom_j == 0) NA_real_ else S12_obs / denom_j
  
  # Classic Sørenson: Lclas = 2*S12 / (S1 + S2)
  denom_l <- (S1_obs + S2_obs)
  Lclas <- if (denom_l == 0) NA_real_ else (2 * S12_obs) / denom_l
  
  # ----------------------------
  # 2) Abundance-based overlap estimators U_est and V_est: Eq. (5)–(6) :contentReference[oaicite:4]{index=4}
  # ----------------------------
  # Notation in paper:
  # - Xi, Yi are abundances of i-th shared OTU in A and B
  # - Summations are over i = 1..D12 (D12 = number of shared OTUs)
  # - I( Yi == 1 ) indicator used in correction term (singletons in the other community)
  #
  # We implement exactly:
  # U_est = sum_{shared} Xi/ntotal + ((mtotal-1)/mtotal) * (f_.1/(2*f_.2)) * sum_{shared} Xi/ntotal * I(Yi==1)
  # V_est = sum_{shared} Yi/mtotal + ((ntotal-1)/ntotal) * (f_1./(2*f_2.)) * sum_{shared} Yi/mtotal * I(Xi==1)
  #
  # where:
  # f_.1 = number of shared OTUs with Yi==1 (singletons in B among shared OTUs)
  # f_.2 = number of shared OTUs with Yi==2 (doubletons in B among shared OTUs)
  # f_1. = number of shared OTUs with Xi==1 (singletons in A among shared OTUs)
  # f_2. = number of shared OTUs with Xi==2 (doubletons in A among shared OTUs)
  
  x_shared <- x[shared]
  y_shared <- y[shared]
  
  # If there are no shared OTUs, the sums are 0 and U_est = V_est = 0
  # (interpretable: no sequences belong to shared OTUs).
  # If ntotal or mtotal are 0, metrics are undefined (no sequences sampled).
  if (ntotal == 0 || mtotal == 0) {
    U_est <- NA_real_
    V_est <- NA_real_
  } else if (length(x_shared) == 0) {
    U_est <- 0
    V_est <- 0
  } else {
    # Base terms: sum Xi/ntotal and sum Yi/mtotal over shared OTUs
    base_U <- sum(x_shared / ntotal)
    base_V <- sum(y_shared / mtotal)
    
    # Singleton/doubleton counts among *shared* OTUs
    f_dot1 <- sum(y_shared == 1)  # f_.1
    f_dot2 <- sum(y_shared == 2)  # f_.2
    f_1dot <- sum(x_shared == 1)  # f_1.
    f_2dot <- sum(x_shared == 2)  # f_2.
    
    # Indicator-weighted sums in the correction terms
    sumU_ind <- sum((x_shared / ntotal) * (y_shared == 1))
    sumV_ind <- sum((y_shared / mtotal) * (x_shared == 1))
    
    # Correction multipliers:
    # ((mtotal - 1)/mtotal) * (f_.1/(2*f_.2))  for U_est
    # ((ntotal - 1)/ntotal) * (f_1./(2*f_2.))  for V_est
    #
    # Edge case: f_.2 == 0 or f_2. == 0 → ratio would blow up.
    # The paper does not spell out a special-case formula here, but the standard
    # practice in Chao-type estimators is: if there are no doubletons, the
    # adjustment cannot be estimated reliably.
    #
    # Conservative choice: set the correction to 0 when f2 == 0.
    # (This avoids infinite values and yields a lower-bound style estimate.)
    corr_U <- 0
    if (f_dot2 > 0 && mtotal > 0) {
      corr_U <- ((mtotal - 1) / mtotal) * (f_dot1 / (2 * f_dot2)) * sumU_ind
    }
    
    corr_V <- 0
    if (f_2dot > 0 && ntotal > 0) {
      corr_V <- ((ntotal - 1) / ntotal) * (f_1dot / (2 * f_2dot)) * sumV_ind
    }
    
    U_est <- base_U + corr_U
    V_est <- base_V + corr_V
    
    # Numerical safety: keep within [0, 1] (should be anyway, but rounding/correction can overshoot)
    U_est <- min(max(U_est, 0), 1)
    V_est <- min(max(V_est, 0), 1)
  }
  
  # Abundance-based Jaccard and Sørenson from Eq. (7)–(8) :contentReference[oaicite:5]{index=5}
  if (is.na(U_est) || is.na(V_est)) {
    Jabund <- NA_real_
    Labund <- NA_real_
  } else {
    # Eq. (7): Jabund = (U*V) / (U + V - U*V)
    denom_jab <- (U_est + V_est - U_est * V_est)
    Jabund <- if (denom_jab == 0) NA_real_ else (U_est * V_est) / denom_jab
    
    # Eq. (8): Labund = (2*U*V) / (U + V)
    denom_lab <- (U_est + V_est)
    Labund <- if (denom_lab == 0) NA_real_ else (2 * U_est * V_est) / denom_lab
  }
  
  # ----------------------------
  # 3) Yue & Clayton community structure similarity: Eq. (9) :contentReference[oaicite:6]{index=6}
  # ----------------------------
  # The paper presents theta as:
  # theta = [ sum_{i in shared} (Xi/ntotal)*(Yi/mtotal) ] /
  #         [ sum_{i in A} (Xi/ntotal)^2 + sum_{i in B} (Yi/mtotal)^2 - sum_{i in shared} (Xi/ntotal)*(Yi/mtotal) ]
  #
  # We compute using full vectors (including zeros) because:
  # - sum over A includes OTUs not in B (their Yi = 0)
  # - sum over B includes OTUs not in A (their Xi = 0)
  # - the cross-term only contributes for shared OTUs anyway
  #
  # That makes the computation clear and robust.
  
  if (ntotal == 0 || mtotal == 0) {
    theta_yc <- NA_real_
  } else {
    pA <- x / ntotal
    pB <- y / mtotal
    
    cross <- sum(pA * pB)   # equivalent to sum over shared (since pA*pB = 0 when not shared)
    denom_theta <- (sum(pA^2) + sum(pB^2) - cross)
    
    theta_yc <- if (denom_theta == 0) NA_real_ else cross / denom_theta
    
    # Keep within [0,1] for numerical stability
    theta_yc <- min(max(theta_yc, 0), 1)
  }
  
  # Commonly used distance is (1 - theta) (paper shows distance = 1 - theta in dendrograms) :contentReference[oaicite:7]{index=7}
  dist_theta_yc <- if (is.na(theta_yc)) NA_real_ else (1 - theta_yc)
  
  # ----------------------------
  # 4) Optional bootstrap confidence intervals (basic, reproducible)
  # ----------------------------
  # The paper mentions bootstrapping for Eq. (5)–(8). :contentReference[oaicite:8]{index=8}
  # We implement a pragmatic multinomial bootstrap:
  # - For each community, resample ntotal sequences among OTUs according to pA,
  #   and mtotal sequences according to pB.
  # - Recompute U_est, V_est, Jabund, Labund (and optionally theta_yc).
  #
  # Caveat: This bootstrap conditions on observed composition; it is not the only
  # possible bootstrap scheme, but is a standard and reproducible option.
  
  boot <- NULL
  if (bootstrap) {
    if (!is.null(seed)) set.seed(seed)
    
    if (ntotal == 0 || mtotal == 0) {
      boot <- list(note = "Bootstrap skipped because ntotal==0 or mtotal==0.")
    } else {
      alpha <- (1 - conf) / 2
      qs <- c(alpha, 1 - alpha)
      
      # Helper to compute only the estimators we care about for each bootstrap replicate
      compute_core <- function(xb, yb) {
        # re-use the same function but with bootstrap disabled
        tmp <- sons_pairwise(xb, yb, bootstrap = FALSE)
        c(U_est = tmp$U_est, V_est = tmp$V_est, Jabund = tmp$Jabund, Labund = tmp$Labund, theta_yc = tmp$theta_yc)
      }
      
      pA <- x / ntotal
      pB <- y / mtotal
      
      # If a community has all mass in a single OTU, multinomial is still fine.
      mat <- replicate(B, {
        xb <- as.numeric(rmultinom(1, size = ntotal, prob = pA))
        yb <- as.numeric(rmultinom(1, size = mtotal, prob = pB))
        compute_core(xb, yb)
      })
      
      # replicate returns a matrix with metrics as rows
      boot_est <- t(mat)
      
      boot <- list(
        B = B,
        conf = conf,
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
    dist_theta_yc = dist_theta_yc
  )
  
  if (bootstrap) out$bootstrap <- boot
  
  out
}


#' Convenience wrapper: compute SONS metrics from an OTU table for two named samples
#'
#' @param otu_mat Matrix/data.frame with OTUs in rows and samples in columns (counts).
#' @param sample_a Column name for community A.
#' @param sample_b Column name for community B.
#' @param ... Passed to sons_pairwise (e.g., bootstrap=TRUE, B=999, seed=1).
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
