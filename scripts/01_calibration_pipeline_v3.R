# =============================================================================
# Low-Cost NIR soil analysis pipeline for TC, TN, pH (H2O), and Available P
# Script 1: Calibration Data Preparation and Modelling
#
# Project: Developing a feasible and low-cost approach through machine learning,
#          near-infrared sensor, and open-access dataset for environmentally
#          sustainable analysis of soils in agriculture
#
# Instruments: ASD FieldSpec (350-2500 nm) and NIRVascan (900-1700 nm)
# Targets:     Total Carbon (%), Total Nitrogen (%), Phosphorus (mg/100g), pH
# Models:      PLSR and Cubist
# Indices:     1=DS, 2=CF global, 3=CF band-wise (LB), 4=PDS,
#              5=CF+DS, 6=RBR (novel - Regularised Band-wise Regression)
#
# Output:      Saved models (.rds), validation plots, metrics tables
# =============================================================================


# =============================================================================
# SECTION 0 - CONFIGURATION (edit this block only)
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
cat("Working directory set to:", getwd(), "\n\n")

cfg <- list(
  dir_root    = getwd(),
  wl_min      = 950,
  wl_max      = 1650,
  wl_step     = 5,
  asd_splice  = c(1000, 1800),
  asd_interp  = 10,
  sg_m        = 0,
  sg_p        = 2,
  sg_w        = 11,
  lb_lambda   = 1400,
  cal_frac    = 0.80,
  clhs_iter   = 500,
  random_seed = 42,
  pls_ncomp   = 15,
  cub_folds   = 10,
  cub_repeats = 10,
  n_cores     = NULL,
  # RBR (Index 6) parameters
  rbr_lambda  = NULL    # NULL = auto-select via LOO-CV per band; or set a fixed value e.g. 0.01
)

cfg$dir_data    <- file.path(cfg$dir_root, "data")
cfg$dir_results <- file.path(cfg$dir_root, "results")
cfg$dir_figs    <- file.path(cfg$dir_root, "results", "figs")
cfg$dir_models  <- file.path(cfg$dir_root, "models")
cfg$dir_tables  <- file.path(cfg$dir_root, "results", "tables")

for (d in c(cfg$dir_data, cfg$dir_results, cfg$dir_figs,
            cfg$dir_models, cfg$dir_tables)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

#load(file.path(cfg$dir_root, "scripts", "calibration_pipeline.RData"))

# =============================================================================
# SECTION 1 - PACKAGES
# =============================================================================

pkgs <- c("dplyr", "tidyr", "readr", "stringr", "tibble",
          "prospectr", "pls", "Cubist", "caret",
          "clhs", "e1071", "cowplot", "gtable",
          "ggplot2", "patchwork",
          "parallel", "doParallel", "purrr")

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

n_cores <- if (is.null(cfg$n_cores)) max(1, parallel::detectCores() - 2) else cfg$n_cores
message(sprintf("Using %d parallel cores.", n_cores))


# =============================================================================
# SECTION 2 - HELPER FUNCTIONS
# =============================================================================

# 2.1  Goodness-of-fit ----------------------------------------------------------
goof <- function(observed, predicted, type = "spec") {
  fit  <- lm(observed ~ predicted)
  MEC  <- summary(fit)$adj.r.squared
  RMSE <- sqrt(mean((observed - predicted)^2))
  bias <- mean(predicted) - mean(observed)
  RPD  <- sd(observed) / RMSE
  IQ   <- as.numeric(quantile(observed, 0.75) - quantile(observed, 0.25))
  RPIQ <- IQ / RMSE
  mx <- mean(observed); my <- mean(predicted)
  s2x <- var(observed); s2y <- var(predicted)
  sxy <- mean((observed - mx) * (predicted - my))
  CCC <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  data.frame(MEC  = round(MEC,  3), CCC  = round(CCC,  3),
             RMSE = round(RMSE, 3), bias = round(bias, 3),
             RPD  = round(RPD,  3), RPIQ = round(RPIQ, 3))
}

# 2.2  Publication theme --------------------------------------------------------
theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
      panel.border      = element_rect(colour = "grey40"),
      axis.title        = element_text(size = base_size, colour = "black"),
      axis.text         = element_text(size = base_size - 1, colour = "black"),
      plot.title        = element_text(size = base_size, face = "bold",
                                       hjust = 0, margin = margin(b = 4)),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key.size   = unit(0.4, "cm"),
      legend.text       = element_text(size = base_size - 1),
      strip.background  = element_rect(fill = "grey95"),
      strip.text        = element_text(size = base_size - 1, face = "bold")
    )
}

# 2.3  Obs vs pred plot ---------------------------------------------------------
obs_pred_plot <- function(pred, obs, title = "", unit = "%", axis_lim = NULL) {
  metrics  <- goof(observed = obs, predicted = pred)
  if (is.null(axis_lim)) {
    rng      <- range(c(obs, pred), na.rm = TRUE)
    axis_lim <- c(floor(rng[1] * 0.95), ceiling(rng[2] * 1.05))
  }
  ax <- axis_lim[1] + 0.05 * diff(axis_lim)
  ay <- axis_lim[2];  dy <- 0.07 * diff(axis_lim)
  ann <- data.frame(
    x = ax, y = ay - (0:4) * dy,
    label = c(sprintf("RMSE = %.3f", metrics$RMSE),
              sprintf("MEC  = %.3f", metrics$MEC),
              sprintf("CCC  = %.3f", metrics$CCC),
              sprintf("bias = %.3f", metrics$bias),
              sprintf("RPIQ = %.3f", metrics$RPIQ))
  )
  ggplot(data.frame(pred = pred, obs = obs), aes(x = pred, y = obs)) +
    geom_point(shape = 16, size = 1.8, colour = "grey30", alpha = 0.75) +
    geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "#D62728", linewidth = 0.7, linetype = "dashed") +
    geom_text(data = ann, aes(x = x, y = y, label = label),
              hjust = 0, vjust = 1, size = 2.8, family = "mono") +
    scale_x_continuous(limits = axis_lim, expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = axis_lim, expand = expansion(mult = 0.01)) +
    labs(title = title,
         x = bquote("Predicted (" * .(unit) * ")"),
         y = bquote("Observed ("  * .(unit) * ")")) +
    theme_pub()
}

# 2.4  Spectral pre-processing --------------------------------------------------
preprocess_asd <- function(df, wl, splice = c(1000, 1800), interp = 10,
                            sg_m = 0, sg_p = 2, sg_w = 11) {
  mat <- as.matrix(df[, -1]);  rownames(mat) <- df[[1]]
  mat <- apply(mat, 2, as.numeric);  rownames(mat) <- df[[1]]
  mat_spl <- prospectr::spliceCorrection(mat, wl, splice = splice,
                                          interpol.bands = interp)
  mat_sg  <- prospectr::savitzkyGolay(mat_spl, m = sg_m, p = sg_p, w = sg_w)
  list(orig = data.frame(id = df[[1]], as.data.frame(mat_spl)),
       sg   = data.frame(id = df[[1]], as.data.frame(mat_sg)))
}

preprocess_nirvascan <- function(df, sg_m = 0, sg_p = 2, sg_w = 11) {
  mat <- apply(as.matrix(df[, -1]), 2, as.numeric)
  rownames(mat) <- df[[1]]
  mat_sg <- prospectr::savitzkyGolay(mat, m = sg_m, p = sg_p, w = sg_w)
  list(orig = df,
       sg   = data.frame(id = df[[1]], as.data.frame(mat_sg)))
}

resample_spectra <- function(df, wl_from, wl_to) {
  mat <- apply(as.matrix(df[, -1]), 2, as.numeric)
  colnames(mat) <- wl_from
  mat_rs <- prospectr::resample(mat, wl_from, wl_to, interpol = "spline")
  out <- data.frame(id = df[[1]], as.data.frame(mat_rs))
  colnames(out)[-1] <- as.character(wl_to)
  out
}

# 2.5  Spectral standardisation functions (Indices 1–6) -----------------------

## Index 1: Direct Standardisation (DS)
perform_ds <- function(target_mat, source_mat, intercept = TRUE) {
  X <- if (intercept) cbind(1, source_mat) else source_mat
  B <- qr.solve(X, target_mat)
  list(B = B, intercept = intercept)
}
apply_ds <- function(source_mat, ds_obj) {
  X <- if (ds_obj$intercept) cbind(1, source_mat) else source_mat
  X %*% ds_obj$B
}

## Index 2: Global CF at lambda_ref
perform_cf_global <- function(target_mat, source_mat, wl, lambda_ref = 1400,
                               eps = 1e-8, robust = TRUE) {
  idx        <- which.min(abs(wl - lambda_ref))
  target_ref <- target_mat[, idx];  source_ref <- source_mat[, idx]
  num <- if (robust) median(target_ref, na.rm = TRUE) else mean(target_ref, na.rm = TRUE)
  den <- if (robust) median(source_ref, na.rm = TRUE) else mean(source_ref, na.rm = TRUE)
  if (is.na(den) || abs(den) < eps) stop("CF denominator too small.")
  list(cf = num / den, lambda_ref = lambda_ref, idx = idx, robust = robust)
}
apply_cf_global <- function(source_mat, cf_obj) source_mat * cf_obj$cf

## Index 3: Band-wise CF using Lucky Bay standard
perform_cf_lb <- function(target_lb_vec, source_lb_vec, eps = 1e-8) {
  list(cf_vec = target_lb_vec / pmax(source_lb_vec, eps))
}
apply_cf_lb <- function(source_mat, cf_obj) sweep(source_mat, 2, cf_obj$cf_vec, "*")

## Index 4: Piecewise DS (PDS)
perform_pds <- function(target, source, MWsize = 2, Ncomp = 2) {
  i <- MWsize;  k <- i - 1
  P <- matrix(0, nrow = ncol(target), ncol = ncol(target) - 2 * i + 2)
  Intercept <- numeric(0)
  while (i <= (ncol(target) - k)) {
    fit    <- pls::plsr(target[, i] ~ as.matrix(source[, (i-k):(i+k)]),
                        ncomp = Ncomp, scale = FALSE, method = "oscorespls")
    coef_r <- as.numeric(coef(fit, ncomp = Ncomp, intercept = TRUE))
    Intercept <- c(Intercept, coef_r[1])
    P[(i-k):(i+k), i-k] <- coef_r[-1]
    i <- i + 1
  }
  cat("\n")
  P_full <- cbind(matrix(0, ncol(target), k), P, matrix(0, ncol(target), k))
  Intercept_full <- c(rep(0, k), Intercept, rep(0, k))
  list(P = P_full, Intercept = Intercept_full)
}
apply_pds <- function(source_mat, pds_obj) {
  res <- source_mat %*% pds_obj$P
  sweep(res, 2, pds_obj$Intercept, "+")
}

## Index 5: CF global + DS
apply_cf_then_ds <- function(source_mat, cf_obj, ds_obj) {
  apply_ds(apply_cf_global(source_mat, cf_obj), ds_obj)
}

## Index 6: Regularised Band-wise Regression (RBR) — novel method
# Rationale: DS fits ~20,000 parameters simultaneously → massive overfitting.
# RBR fits one Ridge regression per spectral band independently (1 intercept +
# 1 slope each), with Ridge shrinkage controlled by lambda. This gives only
# 282 effective parameters, is sample-specific unlike CF methods, and
# generalises far better than full DS on small paired calibration sets.
# Lambda is selected automatically via leave-one-out CV if not provided.
perform_rbr <- function(target_mat, source_mat, lambda = NULL) {
  
  n_bands <- ncol(target_mat)
  coefs   <- matrix(NA, nrow = 2, ncol = n_bands)  # row 1 = intercept, row 2 = slope
  lambda_used <- numeric(n_bands)
  
  # For a single predictor, Ridge solution is analytic:
  # slope     = sum((x - x_mean)(y - y_mean)) / (sum((x - x_mean)^2) + lambda)
  # intercept = y_mean - slope * x_mean
  # lambda = 0 recovers OLS; larger lambda shrinks slope toward 0.
  
  # Auto-select lambda via LOO-CV on first band if not provided
  if (is.null(lambda)) {
    message("  RBR: selecting lambda via LOO-CV on band 1...")
    y_cv <- target_mat[, 1]
    x_cv <- source_mat[, 1]
    lambdas  <- 10^seq(-4, 1, length.out = 30)
    rmse_cv  <- sapply(lambdas, function(lam) {
      n     <- length(y_cv)
      xc    <- x_cv - mean(x_cv);  yc <- y_cv - mean(y_cv)
      denom <- sum(xc^2) + lam
      # LOO residuals via hat matrix shortcut
      h     <- xc^2 / denom
      slope <- sum(xc * yc) / denom
      resid <- yc - slope * xc
      sqrt(mean((resid / (1 - h))^2))
    })
    lambda <- lambdas[which.min(rmse_cv)]
    message(sprintf("  RBR: optimal lambda = %.5f", lambda))
  }
  
  for (j in seq_len(n_bands)) {
    y      <- target_mat[, j]
    x      <- source_mat[, j]
    xc     <- x - mean(x);  yc <- y - mean(y)
    slope  <- sum(xc * yc) / (sum(xc^2) + lambda)
    interc <- mean(y) - slope * mean(x)
    coefs[1, j]    <- interc
    coefs[2, j]    <- slope
    lambda_used[j] <- lambda
    if (j %% 20 == 0 || j == n_bands)
      cat("\r  RBR progress:", round(j / n_bands * 100), "%")
  }
  cat("\n")
  
  list(coefs       = coefs,
       lambda_used = lambda_used,
       lambda_mean = mean(lambda_used))
}

apply_rbr <- function(source_mat, rbr_obj) {
  # Apply band-wise: corrected[i,j] = intercept[j] + slope[j] * source[i,j]
  out <- sweep(source_mat, 2, rbr_obj$coefs[2, ], "*")   # multiply by slopes
  out <- sweep(out,        2, rbr_obj$coefs[1, ], "+")   # add intercepts
  out
}

## RBR diagnostic plot: fitted slopes and intercepts per band
plot_rbr_coefs <- function(rbr_obj, wl) {
  data.frame(
    wavelength = wl,
    Intercept  = rbr_obj$coefs[1, ],
    Slope      = rbr_obj$coefs[2, ],
    Lambda     = rbr_obj$lambda_used
  ) %>%
    tidyr::pivot_longer(-wavelength, names_to = "parameter", values_to = "value") %>%
    mutate(parameter = factor(parameter, levels = c("Slope", "Intercept", "Lambda"))) %>%
    ggplot(aes(x = wavelength, y = value)) +
    geom_line(colour = "#C0392B", linewidth = 0.6) +
    geom_hline(data = data.frame(parameter = factor(c("Slope","Intercept","Lambda"),
                                                     levels = c("Slope","Intercept","Lambda")),
                                  yint = c(1, 0, NA)),
               aes(yintercept = yint), colour = "grey50",
               linetype = "dotted", linewidth = 0.4, na.rm = TRUE) +
    facet_wrap(~parameter, scales = "free_y", ncol = 1) +
    scale_x_continuous(breaks = seq(950, 1650, 100),
                       expand = expansion(mult = 0.01)) +
    labs(x = "Wavelength (nm)", y = "Value",
         title    = "Index 6 – RBR: band-wise Ridge regression coefficients",
         subtitle = sprintf("Mean lambda = %.4f | Slope ~1 and Intercept ~0 = good alignment",
                            rbr_obj$lambda_mean)) +
    theme_pub(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# 2.6  Model fitting ------------------------------------------------------------
fit_plsr <- function(x_cal, y_cal, ncomp = 15, n_cores = 1) {
  df  <- data.frame(y = y_cal, x_cal)
  cl  <- parallel::makeCluster(n_cores, type = "PSOCK")
  pls::pls.options(parallel = cl)
  m   <- pls::plsr(y ~ ., ncomp = ncomp, data = df,
                   validation = "CV", scale = FALSE)
  parallel::stopCluster(cl)
  press     <- m$validation$PRESS[1, ]
  ncomp_opt <- which.min(press)
  list(model = m, ncomp_opt = ncomp_opt)
}
predict_plsr <- function(fit_obj, x_new) {
  as.numeric(stats::predict(fit_obj$model, newdata = x_new,
                            ncomp = fit_obj$ncomp_opt)[, 1, 1])
}
fit_cubist <- function(x_cal, y_cal, folds = 10, repeats = 10) {
  ctrl <- caret::trainControl(method = "repeatedcv",
                              number = folds, repeats = repeats)
  grid <- expand.grid(committees = c(1, 5, 10, 20, 50, 100),
                      neighbors  = c(0, 1, 5, 9))
  caret::train(x = x_cal, y = y_cal, method = "cubist",
               trControl = ctrl, tuneGrid = grid)
}
predict_cubist <- function(fit_obj, x_new) {
  as.numeric(stats::predict(fit_obj, newdata = x_new))
}


# =============================================================================
# SECTION 3 - DATA IMPORT
# =============================================================================

soil <- read.csv(file.path(cfg$dir_data, "db_soil_analysis_all.csv"),
                 stringsAsFactors = FALSE) %>%
  dplyr::mutate(across(c(country, landuse), as.factor))
message(sprintf("  Soil data: %d observations, %d variables", nrow(soil), ncol(soil)))

asd_raw <- read.csv(file.path(cfg$dir_data, "asd_vnir_raw_all.csv"),
                    stringsAsFactors = FALSE)
colnames(asd_raw)[1] <- "id"
colnames(asd_raw)     <- gsub("X", "", colnames(asd_raw))
asd_raw$id            <- as.numeric(gsub("sample_|\\.asd", "", asd_raw$id))
asd_wl                <- as.numeric(colnames(asd_raw)[-1])
message(sprintf("  ASD: %d samples, %d bands (%.0f–%.0f nm)",
                nrow(asd_raw), length(asd_wl), min(asd_wl), max(asd_wl)))

nirv_raw <- read.csv(file.path(cfg$dir_data, "isc_nir_raw_all.csv"),
                     stringsAsFactors = FALSE)
colnames(nirv_raw) <- gsub("X", "", colnames(nirv_raw))
nirv_wl            <- as.numeric(colnames(nirv_raw)[-1])
message(sprintf("  NIRVascan: %d samples, %d bands (%.0f–%.0f nm)",
                nrow(nirv_raw), length(nirv_wl), min(nirv_wl), max(nirv_wl)))


# =============================================================================
# SECTION 4 - SPECTRAL PRE-PROCESSING
# =============================================================================

wl_target <- seq(cfg$wl_min, cfg$wl_max, by = cfg$wl_step)

asd_pp    <- preprocess_asd(asd_raw, asd_wl, splice = cfg$asd_splice,
                             interp = cfg$asd_interp,
                             sg_m = cfg$sg_m, sg_p = cfg$sg_p, sg_w = cfg$sg_w)

asd_nir_orig <- resample_spectra(asd_pp$orig, asd_wl, wl_target)
wl_sg_asd    <- as.numeric(gsub("X", "", colnames(asd_pp$sg)[-1]))
asd_nir_sg   <- resample_spectra(asd_pp$sg, wl_sg_asd, wl_target)

nirv_pp      <- preprocess_nirvascan(nirv_raw, sg_m = cfg$sg_m,
                                      sg_p = cfg$sg_p, sg_w = cfg$sg_w)
nirv_nir_orig <- resample_spectra(nirv_pp$orig, nirv_wl, wl_target)
wl_sg_nirv    <- as.numeric(gsub("X", "", colnames(nirv_pp$sg)[-1]))
nirv_nir_sg   <- resample_spectra(nirv_pp$sg, wl_sg_nirv, wl_target)

message("  Pre-processing complete.")

# --- Spectral plots (Section 4) ---
ggsave(file.path(cfg$dir_figs, "asd_spectra.png"),
  bind_rows(
    asd_pp$orig %>% filter(!is.na(id)) %>% mutate(pp = "Original (splice corrected)"),
    asd_pp$sg   %>% filter(!is.na(id)) %>% mutate(pp = "Savitzky-Golay smoothed")
  ) %>%
    pivot_longer(-c(id, pp), names_to = "wl", values_to = "ref") %>%
    mutate(wl = as.numeric(gsub("X", "", wl))) %>% filter(!is.na(wl)) %>%
    group_by(pp, wl) %>%
    summarise(mn = mean(ref, na.rm=T), lo = min(ref, na.rm=T),
              hi = max(ref, na.rm=T), .groups="drop") %>%
    mutate(pp = factor(pp, levels = c("Original (splice corrected)",
                                      "Savitzky-Golay smoothed"))) %>%
    ggplot(aes(x = wl)) +
    geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey60", alpha=0.25) +
    geom_line(aes(y=mn), colour="black", linewidth=0.55) +
    facet_wrap(~pp, ncol=2, scales="free_y") +
    scale_x_continuous(breaks=seq(400,2500,200), expand=expansion(mult=0.01)) +
    labs(x="Wavelength (nm)", y="Reflectance Factor",
         title="ASD vis-NIR spectra (400\u20132450 nm)",
         subtitle="Mean (solid line) \u00b1 full range (shaded ribbon); n = 253 soil samples") +
    theme_pub(11) + theme(axis.text.x=element_text(angle=45, hjust=1)),
  width=10, height=5, dpi=300, units="in")

ggsave(file.path(cfg$dir_figs, "asd_vs_nirvascan_spectra.png"),
  bind_rows(
    asd_nir_orig  %>% filter(!is.na(id)) %>% mutate(instrument="ASD",       pp="No pre-processing"),
    nirv_nir_orig %>% filter(!is.na(id)) %>% mutate(instrument="NIRVascan", pp="No pre-processing"),
    asd_nir_sg    %>% filter(!is.na(id)) %>% mutate(instrument="ASD",       pp="Savitzky-Golay"),
    nirv_nir_sg   %>% filter(!is.na(id)) %>% mutate(instrument="NIRVascan", pp="Savitzky-Golay")
  ) %>%
    pivot_longer(-c(id,instrument,pp), names_to="wl", values_to="ref") %>%
    mutate(wl=as.numeric(gsub("X","",wl))) %>% filter(!is.na(wl)) %>%
    mutate(instrument=factor(instrument,levels=c("ASD","NIRVascan")),
           pp=factor(pp,levels=c("No pre-processing","Savitzky-Golay"))) %>%
    group_by(instrument,pp,wl) %>%
    summarise(mn=mean(ref,na.rm=T), lo=min(ref,na.rm=T),
              hi=max(ref,na.rm=T), .groups="drop") %>%
    ggplot(aes(x=wl)) +
    geom_ribbon(aes(ymin=lo,ymax=hi), fill="grey60", alpha=0.25) +
    geom_line(aes(y=mn), colour="black", linewidth=0.55) +
    facet_grid(pp~instrument, scales="free_y") +
    scale_x_continuous(breaks=seq(950,1650,100), expand=expansion(mult=0.01)) +
    labs(x="Wavelength (nm)", y="Reflectance Factor",
         title="Resampled NIR spectra: ASD vs. NIRVascan (950\u20131650 nm, 5 nm)",
         subtitle="Mean (solid line) \u00b1 full range (shaded ribbon); n = 254") +
    theme_pub(11) + theme(axis.text.x=element_text(angle=45, hjust=1)),
  width=10, height=5, dpi=300, units="in")

message("  Section 4 plots saved.")

# ---- Helper: summarise — filter NAs that cause ribbon gaps -----------------
# (a) asd_pp$orig    → still has 350-399 & 2451-2500 nm noise → ADD wl filter
# (b) asd_pp$sg      → same issue → ADD wl filter
# (c) nirv_pp$orig   → still has 900-949 & 1651-1700 nm noise → use nirv_nir_orig
# (d) nirv_pp$sg     → same issue → use nirv_nir_sg
# (e) asd_nir_orig   → already 950-1650 nm, 5 nm ✓
# (f) asd_nir_sg     → already 950-1650 nm, 5 nm ✓

# Updated spec_to_long with optional wavelength range filter
spec_to_long <- function(df, wl_min = -Inf, wl_max = Inf) {
  df %>%
    dplyr::filter(!is.na(id), id != lb_id) %>%
    tidyr::pivot_longer(-id, names_to = "wl", values_to = "ref") %>%
    dplyr::mutate(wl = as.numeric(gsub("X", "", wl))) %>%
    dplyr::filter(!is.na(wl), !is.na(ref),
                  wl >= wl_min, wl <= wl_max) %>%   # <-- trim noisy edges
    dplyr::group_by(wl) %>%
    dplyr::summarise(mn = mean(ref, na.rm = TRUE),
                     lo = min(ref,  na.rm = TRUE),
                     hi = max(ref,  na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::filter(is.finite(mn), is.finite(lo), is.finite(hi))
}

# (a) ASD splice corrected — trim to 400-2450 nm
s_asd_orig   <- spec_to_long(asd_pp$orig,   wl_min = 400, wl_max = 2450)

# (b) ASD Savitzky-Golay — trim to 400-2450 nm
s_asd_sg     <- spec_to_long(asd_pp$sg,     wl_min = 400, wl_max = 2450)

# (c) NIRVascan original — use resampled object (already 950-1650 nm, 5 nm)
s_nirv_orig  <- spec_to_long(nirv_nir_orig, wl_min = 950, wl_max = 1650)

# (d) NIRVascan Savitzky-Golay — use resampled SG object (already 950-1650 nm, 5 nm)
s_nirv_sg    <- spec_to_long(nirv_nir_sg,   wl_min = 950, wl_max = 1650)

# (e) ASD NIR splice corrected — already 950-1650 nm, 5 nm
s_asd_5nm    <- spec_to_long(asd_nir_orig,  wl_min = 950, wl_max = 1650)

# (f) ASD NIR Savitzky-Golay — already 950-1650 nm, 5 nm
s_asd_sg_5nm <- spec_to_long(asd_nir_sg,    wl_min = 950, wl_max = 1650)

# Verify ranges
message(sprintf("(a) ASD orig:     %.0f-%.0f nm", min(s_asd_orig$wl),   max(s_asd_orig$wl)))
message(sprintf("(b) ASD SG:       %.0f-%.0f nm", min(s_asd_sg$wl),     max(s_asd_sg$wl)))
message(sprintf("(c) NIRVascan:    %.0f-%.0f nm", min(s_nirv_orig$wl),  max(s_nirv_orig$wl)))
message(sprintf("(d) NIRVascan SG: %.0f-%.0f nm", min(s_nirv_sg$wl),    max(s_nirv_sg$wl)))
message(sprintf("(e) ASD 5nm:      %.0f-%.0f nm", min(s_asd_5nm$wl),    max(s_asd_5nm$wl)))
message(sprintf("(f) ASD SG 5nm:   %.0f-%.0f nm", min(s_asd_sg_5nm$wl), max(s_asd_sg_5nm$wl)))

# ---- Generic panel builder -------------------------------------------------
make_spec_panel <- function(df, x_breaks, tag, title, y_lim = c(0, 1)) {
  ggplot(df, aes(x = wl)) +
    geom_ribbon(aes(ymin = lo, ymax = hi),
                fill  = "grey60",
                alpha = 0.25,
                na.rm = TRUE) +          # <-- add this
    geom_line(aes(y = mn),
              colour    = "black",
              linewidth = 0.55,
              na.rm     = TRUE) +        # <-- add this
    scale_x_continuous(breaks = x_breaks,
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = y_lim,
                       expand = expansion(mult = 0.02)) +
    labs(x = "Wavelength (nm)",
         y = "Reflectance Factor",
         tag   = tag,
         title = title) +
    theme_pub(10) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title   = element_text(size = 9),
      plot.tag     = element_text(face = "bold", size = 12),
      plot.title   = element_text(face = "bold", size = 9, hjust = 0.5),
      plot.margin  = margin(6, 10, 6, 6)
    )
}


# ---- ASD panels: full vis-NIR range (400-2450 nm) --------------------------
brk_full <- seq(400, 2500, 200)

pa <- make_spec_panel(s_asd_orig, brk_full, "(a)",
                      "ASD splice corrected (400\u20132450 nm, 1 nm)",
                      y_lim = c(0, 1.2))
pb <- make_spec_panel(s_asd_sg,   brk_full, "(b)",
                      "ASD Savitzky-Golay (400\u20132450 nm, 1 nm)",
                      y_lim = c(0, 1.2))

# ---- NIRVascan + resampled panels: NIR range (950-1650 nm) -----------------
brk_nir  <- seq(950, 1650, 100)
ylim_nir <- c(0, 1)

pc <- make_spec_panel(s_nirv_orig,  brk_nir, "(c)",
                      "NIRVascan original (950\u20131650 nm, 5 nm)")
pd <- make_spec_panel(s_nirv_sg,    brk_nir, "(d)",
                      "NIRVascan Savitzky-Golay (950\u20131650 nm, 5 nm)")
pe <- make_spec_panel(s_asd_5nm,    brk_nir, "(e)",
                      "ASD NIR splice corrected (950\u20131650 nm, 5 nm)")
pf <- make_spec_panel(s_asd_sg_5nm, brk_nir, "(f)",
                      "ASD NIR Savitzky-Golay (950\u20131650 nm, 5 nm)")
#  ^^ was NIRVascan resampled — now correctly ASD NIR SG-filtered
# ---- Common subtitle annotation --------------------------------------------
subtitle_grob <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Mean (solid line) \u00b1 full range (shaded ribbon); n = 253 soil samples",
    fontface = "plain", size = 8.5, colour = "grey40", hjust = 0, x = 0.01
  )

# ---- Assemble 3-row layout -------------------------------------------------
# Row 1: (a) and (b) — full vis-NIR ASD
# Row 2: (c) and (d) — NIRVascan original + SG
# Row 3: (e) and (f) — ASD and NIRVascan both at 5 nm

row1 <- cowplot::plot_grid(pa, pb, ncol = 2, align = "hv", axis = "tblr")
row2 <- cowplot::plot_grid(pc, pd, ncol = 2, align = "hv", axis = "tblr")
row3 <- cowplot::plot_grid(pe, pf, ncol = 2, align = "hv", axis = "tblr")

# Row labels
row1_lab <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("ASD (400\u20132450 nm)",
                                          fontface = "bold", size = 9, angle = 90),
  row1, ncol = 2, rel_widths = c(0.03, 1))

row2_lab <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("NIRVascan (950\u20131650 nm)",
                                          fontface = "bold", size = 9, angle = 90),
  row2, ncol = 2, rel_widths = c(0.03, 1))

row3_lab <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("Resampled NIR (950\u20131650 nm, 5 nm)",
                                          fontface = "bold", size = 9, angle = 90),
  row3, ncol = 2, rel_widths = c(0.03, 1))

fig_spec <- cowplot::plot_grid(
  pa, pb,
  pc, pd,
  pe, pf,
  ncol  = 2,
  align = "hv",
  axis  = "tblr"
)

fig_spec_final <- cowplot::plot_grid(
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Mean (solid line) \u00b1 full range (shaded ribbon); n = 253 soil samples",
      fontface = "plain", colour = "grey40", size = 9, hjust = 0, x = 0.01
    ),
  fig_spec,
  ncol        = 1,
  rel_heights = c(0.03, 1)
)

ggsave(file.path(cfg$dir_figs, "Fig2_spectral_preprocessing.jpeg"),
       fig_spec, width = 11, height = 13, dpi = 300, units = "in")
message("  Fig2 spectral pre-processing saved.")

# =============================================================================
# SECTION 5 - CALIBRATION / VALIDATION SPLIT (cLHS on NIRVascan PCA scores)
# =============================================================================

# Define Lucky Bay ID once — used here and in Section 6
lb_id <- min(nirv_nir_orig$id)   # id = 0 (standard, no lab measurements)

# ---- cLHS on soil samples only (exclude Lucky Bay) -------------------------
set.seed(cfg$random_seed)
nirv_for_clhs <- nirv_nir_orig %>% filter(id != lb_id)   # 253 soil samples
head(nirv_for_clhs[1:6], 5)
nirv_mat   <- apply(as.matrix(nirv_for_clhs[, -1]), 2, as.numeric)
pca_scores <- prcomp(nirv_mat, scale. = TRUE, center = TRUE)$x[, 1:10]
n_cal      <- round(nrow(nirv_for_clhs) * cfg$cal_frac)
idx        <- clhs::clhs(as.data.frame(pca_scores), size = n_cal,
                         iter = cfg$clhs_iter, simple = TRUE, progress = FALSE)
id_cal <- nirv_for_clhs$id[idx]
id_val <- nirv_for_clhs$id[-idx]
message(sprintf("  Calibration: %d | Validation: %d", length(id_cal), length(id_val)))

# Variance explained — computed on same 253-sample matrix
pca_var <- prcomp(nirv_mat, scale. = TRUE, center = TRUE)
var_exp <- round(summary(pca_var)$importance[2, 1:3] * 100, 1)

# PCA plot data — 253 soil samples only
pca_df  <- as.data.frame(pca_scores) %>%
  mutate(id  = nirv_for_clhs$id,
         set = if_else(id %in% id_cal, "Calibration", "Validation"))

pca_plot <- function(df, x, y, var_exp) {
  xi <- as.integer(gsub("PC","",x)); yi <- as.integer(gsub("PC","",y))
  ggplot(df, aes(.data[[x]], .data[[y]], colour=set, shape=set)) +
    geom_point(size=2, alpha=0.8) +
    stat_ellipse(aes(group=set), linewidth=0.5, linetype="dashed") +
    scale_colour_manual(values=c(Calibration="#1B4F72", Validation="#C0392B"), name=NULL) +
    scale_shape_manual(values=c(Calibration=16, Validation=17), name=NULL) +
    labs(x=sprintf("PC%d (%.1f%%)",xi,var_exp[xi]),
         y=sprintf("PC%d (%.1f%%)",yi,var_exp[yi])) +
    theme_pub(11) + theme(legend.position="bottom", legend.key.size=unit(0.4,"cm"))
}

ggsave(file.path(cfg$dir_figs, "pcs_clhs_selection.png"),
       (pca_plot(pca_df,"PC1","PC2",var_exp) + pca_plot(pca_df,"PC1","PC3",var_exp)) +
         plot_layout(guides="collect") +
         plot_annotation(
           title    = "cLHS calibration/validation split",
           subtitle = sprintf("Calibration: n = %d | Validation: n = %d (Lucky Bay standard excluded from split)",
                              length(id_cal), length(id_val)),
           theme    = theme(plot.title    = element_text(face="bold", size=11),
                            plot.subtitle = element_text(colour="grey40", size=9))) &
         theme(legend.position="bottom"),
       width=10, height=5, dpi=300, units="in")

# ---- Partition all spectral datasets (Lucky Bay excluded via id_cal/id_val) -
partition_spec <- function(df, id_cal, id_val) {
  list(cal = df[df$id %in% id_cal, ], val = df[df$id %in% id_val, ])
}
spec <- list(
  asd_orig  = partition_spec(asd_nir_orig,  id_cal, id_val),
  asd_sg    = partition_spec(asd_nir_sg,    id_cal, id_val),
  nirv_orig = partition_spec(nirv_nir_orig, id_cal, id_val),
  nirv_sg   = partition_spec(nirv_nir_sg,   id_cal, id_val)
)


# =============================================================================
# SECTION 6 - SPECTRAL STANDARDISATION (NIRVascan → ASD NIR space)
# =============================================================================

message("\n[6] Fitting spectral standardisation models (Indices 1–6)...")

# ---- Step 1: Paired calibration matrices -----------------------------------
paired_ids <- intersect(spec$asd_orig$cal$id, spec$nirv_orig$cal$id)
get_mat    <- function(df, ids) {
  m <- apply(as.matrix(df[df$id %in% ids, -1]), 2, as.numeric)
  colnames(m) <- as.character(wl_target);  m
}
target_orig_cal <- get_mat(spec$asd_orig$cal,  paired_ids)
source_orig_cal <- get_mat(spec$nirv_orig$cal, paired_ids)
target_sg_cal   <- get_mat(spec$asd_sg$cal,    paired_ids)
source_sg_cal   <- get_mat(spec$nirv_sg$cal,   paired_ids)

# ---- Step 2: Lucky Bay reference vectors -----------------------------------
lb_idx_target_orig <- as.numeric(asd_nir_orig[asd_nir_orig$id   == min(asd_nir_orig$id),   -1])
lb_idx_source_orig <- as.numeric(nirv_nir_orig[nirv_nir_orig$id == min(nirv_nir_orig$id), -1])
lb_idx_target_sg   <- as.numeric(asd_nir_sg[asd_nir_sg$id       == min(asd_nir_sg$id),     -1])
lb_idx_source_sg   <- as.numeric(nirv_nir_sg[nirv_nir_sg$id     == min(nirv_nir_sg$id),   -1])

# ---- Step 3: Lucky Bay offset plot -----------------------------------------
lb_ribbon <- data.frame(wavelength = wl_target,
                         asd = lb_idx_target_orig, nirv = lb_idx_source_orig,
                         asd_sg = lb_idx_target_sg, nirv_sg = lb_idx_source_sg)

ggsave(file.path(cfg$dir_figs, "06a_lucky_bay_ASD_vs_NIRvascan.png"),
  ggplot(lb_ribbon, aes(x = wavelength)) +
    geom_ribbon(aes(ymin=nirv, ymax=asd, fill="Instrument offset"), alpha=0.25) +
    geom_line(aes(y=asd,     colour="ASD",      linetype="No pre-processing"), linewidth=0.7) +
    geom_line(aes(y=asd_sg,  colour="ASD",      linetype="Savitzky-Golay"),    linewidth=0.5) +
    geom_line(aes(y=nirv,    colour="NIRVascan", linetype="No pre-processing"), linewidth=0.7) +
    geom_line(aes(y=nirv_sg, colour="NIRVascan", linetype="Savitzky-Golay"),    linewidth=0.5) +
    geom_hline(yintercept=1, colour="grey50", linetype="dotted", linewidth=0.4) +
    annotate("text", x=max(wl_target)-30, y=1.002,
             label="Spectralon ref.", colour="grey50", hjust=1, size=3) +
    annotate("segment", x=1300, xend=1300,
             y    = lb_ribbon$nirv[lb_ribbon$wavelength==1300],
             yend = lb_ribbon$asd[lb_ribbon$wavelength==1300],
             arrow=arrow(ends="both", length=unit(0.2,"cm"), type="closed"),
             colour="grey20", linewidth=0.7) +
    annotate("label", x=1310,
             y=mean(c(lb_ribbon$nirv[lb_ribbon$wavelength==1300],
                      lb_ribbon$asd[lb_ribbon$wavelength==1300])),
             label=sprintf("\u0394 = %.3f \u00b1 %.4f",
                           mean(lb_ribbon$asd - lb_ribbon$nirv),
                           sd(lb_ribbon$asd   - lb_ribbon$nirv)),
             hjust=0, size=4.5, fontface="bold", colour="grey20", fill="white",
             label.size=0.3, label.padding=unit(0.35,"lines")) +
    scale_colour_manual(values=c(ASD="#1B4F72", NIRVascan="#C0392B"), name="Instrument") +
    scale_fill_manual(values=c("Instrument offset"="grey60"), name=NULL) +
    scale_linetype_manual(values=c("No pre-processing"="solid","Savitzky-Golay"="dashed"),
                          name="Pre-processing") +
    scale_x_continuous(breaks=seq(950,1650,100), expand=expansion(mult=0.01)) +
    labs(x="Wavelength (nm)", y="Reflectance Factor",
         title="Lucky Bay soil spectral standard",
         subtitle="Shaded area = systematic offset | \u0394 = mean \u00b1 sd across 950\u20131650 nm") +
    theme_pub(11) + theme(axis.text.x=element_text(angle=45,hjust=1),
                          legend.position="bottom", legend.key.width=unit(1.2,"cm")),
  width=10, height=5, dpi=300, units="in")
message("  Lucky Bay plot saved.")

# ---- Step 4: Fit all 6 standardisation models ------------------------------

message("    [1/6] Direct Standardisation (DS) - orig")
idx1_orig <- perform_ds(target_orig_cal, source_orig_cal)
message("    [1/6] Direct Standardisation (DS) - sg")
idx1_sg   <- perform_ds(target_sg_cal,   source_sg_cal)

message("    [2/6] Global CF at lambda=1400 - orig")
idx2_orig <- perform_cf_global(target_orig_cal, source_orig_cal, wl_target, cfg$lb_lambda)
message("    [2/6] Global CF at lambda=1400 - sg")
idx2_sg   <- perform_cf_global(target_sg_cal,   source_sg_cal,   wl_target, cfg$lb_lambda)

message("    [3/6] Band-wise CF (Lucky Bay) - orig")
idx3_orig <- perform_cf_lb(lb_idx_target_orig, lb_idx_source_orig)
message("    [3/6] Band-wise CF (Lucky Bay) - sg")
idx3_sg   <- perform_cf_lb(lb_idx_target_sg,   lb_idx_source_sg)

message("    [4/6] Piecewise DS (PDS) - orig  [slowest step]")
idx4_orig <- perform_pds(target_orig_cal, source_orig_cal, MWsize = 2, Ncomp = 2)
message("    [4/6] Piecewise DS (PDS) - sg")
idx4_sg   <- perform_pds(target_sg_cal,   source_sg_cal,   MWsize = 2, Ncomp = 2)

message("    [5/6] CF global + DS - orig")
idx5_cf_orig <- perform_cf_global(target_orig_cal, source_orig_cal, wl_target, cfg$lb_lambda)
idx5_orig    <- list(cf = idx5_cf_orig,
                     ds = perform_ds(target_orig_cal,
                                     apply_cf_global(source_orig_cal, idx5_cf_orig)))
message("    [5/6] CF global + DS - sg")
idx5_cf_sg <- perform_cf_global(target_sg_cal, source_sg_cal, wl_target, cfg$lb_lambda)
idx5_sg    <- list(cf = idx5_cf_sg,
                   ds = perform_ds(target_sg_cal,
                                   apply_cf_global(source_sg_cal, idx5_cf_sg)))

message("    [6/6] Regularised Band-wise Regression (RBR) - orig  [slow: one Ridge CV per band]")
idx6_orig <- perform_rbr(target_orig_cal, source_orig_cal, lambda = cfg$rbr_lambda)
message(sprintf("      orig complete — mean lambda = %.4f", idx6_orig$lambda_mean))

message("    [6/6] Regularised Band-wise Regression (RBR) - sg")
idx6_sg   <- perform_rbr(target_sg_cal, source_sg_cal, lambda = cfg$rbr_lambda)
message(sprintf("      sg complete — mean lambda = %.4f", idx6_sg$lambda_mean))

# Save RBR coefficient diagnostic plot
ggsave(file.path(cfg$dir_figs, "06d_rbr_coefficients.png"),
       plot_rbr_coefs(idx6_orig, wl_target),
       width = 8, height = 9, dpi = 300, units = "in")
message("  RBR coefficient plot saved.")

std_models <- list(
  idx1_orig = idx1_orig, idx1_sg = idx1_sg,
  idx2_orig = idx2_orig, idx2_sg = idx2_sg,
  idx3_orig = idx3_orig, idx3_sg = idx3_sg,
  idx4_orig = idx4_orig, idx4_sg = idx4_sg,
  idx5_orig = idx5_orig, idx5_sg = idx5_sg,
  idx6_orig = idx6_orig, idx6_sg = idx6_sg   # NEW
)

saveRDS(std_models, file.path(cfg$dir_models, "spectral_std_models.rds"))
saveRDS(wl_target,  file.path(cfg$dir_models, "wavelength_grid.rds"))
message("  All 6 standardisation models fitted and saved.")

# ---- Step 5: Apply models to calibration set -------------------------------

corrected_list <- list(
  idx1_orig = apply_ds(source_orig_cal,         std_models$idx1_orig),
  idx1_sg   = apply_ds(source_sg_cal,           std_models$idx1_sg),
  idx2_orig = apply_cf_global(source_orig_cal,  std_models$idx2_orig),
  idx2_sg   = apply_cf_global(source_sg_cal,    std_models$idx2_sg),
  idx3_orig = apply_cf_lb(source_orig_cal,      std_models$idx3_orig),
  idx3_sg   = apply_cf_lb(source_sg_cal,        std_models$idx3_sg),
  idx4_orig = apply_pds(source_orig_cal,        std_models$idx4_orig),
  idx4_sg   = apply_pds(source_sg_cal,          std_models$idx4_sg),
  idx5_orig = apply_cf_then_ds(source_orig_cal, std_models$idx5_orig$cf, std_models$idx5_orig$ds),
  idx5_sg   = apply_cf_then_ds(source_sg_cal,   std_models$idx5_sg$cf,   std_models$idx5_sg$ds),
  idx6_orig = apply_rbr(source_orig_cal, std_models$idx6_orig),  # RBR
  idx6_sg   = apply_rbr(source_sg_cal,   std_models$idx6_sg)     # RBR
)

corrected_list <- lapply(corrected_list, function(m) {
  colnames(m) <- as.character(wl_target); m
})

# ---- Step 6: Build panel_lines for calibration (Indices 1–6) ---------------

n_idx   <- 6   # total number of indices
methods <- c("DS", "CF global", "CF band-wise (LB)",
             "Piecewise DS", "CF global + DS", "RBR (novel)")

panel_lines <- purrr::map_dfr(c("orig", "sg"), function(pp) {
  pp_label <- if (pp == "orig") "No pre-processing" else "Savitzky-Golay"
  target   <- if (pp == "orig") target_orig_cal else target_sg_cal
  source   <- if (pp == "orig") source_orig_cal else source_sg_cal

  purrr::map_dfr(1:n_idx, function(i) {
    key     <- paste0("idx", i, "_", pp)
    cor     <- corrected_list[[key]]
    wl_use  <- if (i == 4) wl_target[4:(length(wl_target) - 3)] else wl_target
    col_use <- as.character(wl_use)

    bind_rows(
      as.data.frame(target[, col_use]) %>% mutate(sample=seq_len(nrow(target)), source="Target (ASD)"),
      as.data.frame(source[, col_use]) %>% mutate(sample=seq_len(nrow(source)), source="Source (NIRVascan)"),
      as.data.frame(cor[,   col_use]) %>% mutate(sample=seq_len(nrow(cor)),    source="Corrected NIRVascan")
    ) %>%
      tidyr::pivot_longer(-c(sample, source), names_to="wavelength", values_to="reflectance") %>%
      mutate(wavelength = as.numeric(wavelength),
             index      = sprintf("Index %d", i),
             pp         = pp_label)
  })
}) %>%
  mutate(source = factor(source, levels = c("Target (ASD)",
                                             "Source (NIRVascan)",
                                             "Corrected NIRVascan")),
         index  = factor(index,  levels = paste("Index", 1:n_idx)),
         pp     = factor(pp,     levels = c("No pre-processing", "Savitzky-Golay")))

# ---- Step 7: Calibration RMSE + ranking ------------------------------------

rmse_df <- purrr::map_dfr(c("orig", "sg"), function(pp) {
  target <- if (pp == "orig") target_orig_cal else target_sg_cal
  source <- if (pp == "orig") source_orig_cal else source_sg_cal
  rmse_baseline <- sqrt(mean((target - source)^2, na.rm = TRUE))

  purrr::map_dfr(1:n_idx, function(i) {
    key       <- paste0("idx", i, "_", pp)
    corrected <- corrected_list[[key]]
    wl_use    <- if (i == 4) wl_target[4:(length(wl_target) - 3)] else wl_target
    col_use   <- as.character(wl_use)
    rmse_corr <- sqrt(mean((target[, col_use] - corrected[, col_use])^2, na.rm = TRUE))
    data.frame(
      Index           = sprintf("Index %d", i),
      Method          = methods[i],
      Pre_processing  = if (pp == "orig") "No pre-processing" else "Savitzky-Golay",
      RMSE_baseline   = round(rmse_baseline, 4),
      RMSE            = round(rmse_corr, 4),
      Improvement_pct = round((1 - rmse_corr / rmse_baseline) * 100, 1)
    )
  })
})

ranked_table <- rmse_df %>%
  group_by(Index, Method) %>%
  summarise(RMSE_baseline   = round(mean(RMSE_baseline), 4),
            RMSE_orig       = RMSE[Pre_processing == "No pre-processing"],
            RMSE_sg         = RMSE[Pre_processing == "Savitzky-Golay"],
            RMSE_mean       = round(mean(RMSE), 4),
            Improvement_pct = round(mean(Improvement_pct), 1),
            .groups = "drop") %>%
  arrange(RMSE_mean) %>%
  mutate(Rank = row_number(),
         Assessment = case_when(Improvement_pct >= 50 ~ "Good", TRUE ~ "Poor")) %>%
  dplyr::select(Rank, Index, Method, RMSE_baseline, RMSE_orig, RMSE_sg,
                RMSE_mean, Improvement_pct, Assessment)

print(ranked_table)
write.csv(ranked_table, file.path(cfg$dir_tables, "06c_spectral_std_ranking.csv"),
          row.names = FALSE)

# ---- Step 8: Build panel_sum with RMSE strip labels ------------------------

rmse_labels <- rmse_df %>%
  rename(index = Index) %>%
  group_by(index) %>%
  summarise(RMSE_mean       = mean(RMSE),
            Improvement_pct = mean(Improvement_pct),
            .groups = "drop") %>%
  mutate(index_label = sprintf("%s\nRMSE = %.4f | \u0394 = %.1f%%",
                               index, RMSE_mean, Improvement_pct))

panel_sum <- panel_lines %>%
  group_by(index, pp, source, wavelength) %>%
  summarise(mn = mean(reflectance, na.rm=T), lo = min(reflectance, na.rm=T),
            hi = max(reflectance, na.rm=T), .groups="drop") %>%
  left_join(rmse_labels, by = "index") %>%
  mutate(index_label = factor(index_label,
                              levels = rmse_labels$index_label[order(rmse_labels$index)]))

# ---- Step 9: Calibration plot (6 panels) -----------------------------------

fig_std <- ggplot(panel_sum, aes(x = wavelength)) +
  geom_ribbon(aes(ymin=lo, ymax=hi, fill=source), alpha=0.08) +
  geom_line(aes(y=mn, colour=source), linewidth=0.35) +
  facet_grid(pp ~ index_label, scales="free_y") +
  scale_colour_manual(values=c("Target (ASD)"="#0072B2",
                                "Source (NIRVascan)"="#E69F00",
                                "Corrected NIRVascan"="#009E73"), name=NULL) +
  scale_fill_manual(values=c("Target (ASD)"="#0072B2",
                              "Source (NIRVascan)"="#E69F00",
                              "Corrected NIRVascan"="#009E73"), name=NULL) +
  scale_x_continuous(breaks=seq(950,1650,200), expand=expansion(mult=0.01)) +
  scale_y_continuous(expand=expansion(mult=0.02)) +
  labs(x="Wavelength (nm)", y="Reflectance Factor",
       title="Spectral standardisation: NIRVascan \u2192 ASD NIR space (calibration set)",
       subtitle="Mean (solid line) \u00b1 full range (ribbon) \u00b7 strip = RMSE and % improvement over uncorrected NIRVascan") +
  theme_bw(base_size=10) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour="grey92", linewidth=0.3),
        panel.border=element_rect(colour="grey40"),
        strip.background=element_rect(fill="grey95"),
        strip.text=element_text(face="bold", size=7),
        strip.text.y=element_text(angle=0),
        axis.text=element_text(colour="black", size=8),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.title=element_text(colour="black", size=9),
        plot.title=element_text(face="bold", size=11),
        plot.subtitle=element_text(colour="grey40", size=8),
        legend.position="bottom", legend.key.width=unit(1.2,"cm")) +
  guides(fill=guide_legend(override.aes=list(alpha=0.4)),
         colour=guide_legend(override.aes=list(linewidth=1.2)))

ggsave(file.path(cfg$dir_figs, "06b_std_all_indices_cal.png"),
       fig_std, width=19, height=7, dpi=300, units="in")
message("  Calibration standardisation plot saved.")

# ---- Step 10: Apply transfer to validation set -----------------------------

apply_all_std <- function(nirv_df, pp_label, std_models, wl) {
  ids <- nirv_df$id
  m   <- apply(as.matrix(nirv_df[, -1]), 2, as.numeric)
  colnames(m) <- as.character(wl)
  sfx <- pp_label

  out <- list(
    std  = nirv_df,
    idx1 = data.frame(id=ids, as.data.frame(apply_ds(m, std_models[[paste0("idx1_",sfx)]]))),
    idx2 = data.frame(id=ids, as.data.frame(apply_cf_global(m, std_models[[paste0("idx2_",sfx)]]))),
    idx3 = data.frame(id=ids, as.data.frame(apply_cf_lb(m, std_models[[paste0("idx3_",sfx)]]))),
    idx4 = data.frame(id=ids, as.data.frame(apply_pds(m, std_models[[paste0("idx4_",sfx)]]))),
    idx5 = data.frame(id=ids, as.data.frame(
      apply_cf_then_ds(m, std_models[[paste0("idx5_",sfx)]]$cf,
                           std_models[[paste0("idx5_",sfx)]]$ds))),
    idx6 = data.frame(id=ids, as.data.frame(   # RBR
      apply_rbr(m, std_models[[paste0("idx6_",sfx)]])))
  )
  for (nm in names(out)) colnames(out[[nm]])[-1] <- as.character(wl)
  out
}

transferred <- list(
  val_orig = apply_all_std(spec$nirv_orig$val, "orig", std_models, wl_target),
  val_sg   = apply_all_std(spec$nirv_sg$val,   "sg",   std_models, wl_target)
)
message("  Spectral standardisation applied to validation set.")

# ---- Step 11: Validation spectral curves -----------------------------------

target_orig_val <- get_mat(spec$asd_orig$val,
                            intersect(spec$asd_orig$val$id, spec$nirv_orig$val$id))
target_sg_val   <- get_mat(spec$asd_sg$val,
                            intersect(spec$asd_sg$val$id,   spec$nirv_sg$val$id))

panel_lines_val <- purrr::map_dfr(c("orig", "sg"), function(pp) {
  pp_label       <- if (pp == "orig") "No pre-processing" else "Savitzky-Golay"
  target         <- if (pp == "orig") target_orig_val      else target_sg_val
  transferred_pp <- if (pp == "orig") transferred$val_orig  else transferred$val_sg
  source_mat     <- apply(as.matrix(transferred_pp$std[,-1]), 2, as.numeric)
  colnames(source_mat) <- as.character(wl_target)

  purrr::map_dfr(1:n_idx, function(i) {
    key     <- paste0("idx", i)
    cor_mat <- apply(as.matrix(transferred_pp[[key]][,-1]), 2, as.numeric)
    colnames(cor_mat) <- as.character(wl_target)
    wl_use  <- if (i == 4) wl_target[4:(length(wl_target)-3)] else wl_target
    col_use <- as.character(wl_use)

    bind_rows(
      as.data.frame(target[,     col_use]) %>% mutate(sample=seq_len(nrow(target)),     source="Target (ASD)"),
      as.data.frame(source_mat[, col_use]) %>% mutate(sample=seq_len(nrow(source_mat)), source="Source (NIRVascan)"),
      as.data.frame(cor_mat[,    col_use]) %>% mutate(sample=seq_len(nrow(cor_mat)),    source="Corrected NIRVascan")
    ) %>%
      tidyr::pivot_longer(-c(sample,source), names_to="wavelength", values_to="reflectance") %>%
      mutate(wavelength=as.numeric(wavelength), index=sprintf("Index %d",i), pp=pp_label)
  })
}) %>%
  mutate(source=factor(source,levels=c("Target (ASD)","Source (NIRVascan)","Corrected NIRVascan")),
         index=factor(index,levels=paste("Index",1:n_idx)),
         pp=factor(pp,levels=c("No pre-processing","Savitzky-Golay")))

rmse_val_df <- purrr::map_dfr(c("orig","sg"), function(pp) {
  target         <- if (pp == "orig") target_orig_val      else target_sg_val
  transferred_pp <- if (pp == "orig") transferred$val_orig  else transferred$val_sg
  source_mat     <- apply(as.matrix(transferred_pp$std[,-1]), 2, as.numeric)
  colnames(source_mat) <- as.character(wl_target)
  rmse_baseline  <- sqrt(mean((target - source_mat)^2, na.rm=TRUE))

  purrr::map_dfr(1:n_idx, function(i) {
    key     <- paste0("idx", i)
    cor_mat <- apply(as.matrix(transferred_pp[[key]][,-1]), 2, as.numeric)
    colnames(cor_mat) <- as.character(wl_target)
    wl_use  <- if (i == 4) wl_target[4:(length(wl_target)-3)] else wl_target
    col_use <- as.character(wl_use)
    rmse_corr       <- sqrt(mean((target[,col_use]-cor_mat[,col_use])^2, na.rm=TRUE))
    pct_improvement <- round((1 - rmse_corr / rmse_baseline) * 100, 1)
    data.frame(index=sprintf("Index %d",i),
               Pre_processing=if(pp=="orig") "No pre-processing" else "Savitzky-Golay",
               RMSE=round(rmse_corr,4), Improvement_pct=pct_improvement)
  })
})

rmse_val_labels <- rmse_val_df %>%
  group_by(index) %>%
  summarise(RMSE_mean=round(mean(RMSE),4),
            Improvement_pct=round(mean(Improvement_pct),1), .groups="drop") %>%
  mutate(index=as.character(index),
         index_label=sprintf("%s\nRMSE = %.4f | \u0394 = %.1f%%",
                             index, RMSE_mean, Improvement_pct))

panel_sum_val <- panel_lines_val %>%
  group_by(index,pp,source,wavelength) %>%
  summarise(mn=mean(reflectance,na.rm=T), lo=min(reflectance,na.rm=T),
            hi=max(reflectance,na.rm=T), .groups="drop") %>%
  mutate(index=as.character(index)) %>%
  left_join(rmse_val_labels, by="index") %>%
  mutate(index=factor(index,levels=paste("Index",1:n_idx)),
         index_label=factor(index_label,
                            levels=rmse_val_labels$index_label[order(rmse_val_labels$index)]))

fig_std_val <- ggplot(panel_sum_val, aes(x = wavelength)) +
  geom_ribbon(aes(ymin=lo, ymax=hi, fill=source), alpha=0.08) +
  geom_line(aes(y=mn, colour=source), linewidth=0.35) +
  facet_grid(pp ~ index_label, scales="free_y") +
  scale_colour_manual(values=c("Target (ASD)"="#0072B2",
                                "Source (NIRVascan)"="#E69F00",
                                "Corrected NIRVascan"="#009E73"), name=NULL) +
  scale_fill_manual(values=c("Target (ASD)"="#0072B2",
                              "Source (NIRVascan)"="#E69F00",
                              "Corrected NIRVascan"="#009E73"), name=NULL) +
  scale_x_continuous(breaks=seq(950,1650,200), expand=expansion(mult=0.01)) +
  scale_y_continuous(expand=expansion(mult=0.02)) +
  labs(x="Wavelength (nm)", y="Reflectance Factor",
       title="Spectral standardisation: NIRVascan \u2192 ASD NIR space (validation set)",
       subtitle="Mean (solid line) \u00b1 full range (ribbon) \u00b7 strip = validation RMSE and % improvement") +
  theme_bw(base_size=10) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour="grey92", linewidth=0.3),
        panel.border=element_rect(colour="grey40"),
        strip.background=element_rect(fill="grey95"),
        strip.text=element_text(face="bold", size=7),
        strip.text.y=element_text(angle=0),
        axis.text=element_text(colour="black", size=8),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.title=element_text(colour="black", size=9),
        plot.title=element_text(face="bold", size=11),
        plot.subtitle=element_text(colour="grey40", size=8),
        legend.position="bottom", legend.key.width=unit(1.2,"cm")) +
  guides(fill=guide_legend(override.aes=list(alpha=0.4)),
         colour=guide_legend(override.aes=list(linewidth=1.2)))

ggsave(file.path(cfg$dir_figs, "06c_std_all_indices_val.png"),
       fig_std_val, width=19, height=7, dpi=300, units="in")
message("  Validation standardisation plot saved.")

# ---- Step 12: Cal vs val comparison table ----------------------------------

comparison_table <- left_join(
  ranked_table %>%
    dplyr::select(Rank, Index, Method, RMSE_mean, Improvement_pct) %>%
    rename(RMSE_cal = RMSE_mean, Delta_cal = Improvement_pct),
  rmse_val_df %>%
    group_by(index) %>%
    summarise(RMSE_val  = round(mean(RMSE),4),
              Delta_val = round(mean(Improvement_pct),1), .groups="drop") %>%
    rename(Index = index),
  by = "Index"
) %>%
  mutate(
    Generalisation = case_when(
      Delta_val >= 50                  ~ "Good",
      Delta_val >= 25 & Delta_val < 50 ~ "Moderate",
      Delta_val >= 0  & Delta_val < 25 ~ "Weak",
      TRUE                             ~ "Overfits"
    ),
    RMSE_increase_pct = round((RMSE_val - RMSE_cal) / RMSE_cal * 100, 1)
  ) %>%
  arrange(desc(Delta_val))

print(comparison_table)
write.csv(comparison_table,
          file.path(cfg$dir_tables, "06d_spectral_std_cal_vs_val.csv"),
          row.names = FALSE)
message("  Cal vs val comparison table saved.")
message("\n  Section 6 complete. Best method based on validation: ",
        comparison_table$Index[1], " (", comparison_table$Method[1], ")",
        " — Delta_val = ", comparison_table$Delta_val[1], "%")

# =============================================================================
# SECTION 7 - MODELLING
# =============================================================================

message("\n[7] Fitting PLSR and Cubist models...")
summary(soil)
# Soil properties to model
targets <- list(
  tc   = list(col = "tc.perc", unit = "%", name = "Total Carbon"),
  tn   = list(col = "tn.perc", unit = "%", name = "Total Nitrogen"),
  ph   = list(col = "ph.water", unit = "adimensional", name = "pH (water)"),
  p_dl = list(col = "p.dl.mg_100g_ok", unit = "mg/100g", name = "Plant-available phosphorus")
)

# Spectral datasets for modelling
# Calibration always uses ASD (the target instrument)
# Validation uses the transferred NIRVascan spectra for index datasets
spec_datasets <- list(
  # --- Standard instrument comparisons ---
  asd_orig  = list(cal = spec$asd_orig$cal,  val = spec$asd_orig$val,
                   label = "ASD NIR (No pre-processing)"),
  asd_sg    = list(cal = spec$asd_sg$cal,    val = spec$asd_sg$val,
                   label = "ASD NIR (Savitzky-Golay)"),
  nirv_orig = list(cal = spec$nirv_orig$cal, val = spec$nirv_orig$val,
                   label = "NIRVascan NIR (No pre-processing)"),
  nirv_sg   = list(cal = spec$nirv_sg$cal,   val = spec$nirv_sg$val,
                   label = "NIRVascan NIR (Savitzky-Golay)"),
  
  # --- Index 1: Direct Standardisation (DS) ---
  idx1_orig = list(cal = spec$asd_orig$cal, val = transferred$val_orig$idx1,
                   label = "Index 1 – DS (No pre-processing)"),
  idx1_sg   = list(cal = spec$asd_sg$cal,   val = transferred$val_sg$idx1,
                   label = "Index 1 – DS (Savitzky-Golay)"),
  
  # --- Index 2: Global CF ---
  idx2_orig = list(cal = spec$asd_orig$cal, val = transferred$val_orig$idx2,
                   label = "Index 2 – CF global (No pre-processing)"),
  idx2_sg   = list(cal = spec$asd_sg$cal,   val = transferred$val_sg$idx2,
                   label = "Index 2 – CF global (Savitzky-Golay)"),
  
  # --- Index 3: Band-wise CF (Lucky Bay) ---
  idx3_orig = list(cal = spec$asd_orig$cal, val = transferred$val_orig$idx3,
                   label = "Index 3 – CF band-wise (No pre-processing)"),
  idx3_sg   = list(cal = spec$asd_sg$cal,   val = transferred$val_sg$idx3,
                   label = "Index 3 – CF band-wise (Savitzky-Golay)"),
  
  # --- Index 4: Piecewise DS (PDS) ---
  idx4_orig = list(cal = spec$asd_orig$cal, val = transferred$val_orig$idx4,
                   label = "Index 4 – PDS (No pre-processing)"),
  idx4_sg   = list(cal = spec$asd_sg$cal,   val = transferred$val_sg$idx4,
                   label = "Index 4 – PDS (Savitzky-Golay)"),
  
  # --- Index 5: CF global + DS ---
  idx5_orig = list(cal = spec$asd_orig$cal, val = transferred$val_orig$idx5,
                   label = "Index 5 – CF+DS (No pre-processing)"),
  idx5_sg   = list(cal = spec$asd_sg$cal,   val = transferred$val_sg$idx5,
                   label = "Index 5 – CF+DS (Savitzky-Golay)"),
  
  # --- Index 6: Regularised Band-wise Regression (RBR) — novel ---
  idx6_orig = list(cal = spec$asd_orig$cal, val = transferred$val_orig$idx6,
                   label = "Index 6 – RBR (No pre-processing)"),
  idx6_sg   = list(cal = spec$asd_sg$cal,   val = transferred$val_sg$idx6,
                   label = "Index 6 – RBR (Savitzky-Golay)")
)

# Storage
all_models  <- list()
all_metrics <- list()
all_preds   <- list()

for (tgt_nm in names(targets)) {
  tgt <- targets[[tgt_nm]]
  message(sprintf("  --- Property: %s ---", tgt$name))
  
  for (spec_nm in names(spec_datasets)) {
    sd    <- spec_datasets[[spec_nm]]
    label <- sd$label
    
    cal_df <- merge(soil[, c("id", tgt$col)], sd$cal, by = "id")
    val_df <- merge(soil[, c("id", tgt$col)], sd$val, by = "id")
    
    if (nrow(cal_df) < 10 || nrow(val_df) < 3) {
      message(sprintf("    Skipping %s [%s]: insufficient data.", tgt_nm, spec_nm))
      next
    }
    
    y_cal <- cal_df[[tgt$col]]
    x_cal <- apply(as.matrix(cal_df[, -(1:2)]), 2, as.numeric)
    
    y_val <- val_df[[tgt$col]]
    x_val <- apply(as.matrix(val_df[, -(1:2)]), 2, as.numeric)
    
    model_key <- paste(tgt_nm, spec_nm, sep = ".")
    
    # -- PLSR --
    message(sprintf("    PLSR: %s | %s", tgt$name, label))
    tryCatch({
      pls_fit  <- fit_plsr(x_cal, y_cal, ncomp = cfg$pls_ncomp, n_cores = n_cores)
      pls_pred <- predict_plsr(pls_fit, x_val)
      all_models[[paste0(model_key, ".pls")]]  <- pls_fit
      all_preds[[paste0(model_key, ".pls")]]   <- data.frame(id = val_df$id,
                                                             obs = y_val,
                                                             pred = pls_pred)
      all_metrics[[paste0(model_key, ".pls")]] <- cbind(
        data.frame(property = tgt$name, spectra = label, model = "PLSR"),
        goof(y_val, pls_pred))
    }, error = function(e) message("      PLSR error: ", e$message))
    
    # -- Cubist --
    message(sprintf("    Cubist: %s | %s", tgt$name, label))
    tryCatch({
      cub_fit  <- fit_cubist(x_cal, y_cal,
                             folds = cfg$cub_folds, repeats = cfg$cub_repeats)
      cub_pred <- predict_cubist(cub_fit, x_val)
      all_models[[paste0(model_key, ".cub")]]  <- cub_fit
      all_preds[[paste0(model_key, ".cub")]]   <- data.frame(id = val_df$id,
                                                             obs = y_val,
                                                             pred = cub_pred)
      all_metrics[[paste0(model_key, ".cub")]] <- cbind(
        data.frame(property = tgt$name, spectra = label, model = "Cubist"),
        goof(y_val, cub_pred))
    }, error = function(e) message("      Cubist error: ", e$message))
  }
}

saveRDS(all_models, file.path(cfg$dir_models, "all_fitted_models.rds"))
message("\n  All models saved to: ", cfg$dir_models)
beepr::beep(6)


# =============================================================================
# SECTION 8 - METRICS TABLE
# =============================================================================

message("\n[8] Compiling validation metrics...")

metrics_df <- dplyr::bind_rows(all_metrics) %>%
  dplyr::arrange(property, model, spectra) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(metrics_df,
          file.path(cfg$dir_tables, "validation_metrics_all.csv"),
          row.names = FALSE)
message("  Metrics table saved.")
print(metrics_df)


# =============================================================================
# SECTION 9 - PUBLICATION-QUALITY FIGURES
# =============================================================================

message("\n[9] Generating figures...")

# ---- 9.1  Obs vs. predicted: standard instruments -------------------------
ax_limits <- list(
  tc   = c(0, 60),
  tn   = c(0, 4),
  ph   = c(3, 9),
  p_dl = c(0, 35)
)

for (tgt_nm in names(targets)) {
  tgt  <- targets[[tgt_nm]]
  xlim <- ax_limits[[tgt_nm]]
  
  for (algo in c("pls", "cub")) {
    algo_label <- if (algo == "pls") "PLSR" else "Cubist"
    
    # --- 9.1a Standard instrument panel (2×2: ASD and NIRVascan × 2 pp) ---
    std_specs  <- c("asd_orig", "asd_sg", "nirv_orig", "nirv_sg")
    std_panels <- lapply(std_specs, function(sp) {
      key <- paste0(tgt_nm, ".", sp, ".", algo)
      if (is.null(all_preds[[key]])) return(NULL)
      prd <- all_preds[[key]]
      obs_pred_plot(pred = prd$pred, obs = prd$obs,
                    title = spec_datasets[[sp]]$label,
                    unit = tgt$unit, axis_lim = xlim)
    }) %>% Filter(Negate(is.null), .)
    
    if (length(std_panels) > 0) {
      fig <- patchwork::wrap_plots(std_panels, ncol = 2) +
        patchwork::plot_annotation(
          title = sprintf("%s – %s validation (standard instruments)", tgt$name, algo_label),
          theme = theme(plot.title = element_text(size = 11, face = "bold")))
      fn <- sprintf("07_%s_std_instruments_%s", tgt_nm, algo)
      ggsave(file.path(cfg$dir_figs, paste0(fn, ".pdf")),
             fig, width = 8, height = 8, dpi = 300)
      ggsave(file.path(cfg$dir_figs, paste0(fn, ".jpeg")),
             fig, width = 8, height = 8, dpi = 300)
    }
    
    # --- 9.1b All 6 indices panel (6 indices × 2 pp = 12 panels) ----------
    idx_specs <- c(
      paste0("idx", 1:6, "_orig"),
      paste0("idx", 1:6, "_sg")
    )
    
    idx_panels <- lapply(idx_specs, function(sp) {
      key <- paste0(tgt_nm, ".", sp, ".", algo)
      if (is.null(all_preds[[key]])) return(NULL)
      prd <- all_preds[[key]]
      obs_pred_plot(pred = prd$pred, obs = prd$obs,
                    title = spec_datasets[[sp]]$label,
                    unit = tgt$unit, axis_lim = xlim)
    }) %>% Filter(Negate(is.null), .)
    
    if (length(idx_panels) > 0) {
      fig_idx <- patchwork::wrap_plots(idx_panels, ncol = 4) +
        patchwork::plot_annotation(
          title = sprintf("%s – %s: all 6 standardisation indices", tgt$name, algo_label),
          theme = theme(plot.title = element_text(size = 11, face = "bold")))
      fn <- sprintf("08_%s_all_indices_%s", tgt_nm, algo)
      ggsave(file.path(cfg$dir_figs, paste0(fn, ".pdf")),
             fig_idx, width = 20, height = 10, dpi = 300)
      ggsave(file.path(cfg$dir_figs, paste0(fn, ".jpeg")),
             fig_idx, width = 20, height = 10, dpi = 300)
    }
  }
  message(sprintf("  Saved obs-vs-pred figures for: %s", tgt$name))
}

# ---- 9.2  Best index per property summary ----------------------------------
# Select best Cubist model per property based on RPIQ (higher = better)
best_summary <- metrics_df %>%
  dplyr::filter(model == "Cubist") %>%
  dplyr::group_by(property) %>%
  dplyr::slice_max(RPIQ, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

message("\n  Best Cubist model per property (by RPIQ):")
print(best_summary[, c("property", "spectra", "RMSE", "CCC", "RPIQ")])

write.csv(best_summary,
          file.path(cfg$dir_tables, "best_model_per_property.csv"),
          row.names = FALSE)

# Best index for TC from the ranking
#  Extract index-only metrics (exclude raw ASD/NIRVascan) --------
metrics_idx_only <- metrics_df %>%
  dplyr::filter(
    model   == "Cubist",
    grepl("^Index", spectra)           # only standardised NIRVascan indices
  ) %>%
  dplyr::mutate(
    index = stringr::str_extract(spectra, "Index \\d+"),
    pp    = dplyr::if_else(grepl("Savitzky", spectra),
                           "Savitzky-Golay", "No pre-processing")
  )

# Select top 2 indices for TC
best_tc_2 <- metrics_idx_only %>%
  dplyr::filter(property == "Total Carbon") %>%
  dplyr::slice_max(RPIQ, n = 2, with_ties = FALSE) %>%
  dplyr::mutate(
    idx_key = paste0("idx", stringr::str_extract(index, "\\d+")),
    pp_key  = if_else(pp == "Savitzky-Golay", "sg", "orig")
  )

message("Top 2 indices for Total Carbon:")
print(best_tc_2[, c("index", "pp", "RMSE", "CCC", "RPIQ", "idx_key", "pp_key")])

# Store individually for easy reference
best_tc_1 <- best_tc_2[1, ]   # rank 1
best_tc_2r <- best_tc_2[2, ]  # rank 2

message(sprintf("  Rank 1: %s | %s | RPIQ = %.3f",
                best_tc_1$index, best_tc_1$pp, best_tc_1$RPIQ))
message(sprintf("  Rank 2: %s | %s | RPIQ = %.3f",
                best_tc_2r$index, best_tc_2r$pp, best_tc_2r$RPIQ))


# ---- 9.3  Metrics heatmap (all indices, Cubist only) ----------------------
# Compute per-metric range for independent colour scaling
metrics_idx <- metrics_df %>%
  dplyr::filter(model == "Cubist") %>%
  tidyr::pivot_longer(cols = c(RMSE, MEC, CCC, RPIQ),
                      names_to = "metric", values_to = "value") %>%
  dplyr::mutate(
    metric  = factor(metric, levels = c("RPIQ", "CCC", "MEC", "RMSE")),
    spectra = stringr::str_wrap(spectra, width = 20),
    # Flip RMSE direction so that darker = better for ALL metrics
    value_plot = if_else(metric == "RMSE", -value, value)
  )

# Build one plot per metric then combine — this is the only reliable way
# to get truly independent colour scales in ggplot2
make_metric_panel <- function(metric_name, df, low_col, high_col) {
  df_m <- df %>% filter(metric == metric_name)
  
  ggplot(df_m, aes(x = spectra, y = property, fill = value_plot)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", value)),
              size = 2.2, colour = "black") +
    scale_fill_gradient(
      low  = low_col,
      high = high_col,
      name = if (metric_name == "RMSE") "−RMSE\n(darker=lower)" else metric_name,
      guide = guide_colourbar(barwidth = 0.6, barheight = 4)
    ) +
    scale_x_discrete(position = "bottom") +
    labs(x = NULL, y = NULL, title = metric_name) +
    theme_pub(base_size = 8) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 6.5),
      axis.text.y     = element_text(size = 7),
      strip.text      = element_text(face = "bold", size = 9),
      plot.title      = element_text(face = "bold", size = 9, hjust = 0.5),
      legend.position = "right",
      legend.title    = element_text(size = 7),
      legend.text     = element_text(size = 6)
    )
}

# Colorblind-safe: Okabe-Ito sequential pairs
# RPIQ/CCC/MEC: white → teal  (higher = better = darker)
# RMSE:         white → orange (lower = better → we plot -RMSE so darker = better)
p_rpiq <- make_metric_panel("RPIQ", metrics_idx, "#FFFFFF", "#009E73")
p_ccc  <- make_metric_panel("CCC",  metrics_idx, "#FFFFFF", "#0072B2")
p_mec  <- make_metric_panel("MEC",  metrics_idx, "#FFFFFF", "#56B4E9")
p_rmse <- make_metric_panel("RMSE", metrics_idx, "#FFFFFF", "#E69F00")

fig_heat <- (p_rpiq + p_ccc) / (p_mec + p_rmse) +
  patchwork::plot_annotation(
    title    = "Validation metrics \u2014 all standardisation indices (Cubist)",
    subtitle = "Darker = better for all panels | RMSE panel shows \u2212RMSE so colour direction is consistent",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(colour = "grey40", size = 8)
    )
  )

ggsave(file.path(cfg$dir_figs, "09_metrics_heatmap_indices.jpeg"),
       fig_heat, width = 16, height = 10, dpi = 300, units = "in")
ggsave(file.path(cfg$dir_figs, "09_metrics_heatmap_indices.pdf"),
       fig_heat, width = 16, height = 10, dpi = 300, units = "in")

# ---- 9.4  Index ranking dot plot per property ------------------------------
metrics_rank <- metrics_df %>%
  dplyr::filter(model == "Cubist",
                grepl("^Index", spectra)) %>%
  dplyr::mutate(
    index = stringr::str_extract(spectra, "Index \\d+"),
    pp    = dplyr::if_else(grepl("Savitzky", spectra),
                           "Savitzky-Golay", "No pre-processing")
  )

p_rank <- ggplot(metrics_rank,
                 aes(x = RPIQ, y = reorder(index, RPIQ),
                     colour = pp, shape = pp)) +
  geom_point(size = 3, alpha = 0.9) +
  facet_wrap(~property, scales = "free_x", ncol = 2) +
  scale_colour_manual(values = c("No pre-processing" = "#0072B2",
                                 "Savitzky-Golay"    = "#E69F00"),
                      name = "Pre-processing") +
  scale_shape_manual(values  = c("No pre-processing" = 16,
                                 "Savitzky-Golay"    = 17),
                     name = "Pre-processing") +
  labs(x = "RPIQ (validation)", y = NULL,
       title    = "Standardisation index ranking by soil property (Cubist)",
       subtitle = "Higher RPIQ = better generalisation to NIRVascan validation set") +
  theme_pub(base_size = 10) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

ggsave(file.path(cfg$dir_figs, "10_index_ranking_rpiq.pdf"),
       p_rank, width = 10, height = 8, dpi = 300)
ggsave(file.path(cfg$dir_figs, "10_index_ranking_rpiq.jpeg"),
       p_rank, width = 10, height = 8, dpi = 300)

message("  All figures saved to: ", cfg$dir_figs)

# =============================================================================
# ADD-ON: PLSR versions of figure 09 (metrics heatmap) and figure 10 (index
#         ranking dot plot).
#
# Paste this block immediately AFTER Section 9.4 (i.e. after the line
# `message("  All figures saved to: ", cfg$dir_figs)` around line 1505).
#
# It is purely additive: it does NOT modify the existing Cubist code. It reuses
# the globally-defined make_metric_panel() helper, so the look of the PLSR
# panels is identical to the Cubist ones (Okabe-Ito palette, -RMSE flip, layout).
#
# Only change vs. the Cubist blocks: dplyr::filter(model == "PLSR").
# =============================================================================
packageVersion("ggplot2")   # expect 3.5.x
packageVersion("patchwork") # need >= 1.2.0 (1.3.x is current)
packageVersion("gtable")    # 0.3.6

# Guard: only proceed if PLSR rows exist in metrics_df
if (any(metrics_df$model == "PLSR")) {
  
# ---- 9.3-PLSR  Metrics heatmap (all spectra, PLSR only) -------------------
  metrics_idx_plsr <- metrics_df %>%
    dplyr::filter(model == "PLSR") %>%
    tidyr::pivot_longer(cols = c(RMSE, MEC, CCC, RPIQ),
                        names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric  = factor(metric, levels = c("RPIQ", "CCC", "MEC", "RMSE")),
      spectra = stringr::str_wrap(spectra, width = 20),
      # Flip RMSE direction so that darker = better for ALL metrics
      value_plot = dplyr::if_else(metric == "RMSE", -value, value)
    )
  
  # Reuse the existing make_metric_panel() helper (same colours as Cubist fig)
  p_rpiq_pls <- make_metric_panel("RPIQ", metrics_idx_plsr, "#FFFFFF", "#009E73")
  p_ccc_pls  <- make_metric_panel("CCC",  metrics_idx_plsr, "#FFFFFF", "#0072B2")
  p_mec_pls  <- make_metric_panel("MEC",  metrics_idx_plsr, "#FFFFFF", "#56B4E9")
  p_rmse_pls <- make_metric_panel("RMSE", metrics_idx_plsr, "#FFFFFF", "#E69F00")
  
  fig_heat_plsr <- (p_rpiq_pls + p_ccc_pls) / (p_mec_pls + p_rmse_pls) +
    patchwork::plot_annotation(
      title    = "Validation metrics \u2014 all standardisation indices (PLSR)",
      subtitle = "Darker = better for all panels | RMSE panel shows \u2212RMSE so colour direction is consistent",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = "grey40", size = 8)
      )
    )
  
  ggsave(file.path(cfg$dir_figs, "09_metrics_heatmap_indices_plsr.jpeg"),
         fig_heat_plsr, width = 16, height = 10, dpi = 300, units = "in")
  ggsave(file.path(cfg$dir_figs, "09_metrics_heatmap_indices_plsr.pdf"),
         fig_heat_plsr, width = 16, height = 10, dpi = 300, units = "in")
  
# ---- 9.4-PLSR  Index ranking dot plot per property (PLSR only) ------------
  metrics_rank_plsr <- metrics_df %>%
    dplyr::filter(model == "PLSR",
                  grepl("^Index", spectra)) %>%
    dplyr::mutate(
      index = stringr::str_extract(spectra, "Index \\d+"),
      pp    = dplyr::if_else(grepl("Savitzky", spectra),
                             "Savitzky-Golay", "No pre-processing")
    )
  
  p_rank_plsr <- ggplot(metrics_rank_plsr,
                        aes(x = RPIQ, y = reorder(index, RPIQ),
                            colour = pp, shape = pp)) +
    geom_point(size = 3, alpha = 0.9) +
    facet_wrap(~property, scales = "free_x", ncol = 2) +
    scale_colour_manual(values = c("No pre-processing" = "#0072B2",
                                   "Savitzky-Golay"    = "#E69F00"),
                        name = "Pre-processing") +
    scale_shape_manual(values  = c("No pre-processing" = 16,
                                   "Savitzky-Golay"    = 17),
                       name = "Pre-processing") +
    labs(x = "RPIQ (validation)", y = NULL,
         title    = "Standardisation index ranking by soil property (PLSR)",
         subtitle = "Higher RPIQ = better generalisation to NIRVascan validation set") +
    theme_pub(base_size = 10) +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
  
  ggsave(file.path(cfg$dir_figs, "10_index_ranking_rpiq_plsr.pdf"),
         p_rank_plsr, width = 10, height = 8, dpi = 300)
  ggsave(file.path(cfg$dir_figs, "10_index_ranking_rpiq_plsr.jpeg"),
         p_rank_plsr, width = 10, height = 8, dpi = 300)
  
  message("  PLSR figures saved: 09_metrics_heatmap_indices_plsr + 10_index_ranking_rpiq_plsr")
  
} else {
  message("  No PLSR rows found in metrics_df \u2014 skipping PLSR figures.")
}

# ---- 9.5 Normalise RPIQ per property (0-1 scale) ----------------------
# This makes TC, TN, pH, P comparable despite different units
metrics_norm <- metrics_idx_only %>%
  dplyr::group_by(property) %>%
  dplyr::mutate(
    RPIQ_norm = (RPIQ - min(RPIQ, na.rm = TRUE)) /
      (max(RPIQ, na.rm = TRUE) - min(RPIQ, na.rm = TRUE)),
    RMSE_norm = 1 - (RMSE - min(RMSE, na.rm = TRUE)) /
      (max(RMSE, na.rm = TRUE) - min(RMSE, na.rm = TRUE)),
    CCC_norm  = (CCC  - min(CCC,  na.rm = TRUE)) /
      (max(CCC,  na.rm = TRUE)  - min(CCC,  na.rm = TRUE))
  ) %>%
  dplyr::ungroup()

# ---- Aggregate score per index + pp --------------------------------
# Equal weight: RPIQ (40%) + CCC (40%) + RMSE (20%)
index_ranking <- metrics_norm %>%
  dplyr::group_by(index, pp) %>%
  dplyr::summarise(
    score_mean  = round(mean(0.4 * RPIQ_norm +
                               0.4 * CCC_norm  +
                               0.2 * RMSE_norm, na.rm = TRUE), 3),
    RPIQ_mean   = round(mean(RPIQ, na.rm = TRUE), 3),
    CCC_mean    = round(mean(CCC,  na.rm = TRUE), 3),
    RMSE_mean   = round(mean(RMSE, na.rm = TRUE), 3),
    # Count how many properties this index is best for
    n_best      = sum(RPIQ == max(RPIQ)),   # rough check
    .groups     = "drop"
  ) %>%
  dplyr::arrange(desc(score_mean)) %>%
  dplyr::mutate(rank = dplyr::row_number())

message("\n  Index ranking (aggregated across all properties):")
print(index_ranking)

# ---- Best overall index --------------------------------------------
best_index     <- index_ranking$index[1]
best_pp        <- index_ranking$pp[1]
best_score     <- index_ranking$score_mean[1]

message(sprintf("\n  >>> BEST OVERALL: %s | %s | score = %.3f <<<",
                best_index, best_pp, best_score))

write.csv(index_ranking,
          file.path(cfg$dir_tables, "index_overall_ranking.csv"),
          row.names = FALSE)

# ---- Visual ranking — dot plot -------------------------------------
p_rank_overall <- ggplot(index_ranking,
                         aes(x = score_mean,
                             y = reorder(paste(index, pp, sep = "\n"), score_mean),
                             colour = pp, shape = pp)) +
  geom_vline(xintercept = max(index_ranking$score_mean),
             colour = "#009E73", linetype = "dashed",
             linewidth = 0.4, alpha = 0.7) +
  geom_segment(aes(xend = 0,
                   yend = reorder(paste(index, pp, sep = "\n"), score_mean)),
               colour = "grey80", linewidth = 0.4) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.3f", score_mean)),
            hjust = -0.35, size = 3, colour = "grey30") +
  scale_colour_manual(
    values = c("No pre-processing" = "#0072B2",
               "Savitzky-Golay"    = "#E69F00"),
    name   = "Pre-processing") +
  scale_shape_manual(
    values = c("No pre-processing" = 16,
               "Savitzky-Golay"    = 17),
    name   = "Pre-processing") +
  scale_x_continuous(limits = c(0, 1.12),
                     expand = expansion(mult = 0.01)) +
  labs(
    x        = "Aggregate score (0\u20131)",
    y        = NULL,
    title    = "Overall standardisation index ranking",
    subtitle = paste0("Score = 0.4 \u00d7 RPIQ\u2099\u2092\u02b3\u1d39 + ",
                      "0.4 \u00d7 CCC\u2099\u2092\u02b3\u1d39 + ",
                      "0.2 \u00d7 (1\u2212RMSE\u2099\u2092\u02b3\u1d39) ",
                      "| averaged across TC, TN, pH, P")
  ) +
  theme_pub(11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    legend.position    = "bottom"
  )

# ---- Find best cell per property -------------------------------------------
best_per_property <- metrics_idx_only %>%
  dplyr::mutate(label = paste(index, pp, sep = "\n")) %>%
  dplyr::group_by(property) %>%
  dplyr::slice_max(RPIQ, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# ---- Rebuild heatmap with red rectangle ------------------------------------
p_heat_rank <- ggplot(
  metrics_idx_only %>%
    dplyr::mutate(label = paste(index, pp, sep = "\n")),
  aes(x = reorder(label, -RPIQ), y = property, fill = RPIQ)) +
  
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", RPIQ)),
            size = 2.8, colour = "black") +
  
  # Red rectangle on best cell per property
  geom_tile(data = best_per_property %>%
              dplyr::mutate(label = paste(index, pp, sep = "\n")),
            aes(x = label, y = property),
            fill        = NA,
            colour      = "#C0392B",
            linewidth   = 1.2,
            inherit.aes = FALSE) +
  
  scale_fill_gradient(low  = "#F7FCF5",
                      high = "#00441B",
                      name = "RPIQ") +
  labs(x = NULL, y = NULL,
       title    = "RPIQ by index and soil property",
       subtitle = "Darker green = better | red box = best index per property") +
  theme_pub(10) +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right")

# ---- Combine and save ----------------------------------------------
fig_ranking <- p_rank_overall / p_heat_rank +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    theme = theme(plot.margin = margin(5, 5, 5, 5))
  )

ggsave(file.path(cfg$dir_figs, "13_best_index_overall_ranking.jpeg"),
       fig_ranking, width = 12, height = 12, dpi = 300, units = "in")

message(sprintf("  Ranking figure saved. Use %s (%s) for the OSSL transfer step.",
                best_index, best_pp))
# =============================================================================
# SECTION 10 - DESCRIPTIVE STATISTICS TABLE
# =============================================================================

message("\n[10] Descriptive statistics...")
colnames(soil)
desc_stats <- soil %>%
  dplyr::select(country, landuse,
                tc.perc, tn.perc, ph.water, p.dl.mg_100g_ok) %>%
  dplyr::mutate(
    tc.perc      = as.numeric(tc.perc),
    tn.perc      = as.numeric(tn.perc),
    ph.water     = as.numeric(ph.water),
    p.dl.mg_100g_ok = as.numeric(p.dl.mg_100g_ok)
  ) %>%
  tidyr::pivot_longer(cols = -c(country, landuse),
                      names_to = "property", values_to = "value") %>%
  dplyr::group_by(country, landuse, property) %>%
  dplyr::summarise(
    n      = dplyr::n(),
    min    = round(min(value,               na.rm = TRUE), 3),
    q1     = round(quantile(value, 0.25,    na.rm = TRUE), 3),
    median = round(median(value,            na.rm = TRUE), 3),
    mean   = round(mean(value,              na.rm = TRUE), 3),
    q3     = round(quantile(value, 0.75,    na.rm = TRUE), 3),
    max    = round(max(value,               na.rm = TRUE), 3),
    sd     = round(sd(value,                na.rm = TRUE), 3),
    cv_pct = round(sd(value, na.rm = TRUE) /
                     mean(value, na.rm = TRUE) * 100, 1),
    skew   = round(e1071::skewness(value,   na.rm = TRUE), 3),
    kurt   = round(e1071::kurtosis(value,   na.rm = TRUE), 3),
    .groups = "drop"
  )

write.csv(desc_stats,
          file.path(cfg$dir_tables, "descriptive_statistics.csv"),
          row.names = FALSE)
message("  Descriptive statistics saved.")

if (!requireNamespace("ggdist", quietly = TRUE)) install.packages("ggdist")
library(ggdist)
if (!requireNamespace("ggridges", quietly = TRUE)) install.packages("ggridges")
if (!requireNamespace("cowplot", quietly = TRUE))  install.packages("cowplot")
library(ggridges)
library(cowplot)

# ---- n labels per group per property ---------------------------------------
n_labels <- soil_long %>%
  dplyr::group_by(property, group) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0("n = ", n))

# ---- Density panel function ------------------------------------------------
make_density_panel <- function(prop_name, unit, tag = NULL) {
  
  df   <- soil_long %>% dplyr::filter(property == prop_name)
  df_n <- n_labels  %>% dplyr::filter(property == prop_name)
  
  ggplot(df, aes(x = value, y = group, fill = group)) +
    
    ggridges::geom_density_ridges(
      alpha     = 0.65,
      scale     = 1.6,
      linewidth = 0.4,
      colour    = "white",
      bandwidth = NULL
    ) +
    
    # n = label inside each ridge at the far right
    geom_text(
      data        = df_n,
      aes(x = Inf, y = group, label = label),
      hjust       = 1.1,
      vjust       = -0.8,
      size        = 3,
      colour      = "grey30",
      fontface    = "italic",
      inherit.aes = FALSE
    ) +
    
    scale_fill_manual(
      values = c("USA\n(grasslands)" = "#56B4E9",
                 "GE\n(peatlands)"  = "#E69F00"),
      name   = "Country (land use)",
      labels = c("USA (grasslands)", "GE (peatlands)")
    ) +
    
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    
    labs(x = unit, y = NULL,
         title = prop_name,
         tag   = tag) +
    
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
      panel.border      = element_rect(colour = "grey40"),
      axis.text         = element_text(colour = "black", size = 9),
      axis.text.y       = element_text(face = "bold", size = 9),
      axis.title.x      = element_text(colour = "black", size = 9),
      plot.title        = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.tag          = element_text(face = "bold", size = 12),
      legend.position   = "none"
    )
}

# ---- Build four panels -----------------------------------------------------
d1 <- make_density_panel("Total Carbon (%)",             "%",            tag = "(a)")
d2 <- make_density_panel("Total Nitrogen (%)",            "%",            tag = "(b)")
d3 <- make_density_panel("pH in water (adimensional)",    "adimensional", tag = "(c)")
d4 <- make_density_panel("Plant-available P (mg/100g)",   "mg/100g",      tag = "(d)")

# ---- Shared legend ---------------------------------------------------------
legend_grob <- cowplot::get_legend(
  make_density_panel("Total Carbon (%)", "%") +
    theme(legend.position  = "bottom",
          legend.key.size  = unit(0.5, "cm"),
          legend.text      = element_text(size = 9),
          legend.title     = element_text(size = 9, face = "bold"))
)

# ---- Combine ---------------------------------------------------------------
fig1 <- cowplot::plot_grid(
  cowplot::plot_grid(d1, d2, d3, d4, ncol = 2,
                     align = "hv", axis = "tblr"),
  legend_grob,
  ncol        = 1,
  rel_heights = c(10, 0.7)
)

ggsave(file.path(cfg$dir_figs, "Fig1_soil_distribution.jpeg"),
       fig1, width = 10, height = 10, dpi = 300, units = "in")


# =============================================================================
# OSSL CASE STUDY — Total Carbon prediction (OPTIMISED)
#
# Speed fixes applied:
#   1. Workers capped at 4 (ossl_cfg$n_workers)
#   2. CV: folds=5, repeats=1 for benchmark; increase for final run
#   3. Benchmark mode: 3000 OSSL samples (ossl_cfg$benchmark = TRUE)
#   4. PCA compression: 20-30 scores replace 2051 collinear bands
#   5. Single numeric matrix per dataset — no repeated apply(rbind())
#   6. Model C fitted ONCE — NIRVascan sets projected into same PCA, no refit
#
# Models:
#   A: OSSL only           (vis-NIR, PCA)  → val: ASD full range
#   B: OSSL + local ASD    (vis-NIR, PCA)  → val: ASD full range
#   C: OSSL + local ASD    (NIR, PCA)      → val: ASD NIR + top 2 NIRVascan
# =============================================================================

# =============================================================================
# CONFIGURATION — edit here only
# =============================================================================
ossl_cfg <- list(
  benchmark   = FALSE,   # TRUE = fast test on n_bench samples only
  n_bench     = 3000,    # samples used in benchmark mode
  n_pca_full  = 30,      # PCA components, full vis-NIR
  n_pca_nir   = 20,      # PCA components, NIR range
  cub_folds   = 10,      # publication setting
  cub_repeats = 5,       # publication setting
  n_workers   = 4        # parallel workers — keep at 4-6 on 16 GB RAM
)

message("Configuration:")
message(sprintf("  Benchmark : %s (n=%d)", ossl_cfg$benchmark, ossl_cfg$n_bench))
message(sprintf("  PCA       : %d (full) / %d (NIR) components",
                ossl_cfg$n_pca_full, ossl_cfg$n_pca_nir))
message(sprintf("  CV        : %d folds x %d repeats",
                ossl_cfg$cub_folds, ossl_cfg$cub_repeats))
message(sprintf("  Workers   : %d", ossl_cfg$n_workers))

# Local Cubist fitter for this section — overrides global cfg
fit_cubist_ossl <- function(x, y) {
  ctrl <- caret::trainControl(
    method        = "repeatedcv",
    number        = ossl_cfg$cub_folds,
    repeats       = ossl_cfg$cub_repeats,
    allowParallel = TRUE
  )
  grid <- expand.grid(committees = c(1, 10, 50, 100),
                      neighbors  = c(0, 1, 5, 9))
  cl <- parallel::makeCluster(ossl_cfg$n_workers, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  caret::train(x = as.data.frame(x), y = y,
               method    = "cubist",
               trControl = ctrl,
               tuneGrid  = grid)
}

# ---- Top 2 indices for TC --------------------------------------------------
best_tc_top2 <- metrics_idx_only %>%
  dplyr::filter(property == "Total Carbon") %>%
  dplyr::slice_max(RPIQ, n = 2, with_ties = FALSE) %>%
  dplyr::mutate(
    idx_key = paste0("idx", stringr::str_extract(index, "\\d+")),
    pp_key  = dplyr::if_else(pp == "Savitzky-Golay", "sg", "orig")
  )

message("Top 2 TC indices:")
print(best_tc_top2[, c("index", "pp", "RMSE", "CCC", "RPIQ", "idx_key", "pp_key")])


# =============================================================================
# STEP 1: Load OSSL (use data.table::fread — much faster than read.csv)
# =============================================================================
message("\n[OSSL 1] Loading OSSL...")

ossl_soil <- data.table::fread(
  file.path(cfg$dir_data, "ossl_soillab_L0_v1.2.csv"), data.table = FALSE)
ossl_spec <- data.table::fread(
  file.path(cfg$dir_data, "ossl_visnir_L0_v1.2.csv"),  data.table = FALSE)

colnames(ossl_spec) <- gsub("scan_visnir\\.", "", colnames(ossl_spec))
colnames(ossl_spec) <- gsub("_ref",           "", colnames(ossl_spec))

ossl_merged <- ossl_soil %>%
  dplyr::select(id.layer_uuid_txt, tc.perc = c.tot_usda.a622_w.pct) %>%
  dplyr::filter(!is.na(tc.perc)) %>%
  dplyr::inner_join(ossl_spec, by = "id.layer_uuid_txt") %>%
  dplyr::filter(complete.cases(.))

if (ossl_cfg$benchmark) {
  set.seed(cfg$random_seed)
  ossl_merged <- dplyr::slice_sample(ossl_merged, n = ossl_cfg$n_bench)
  message(sprintf("  BENCHMARK: %d OSSL samples", nrow(ossl_merged)))
} else {
  message(sprintf("  OSSL full: %d samples", nrow(ossl_merged)))
}

spec_cols <- colnames(ossl_merged)[
  suppressWarnings(!is.na(as.numeric(colnames(ossl_merged))))
]
wl_ossl <- as.numeric(spec_cols)
y_ossl  <- ossl_merged$tc.perc

message(sprintf("  Bands: %.0f-%.0f nm (%d total)",
                min(wl_ossl), max(wl_ossl), length(wl_ossl)))


# =============================================================================
# STEP 2: Resample OSSL — ONE raw matrix extraction, two resamples
# =============================================================================
ossl_raw <- as.matrix(ossl_merged[, spec_cols])   # extract once

wl_asd_full <- 400:2450

ossl_full_mat <- prospectr::resample(
  ossl_raw[, wl_ossl >= 400 & wl_ossl <= 2450],
  wl_ossl[wl_ossl >= 400 & wl_ossl <= 2450],
  wl_asd_full, interpol = "spline")
colnames(ossl_full_mat) <- as.character(wl_asd_full)

ossl_nir_mat <- prospectr::resample(
  ossl_raw[, wl_ossl >= 950 & wl_ossl <= 1650],
  wl_ossl[wl_ossl >= 950 & wl_ossl <= 1650],
  wl_target, interpol = "spline")
colnames(ossl_nir_mat) <- as.character(wl_target)

rm(ossl_raw); gc()
message(sprintf("  Full: %dx%d | NIR: %dx%d",
                nrow(ossl_full_mat), ncol(ossl_full_mat),
                nrow(ossl_nir_mat),  ncol(ossl_nir_mat)))


# =============================================================================
# STEP 3: Local ASD datasets — extract ONCE, store as numeric matrix
# =============================================================================
# ---- 3a: Full vis-NIR (400-2450 nm) ----------------------------------------
wl_cols_asd <- colnames(asd_pp$orig)[
  suppressWarnings(!is.na(as.numeric(gsub("X", "", colnames(asd_pp$orig)))))]
wl_asd_num  <- as.numeric(gsub("X", "", wl_cols_asd))
wl_keep     <- wl_cols_asd[wl_asd_num >= 400 & wl_asd_num <= 2450]
wl_keep_nm  <- as.character(wl_asd_num[wl_asd_num >= 400 & wl_asd_num <= 2450])

asd_full_df <- asd_pp$orig %>%
  dplyr::filter(!is.na(id), id != lb_id) %>%
  dplyr::select(id, dplyr::all_of(wl_keep))
colnames(asd_full_df)[-1] <- wl_keep_nm

mk_num_mat <- function(df, cols) {
  matrix(as.numeric(as.matrix(df[, cols])),
         nrow = nrow(df), dimnames = list(NULL, cols))
}

local_full_cal_df <- merge(soil[, c("id","tc.perc")],
                           asd_full_df[asd_full_df$id %in% id_cal, ],
                           by = "id") %>% dplyr::filter(!is.na(tc.perc))
local_full_val_df <- merge(soil[, c("id","tc.perc")],
                           asd_full_df[asd_full_df$id %in% id_val, ],
                           by = "id") %>% dplyr::filter(!is.na(tc.perc))
y_local_full_cal  <- local_full_cal_df$tc.perc
y_local_full_val  <- local_full_val_df$tc.perc

shared_full <- intersect(colnames(ossl_full_mat),
                         colnames(local_full_cal_df)[-(1:2)])
x_ossl_full      <- ossl_full_mat[, shared_full]
x_local_full_cal <- mk_num_mat(local_full_cal_df, shared_full)
x_local_full_val <- mk_num_mat(local_full_val_df, shared_full)

message(sprintf("  Full — cal: %d | val: %d | bands: %d",
                nrow(x_local_full_cal), nrow(x_local_full_val),
                length(shared_full)))

# ---- 3b: NIR (950-1650 nm, 5 nm) -------------------------------------------
local_nir_cal_df <- merge(soil[, c("id","tc.perc")],
                          spec$asd_orig$cal, by = "id") %>%
  dplyr::filter(!is.na(tc.perc))
local_nir_val_df <- merge(soil[, c("id","tc.perc")],
                          spec$asd_orig$val, by = "id") %>%
  dplyr::filter(!is.na(tc.perc))
y_local_nir_cal  <- local_nir_cal_df$tc.perc
y_local_nir_val  <- local_nir_val_df$tc.perc

shared_nir <- intersect(colnames(ossl_nir_mat),
                        colnames(local_nir_cal_df)[-(1:2)])
x_ossl_nir      <- ossl_nir_mat[, shared_nir]
x_local_nir_cal <- mk_num_mat(local_nir_cal_df, shared_nir)
x_local_nir_val <- mk_num_mat(local_nir_val_df, shared_nir)

message(sprintf("  NIR  — cal: %d | val: %d | bands: %d",
                nrow(x_local_nir_cal), nrow(x_local_nir_val),
                length(shared_nir)))


# =============================================================================
# STEP 4: PCA compression — biggest single speed improvement
#         Replaces 2051 collinear bands with 20-30 orthogonal scores
# =============================================================================
pca_var_kept <- function(pca, n) {
  round(cumsum(pca$sdev^2)[n] / sum(pca$sdev^2) * 100, 1)
}

# ---- Full PCA (fit on combined training data) ------------------------------
x_AB_combined <- rbind(x_ossl_full, x_local_full_cal)
y_AB_combined <- c(y_ossl, y_local_full_cal)

pca_full  <- prcomp(x_AB_combined, center = TRUE, scale. = FALSE)
n_pc_full <- min(ossl_cfg$n_pca_full,
                 which(cumsum(pca_full$sdev^2)/sum(pca_full$sdev^2) >= 0.999)[1])

sc_ossl_full_A  <- predict(pca_full, x_ossl_full)[,      1:n_pc_full]
sc_local_full_cal <- predict(pca_full, x_local_full_cal)[, 1:n_pc_full]
sc_local_full_val <- predict(pca_full, x_local_full_val)[, 1:n_pc_full]
sc_AB_cal         <- rbind(sc_ossl_full_A, sc_local_full_cal)

message(sprintf("  Full PCA: %d PCs = %.1f%% variance",
                n_pc_full, pca_var_kept(pca_full, n_pc_full)))

rm(x_ossl_full, x_local_full_cal, x_AB_combined); gc()

# ---- NIR PCA ---------------------------------------------------------------
x_C_combined <- rbind(x_ossl_nir, x_local_nir_cal)
y_C_combined <- c(y_ossl, y_local_nir_cal)

pca_nir  <- prcomp(x_C_combined, center = TRUE, scale. = FALSE)
n_pc_nir <- min(ossl_cfg$n_pca_nir,
                which(cumsum(pca_nir$sdev^2)/sum(pca_nir$sdev^2) >= 0.999)[1])

sc_ossl_nir_C   <- predict(pca_nir, x_ossl_nir)[,      1:n_pc_nir]
sc_local_nir_cal <- predict(pca_nir, x_local_nir_cal)[, 1:n_pc_nir]
sc_local_nir_val <- predict(pca_nir, x_local_nir_val)[, 1:n_pc_nir]
sc_C_cal         <- rbind(sc_ossl_nir_C, sc_local_nir_cal)
y_C_cal          <- y_C_combined

message(sprintf("  NIR  PCA: %d PCs = %.1f%% variance",
                n_pc_nir, pca_var_kept(pca_nir, n_pc_nir)))
message(sprintf("  Training: A=%d x %d | B=%d x %d | C=%d x %d",
                nrow(sc_ossl_full_A), ncol(sc_ossl_full_A),
                nrow(sc_AB_cal),      ncol(sc_AB_cal),
                nrow(sc_C_cal),       ncol(sc_C_cal)))

rm(x_ossl_nir, x_local_nir_cal, x_C_combined); gc()


# =============================================================================
# STEP 5: MODEL A — OSSL only (full vis-NIR, PCA) → ASD full val
# =============================================================================
t_A <- system.time({
  cub_A  <- fit_cubist_ossl(sc_ossl_full_A, y_ossl)
  pred_A <- predict(cub_A, as.data.frame(sc_local_full_val))
  met_A  <- goof(observed = y_local_full_val, predicted = pred_A)
})
message(sprintf("  RMSE=%.3f CCC=%.3f RPIQ=%.3f [%.1f min]",
                met_A$RMSE, met_A$CCC, met_A$RPIQ, t_A["elapsed"]/60))


# =============================================================================
# STEP 6: MODEL B — OSSL + local ASD (full vis-NIR, PCA) → ASD full val
# =============================================================================
t_B <- system.time({
  cub_B  <- fit_cubist_ossl(sc_AB_cal, y_AB_combined)
  pred_B <- predict(cub_B, as.data.frame(sc_local_full_val))
  met_B  <- goof(observed = y_local_full_val, predicted = pred_B)
})
message(sprintf("  RMSE=%.3f CCC=%.3f RPIQ=%.3f [%.1f min]",
                met_B$RMSE, met_B$CCC, met_B$RPIQ, t_B["elapsed"]/60))

rm(sc_ossl_full_A, sc_AB_cal); gc()


# =============================================================================
# STEP 7: MODEL C — OSSL + local ASD (NIR, PCA) — fit ONCE
#         Validate on multiple sets by projecting into pca_nir — no refit
# =============================================================================
t_C <- system.time({
  cub_C <- fit_cubist_ossl(sc_C_cal, y_C_cal)
})
message(sprintf("  Fitted [%.1f min]", t_C["elapsed"]/60))

# ---- C-i: ASD NIR val ------------------------------------------------------
pred_C_asd <- predict(cub_C, as.data.frame(sc_local_nir_val))
met_C_asd  <- goof(observed = y_local_nir_val, predicted = pred_C_asd)
message(sprintf("  -> ASD NIR val  RMSE=%.3f CCC=%.3f RPIQ=%.3f",
                met_C_asd$RMSE, met_C_asd$CCC, met_C_asd$RPIQ))

# ---- C-ii: top 2 transferred NIRVascan — project into pca_nir, no refit ---
nirv_results <- purrr::map(seq_len(nrow(best_tc_top2)), function(i) {
  
  idx_key <- best_tc_top2$idx_key[i]
  pp_key  <- best_tc_top2$pp_key[i]
  label   <- sprintf("%s (%s)", best_tc_top2$index[i], best_tc_top2$pp[i])
  
  trans_val <- transferred[[paste0("val_", pp_key)]][[idx_key]]
  
  val_tc <- merge(soil[, c("id","tc.perc")], trans_val, by = "id") %>%
    dplyr::filter(!is.na(tc.perc))
  
  x_val_aligned <- mk_num_mat(val_tc, shared_nir)
  sc_nirv        <- predict(pca_nir, x_val_aligned)[, 1:n_pc_nir]
  
  pred <- predict(cub_C, as.data.frame(sc_nirv))  # no refit
  met  <- goof(observed = val_tc$tc.perc, predicted = pred)
  
  message(sprintf("  -> %s  RMSE=%.3f CCC=%.3f RPIQ=%.3f",
                  label, met$RMSE, met$CCC, met$RPIQ))
  
  list(label = label, idx_key = idx_key, pp_key = pp_key,
       pred = pred, obs = val_tc$tc.perc, met = met, n_val = nrow(val_tc))
})


# =============================================================================
# STEP 8: Compile metrics
# =============================================================================
metrics_all <- dplyr::bind_rows(
  data.frame(model   = "OSSL only",
             val_set = "ASD full range",
             n_cal   = length(y_ossl),
             n_val   = length(y_local_full_val), met_A),
  data.frame(model   = "OSSL + local ASD",
             val_set = "ASD full range",
             n_cal   = length(y_AB_combined),
             n_val   = length(y_local_full_val), met_B),
  data.frame(model   = "OSSL + local ASD",
             val_set = "ASD NIR range",
             n_cal   = length(y_C_cal),
             n_val   = length(y_local_nir_val), met_C_asd),
  purrr::map_dfr(nirv_results, function(r)
    data.frame(model   = "OSSL + local ASD",
               val_set = sprintf("NIRVascan into ASD | %s", r$label),
               n_cal   = length(y_C_cal),
               n_val   = r$n_val, r$met))
) %>%
  dplyr::select(model, val_set, n_cal, n_val, RMSE, MEC, CCC, RPIQ, bias) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

print(metrics_all)
write.csv(metrics_all,
          file.path(cfg$dir_tables, "15c_ossl_tc_all_metrics.csv"),
          row.names = FALSE)


# =============================================================================
# STEP 9: Figures
# =============================================================================
# ---- 15a: OSSL NIR spectral overview ----------------------------------------
ossl_spec_summ <- as.data.frame(ossl_nir_mat) %>%
  dplyr::slice_sample(n = min(500, nrow(ossl_nir_mat))) %>%
  dplyr::mutate(rid = dplyr::row_number()) %>%
  tidyr::pivot_longer(-rid, names_to = "wl", values_to = "ref") %>%
  dplyr::mutate(wl = as.numeric(wl)) %>%
  dplyr::filter(!is.na(ref)) %>%
  dplyr::group_by(wl) %>%
  dplyr::summarise(mn = mean(ref, na.rm = TRUE),
                   lo = min(ref,  na.rm = TRUE),
                   hi = max(ref,  na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(is.finite(mn))

p15a <- ggplot(ossl_spec_summ, aes(x = wl)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey60",
              alpha = 0.2, na.rm = TRUE) +
  geom_line(aes(y = mn), colour = "#0072B2",
            linewidth = 0.6, na.rm = TRUE) +
  scale_x_continuous(breaks = seq(950, 1650, 100),
                     expand = expansion(mult = 0.01)) +
  labs(x = "Wavelength (nm)", y = "Reflectance Factor",
       title    = "OSSL spectra resampled to NIR range (950\u20131650 nm, 5 nm)",
       subtitle = sprintf("Mean \u00b1 full range | n = %d samples",
                          nrow(ossl_nir_mat))) +
  theme_pub(11)

ggsave(file.path(cfg$dir_figs, "15a_ossl_spectra_overview.jpeg"),
       p15a, width = 8, height = 4, dpi = 300)

# ---- 15b: Models A and B obs vs pred ----------------------------------------
p_A <- obs_pred_plot(pred_A, y_local_full_val,
                     title = sprintf("A: OSSL only (PCA)\ncal = %d | RPIQ = %.2f",
                                     length(y_ossl), met_A$RPIQ),
                     unit = "%", axis_lim = c(0, 60))
p_B <- obs_pred_plot(pred_B, y_local_full_val,
                     title = sprintf("B: OSSL + local ASD (PCA)\ncal = %d | RPIQ = %.2f",
                                     length(y_AB_combined), met_B$RPIQ),
                     unit = "%", axis_lim = c(0, 60))

fig_15b <- p_A + p_B +
  patchwork::plot_annotation(
    title    = "Total Carbon \u2014 Models A and B | full vis-NIR PCA | ASD validation",
    subtitle = sprintf("Validation: local ASD holdout set | val. = %d",
                       length(y_local_full_val)),
    theme    = theme(plot.title    = element_text(face = "bold", size = 11),
                     plot.subtitle = element_text(colour = "grey40", size = 9)))

ggsave(file.path(cfg$dir_figs, "15b_ossl_models_AB.jpeg"),
       fig_15b, width = 10, height = 5, dpi = 300)

# ---- 15c: Model C all validation sets ---------------------------------------
p_C_asd <- obs_pred_plot(pred_C_asd, y_local_nir_val,
                         title    = "ASD NIR val",
                         unit     = "%",
                         axis_lim = c(0, 60))

panels_nirv <- purrr::map(nirv_results, function(r)
  obs_pred_plot(r$pred, r$obs,
                title    = r$label,       # label only — no RPIQ, no n_val
                unit     = "%",
                axis_lim = c(0, 60)))

fig_15c <- patchwork::wrap_plots(c(list(p_C_asd), panels_nirv), ncol = 3) +
  patchwork::plot_annotation(
    title    = sprintf("Total Carbon \u2014 Model C | OSSL + local ASD (NIR PCA) | cal = %d",
                       length(y_C_cal)),
    subtitle = "Left: ASD NIR val | Centre/Right: NIRVascan transferred val (same fitted model, no refit)",
    theme    = theme(plot.title    = element_text(face = "bold", size = 11),
                     plot.subtitle = element_text(colour = "grey40", size = 9)))

ggsave(file.path(cfg$dir_figs, "15c_ossl_model_C.jpeg"),
       fig_15c, width = 14, height = 5, dpi = 300)

# ---- 15d: All models RPIQ/CCC/RMSE comparison bar chart --------------------
fig_15d <- metrics_all %>%
  dplyr::mutate(run = stringr::str_wrap(paste0(model, ": ", val_set), 35)) %>%
  tidyr::pivot_longer(cols = c(RPIQ, CCC, RMSE),
                      names_to = "metric", values_to = "value") %>%
  dplyr::mutate(metric = factor(metric, levels = c("RPIQ", "CCC", "RMSE"))) %>%
  ggplot(aes(x = reorder(run, value), y = value, fill = metric)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f", value)),
            hjust = -0.15, size = 2.8, colour = "grey30") +
  facet_wrap(~metric, scales = "free_x", ncol = 3) +
  scale_fill_manual(
    values = c(RPIQ = "#009E73", CCC = "#0072B2", RMSE = "#E69F00"),
    guide  = "none") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = NULL,
       title    = "Total Carbon: all OSSL model configurations",
       subtitle = "RPIQ and CCC: higher = better | RMSE: lower = better") +
  theme_pub(10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        strip.text         = element_text(face = "bold"),
        axis.text.y        = element_text(size = 8))

ggsave(file.path(cfg$dir_figs, "15d_ossl_all_models_comparison.jpeg"),
       fig_15d, width = 14, height = 7, dpi = 300)

# ---- 15e: PCA coverage — OSSL (calibration) vs local ASD val ---------------
# Project both sets into the shared pca_nir space already fitted
scores_ossl_plot <- as.data.frame(sc_ossl_nir_C[, 1:3]) %>%
  dplyr::mutate(set = "OSSL (calibration)")

scores_val_plot <- as.data.frame(sc_local_nir_val[, 1:3]) %>%
  dplyr::mutate(set = "Local ASD (validation)")

pca_df_ossl <- dplyr::bind_rows(scores_ossl_plot, scores_val_plot) %>%
  dplyr::mutate(set = factor(set, levels = c("OSSL (calibration)",
                                             "Local ASD (validation)")))

var_exp_nir <- round(summary(pca_nir)$importance[2, 1:3] * 100, 1)

pca_ossl_panel <- function(df, x, y, var_exp) {
  xi <- as.integer(gsub("PC", "", x))
  yi <- as.integer(gsub("PC", "", y))
  
  df_ossl <- df %>% dplyr::filter(set == "OSSL (calibration)")
  df_val  <- df %>% dplyr::filter(set == "Local ASD (validation)")
  
  ggplot(df, aes(.data[[x]], .data[[y]])) +
    
    # OSSL: small, transparent, background
    geom_point(data    = df_ossl,
               aes(colour = set, shape = set),
               size = 0.9, alpha = 0.2) +
    stat_ellipse(data      = df_ossl,
                 aes(colour = set, group = set),
                 linewidth = 0.5, linetype = "dashed") +
    
    # Local ASD val: larger, solid, foreground
    geom_point(data    = df_val,
               aes(colour = set, shape = set),
               size = 2.2, alpha = 0.9) +
    stat_ellipse(data      = df_val,
                 aes(colour = set, group = set),
                 linewidth = 0.7, linetype = "dashed") +
    
    scale_colour_manual(
      values = c("OSSL (calibration)"     = "#0072B2",
                 "Local ASD (validation)" = "#C0392B"),
      name = NULL) +
    scale_shape_manual(
      values = c("OSSL (calibration)"     = 16,
                 "Local ASD (validation)" = 17),
      name = NULL) +
    labs(x = sprintf("PC%d (%.1f%%)", xi, var_exp[xi]),
         y = sprintf("PC%d (%.1f%%)", yi, var_exp[yi])) +
    theme_pub(11) +
    theme(legend.position  = "bottom",
          legend.key.size  = unit(0.4, "cm"),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))
}

p_pc12 <- pca_ossl_panel(pca_df_ossl, "PC1", "PC2", var_exp_nir)
p_pc13 <- pca_ossl_panel(pca_df_ossl, "PC1", "PC3", var_exp_nir)

fig_15e <- (p_pc12 + p_pc13) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = "OSSL spectral coverage vs. local ASD validation set (NIR 950\u20131650 nm, PCA)",
    subtitle = sprintf(
      "OSSL: n = %d | Local ASD val: n = %d | Ellipse overlap = OSSL covers local spectral range",
      nrow(scores_ossl_plot), nrow(scores_val_plot)),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(colour = "grey40", size = 9))
  ) &
  theme(legend.position = "bottom")

ggsave(file.path(cfg$dir_figs, "15e_ossl_pca_coverage.jpeg"),
       fig_15e, width = 11, height = 5.5, dpi = 300)


# =============================================================================
# STEP 10: Save models (include PCA object for future prediction)
# =============================================================================
saveRDS(list(model=cub_A, pca=pca_full, n_pc=n_pc_full, shared_bands=shared_full),
        file.path(cfg$dir_models, "tc_modelA_ossl_only.rds"))
saveRDS(list(model=cub_B, pca=pca_full, n_pc=n_pc_full, shared_bands=shared_full),
        file.path(cfg$dir_models, "tc_modelB_ossl_local_full.rds"))
saveRDS(list(model=cub_C, pca=pca_nir,  n_pc=n_pc_nir,  shared_bands=shared_nir),
        file.path(cfg$dir_models, "tc_modelC_ossl_local_nir.rds"))

message(sprintf("  Model A  RPIQ=%.3f  [OSSL only, full PCA, n_cal=%d]",
                met_A$RPIQ, length(y_ossl)))
message(sprintf("  Model B  RPIQ=%.3f  [OSSL+local ASD, full PCA, n_cal=%d]",
                met_B$RPIQ, length(y_AB_combined)))
message(sprintf("  Model C  RPIQ=%.3f  [OSSL+local ASD, NIR PCA, ASD val]",
                met_C_asd$RPIQ))
purrr::walk(nirv_results, function(r)
  message(sprintf("  Model C  RPIQ=%.3f  [NIRVascan \u2192 %s]",
                  r$met$RPIQ, r$label)))
message(sprintf("\n  CV: %d folds x %d repeats (publication settings)",
                ossl_cfg$cub_folds, ossl_cfg$cub_repeats))
message("  Figures: 15a-15e saved to ", cfg$dir_figs)
message("  Tables:  15c_ossl_tc_all_metrics.csv saved to ", cfg$dir_tables)
message("  Models:  tc_modelA/B/C_*.rds saved to ", cfg$dir_models)

# =============================================================================
# SECTION 11 - SAVE WORKSPACE SNAPSHOT
# =============================================================================
save.image(file.path(cfg$dir_root, "scripts", "calibration_pipeline.RData"))
message("\n Pipeline complete. All outputs in: ", cfg$dir_root)
