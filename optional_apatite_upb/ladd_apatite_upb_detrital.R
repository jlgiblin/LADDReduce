#!/usr/bin/env Rscript

# Grain-by-grain detrital apatite U-Pb common-Pb correction followed by
# error-aware Gaussian mixture modelling of each sample's corrected ages.

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (i == length(x)) stop("Missing value for --", key)
    out[[key]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

required_arg <- function(args, key) {
  value <- args[[key]]
  if (is.null(value) || !nzchar(value)) stop("Missing --", key)
  value
}

as_bool <- function(x) tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
safe_name <- function(x) gsub("_+", "_", gsub("[^A-Za-z0-9._-]+", "_", trimws(as.character(x))))

logsumexp_rows <- function(x) {
  m <- apply(x, 1L, max)
  m + log(rowSums(exp(x - m)))
}

mixture_loglik <- function(y, sy, mu, tau, weight) {
  k <- length(mu)
  logd <- vapply(seq_len(k), function(j) {
    v <- sy^2 + tau[[j]]^2
    log(weight[[j]]) - 0.5 * (log(2 * pi * v) + (y - mu[[j]])^2 / v)
  }, numeric(length(y)))
  if (k == 1L) logd <- matrix(logd, ncol = 1L)
  sum(logsumexp_rows(logd))
}

run_em <- function(y, sy, initial_mu, max_iter = 300L, tolerance = 1e-7) {
  k <- length(initial_mu)
  n <- length(y)
  mu <- sort(initial_mu)
  spread <- max(stats::sd(y), diff(range(y)) / 6, stats::median(sy), 0.1)
  tau <- rep(max(spread / max(k, 1), 0.1), k)
  weight <- rep(1 / k, k)
  previous <- -Inf
  upper_tau <- max(diff(range(y)), 4 * stats::sd(y), 10 * stats::median(sy), 1)

  for (iteration in seq_len(max_iter)) {
    logd <- vapply(seq_len(k), function(j) {
      v <- sy^2 + tau[[j]]^2
      log(max(weight[[j]], 1e-12)) - 0.5 * (log(2 * pi * v) + (y - mu[[j]])^2 / v)
    }, numeric(n))
    if (k == 1L) logd <- matrix(logd, ncol = 1L)
    denom <- logsumexp_rows(logd)
    responsibility <- exp(logd - denom)
    weight <- pmax(colMeans(responsibility), 1e-8)
    weight <- weight / sum(weight)

    for (j in seq_len(k)) {
      r <- responsibility[, j]
      if (sum(r) <= 1e-8) next
      for (inner in 1:2) {
        v <- sy^2 + tau[[j]]^2
        mu[[j]] <- sum(r * y / v) / sum(r / v)
        objective <- function(candidate_tau) {
          vv <- sy^2 + candidate_tau^2
          0.5 * sum(r * (log(vv) + (y - mu[[j]])^2 / vv))
        }
        tau[[j]] <- stats::optimize(objective, interval = c(0, upper_tau), tol = 1e-8)$minimum
      }
    }

    order_mu <- order(mu)
    mu <- mu[order_mu]
    tau <- tau[order_mu]
    weight <- weight[order_mu]
    responsibility <- responsibility[, order_mu, drop = FALSE]
    current <- mixture_loglik(y, sy, mu, tau, weight)
    if (is.finite(previous) && abs(current - previous) <= tolerance * (1 + abs(previous))) break
    previous <- current
  }

  logd <- vapply(seq_len(k), function(j) {
    v <- sy^2 + tau[[j]]^2
    log(max(weight[[j]], 1e-12)) - 0.5 * (log(2 * pi * v) + (y - mu[[j]])^2 / v)
  }, numeric(n))
  if (k == 1L) logd <- matrix(logd, ncol = 1L)
  responsibility <- exp(logd - logsumexp_rows(logd))
  loglik <- mixture_loglik(y, sy, mu, tau, weight)
  effective_n <- colSums(responsibility)
  mean_se <- vapply(seq_len(k), function(j) {
    sqrt(1 / sum(responsibility[, j] / (sy^2 + tau[[j]]^2)))
  }, numeric(1))
  list(mu = mu, tau = tau, weight = weight, responsibility = responsibility,
       effective_n = effective_n, mean_se = mean_se, loglik = loglik,
       iterations = iteration)
}

fit_mixture <- function(y, sy, k, starts = 6L) {
  n <- length(y)
  if (k == 1L) {
    candidates <- list(stats::median(y))
  } else {
    probs <- seq(0, 1, length.out = k + 2L)[2:(k + 1L)]
    candidates <- list(as.numeric(stats::quantile(y, probs = probs, names = FALSE)))
    km <- tryCatch(stats::kmeans(y, centers = k, nstart = 25), error = function(e) NULL)
    if (!is.null(km)) candidates[[length(candidates) + 1L]] <- sort(as.numeric(km$centers))
    set.seed(7300 + k)
    while (length(candidates) < starts) candidates[[length(candidates) + 1L]] <- sort(sample(y, k))
  }
  fits <- lapply(candidates, function(initial) tryCatch(run_em(y, sy, initial), error = function(e) NULL))
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (!length(fits)) return(NULL)
  fit <- fits[[which.max(vapply(fits, function(z) z$loglik, numeric(1)))]]
  parameters <- 3 * k - 1
  fit$K <- k
  fit$AIC <- -2 * fit$loglik + 2 * parameters
  fit$BIC <- -2 * fit$loglik + log(n) * parameters
  fit$valid <- all(fit$effective_n >= 5) && all(fit$weight >= 0.03) &&
    (k == 1L || all(diff(fit$mu) >= 1))
  fit
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_file <- normalizePath(required_arg(args, "input"), mustWork = TRUE)
output_dir <- required_arg(args, "output")
isochron_summary_file <- if (is.null(args$isochron_summary)) NA_character_ else args$isochron_summary
max_components <- as.integer(if (is.null(args$max_components)) "4" else args$max_components)
if (is.na(max_components) || max_components < 1L || max_components > 6L) stop("max_components must be 1-6")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
local_library <- file.path(script_dir, ".r-lib")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
if (!requireNamespace("IsoplotR", quietly = TRUE)) stop("IsoplotR is not installed in .r-lib")

dat <- utils::read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("", "NA", "NaN", "nan"))
required_columns <- c("InputRow", "Sample", "Analysis", "Pb207U235", "SE_Pb207U235",
                      "Pb206U238", "SE_Pb206U238", "Rho", "Include")
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns)) stop("Input is missing: ", paste(missing_columns, collapse = ", "))
for (nm in c("InputRow", "Pb207U235", "SE_Pb207U235", "Pb206U238", "SE_Pb206U238", "Rho")) {
  dat[[nm]] <- safe_num(dat[[nm]])
}
dat$Include <- as_bool(dat$Include)
dat$Sample <- trimws(as.character(dat$Sample))
dat$Analysis <- trimws(as.character(dat$Analysis))
valid <- dat$Include & is.finite(dat$Pb207U235) & dat$Pb207U235 > 0 &
  is.finite(dat$SE_Pb207U235) & dat$SE_Pb207U235 > 0 &
  is.finite(dat$Pb206U238) & dat$Pb206U238 > 0 &
  is.finite(dat$SE_Pb206U238) & dat$SE_Pb206U238 > 0 &
  is.finite(dat$Rho) & abs(dat$Rho) < 1 & nzchar(dat$Sample)
g <- dat[valid, , drop = FALSE]

measurement_matrix <- as.matrix(g[, c("Pb207U235", "SE_Pb207U235", "Pb206U238", "SE_Pb206U238", "Rho")])
upb <- IsoplotR::read.data(measurement_matrix, method = "U-Pb", format = 1, ierr = 1)
corrected <- suppressWarnings(IsoplotR::Pb0corr(upb, option = 3))
ages <- suppressWarnings(IsoplotR::age(corrected, type = 1, exterr = FALSE))
lambda238 <- IsoplotR::settings("lambda", "U238")[1]

raw_age <- log1p(g$Pb206U238) / lambda238
raw_age_se <- g$SE_Pb206U238 / ((1 + g$Pb206U238) * lambda238)
corrected_age <- ages[, "t.conc"]
corrected_age_se <- ages[, "err[t.conc]"]
correction_ok <- is.finite(corrected_age) & corrected_age > 0 &
  is.finite(corrected_age_se) & corrected_age_se > 0
common_fraction <- 1 - corrected$x[, "Pb206U238"] / g$Pb206U238

grain_table <- data.frame(
  InputRow = g$InputRow, Sample = g$Sample, Analysis = g$Analysis,
  Raw_Pb207U235 = g$Pb207U235, Raw_SE_Pb207U235 = g$SE_Pb207U235,
  Raw_Pb206U238 = g$Pb206U238, Raw_SE_Pb206U238 = g$SE_Pb206U238, Raw_Rho = g$Rho,
  Raw_apparent_206Pb238U_age_Ma = raw_age,
  Raw_apparent_206Pb238U_age_1SE_Ma = raw_age_se,
  SK_corrected_Pb207U235 = corrected$x[, "Pb207U235"],
  SK_corrected_SE_Pb207U235 = corrected$x[, "errPb207U235"],
  SK_corrected_Pb206U238 = corrected$x[, "Pb206U238"],
  SK_corrected_SE_Pb206U238 = corrected$x[, "errPb206U238"],
  SK_corrected_Rho = corrected$x[, "rXY"],
  SK_corrected_concordia_age_Ma = corrected_age,
  SK_corrected_concordia_age_1SE_Ma = corrected_age_se,
  Approx_common_206Pb_fraction = common_fraction,
  Approx_common_206Pb_percent = 100 * common_fraction,
  CorrectionStatus = ifelse(correction_ok, "CORRECTED", "NO_POSITIVE_SOLUTION"),
  stringsAsFactors = FALSE
)
grain_table$Relative_age_1SE_pct <- 100 * grain_table$SK_corrected_concordia_age_1SE_Ma /
  grain_table$SK_corrected_concordia_age_Ma
grain_table$MixtureComponent <- NA_integer_
grain_table$MixtureMembershipProbability <- NA_real_

old_isochron <- data.frame()
if (!is.na(isochron_summary_file) && file.exists(isochron_summary_file)) {
  old_isochron <- utils::read.csv(isochron_summary_file, check.names = FALSE, stringsAsFactors = FALSE)
}

summary_rows <- list()
selection_rows <- list()
component_rows <- list()
samples <- unique(g$Sample)

for (sample_name in samples) {
  in_sample <- which(g$Sample == sample_name)
  use <- in_sample[correction_ok[in_sample]]
  y <- corrected_age[use]
  sy <- corrected_age_se[use]
  raw_y <- raw_age[in_sample]
  n <- length(y)
  fits <- lapply(seq_len(min(max_components, max(1L, floor(n / 10L)))), function(k) fit_mixture(y, sy, k))
  fits <- fits[!vapply(fits, is.null, logical(1))]
  valid_fits <- fits[vapply(fits, function(z) isTRUE(z$valid), logical(1))]
  if (!length(valid_fits)) valid_fits <- fits
  best_bic <- min(vapply(valid_fits, function(z) z$BIC, numeric(1)))
  eligible <- valid_fits[vapply(valid_fits, function(z) z$BIC <= best_bic + 6, logical(1))]
  selected <- eligible[[which.min(vapply(eligible, function(z) z$K, integer(1)))]]
  sorted_bic <- sort(vapply(valid_fits, function(z) z$BIC, numeric(1)))
  runner_delta <- if (length(sorted_bic) >= 2) sorted_bic[[2]] - sorted_bic[[1]] else NA_real_

  for (fit in fits) {
    selection_rows[[length(selection_rows) + 1L]] <- data.frame(
      Sample = sample_name, K = fit$K, N = n, LogLikelihood = fit$loglik,
      AIC = fit$AIC, BIC = fit$BIC,
      DeltaBIC = if (isTRUE(fit$valid)) fit$BIC - best_bic else NA_real_,
      ValidComponentSizes = fit$valid, Selected = fit$K == selected$K,
      Iterations = fit$iterations, stringsAsFactors = FALSE
    )
  }

  assigned <- max.col(selected$responsibility, ties.method = "first")
  probability <- apply(selected$responsibility, 1L, max)
  component_role <- ifelse(selected$effective_n < 10 | selected$weight < 0.10,
    "MINOR_TAIL_OR_OUTLIERS",
    ifelse(selected$tau / selected$mu > 0.25, "BROAD_COMPONENT", "CANDIDATE_TIGHT_MODE"))
  n_tight_modes <- sum(component_role == "CANDIDATE_TIGHT_MODE")
  sample_interpretation <- if (n_tight_modes >= 2) {
    "MULTIPLE_CANDIDATE_TIGHT_MODES"
  } else if (n_tight_modes == 1 && selected$K > 1) {
    "ONE_TIGHT_MODE_PLUS_BROAD_OR_MINOR_TAIL"
  } else if (n_tight_modes == 1) {
    "ONE_CANDIDATE_TIGHT_MODE"
  } else {
    "BROAD_DISTRIBUTION_NO_TIGHT_MODE"
  }
  grain_table$MixtureComponent[use] <- assigned
  grain_table$MixtureMembershipProbability[use] <- probability
  for (k in seq_len(selected$K)) {
    component_rows[[length(component_rows) + 1L]] <- data.frame(
      Sample = sample_name, Recommended_K = selected$K, Component = k,
      Mode_age_Ma = selected$mu[[k]], Approx_mode_age_1SE_Ma = selected$mean_se[[k]],
      Approx_mode_age_95CI_Ma = 1.96 * selected$mean_se[[k]],
      Intrinsic_dispersion_1SD_Ma = selected$tau[[k]],
      Relative_intrinsic_dispersion = selected$tau[[k]] / selected$mu[[k]],
      Proportion = selected$weight[[k]], Effective_N = selected$effective_n[[k]],
      Hard_assigned_N = sum(assigned == k),
      Median_membership_probability = stats::median(probability[assigned == k]),
      Component_role = component_role[[k]],
      stringsAsFactors = FALSE
    )
  }

  old_row <- old_isochron[old_isochron$Sample == sample_name, , drop = FALSE]
  old_age <- if (nrow(old_row) && "Age_Ma" %in% names(old_row)) safe_num(old_row$Age_Ma[[1]]) else NA_real_
  old_ci <- if (nrow(old_row) && "Age_95CI_preferred_Ma" %in% names(old_row)) safe_num(old_row$Age_95CI_preferred_Ma[[1]]) else NA_real_
  old_mswd <- if (nrow(old_row) && "MSWD" %in% names(old_row)) safe_num(old_row$MSWD[[1]]) else NA_real_

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    Sample = sample_name, N_input_valid = length(in_sample), N_corrected = n,
    N_no_positive_solution = length(in_sample) - n,
    Raw_apparent_mean_Ma = mean(raw_y), Raw_apparent_median_Ma = stats::median(raw_y),
    Raw_apparent_SD_Ma = stats::sd(raw_y), Raw_apparent_IQR_Ma = stats::IQR(raw_y),
    Corrected_mean_Ma = mean(y), Corrected_median_Ma = stats::median(y),
    Corrected_SD_Ma = stats::sd(y), Corrected_IQR_Ma = stats::IQR(y),
    Corrected_MAD_Ma = stats::mad(y), Corrected_P05_Ma = unname(stats::quantile(y, 0.05)),
    Corrected_P16_Ma = unname(stats::quantile(y, 0.16)),
    Corrected_P84_Ma = unname(stats::quantile(y, 0.84)),
    Corrected_P95_Ma = unname(stats::quantile(y, 0.95)),
    Median_shift_after_correction_Ma = stats::median(y) - stats::median(raw_y),
    IQR_retained_fraction = stats::IQR(y) / stats::IQR(raw_y),
    SD_retained_fraction = stats::sd(y) / stats::sd(raw_y),
    Median_common_206Pb_percent = stats::median(100 * common_fraction[use], na.rm = TRUE),
    Recommended_K = selected$K, Selected_BIC = selected$BIC,
    Runner_up_DeltaBIC = runner_delta,
    Population_model_strength = ifelse(is.na(runner_delta), "NO_COMPARISON",
      ifelse(runner_delta >= 10, "STRONG", ifelse(runner_delta >= 6, "MODERATE", "AMBIGUOUS"))),
    Candidate_tight_modes = n_tight_modes,
    Distribution_interpretation = sample_interpretation,
    Whole_sample_isochron_age_Ma_DIAGNOSTIC = old_age,
    Whole_sample_isochron_95CI_Ma_DIAGNOSTIC = old_ci,
    Whole_sample_isochron_MSWD_DIAGNOSTIC = old_mswd,
    stringsAsFactors = FALSE
  )

  plot_file <- file.path(output_dir, paste0(safe_name(sample_name), "_detrital_age_distribution.png"))
  grDevices::png(plot_file, width = 2200, height = 1050, res = 180)
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.8, 1.2), oma = c(0, 0, 1.5, 0))
  x_range <- range(c(raw_y, y), finite = TRUE)
  raw_density <- stats::density(raw_y, from = x_range[[1]], to = x_range[[2]], n = 1024)
  corrected_density <- stats::density(y, from = x_range[[1]], to = x_range[[2]], n = 1024)
  graphics::plot(raw_density, col = "#B24A3A", lwd = 3, xlab = "Age (Ma)", ylab = "Density",
                 main = "Before vs Stacey–Kramers correction", xlim = x_range)
  graphics::lines(corrected_density, col = "#1F6E8C", lwd = 3)
  graphics::abline(v = stats::median(raw_y), col = "#B24A3A", lty = 2)
  graphics::abline(v = stats::median(y), col = "#1F6E8C", lty = 2)
  graphics::legend("topright", legend = c("Raw apparent 206Pb/238U", "SK-corrected concordia"),
                   col = c("#B24A3A", "#1F6E8C"), lwd = 3, bty = "n")

  x2 <- seq(min(y), max(y), length.out = 1200)
  graphics::hist(y, breaks = "FD", probability = TRUE, col = "#DCEAF2", border = "white",
                 xlab = "SK-corrected concordia age (Ma)", ylab = "Density",
                 main = paste0("Error-aware Gaussian mixture (K=", selected$K, ")"))
  total_density <- rep(0, length(x2))
  colours <- grDevices::hcl.colors(max(selected$K, 3), "Dark 3")
  for (k in seq_len(selected$K)) {
    component_density <- selected$weight[[k]] * stats::dnorm(x2, selected$mu[[k]],
      sqrt(selected$tau[[k]]^2 + stats::median(sy)^2))
    total_density <- total_density + component_density
    graphics::lines(x2, component_density, col = colours[[k]], lwd = 2, lty = 2)
    graphics::abline(v = selected$mu[[k]], col = colours[[k]], lwd = 2)
  }
  graphics::lines(x2, total_density, col = "#222222", lwd = 3)
  graphics::legend("topright",
    legend = paste0("C", seq_len(selected$K), ": ", sprintf("%.1f", selected$mu),
                    " Ma (", sprintf("%.0f", 100 * selected$weight), "%)"),
    col = colours[seq_len(selected$K)], lwd = 2, bty = "n", cex = 0.82)
  graphics::mtext(paste0(sample_name, ": corrected median = ", sprintf("%.1f", stats::median(y)),
                         " Ma; n = ", n), outer = TRUE, cex = 1.25, font = 2)
  grDevices::dev.off()
}

summary_table <- do.call(rbind, summary_rows)
selection_table <- do.call(rbind, selection_rows)
component_table <- do.call(rbind, component_rows)

utils::write.csv(summary_table, file.path(output_dir, "apatite_upb_detrital_distribution_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(component_table, file.path(output_dir, "apatite_upb_detrital_population_modes.csv"), row.names = FALSE, na = "")
utils::write.csv(selection_table, file.path(output_dir, "apatite_upb_detrital_model_selection.csv"), row.names = FALSE, na = "")
utils::write.csv(grain_table, file.path(output_dir, "apatite_upb_detrital_corrected_grains.csv"), row.names = FALSE, na = "")

methods <- c(
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("IsoplotR version: ", as.character(utils::packageVersion("IsoplotR"))),
  "Purpose: detrital apatite U-Pb correction and age-distribution diagnostics.",
  "Common-Pb correction: IsoplotR Pb0corr option 3, using the age-dependent two-stage Stacey-Kramers terrestrial Pb evolution model for each grain independently.",
  "The correction uses measured 207Pb/235U and 206Pb/238U ratios, their absolute 1SE uncertainties, and rho.",
  "The previous whole-sample isochron ages are retained only as diagnostics because a detrital sample need not share one age or common-Pb composition.",
  "Raw apparent ages are uncorrected 206Pb/238U ages and are not interpreted as crystallisation ages.",
  "Corrected means and medians are unweighted distribution summaries, not weighted-mean crystallisation ages.",
  "Population model: finite Gaussian mixture with observed variance = analytical age variance + component intrinsic variance.",
  paste0("Candidate component counts: 1 through ", max_components, ". Models require effective N >= 5 and component proportion >= 0.03."),
  "Selection: lowest-BIC family, with the simplest model within 6 BIC units of the minimum preferred to avoid unnecessary splitting.",
  "Mode-age uncertainties are local likelihood approximations and do not propagate uncertainty in the Stacey-Kramers evolution model.",
  "A candidate tight mode has effective N >= 10, proportion >= 0.10, and intrinsic 1SD <= 25% of its mode age. Other components are labelled broad or minor-tail/outlier components rather than firm populations.",
  "Mixture components are exploratory provenance/thermal populations and require geological validation; they are not automatic rejection groups.",
  "No grain was removed for being old, young, or inconvenient. Failed positive-age solutions are retained and flagged."
)
writeLines(methods, file.path(output_dir, "apatite_upb_detrital_methods.txt"))

cat("Detrital apatite U-Pb correction complete.\n")
cat("Samples:", length(samples), "\n")
cat("Corrected grains:", sum(correction_ok), "of", nrow(g), "valid inputs\n")
cat("Output:", output_dir, "\n")
