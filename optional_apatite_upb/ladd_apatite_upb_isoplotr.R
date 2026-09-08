#!/usr/bin/env Rscript

# Batch apatite U-Pb semitotal-Pb/U isochron reduction using IsoplotR.
# The input is a normalized CSV written by ladd_apatite_upb_isochron.m.

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

as_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_name <- function(x) {
  y <- gsub("[^A-Za-z0-9._-]+", "_", trimws(as.character(x)))
  y <- gsub("_+", "_", y)
  ifelse(nzchar(y), y, "sample")
}

capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(list(error = conditionMessage(e)), class = "fit_error")),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

is_fit_error <- function(x) inherits(x, "fit_error")

fit_fields <- function(fit, alpha) {
  age <- unname(fit$par["t"])
  se <- sqrt(unname(fit$cov["t", "t"]))
  mswd <- if (is.null(fit$mswd)) NA_real_ else unname(fit$mswd)
  p <- if (is.null(fit$p.value)) NA_real_ else unname(fit$p.value)
  df <- if (is.null(fit$df)) NA_real_ else unname(fit$df)
  common <- unname(fit$par["a0"])
  common_se <- sqrt(unname(fit$cov["a0", "a0"]))
  z <- stats::qnorm(1 - alpha / 2)
  ci_internal <- z * se
  ci_common <- z * common_se
  inflated <- is.finite(p) && p < alpha && is.finite(mswd) && mswd > 0
  if (inflated) {
    se_preferred <- se * sqrt(mswd)
    common_se_preferred <- common_se * sqrt(mswd)
    t_factor <- if (is.finite(df) && df > 0) stats::qt(1 - alpha / 2, df) else z
    ci_preferred <- t_factor * se_preferred
    ci_common_preferred <- t_factor * common_se_preferred
  } else {
    se_preferred <- se
    common_se_preferred <- common_se
    ci_preferred <- ci_internal
    ci_common_preferred <- ci_common
  }
  list(
    age = age,
    se = se,
    ci = ci_internal,
    se_preferred = se_preferred,
    ci_preferred = ci_preferred,
    common = common,
    common_se = common_se,
    common_ci = ci_common,
    common_se_preferred = common_se_preferred,
    common_ci_preferred = ci_common_preferred,
    mswd = mswd,
    p = p,
    df = df,
    inflated = inflated,
    convergence = if (is.null(fit$convergence)) NA_integer_ else fit$convergence,
    message = if (is.null(fit$message)) "" else fit$message
  )
}

make_diseq <- function(mode, ratio, ratio_se) {
  if (identical(mode, "none")) return(IsoplotR::diseq())
  if (identical(mode, "fixed_th230_u238")) {
    return(IsoplotR::diseq(ThU = list(x = ratio, sx = ratio_se, option = 1)))
  }
  stop("Unsupported disequilibrium mode: ", mode,
       ". Use 'none' or 'fixed_th230_u238'.")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_file <- normalizePath(required_arg(args, "input"), mustWork = TRUE)
output_dir <- required_arg(args, "output")
alpha <- safe_num(if (is.null(args$alpha)) "0.05" else args$alpha)
minimum_n <- as.integer(if (is.null(args$minimum_n)) "4" else args$minimum_n)
diseq_mode <- if (is.null(args$diseq_mode)) "none" else tolower(args$diseq_mode)
th230_u238 <- safe_num(if (is.null(args$th230_u238)) "1" else args$th230_u238)
th230_u238_se <- safe_num(if (is.null(args$th230_u238_se)) "0" else args$th230_u238_se)

if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
if (is.na(minimum_n) || minimum_n < 3) stop("minimum_n must be at least 3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
local_library <- file.path(script_dir, ".r-lib")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
if (!requireNamespace("IsoplotR", quietly = TRUE)) {
  stop("IsoplotR is not installed. From this folder run: ",
       "Rscript -e \"dir.create('.r-lib',showWarnings=FALSE); ",
       "install.packages('IsoplotR',repos='https://cloud.r-project.org',lib='.r-lib')\"")
}

IsoplotR::settings("alpha", alpha)
dat <- utils::read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("", "NA", "NaN", "nan"))
required_columns <- c("InputRow", "Sample", "Analysis", "Pb207U235", "SE_Pb207U235",
                      "Pb206U238", "SE_Pb206U238", "Rho", "Include", "RhoWasAssumed")
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns)) stop("Normalized input is missing: ", paste(missing_columns, collapse = ", "))

numeric_columns <- c("InputRow", "Pb207U235", "SE_Pb207U235", "Pb206U238", "SE_Pb206U238", "Rho")
for (nm in numeric_columns) dat[[nm]] <- safe_num(dat[[nm]])
dat$Include <- as_bool(dat$Include)
dat$RhoWasAssumed <- as_bool(dat$RhoWasAssumed)
dat$Sample <- trimws(as.character(dat$Sample))
dat$Analysis <- trimws(as.character(dat$Analysis))

valid_measurement <- is.finite(dat$Pb207U235) & dat$Pb207U235 > 0 &
  is.finite(dat$SE_Pb207U235) & dat$SE_Pb207U235 > 0 &
  is.finite(dat$Pb206U238) & dat$Pb206U238 > 0 &
  is.finite(dat$SE_Pb206U238) & dat$SE_Pb206U238 > 0 &
  is.finite(dat$Rho) & abs(dat$Rho) < 1 & nzchar(dat$Sample)
dat$FitEligible <- dat$Include & valid_measurement
dat$InputStatus <- ifelse(!dat$Include, "EXCLUDED_BY_INPUT",
                          ifelse(!valid_measurement, "INVALID_MEASUREMENT", "FIT_ELIGIBLE"))

diseq_object <- make_diseq(diseq_mode, th230_u238, th230_u238_se)
summary_rows <- list()
corrected_rows <- list()
samples <- unique(dat$Sample[nzchar(dat$Sample)])
u_ratio <- IsoplotR::settings("iratio", "U238U235")[1]
lambda238 <- IsoplotR::settings("lambda", "U238")[1]
lambda235 <- IsoplotR::settings("lambda", "U235")[1]

for (sample_name in samples) {
  take <- which(dat$Sample == sample_name & dat$FitEligible)
  n <- length(take)
  if (n < minimum_n) {
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      Sample = sample_name, N = n, Status = "LOW_N", Notes = paste0("Need at least ", minimum_n, " valid analyses"),
      Age_Ma = NA, Age_1SE_internal_Ma = NA, Age_95CI_internal_Ma = NA,
      Age_1SE_preferred_Ma = NA, Age_95CI_preferred_Ma = NA,
      Age_1SE_external_Ma = NA, Age_95CI_external_Ma = NA,
      CommonPb207Pb206 = NA, CommonPb207Pb206_1SE = NA, CommonPb207Pb206_95CI = NA,
      MSWD = NA, P_value = NA, DF = NA, ErrorInflatedForMSWD = FALSE,
      U238Pb206_range = NA, U238Pb206_relative_range = NA,
      Rho_assumed_zero = any(dat$RhoWasAssumed[take]),
      Disequilibrium_mode = diseq_mode, Th230U238_initial = ifelse(diseq_mode == "none", NA, th230_u238),
      Model3_age_Ma = NA, Model3_age_1SE_Ma = NA, Model3_overdispersion_Ma = NA,
      Model3_convergence = NA, stringsAsFactors = FALSE
    )
    next
  }

  g <- dat[take, , drop = FALSE]
  measurement_matrix <- as.matrix(g[, c("Pb207U235", "SE_Pb207U235", "Pb206U238", "SE_Pb206U238", "Rho")])
  upb <- IsoplotR::read.data(measurement_matrix, method = "U-Pb", format = 1,
                             ierr = 1, d = diseq_object)
  internal_attempt <- capture_warnings(IsoplotR::isochron(upb, plot = FALSE, model = 1, exterr = FALSE))
  if (is_fit_error(internal_attempt$value)) {
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      Sample = sample_name, N = n, Status = "FIT_FAILED", Notes = internal_attempt$value$error,
      Age_Ma = NA, Age_1SE_internal_Ma = NA, Age_95CI_internal_Ma = NA,
      Age_1SE_preferred_Ma = NA, Age_95CI_preferred_Ma = NA,
      Age_1SE_external_Ma = NA, Age_95CI_external_Ma = NA,
      CommonPb207Pb206 = NA, CommonPb207Pb206_1SE = NA, CommonPb207Pb206_95CI = NA,
      MSWD = NA, P_value = NA, DF = NA, ErrorInflatedForMSWD = FALSE,
      U238Pb206_range = diff(range(1 / g$Pb206U238)),
      U238Pb206_relative_range = diff(range(1 / g$Pb206U238)) / stats::median(1 / g$Pb206U238),
      Rho_assumed_zero = any(g$RhoWasAssumed), Disequilibrium_mode = diseq_mode,
      Th230U238_initial = ifelse(diseq_mode == "none", NA, th230_u238),
      Model3_age_Ma = NA, Model3_age_1SE_Ma = NA, Model3_overdispersion_Ma = NA,
      Model3_convergence = NA, stringsAsFactors = FALSE
    )
    next
  }

  fit_internal <- internal_attempt$value
  fields <- fit_fields(fit_internal, alpha)
  external_attempt <- capture_warnings(IsoplotR::isochron(upb, plot = FALSE, model = 1, exterr = TRUE))
  if (!is_fit_error(external_attempt$value)) {
    external_fields <- fit_fields(external_attempt$value, alpha)
    external_se <- external_fields$se_preferred
    external_ci <- external_fields$ci_preferred
  } else {
    external_se <- NA_real_
    external_ci <- NA_real_
  }

  model3_age <- model3_se <- model3_disp <- model3_convergence <- NA_real_
  model3_notes <- character()
  if (is.finite(fields$p) && fields$p < alpha && n >= 4) {
    model3_attempt <- capture_warnings(IsoplotR::isochron(upb, plot = FALSE, model = 3, exterr = FALSE))
    model3_notes <- model3_attempt$warnings
    if (!is_fit_error(model3_attempt$value)) {
      model3_fit <- model3_attempt$value
      model3_age <- unname(model3_fit$par["t"])
      model3_se <- sqrt(unname(model3_fit$cov["t", "t"]))
      model3_disp <- if (is.null(model3_fit$disp)) NA_real_ else unname(model3_fit$disp["w"])
      model3_convergence <- if (is.null(model3_fit$convergence)) NA_real_ else model3_fit$convergence
    } else {
      model3_notes <- c(model3_notes, model3_attempt$value$error)
    }
  }

  status <- "OK"
  status_notes <- c(internal_attempt$warnings, external_attempt$warnings, model3_notes)
  if (is.finite(fields$p) && fields$p < alpha) status <- "EXCESS_SCATTER"
  if (!identical(fields$convergence, 0L) || !is.finite(fields$age) || fields$age <= 0) status <- "FIT_WARNING"
  if (any(g$RhoWasAssumed)) status_notes <- c(status_notes, "Ratio correlation was unavailable and set to zero")
  if (fields$common < 0 || fields$common > 2) {
    status <- "FIT_WARNING"
    status_notes <- c(status_notes, "Fitted common 207Pb/206Pb is outside 0-2")
  }
  status_notes <- paste(unique(status_notes[nzchar(status_notes)]), collapse = " | ")

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    Sample = sample_name, N = n, Status = status, Notes = status_notes,
    Age_Ma = fields$age,
    Age_1SE_internal_Ma = fields$se,
    Age_95CI_internal_Ma = fields$ci,
    Age_1SE_preferred_Ma = fields$se_preferred,
    Age_95CI_preferred_Ma = fields$ci_preferred,
    Age_1SE_external_Ma = external_se,
    Age_95CI_external_Ma = external_ci,
    CommonPb207Pb206 = fields$common,
    CommonPb207Pb206_1SE = fields$common_se_preferred,
    CommonPb207Pb206_95CI = fields$common_ci_preferred,
    MSWD = fields$mswd, P_value = fields$p, DF = fields$df,
    ErrorInflatedForMSWD = fields$inflated,
    U238Pb206_range = diff(range(1 / g$Pb206U238)),
    U238Pb206_relative_range = diff(range(1 / g$Pb206U238)) / stats::median(1 / g$Pb206U238),
    Rho_assumed_zero = any(g$RhoWasAssumed), Disequilibrium_mode = diseq_mode,
    Th230U238_initial = ifelse(diseq_mode == "none", NA, th230_u238),
    Model3_age_Ma = model3_age, Model3_age_1SE_Ma = model3_se,
    Model3_overdispersion_Ma = model3_disp, Model3_convergence = model3_convergence,
    stringsAsFactors = FALSE
  )

  correction_attempt <- capture_warnings(IsoplotR::Pb0corr(upb, option = 2))
  if (!is_fit_error(correction_attempt$value)) {
    corrected <- correction_attempt$value
    ages_attempt <- capture_warnings(IsoplotR::age(corrected, type = 1, exterr = FALSE))
    corrected_matrix <- corrected$x
    if (!is_fit_error(ages_attempt$value)) {
      ages <- ages_attempt$value
      age68 <- ages[, "t.68"]
      se68 <- ages[, "err[t.68]"]
      age75 <- ages[, "t.75"]
      se75 <- ages[, "err[t.75]"]
    } else {
      age68 <- se68 <- age75 <- se75 <- rep(NA_real_, n)
    }
    y_measured <- g$Pb207U235 / (g$Pb206U238 * u_ratio)
    y_radiogenic <- expm1(lambda235 * fields$age) / (expm1(lambda238 * fields$age) * u_ratio)
    common_fraction_206 <- (y_measured - y_radiogenic) / (fields$common - y_radiogenic)
    corrected_rows[[length(corrected_rows) + 1L]] <- data.frame(
      InputRow = g$InputRow, Sample = sample_name, Analysis = g$Analysis,
      Raw_Pb207U235 = g$Pb207U235, Raw_SE_Pb207U235 = g$SE_Pb207U235,
      Raw_Pb206U238 = g$Pb206U238, Raw_SE_Pb206U238 = g$SE_Pb206U238,
      Raw_Rho = g$Rho,
      Corrected_Pb207U235 = corrected_matrix[, "Pb207U235"],
      Corrected_SE_Pb207U235 = corrected_matrix[, "errPb207U235"],
      Corrected_Pb206U238 = corrected_matrix[, "Pb206U238"],
      Corrected_SE_Pb206U238 = corrected_matrix[, "errPb206U238"],
      Corrected_Rho = corrected_matrix[, "rXY"],
      Corrected_206Pb238U_age_Ma = age68,
      Corrected_206Pb238U_1SE_analytical_Ma = se68,
      Corrected_207Pb235U_age_Ma = age75,
      Corrected_207Pb235U_1SE_analytical_Ma = se75,
      Approx_common_206Pb_fraction = common_fraction_206,
      Approx_common_206Pb_percent = 100 * common_fraction_206,
      stringsAsFactors = FALSE
    )
  }

  plot_file <- file.path(output_dir, paste0(safe_name(sample_name), "_TeraWasserburg_isochron.png"))
  plot_attempt <- capture_warnings({
    grDevices::png(plot_file, width = 1800, height = 1400, res = 200)
    on.exit(grDevices::dev.off(), add = TRUE)
    IsoplotR::isochron(upb, model = 1, type = 2, plot = TRUE, oerr = 3,
                      sigdig = 2, show.numbers = TRUE, title = TRUE)
    grDevices::dev.off()
    on.exit(NULL, add = FALSE)
  })
}

summary_table <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
corrected_table <- if (length(corrected_rows)) do.call(rbind, corrected_rows) else data.frame()
excluded_table <- dat[!dat$FitEligible, c("InputRow", "Sample", "Analysis", "InputStatus"), drop = FALSE]

utils::write.csv(summary_table, file.path(output_dir, "apatite_upb_isochron_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(corrected_table, file.path(output_dir, "apatite_upb_commonPb_corrected_analyses.csv"), row.names = FALSE, na = "")
utils::write.csv(excluded_table, file.path(output_dir, "apatite_upb_excluded_or_invalid.csv"), row.names = FALSE, na = "")

methods_lines <- c(
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("IsoplotR version: ", as.character(utils::packageVersion("IsoplotR"))),
  "Regression: Ludwig (1998) semitotal-Pb/U maximum likelihood isochron.",
  "Input format: 207Pb/235U, 206Pb/238U, 1SE absolute uncertainties, and their correlation.",
  "A finite ratio-error correlation strictly between -1 and 1 is required; exact +/-1 gives a singular covariance matrix and is excluded as invalid.",
  paste0("Probability cutoff alpha: ", alpha),
  "Model 1 uncertainties are reported both internally and, when p < alpha, expanded by sqrt(MSWD) with a t-based confidence interval.",
  "External uncertainty includes IsoplotR decay-constant and 238U/235U uncertainties; it does not invent missing iolite session-level systematic uncertainty.",
  paste0("Disequilibrium mode: ", diseq_mode),
  if (diseq_mode == "none") "No U-series disequilibrium correction applied." else paste0("Fixed initial 230Th/238U activity ratio: ", th230_u238, " +/- ", th230_u238_se, " (1SE)."),
  "Per-analysis corrected dates are diagnostics; their listed analytical errors do not represent the shared uncertainty of the fitted common-Pb composition.",
  "No row was automatically rejected for statistical discordance or excess scatter."
)
writeLines(methods_lines, file.path(output_dir, "apatite_upb_methods.txt"))

cat("Apatite U-Pb isochron reduction complete.\n")
cat("Samples attempted:", length(samples), "\n")
cat("Summary:", file.path(output_dir, "apatite_upb_isochron_summary.csv"), "\n")
