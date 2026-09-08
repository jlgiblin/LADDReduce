#!/usr/bin/env Rscript

# Normalize a folder of iolite apatite U-Pb CSV files and run the project
# IsoplotR semitotal-Pb/U isochron backend in one batch. Raw CSVs are not edited.

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

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

as_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c(
    "1", "true", "t", "yes", "y", "include", "included", "keep"
  )
}

normalized_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

resolve_column <- function(names, role, required = TRUE) {
  n <- normalized_name(names)
  is_error <- grepl("err|error|unc|sigma|se|sd|2s|1s", n)
  is_age <- grepl("age", n)

  score <- switch(
    role,
    sample = 100 * n %in% c("sample", "samplename", "sampleid", "group", "groupname") +
      10 * grepl("sample", n) + 3 * grepl("group", n),
    analysis = 100 * n %in% c("analysis", "analysisname", "name", "grainid", "spot", "spotname", "label") +
      10 * grepl("analysis", n) + 8 * grepl("grain", n) + 6 * grepl("spot", n),
    ratio75 = 20 * (grepl("207", n) & grepl("235", n)) - 30 * is_error - 30 * is_age + 5 * grepl("mean", n),
    error75 = 20 * (grepl("207", n) & grepl("235", n)) + 20 * is_error - 30 * is_age,
    ratio68 = 20 * (grepl("206", n) & grepl("238", n)) - 30 * is_error - 30 * is_age + 5 * grepl("mean", n),
    error68 = 20 * (grepl("206", n) & grepl("238", n)) + 20 * is_error - 30 * is_age,
    rho = 100 * n %in% c("rho", "rxy", "errorcorrelation", "correlation") +
      20 * grepl("rho", n) + 10 * grepl("correlation", n) +
      3 * (grepl("207", n) & grepl("235", n) & grepl("206", n) & grepl("238", n)),
    include = 100 * n %in% c("include", "use", "accepted", "accept", "keep"),
    stop("Unknown column role: ", role)
  )

  threshold <- if (role %in% c("ratio75", "error75", "ratio68", "error68")) 20 else 1
  best <- max(score)
  if (!is.finite(best) || best < threshold) {
    if (required) stop("Could not identify ", role, " column from: ", paste(names, collapse = ", "))
    return(NA_character_)
  }
  ties <- which(score == best)
  if (required && length(ties) > 1L) {
    stop("More than one possible ", role, " column: ", paste(names[ties], collapse = ", "))
  }
  names[ties[[1L]]]
}

convert_errors <- function(r75, e75, r68, e68, mode) {
  switch(
    tolower(mode),
    "1se_abs" = list(e75 = e75, e68 = e68),
    "2se_abs" = list(e75 = e75 / 2, e68 = e68 / 2),
    "1se_pct" = list(e75 = abs(r75) * e75 / 100, e68 = abs(r68) * e68 / 100),
    "2se_pct" = list(e75 = abs(r75) * e75 / 200, e68 = abs(r68) * e68 / 200),
    stop("input_errors must be 1se_abs, 2se_abs, 1se_pct, or 2se_pct")
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- normalizePath(required_arg(args, "input_dir"), mustWork = TRUE)
output_dir <- if (is.null(args$output_dir)) file.path(input_dir, "isochron_results") else args$output_dir
input_errors <- if (is.null(args$input_errors)) "2se_abs" else tolower(args$input_errors)
column_overrides_file <- if (is.null(args$column_overrides)) NA_character_ else args$column_overrides
alpha <- if (is.null(args$alpha)) "0.05" else args$alpha
minimum_n <- if (is.null(args$minimum_n)) "4" else args$minimum_n
diseq_mode <- if (is.null(args$diseq_mode)) "none" else args$diseq_mode
th230_u238 <- if (is.null(args$th230_u238)) "1" else args$th230_u238
th230_u238_se <- if (is.null(args$th230_u238_se)) "0" else args$th230_u238_se

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

files <- sort(list.files(input_dir, pattern = "[.]csv$", full.names = TRUE, recursive = FALSE))
if (!length(files)) stop("No CSV files found in ", input_dir)

column_overrides <- data.frame()
if (!is.na(column_overrides_file)) {
  column_overrides_file <- normalizePath(column_overrides_file, mustWork = TRUE)
  column_overrides <- utils::read.csv(column_overrides_file, check.names = FALSE,
                                      stringsAsFactors = FALSE)
  override_fields <- c("SourceFile", "NormalizedField", "SourceColumn")
  missing_override_fields <- setdiff(override_fields, names(column_overrides))
  if (length(missing_override_fields)) {
    stop("Column override file is missing: ", paste(missing_override_fields, collapse = ", "))
  }
}

normalized_rows <- list()
mapping_rows <- list()
inventory_rows <- list()
row_map_rows <- list()
next_input_row <- 1L

for (file in files) {
  raw <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("", "NA", "NaN", "nan"))
  names_raw <- names(raw)
  cols <- c(
    Sample = resolve_column(names_raw, "sample", FALSE),
    Analysis = resolve_column(names_raw, "analysis", FALSE),
    Pb207U235 = resolve_column(names_raw, "ratio75", TRUE),
    SE_Pb207U235 = resolve_column(names_raw, "error75", TRUE),
    Pb206U238 = resolve_column(names_raw, "ratio68", TRUE),
    SE_Pb206U238 = resolve_column(names_raw, "error68", TRUE),
    Rho = resolve_column(names_raw, "rho", FALSE),
    Include = resolve_column(names_raw, "include", FALSE)
  )
  if (nrow(column_overrides)) {
    use_overrides <- column_overrides[column_overrides$SourceFile == basename(file), , drop = FALSE]
    for (j in seq_len(nrow(use_overrides))) {
      field <- use_overrides$NormalizedField[[j]]
      source <- use_overrides$SourceColumn[[j]]
      if (!field %in% names(cols)) stop("Unknown override field for ", basename(file), ": ", field)
      if (!source %in% names_raw) stop("Override source column not found in ", basename(file), ": ", source)
      cols[[field]] <- source
    }
  }

  n <- nrow(raw)
  sample <- if (is.na(cols[["Sample"]])) rep(tools::file_path_sans_ext(basename(file)), n) else trimws(as.character(raw[[cols[["Sample"]]]]))
  analysis <- if (is.na(cols[["Analysis"]])) paste0("row_", seq_len(n)) else trimws(as.character(raw[[cols[["Analysis"]]]]))
  r75 <- safe_num(raw[[cols[["Pb207U235"]]]])
  e75 <- safe_num(raw[[cols[["SE_Pb207U235"]]]])
  r68 <- safe_num(raw[[cols[["Pb206U238"]]]])
  e68 <- safe_num(raw[[cols[["SE_Pb206U238"]]]])
  converted <- convert_errors(r75, e75, r68, e68, input_errors)

  rho_missing <- is.na(cols[["Rho"]])
  rho <- if (rho_missing) rep(0, n) else safe_num(raw[[cols[["Rho"]]]])
  rho_assumed <- rep(rho_missing, n)
  bad_rho <- !is.finite(rho) | abs(rho) > 1
  rho[bad_rho] <- 0
  rho_assumed[bad_rho] <- TRUE
  include <- if (is.na(cols[["Include"]])) rep(TRUE, n) else as_bool(raw[[cols[["Include"]]]])

  input_rows <- seq.int(next_input_row, length.out = n)
  next_input_row <- next_input_row + n
  normalized_rows[[length(normalized_rows) + 1L]] <- data.frame(
    InputRow = input_rows, Sample = sample, Analysis = analysis,
    Pb207U235 = r75, SE_Pb207U235 = converted$e75,
    Pb206U238 = r68, SE_Pb206U238 = converted$e68,
    Rho = rho, Include = include, RhoWasAssumed = rho_assumed,
    stringsAsFactors = FALSE
  )
  row_map_rows[[length(row_map_rows) + 1L]] <- data.frame(
    InputRow = input_rows, SourceFile = basename(file), SourceRow = seq_len(n) + 1L,
    Sample = sample, Analysis = analysis, stringsAsFactors = FALSE
  )
  mapping_rows[[length(mapping_rows) + 1L]] <- data.frame(
    SourceFile = basename(file), NormalizedField = names(cols), SourceColumn = unname(cols),
    InputErrorMode = input_errors, stringsAsFactors = FALSE
  )

  valid <- is.finite(r75) & r75 > 0 & is.finite(converted$e75) & converted$e75 > 0 &
    is.finite(r68) & r68 > 0 & is.finite(converted$e68) & converted$e68 > 0 &
    is.finite(rho) & abs(rho) < 1 & nzchar(sample)
  inventory_rows[[length(inventory_rows) + 1L]] <- data.frame(
    SourceFile = basename(file), Rows = n, Samples = paste(unique(sample[nzchar(sample)]), collapse = "|"),
    IncludedRows = sum(include), ValidIncludedRows = sum(include & valid),
    InvalidOrExcludedRows = sum(!include | !valid), RhoAssumedZeroRows = sum(rho_assumed),
    stringsAsFactors = FALSE
  )
}

normalized <- do.call(rbind, normalized_rows)
row_map <- do.call(rbind, row_map_rows)
mapping <- do.call(rbind, mapping_rows)
inventory <- do.call(rbind, inventory_rows)

normalized_file <- file.path(output_dir, "apatite_upb_normalized_input.csv")
utils::write.csv(normalized, normalized_file, row.names = FALSE, na = "")
utils::write.csv(row_map, file.path(output_dir, "apatite_upb_input_row_mapping.csv"), row.names = FALSE, na = "")
utils::write.csv(mapping, file.path(output_dir, "apatite_upb_column_mapping.csv"), row.names = FALSE, na = "")
utils::write.csv(inventory, file.path(output_dir, "apatite_upb_input_inventory.csv"), row.names = FALSE, na = "")
if (nrow(column_overrides)) {
  utils::write.csv(column_overrides, file.path(output_dir, "apatite_upb_applied_column_overrides.csv"),
                   row.names = FALSE, na = "")
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
backend <- file.path(script_dir, "ladd_apatite_upb_isoplotr.R")
if (!file.exists(backend)) stop("Backend not found: ", backend)

status <- system2(
  "Rscript",
  c(
    shQuote(backend), "--input", shQuote(normalized_file),
    "--output", shQuote(output_dir), "--alpha", alpha,
    "--minimum_n", minimum_n, "--diseq_mode", diseq_mode,
    "--th230_u238", th230_u238, "--th230_u238_se", th230_u238_se
  )
)
if (!identical(status, 0L)) stop("Isochron backend failed with status ", status)

cat("Batch reduction complete.\n")
cat("Input files:", length(files), "\n")
cat("Input rows:", nrow(normalized), "\n")
cat("Output:", output_dir, "\n")
