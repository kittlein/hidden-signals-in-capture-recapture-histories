# ==================================================================================================
# Evaluate generated rodent capture histories with the five POPAN models exported by their generator.
# Run the entire script after editing the settings below, using one species and sampling design per run.
# ==================================================================================================

# 1. Select the generated dataset and the histories to evaluate, then check the required input files.
dataset_directory <- "akodon_R_10_sessions_Nx1" # Select the experiment root, not its images folder.
evaluation_source <- "cnn_validation" # Use cnn_validation for paired comparison, or generated_sample before CNN training.
cnn_run_directory <- "" # Set a CNN run folder, or leave empty to select the latest run with validation predictions.
histories_per_class <- 100L # Evaluate at most this many histories per class, or use zero for all eligible histories.
evaluation_seed <- 123L # Reproduce the sample selected within each generating class without changing the CNN split.
evaluate_observed <- TRUE # Also fit observed_capture_histories.csv.gz when the generator exported that file.
require_all_models <- TRUE # Leave a history unclassified when any of the five candidates lacks a usable AICc.
adjust_parameter_count <- TRUE # Use RMark's design-based parameter count and also retain the unadjusted MARK values.
use_initial_model <- TRUE # Fit a common simple model to initialize candidates using only the current capture history.
mark_threads <- 1L # Control the MARK threads used by each fit, without evaluating different histories in parallel.
mark_retry <- 0L # Set the number of additional RMark attempts for models reporting singular beta parameters.
mark_directory <- "" # Optional directory containing the installed MARK executable, for example C:/Program Files/Mark.
mark_work_directory <- "" # Optional temporary parent directory, otherwise use the current R temporary directory.
classes <- c("p1", "p2", "p3", "p4", "p5") # Read the meaning of these labels from the selected generator metadata.
parameter_names <- c("Phi", "p", "pent", "N") # Fit survival, capture, entry and superpopulation size in POPAN.

if (!dir.exists(dataset_directory)) stop("Dataset directory does not exist: ", dataset_directory)
dataset_directory <- normalizePath(dataset_directory, winslash = "/", mustWork = TRUE) # Resolve the experiment root.
if (!evaluation_source %in% c("cnn_validation", "generated_sample")) stop("Use cnn_validation or generated_sample.")
if (length(histories_per_class) != 1L || !is.finite(histories_per_class) || histories_per_class < 0 || histories_per_class != floor(histories_per_class)) stop("histories_per_class must be a nonnegative integer.")
if (length(evaluation_seed) != 1L || !is.finite(evaluation_seed) || evaluation_seed < 0 || evaluation_seed > .Machine$integer.max || evaluation_seed != floor(evaluation_seed)) stop("evaluation_seed must fit R's nonnegative integer range.")
if (length(mark_threads) != 1L || !is.finite(mark_threads) || mark_threads == 0 || mark_threads != floor(mark_threads)) stop("mark_threads must be a nonzero integer accepted by MARK.")
if (length(mark_retry) != 1L || !is.finite(mark_retry) || mark_retry < 0 || mark_retry != floor(mark_retry)) stop("mark_retry must be a nonnegative integer.")
logical_settings <- c(evaluate_observed, require_all_models, adjust_parameter_count, use_initial_model)
if (length(logical_settings) != 4L || !is.logical(logical_settings) || anyNA(logical_settings)) stop("The four logical settings must each be TRUE or FALSE.")
required_files <- file.path(dataset_directory, c("capture_history_index.csv", "class_definitions.csv", "sessions.csv"))
if (any(!file.exists(required_files))) stop("Missing generator files: ", paste(required_files[!file.exists(required_files)], collapse = ", "))
pair_index <- read.csv(required_files[1L], colClasses = "character", stringsAsFactors = FALSE, check.names = FALSE)
model_definitions <- read.csv(required_files[2L], colClasses = "character", stringsAsFactors = FALSE, check.names = FALSE)
sessions <- read.csv(required_files[3L], stringsAsFactors = FALSE, check.names = FALSE) # Preserve the exported calendar.
if (!all(c("scenario", "replicate", "image_file", "history_file") %in% names(pair_index))) stop("capture_history_index.csv must contain scenario, replicate, image_file and history_file.")
if (!all(c("scenario", parameter_names) %in% names(model_definitions))) stop("class_definitions.csv must contain scenario, Phi, p, pent and N.")
if (anyDuplicated(model_definitions$scenario) || !setequal(model_definitions$scenario, classes)) stop("Provide exactly one model definition for each class p1 to p5.")
model_definitions <- model_definitions[match(classes, model_definitions$scenario), , drop = FALSE]
if (!nrow(pair_index) || anyNA(pair_index[, c("scenario", "replicate", "image_file", "history_file")]) || any(!pair_index$scenario %in% classes)) stop("The paired-file index contains missing values or unsupported class labels.")
for (path_column in c("image_file", "history_file")) { # Resolve the actual generator pairs rather than guessing filenames.
    index_paths <- gsub("\\\\", "/", pair_index[[path_column]]) # Accept paths exported on either Windows or Linux.
    absolute_paths <- grepl("^(/|[A-Za-z]:/)", index_paths) # Keep any absolute path already present in the index.
    if (any(!nzchar(index_paths))) stop("The paired-file index contains an empty ", path_column, " path.")
    relative_paths <- index_paths[!absolute_paths] # Generators store paths relative to their images or histories folder.
    component_directory <- if (path_column == "image_file") "images" else "capture_histories" # Respect the index schema.
    root_candidates <- file.path(dataset_directory, relative_paths) # Also accept index paths already including this folder.
    component_candidates <- file.path(dataset_directory, component_directory, relative_paths) # Resolve the current schema.
    if (any(file.exists(root_candidates) & file.exists(component_candidates))) stop("Ambiguous indexed paths exist in both the dataset root and its component folder.")
    index_paths[!absolute_paths] <- ifelse(file.exists(root_candidates), root_candidates, component_candidates)
    if (any(!file.exists(index_paths))) stop("Indexed files do not exist: ", paste(head(index_paths[!file.exists(index_paths)], 5L), collapse = ", "))
    if (any(file.info(index_paths)$isdir)) stop("An indexed ", path_column, " path is a directory instead of a file.")
    pair_index[[path_column]] <- normalizePath(index_paths, winslash = "/", mustWork = TRUE)
}
if (anyDuplicated(pair_index$image_file) || anyDuplicated(pair_index$history_file)) stop("Each image and history must appear in exactly one index row.")

# 2. Read the calendar in months and prepare the covariates required by the exported candidate formulas.
n_occasions <- nrow(sessions) # Use ten, thirty or 370 sessions according to the selected generated experiment.
if (n_occasions < 2L || !"session" %in% names(sessions) || anyNA(sessions$session) || !identical(as.integer(sessions$session), seq_len(n_occasions))) stop("sessions.csv must be ordered and numbered consecutively from one.")
interval_column <- intersect(c("interval_months", "time_interval_months"), names(sessions))
if (length(interval_column) != 1L) stop("sessions.csv must contain one interval_months or time_interval_months column.")
time_intervals <- as.numeric(sessions[[interval_column]][-1L]) # MARK survival intervals use the generator's month units.
if (any(!is.finite(time_intervals)) || any(time_intervals <= 0)) stop("All intervals after the first session must be positive months.")
occasion_times <- c(1, 1 + cumsum(time_intervals)) # Match RMark design times using elapsed months, not session numbers.
calendar_covariates <- list() # Keep factor levels consistent across both sexes and all five competing models.
if ("season" %in% names(sessions)) calendar_covariates$season <- factor(sessions$season, levels = c("cold", "warm"))
if ("year" %in% names(sessions)) calendar_covariates$yr <- factor(sessions$year, levels = sort(unique(sessions$year)))
if ("season2" %in% names(sessions)) calendar_covariates$season2 <- factor(sessions$season2, levels = c("N", "B"))
if ("season_rain2" %in% names(sessions)) calendar_covariates$season_rain2 <- factor(sessions$season_rain2, levels = c("D", "W"))
if ("wet_dry2" %in% names(sessions)) calendar_covariates$wet_dry2 <- factor(sessions$wet_dry2, levels = c("D", "W"))
model_parameters <- vector("list", length(classes)) # Translate the saved formulas without assuming a species.
for (class_index in seq_along(classes)) { # Build the same four-parameter candidate specification for each class.
    candidate_parameters <- setNames(vector("list", length(parameter_names)), parameter_names)
    for (parameter_name in parameter_names) { # Only known calendar factors, time and sex may define these models.
        parameter_formula <- as.formula(model_definitions[class_index, parameter_name]) # Preserve the exported formula.
        if (length(parameter_formula) != 2L) stop("Candidate formulas must be one-sided formulas beginning with ~.")
        allowed_variables <- c("sex", "time", names(calendar_covariates)) # Exclude simulated truth and image labels.
        if (any(!all.vars(parameter_formula) %in% allowed_variables)) stop("Unsupported variable in ", classes[class_index], ": ", parameter_name)
        if (parameter_name == "N" && any(!all.vars(parameter_formula) %in% c("sex"))) stop("Superpopulation size may use ~1 or sex, but not occasion covariates.")
        candidate_parameters[[parameter_name]] <- list(formula = parameter_formula) # Let RMark use its POPAN links.
    }
    model_parameters[[class_index]] <- candidate_parameters # Retain candidate order p1 to p5 for reporting ties.
}
names(model_parameters) <- classes # The generating class is used only to assess classification after fitting.
if (length(calendar_covariates) && any(vapply(calendar_covariates, anyNA, logical(1L)))) stop("Calendar factors contain missing or unrecognized values.")
fixed_zero_sessions <- integer(0L) # Akodon and Oligoryzomys have no fixed-zero capture occasions in this metadata.
if ("p_fixed_zero" %in% names(sessions)) { # Phyllotis records the same fixed-zero occasions used during generation.
    fixed_zero_flags <- tolower(as.character(sessions$p_fixed_zero)) # Accept exported TRUE/FALSE or one/zero flags.
    if (anyNA(fixed_zero_flags) || any(!fixed_zero_flags %in% c("true", "false", "1", "0"))) stop("p_fixed_zero must contain TRUE/FALSE or 1/0.")
    fixed_zero_sessions <- which(fixed_zero_flags %in% c("true", "1")) # Fix capture at zero for both sex groups.
}

# 3. Select validation pairs from one CNN run, or sample generated histories without requiring a CNN.
cnn_predictions <- NULL # This remains empty when evaluating a standalone sample of generated capture histories.
eligible_pairs <- pair_index # The standalone evaluation can sample from all classes in the generator index.
if (evaluation_source == "cnn_validation") { # Restrict the comparison to the held-out images from one trained network.
    if (!nzchar(cnn_run_directory)) { # Find the latest completed prediction export rather than a partially trained run.
        cnn_parent <- file.path(dataset_directory, "cnn_results") # Generic CNN training stores each run separately.
        run_candidates <- unique(c(cnn_parent, list.dirs(cnn_parent, recursive = FALSE, full.names = TRUE)))
        completed_runs <- run_candidates[file.exists(file.path(run_candidates, "split_manifest.csv")) & file.exists(file.path(run_candidates, "validation_predictions.csv"))]
        if (!length(completed_runs)) stop("No CNN validation outputs were found. Set cnn_run_directory or use evaluation_source = 'generated_sample'.")
        run_times <- file.info(file.path(completed_runs, "validation_predictions.csv"))$mtime # Order completed exports.
        cnn_run_directory <- completed_runs[order(run_times, decreasing = TRUE)[1L]] # Select one complete CNN run.
    }
    cnn_run_directory <- normalizePath(cnn_run_directory, winslash = "/", mustWork = TRUE)
    cnn_manifest <- read.csv(file.path(cnn_run_directory, "split_manifest.csv"), colClasses = "character", stringsAsFactors = FALSE)
    cnn_predictions <- read.csv(file.path(cnn_run_directory, "validation_predictions.csv"), stringsAsFactors = FALSE)
    if (!all(c("image_file", "class", "subset") %in% names(cnn_manifest))) stop("The CNN manifest must contain image_file, class and subset.")
    if (!all(c("image_file", "true_class", "predicted_class") %in% names(cnn_predictions))) stop("The CNN predictions must contain image_file, true_class and predicted_class.")
    cnn_manifest$image_file <- normalizePath(cnn_manifest$image_file, winslash = "/", mustWork = TRUE)
    cnn_predictions$image_file <- normalizePath(cnn_predictions$image_file, winslash = "/", mustWork = TRUE)
    if (anyDuplicated(cnn_manifest$image_file) || anyDuplicated(cnn_predictions$image_file)) stop("CNN outputs contain repeated image paths.")
    validation_rows <- cnn_manifest[cnn_manifest$subset == "validation", , drop = FALSE] # Exclude every training image.
    index_positions <- match(validation_rows$image_file, pair_index$image_file) # Join the exact image/history pairs.
    prediction_positions <- match(validation_rows$image_file, cnn_predictions$image_file) # Keep the same CNN run.
    if (anyNA(index_positions) || anyNA(prediction_positions)) stop("Every CNN validation image must match the generator index and the prediction export.")
    eligible_pairs <- pair_index[index_positions, , drop = FALSE] # Keep only the paired validation histories.
    cnn_predictions <- cnn_predictions[prediction_positions, , drop = FALSE] # Align predictions with the same pairs.
    if (anyNA(validation_rows$class) || anyNA(cnn_predictions$true_class) || anyNA(cnn_predictions$predicted_class) || any(validation_rows$class != eligible_pairs$scenario) || any(cnn_predictions$true_class != eligible_pairs$scenario) || any(!cnn_predictions$predicted_class %in% classes)) stop("CNN labels disagree with the generator index or contain unsupported classes.")
}
set.seed(as.integer(evaluation_seed)) # Sample within classes independently of the already completed CNN training.
selected_rows <- vector("list", length(classes)) # Preserve balanced sampling where sufficient histories are available.
for (class_index in seq_along(classes)) { # Select up to the configured number from each eligible generating class.
    class_positions <- which(eligible_pairs$scenario == classes[class_index]) # The class label is sampling metadata.
    if (!length(class_positions)) stop("No eligible histories were found for class ", classes[class_index])
    sample_size <- if (histories_per_class == 0L) length(class_positions) else min(length(class_positions), histories_per_class)
    selected_positions <- class_positions[sample.int(length(class_positions), size = sample_size)]
    class_selection <- eligible_pairs[selected_positions, , drop = FALSE] # Preserve each image and history association.
    class_selection$within_class_order <- seq_len(nrow(class_selection)) # Interleave classes during the fitting loop.
    selected_rows[[class_index]] <- class_selection # Keep the chosen replicate metadata without passing it to MARK.
}
evaluation_pairs <- do.call(rbind, selected_rows) # Combine the selected generator rows for one evaluation experiment.
evaluation_pairs <- evaluation_pairs[order(evaluation_pairs$within_class_order, match(evaluation_pairs$scenario, classes)), , drop = FALSE]
evaluation_pairs$evaluation_id <- sprintf("history_%05d", seq_len(nrow(evaluation_pairs))) # Name diagnostic files safely.
if (!is.null(cnn_predictions)) cnn_predictions <- cnn_predictions[match(evaluation_pairs$image_file, cnn_predictions$image_file), , drop = FALSE]

# 4. Check RMark and MARK before fitting, then create a fresh results directory for this evaluation.
if (nzchar(mark_directory)) { # Add the explicitly configured executable directory before loading RMark.
    if (!dir.exists(mark_directory)) stop("The configured mark_directory does not exist: ", mark_directory)
    mark_directory <- normalizePath(mark_directory, winslash = "/", mustWork = TRUE) # Use one executable directory.
    Sys.setenv(PATH = paste(mark_directory, Sys.getenv("PATH"), sep = .Platform$path.sep)) # Linux also uses this PATH.
    assign("MarkPath", paste0(mark_directory, "/"), envir = .GlobalEnv) # RMark uses MarkPath on Windows installations.
}
if (!requireNamespace("RMark", quietly = TRUE)) stop("Install the RMark package with install.packages('RMark').")
if (.Platform$OS.type == "windows") { # RMark may detect a standard MARK installation outside the current PATH.
    detected_mark_path <- get0("MarkPath", envir = .GlobalEnv, ifnotfound = "") # RMark sets this during package loading.
    executable_names <- c("mark.exe", "mark64.exe", "mark32.exe") # Accept the standard Windows MARK executable names.
    path_executables <- unname(Sys.which(executable_names)) # A configured or standard PATH may provide MARK directly.
    search_directories <- unique(c(mark_directory, detected_mark_path, "C:/Program Files/Mark", "C:/Program Files (x86)/Mark"))
    directory_executables <- as.vector(outer(search_directories[nzchar(search_directories)], executable_names, file.path))
    if (!any(nzchar(path_executables)) && !any(file.exists(directory_executables))) stop("Install Program MARK and set mark_directory to its executable directory.")
    if (any(file.exists(directory_executables))) { # Honor a standard installation even when the package is not attached.
        mark_directory <- dirname(directory_executables[which(file.exists(directory_executables))[1L]])
        assign("MarkPath", paste0(mark_directory, "/"), envir = .GlobalEnv) # Retain compatibility with older RMark releases.
        Sys.setenv(PATH = paste(mark_directory, Sys.getenv("PATH"), sep = .Platform$path.sep)) # Make the executable visible.
    }
} else {
    if (!nzchar(Sys.which("mark"))) stop("Install a compatible Program MARK executable named mark and add its directory to PATH.")
}
if (!nzchar(mark_work_directory)) mark_work_directory <- tempdir() # Each history receives its own temporary subfolder.
if (!dir.exists(mark_work_directory)) dir.create(mark_work_directory, recursive = TRUE) # Allow a chosen local parent folder.
mark_work_directory <- normalizePath(mark_work_directory, winslash = "/", mustWork = TRUE) # Resolve temporary storage.
results_parent <- file.path(dataset_directory, "rmark_results") # Keep evaluations separate from CNN training outputs.
dir.create(results_parent, recursive = TRUE, showWarnings = FALSE) # Preserve all previously completed evaluation runs.
results_directory <- file.path(results_parent, format(Sys.time(), "run_%Y%m%d_%H%M%S")) # Name a new evaluation run.
if (dir.exists(results_directory)) results_directory <- tempfile(pattern = "run_", tmpdir = results_parent)
dir.create(results_directory, recursive = TRUE) # Create the new run only after validating all required inputs.
logs_directory <- file.path(results_directory, "fit_logs") # Keep MARK messages and per-history summaries for diagnosis.
dir.create(logs_directory) # The remaining temporary MARK files are removed after each history has been processed.
write.csv(evaluation_pairs, file.path(results_directory, "evaluation_manifest.csv"), row.names = FALSE)
write.csv(model_definitions, file.path(results_directory, "class_definitions_used.csv"), row.names = FALSE)
write.csv(sessions, file.path(results_directory, "sessions_used.csv"), row.names = FALSE) # Record the actual fitted calendar.
configuration <- list(dataset_directory = dataset_directory, evaluation_source = evaluation_source, cnn_run_directory = cnn_run_directory, histories_per_class = histories_per_class, evaluation_seed = evaluation_seed, evaluate_observed = evaluate_observed, require_all_models = require_all_models, adjust_parameter_count = adjust_parameter_count, use_initial_model = use_initial_model, mark_threads = mark_threads, mark_retry = mark_retry, mark_directory = mark_directory, mark_work_directory = mark_work_directory, n_occasions = n_occasions, time_intervals = time_intervals, fixed_zero_sessions = fixed_zero_sessions)
saveRDS(configuration, file.path(results_directory, "evaluation_configuration.rds")) # Save all settings used in this run.
writeLines(capture.output(sessionInfo()), file.path(results_directory, "runtime_information.txt")) # Record RMark's version.
if (!is.null(cnn_predictions)) write.csv(cnn_predictions, file.path(results_directory, "cnn_predictions_selected.csv"), row.names = FALSE)
cat(sprintf("\nEvaluating %d histories across %d classes with %d sessions.\nResults: %s\n", nrow(evaluation_pairs), length(classes), n_occasions, results_directory))

# 5. Fit one history table with the five candidate models, keeping all preparation and fitting steps together.
fit_capture_history <- function(history_file, evaluation_id) { # Reuse this fitting sequence for simulated and observed data.
    history_connection <- if (grepl("\\.gz$", history_file, ignore.case = TRUE)) gzfile(history_file, open = "rt") else file(history_file, open = "rt")
    histories <- tryCatch(read.csv(history_connection, colClasses = "character", stringsAsFactors = FALSE), finally = close(history_connection))
    if (!all(c("ch", "freq", "sex") %in% names(histories)) || !nrow(histories)) stop("The history table must contain nonempty ch, freq and sex columns.")
    if (anyNA(histories$ch) || any(nchar(histories$ch) != n_occasions) || any(!grepl("^[01]+$", histories$ch)) || any(!grepl("1", histories$ch, fixed = TRUE))) stop("Histories must be binary text of the calendar length and contain at least one capture.")
    frequencies <- suppressWarnings(as.numeric(histories$freq)) # Preserve negative observed frequencies for trap deaths.
    if (any(!is.finite(frequencies)) || any(frequencies == 0) || any(frequencies != floor(frequencies))) stop("freq must contain nonzero integer frequencies.")
    if (anyNA(histories$sex) || any(!histories$sex %in% c("F", "M")) || !all(c("F", "M") %in% histories$sex)) stop("Every fit needs female and male records labeled F and M.")
    for (session_number in fixed_zero_sessions) { # Reject histories contradicting the supplied capture-zero calendar.
        if (any(substr(histories$ch, session_number, session_number) != "0")) stop("A capture occurs in a fixed-zero session: ", session_number)
    }
    capture_data <- data.frame(ch = histories$ch, freq = frequencies, sex = factor(histories$sex, levels = c("F", "M")))
    processed_data <- RMark::process.data(capture_data, model = "POPAN", groups = "sex", begin.time = 1, time.intervals = time_intervals)
    design_data <- RMark::make.design.data(processed_data) # Build the correct POPAN time rows for both sex groups.
    for (parameter_name in c("Phi", "p", "pent")) { # Survival uses interval starts and entry uses subsequent occasions.
        design_times <- as.numeric(as.character(design_data[[parameter_name]]$time)) # Time is the actual factor time.
        calendar_positions <- match(round(design_times, 8L), round(occasion_times, 8L)) # Time is not the centered Time column.
        if (anyNA(calendar_positions)) stop("RMark design times do not match the exported monthly calendar.")
        design_data[[parameter_name]]$session <- calendar_positions # Record the matched occasion for fixing capture.
        for (covariate_name in names(calendar_covariates)) { # Copy factors with their original full calendar levels.
            design_data[[parameter_name]][[covariate_name]] <- calendar_covariates[[covariate_name]][calendar_positions]
        }
    }
    if (length(fixed_zero_sessions)) { # Add a fix column only when the calendar actually contains constrained captures.
        design_data$p$fix <- NA_real_ # A missing fix leaves a capture parameter free under its candidate formula.
        design_data$p$fix[design_data$p$session %in% fixed_zero_sessions] <- 0 # Apply the same zeros to all five models.
    }
    previous_directory <- getwd() # Restore the caller's directory even if MARK fails to produce valid output.
    history_work_directory <- tempfile(pattern = paste0(evaluation_id, "_"), tmpdir = mark_work_directory)
    dir.create(history_work_directory) # Fresh temporary folders prevent prompts about reusing previous MARK files.
    on.exit({ # Always restore the original working directory and remove this history's own temporary folder.
        setwd(previous_directory) # Restore the directory used before starting this individual history fit.
        unlink(history_work_directory, recursive = TRUE) # Delete only the temporary directory owned by this fit.
    }, add = TRUE) # Remove only this history's temporary files after its logs have been copied to the results folder.
    setwd(history_work_directory) # Relative MARK prefixes do not include spaces from the parent directory path.
    initial_model <- NULL # Candidate initialization never reads generating parameter values or true entry/exit sessions.
    if (use_initial_model) { # Start with constant probabilities and the same first candidate superpopulation formula.
        initial_parameters <- list(Phi = list(formula = ~1), p = list(formula = ~1), pent = list(formula = ~1), N = model_parameters[[1L]]$N)
        initial_output <- capture.output(initial_result <- tryCatch(RMark::mark(processed_data, ddl = design_data, model.parameters = initial_parameters, adjust = adjust_parameter_count, prefix = "initial_", threads = mark_threads, retry = mark_retry, realvcv = FALSE, delete = TRUE, output = FALSE, silent = TRUE, invisible = TRUE), error = identity))
        if (inherits(initial_result, "error")) initial_output <- c(initial_output, conditionMessage(initial_result))
        if (!inherits(initial_result, "error") && is.list(initial_result) && !is.null(initial_result$results$beta) && all(is.finite(initial_result$results$beta$estimate))) initial_model <- initial_result
        writeLines(initial_output, file.path(logs_directory, paste0(evaluation_id, "_initial.txt"))) # Retain starter diagnostics.
    }
    candidate_rows <- vector("list", length(classes)) # Keep five rows even when a candidate fails to fit.
    for (class_index in seq_along(classes)) { # Fit every candidate independently of the dataset's generating class.
        warning_messages <- character(0L) # Store warnings alongside the criterion instead of treating them as convergence.
        fit_start <- proc.time()[["elapsed"]] # Measure the cost of the current MARK candidate fit.
        fit_output <- capture.output(fitted_model <- tryCatch(withCallingHandlers(RMark::mark(processed_data, ddl = design_data, model.parameters = model_parameters[[class_index]], initial = initial_model, adjust = adjust_parameter_count, prefix = paste0(classes[class_index], "_"), threads = mark_threads, retry = mark_retry, realvcv = FALSE, delete = TRUE, output = FALSE, silent = TRUE, invisible = TRUE), warning = function(warning_condition) {
            warning_messages <<- c(warning_messages, conditionMessage(warning_condition)) # Record backend warnings for review.
            invokeRestart("muffleWarning") # The same warning is written once in the exported candidate diagnostics.
        }), error = identity)) # An unsuccessful candidate must not prevent the remaining models from being evaluated.
        candidate_row <- data.frame(model = classes[class_index], AICc = NA_real_, AICc_unadjusted = NA_real_, npar = NA_real_, npar_unadjusted = NA_real_, effective_sample_size = NA_real_, deviance = NA_real_, singular_count = NA_integer_, elapsed_seconds = proc.time()[["elapsed"]] - fit_start, usable_AICc = FALSE, warning = paste(unique(warning_messages), collapse = " | "), error = "", stringsAsFactors = FALSE)
        if (inherits(fitted_model, "error")) { # Errors retain their text and have no admissible selection criterion.
            candidate_row$error <- conditionMessage(fitted_model) # Preserve the actual backend failure in the results table.
        } else if (!is.list(fitted_model) || is.null(fitted_model$results)) { # Some MARK errors return a character result.
            candidate_row$error <- paste("MARK returned no model results:", paste(as.character(fitted_model), collapse = " "))
        } else {
            result_fields <- c(AICc = "AICc", AICc_unadjusted = "AICc.unadjusted", npar = "npar", npar_unadjusted = "npar.unadjusted", effective_sample_size = "n", deviance = "deviance")
            for (field_name in names(result_fields)) { # Safely extract scalar model diagnostics from each MARK result.
                field_value <- fitted_model$results[[result_fields[[field_name]]]] # An unavailable field remains missing.
                if (is.numeric(field_value) && length(field_value) == 1L) candidate_row[[field_name]] <- field_value
            }
            if (is.na(candidate_row$AICc_unadjusted)) candidate_row$AICc_unadjusted <- candidate_row$AICc # No adjustment needed.
            if (is.na(candidate_row$npar_unadjusted)) candidate_row$npar_unadjusted <- candidate_row$npar # Retain equal counts.
            candidate_row$singular_count <- length(fitted_model$results$singular) # Boundary/nonestimable beta diagnostics.
            candidate_row$usable_AICc <- is.finite(candidate_row$AICc) && is.finite(candidate_row$npar) && is.finite(candidate_row$effective_sample_size) && candidate_row$effective_sample_size > candidate_row$npar + 1
            if (!candidate_row$usable_AICc) candidate_row$error <- "AICc is missing, nonfinite or has insufficient effective sample size."
        }
        candidate_rows[[class_index]] <- candidate_row # Preserve failures as rows instead of silently dropping them.
        writeLines(c(fit_output, warning_messages, candidate_row$error), file.path(logs_directory, paste0(evaluation_id, "_", classes[class_index], ".txt")))
    }
    return(list(candidates = do.call(rbind, candidate_rows), observed_individuals = sum(abs(frequencies)), history_rows = nrow(histories)))
}

# 6. Evaluate the selected histories sequentially and export progress after each complete five-model comparison.
evaluation_tasks <- evaluation_pairs # The simulated tasks retain their true class only for the final assessment.
evaluation_tasks$data_type <- "simulated" # Observed histories have no known generating label for a confusion matrix.
observed_history_file <- file.path(dataset_directory, "observed_capture_histories.csv.gz") # Standard generator output.
if (evaluate_observed && file.exists(observed_history_file)) { # An expanded thirty-session design has no observed task.
    observed_task <- evaluation_tasks[1L, , drop = FALSE] # Use the same table schema without assigning a true model.
    observed_task[1L, ] <- NA # Clear all simulated seed, replicate and population metadata from the observed task.
    observed_task$evaluation_id <- "observed" # Keep observed logs separate from simulated validation diagnostics.
    observed_task$history_file <- normalizePath(observed_history_file, winslash = "/", mustWork = TRUE)
    observed_task$image_file <- file.path(dataset_directory, "observed.png") # Image rendering is not an RMark input.
    observed_task$data_type <- "observed" # Assess empirical histories separately from simulation classification accuracy.
    evaluation_tasks <- rbind(evaluation_tasks, observed_task) # Fit the same five model formulas to the observed data.
}
if (evaluate_observed && !file.exists(observed_history_file)) message("No observed capture-history export is available. Evaluating simulations only.")
summary_rows <- vector("list", nrow(evaluation_tasks)) # Reserve one progress record per simulated or observed task.
model_rows <- vector("list", nrow(evaluation_tasks)) # Preserve the full five-model ranking for every fitted history.
for (task_index in seq_len(nrow(evaluation_tasks))) { # Follow one visible sequence from history input to model selection.
    current_task <- evaluation_tasks[task_index, , drop = FALSE] # Class, seed and replicate remain reporting metadata.
    cat(sprintf("[%d/%d] %s | class %s | history %s\n", task_index, nrow(evaluation_tasks), current_task$evaluation_id, current_task$scenario, basename(current_task$history_file)))
    task_start <- proc.time()[["elapsed"]] # Record total preparation, initialization and candidate-fitting time.
    fit_result <- tryCatch(fit_capture_history(current_task$history_file, current_task$evaluation_id), error = identity)
    preparation_error <- "" # A preparation error is retained and leaves all five candidates without a criterion.
    if (inherits(fit_result, "error")) { # Continue with other histories while keeping the failed task in every denominator.
        preparation_error <- conditionMessage(fit_result) # Report invalid input or backend preparation errors explicitly.
        candidates <- data.frame(model = classes, AICc = NA_real_, AICc_unadjusted = NA_real_, npar = NA_real_, npar_unadjusted = NA_real_, effective_sample_size = NA_real_, deviance = NA_real_, singular_count = NA_integer_, elapsed_seconds = NA_real_, usable_AICc = FALSE, warning = "", error = preparation_error, stringsAsFactors = FALSE)
        observed_individuals <- NA_real_ # An invalid history cannot provide a validated individual count.
        history_rows <- NA_integer_ # Keep the failure distinct from a valid capture history with zero records.
    } else {
        candidates <- fit_result$candidates # Work with the five candidates produced by this history's own fit sequence.
        observed_individuals <- fit_result$observed_individuals # Only captured individuals are part of the input table.
        history_rows <- fit_result$history_rows # Repeated frequency counts can represent more individuals than rows.
    }
    usable_positions <- which(candidates$usable_AICc) # Select models by AICc, not by their generating parameter values.
    complete_candidate_set <- length(usable_positions) == length(classes) # Verify comparable criteria for all five models.
    candidates$delta_AICc <- NA_real_ # An unavailable or invalid criterion cannot receive an Akaike weight.
    candidates$Akaike_weight <- NA_real_ # Partial-set weights are flagged by complete_candidate_set in the export.
    if (length(usable_positions)) { # Calculate rankings among the usable models while retaining missing candidates.
        candidates$delta_AICc[usable_positions] <- candidates$AICc[usable_positions] - min(candidates$AICc[usable_positions])
        relative_support <- exp(-0.5 * candidates$delta_AICc[usable_positions]) # Stabilize weights by subtracting the minimum.
        candidates$Akaike_weight[usable_positions] <- relative_support / sum(relative_support) # Sum to one within usable fits.
    }
    selected_class <- "no_selection" # This category remains in the classification denominator when selection fails.
    selected_AICc <- NA_real_ # Do not report an unavailable candidate criterion as a selected model.
    selected_weight <- NA_real_ # No selection also has no selected-model Akaike support value.
    exact_tie <- FALSE # A numerical tie is resolved by p1 to p5 order and remains explicitly flagged.
    if (length(usable_positions) && (!require_all_models || complete_candidate_set)) { # Default selection needs five criteria.
        minimum_AICc <- min(candidates$AICc[usable_positions]) # Rank candidates only within the current history table.
        tied_positions <- usable_positions[abs(candidates$AICc[usable_positions] - minimum_AICc) <= 1e-8]
        best_position <- tied_positions[1L] # Use deterministic candidate order when the AICc minimum is tied.
        selected_class <- candidates$model[best_position] # The label identifies the model with the lowest AICc.
        selected_AICc <- candidates$AICc[best_position] # Save the selected criterion alongside all model-level diagnostics.
        selected_weight <- candidates$Akaike_weight[best_position] # This weight is conditional on the usable candidate set.
        exact_tie <- length(tied_positions) > 1L # Retain ties even when one deterministic prediction is exported.
    }
    ranking_order <- usable_positions[order(candidates$AICc[usable_positions])] # Keep the smallest two usable criteria.
    second_model_delta <- if (length(ranking_order) >= 2L) candidates$delta_AICc[ranking_order[2L]] else NA_real_
    candidates$evaluation_id <- current_task$evaluation_id # Join each candidate back to its selected image/history pair.
    candidates$data_type <- current_task$data_type # Observed and simulated rows are exported separately below.
    candidates$true_class <- current_task$scenario # Truth labels are added only after all likelihood fits have finished.
    candidates$complete_candidate_set <- complete_candidate_set # Identify weights computed from an incomplete model set.
    candidates$selected <- candidates$model == selected_class # A no_selection result leaves all five flags false.
    model_rows[[task_index]] <- candidates # Retain all criterion values, parameter counts, warnings and errors.
    summary_rows[[task_index]] <- data.frame(evaluation_id = current_task$evaluation_id, data_type = current_task$data_type, true_class = current_task$scenario, replicate = current_task$replicate, image_file = current_task$image_file, history_file = current_task$history_file, selected_class = selected_class, selected_AICc = selected_AICc, selected_weight = selected_weight, second_model_delta = second_model_delta, exact_AICc_tie = exact_tie, usable_models = length(usable_positions), complete_candidate_set = complete_candidate_set, singular_models = sum(candidates$singular_count > 0, na.rm = TRUE), observed_individuals = observed_individuals, history_rows = history_rows, elapsed_seconds = proc.time()[["elapsed"]] - task_start, preparation_error = preparation_error, stringsAsFactors = FALSE)
    completed_summary <- do.call(rbind, summary_rows[seq_len(task_index)]) # Export progress without waiting for the last task.
    completed_models <- do.call(rbind, model_rows[seq_len(task_index)]) # Preserve every completed candidate table.
    write.csv(completed_summary, file.path(results_directory, "selection_summary.csv"), row.names = FALSE, na = "")
    write.csv(completed_models, file.path(results_directory, "model_selection.csv"), row.names = FALSE, na = "")
    saveRDS(list(summary = summary_rows[[task_index]], candidates = candidates), file.path(logs_directory, paste0(current_task$evaluation_id, "_summary.rds")))
}

# 7. Summarize RMark classifications, keeping unavailable selections in the full evaluation denominator.
simulation_results <- completed_summary[completed_summary$data_type == "simulated", , drop = FALSE]
simulation_results$correct <- simulation_results$selected_class == simulation_results$true_class # Failures count as incorrect.
write.csv(simulation_results, file.path(results_directory, "simulation_predictions.csv"), row.names = FALSE, na = "")
methods <- list(rmark = data.frame(true_class = simulation_results$true_class, predicted_class = simulation_results$selected_class))
if (!is.null(cnn_predictions)) { # Calculate CNN accuracy on exactly the histories selected for RMark evaluation.
    paired_comparison <- cnn_predictions # Preserve the original softmax scores along with the paired model predictions.
    paired_comparison$history_file <- simulation_results$history_file # Row order is the same selected evaluation order.
    paired_comparison$rmark_selected_class <- simulation_results$selected_class # Pair predictions by exported image path.
    paired_comparison$rmark_correct <- simulation_results$correct # Include incomplete-fit tasks in the paired table.
    paired_comparison$rmark_complete_candidate_set <- simulation_results$complete_candidate_set # Expose fit completeness.
    paired_comparison$rmark_singular_models <- simulation_results$singular_models # Expose boundary/nonestimable diagnostics.
    paired_comparison$rmark_selected_weight <- simulation_results$selected_weight # AICc support is not a softmax probability.
    paired_comparison$rmark_second_model_delta <- simulation_results$second_model_delta # Record close model alternatives.
    write.csv(paired_comparison, file.path(results_directory, "cnn_rmark_paired_predictions.csv"), row.names = FALSE, na = "")
    methods$cnn <- data.frame(true_class = cnn_predictions$true_class, predicted_class = cnn_predictions$predicted_class)
}
prediction_levels <- c(classes, "no_selection") # Use a sixth column to make unavailable RMark selections visible.
metrics_rows <- vector("list", length(methods) * length(classes)) # Report both full-denominator and selected-only accuracy.
matrix_tables <- list() # Retain matrices for one comparison figure after all numeric exports have been written.
metric_position <- 0L # Advance one metric row for each method and generating class in the evaluation sample.
for (method_name in names(methods)) { # The same denominator and label order apply to both classification methods.
    method_predictions <- methods[[method_name]] # CNN predictions use the same paired simulations as the RMark results.
    confusion_counts <- table(factor(method_predictions$true_class, levels = classes), factor(method_predictions$predicted_class, levels = prediction_levels), dnn = c("true_class", "predicted_class"))
    confusion_proportions <- sweep(confusion_counts, 1L, rowSums(confusion_counts), "/") # Normalize within each true class.
    counts_export <- data.frame(true_class = classes, unclass(confusion_counts), check.names = FALSE) # Preserve labels.
    proportions_export <- data.frame(true_class = classes, unclass(confusion_proportions), check.names = FALSE)
    write.csv(counts_export, file.path(results_directory, paste0(method_name, "_confusion_counts.csv")), row.names = FALSE)
    write.csv(proportions_export, file.path(results_directory, paste0(method_name, "_confusion_proportions.csv")), row.names = FALSE)
    matrix_tables[[method_name]] <- confusion_proportions # Store row-normalized values for the paired comparison figure.
    for (class_name in classes) { # Calculate classification accuracy using all selected histories in this class.
        class_predictions <- method_predictions$predicted_class[method_predictions$true_class == class_name]
        attempted_count <- length(class_predictions) # All attempted histories remain in this evaluation denominator.
        selected_count <- sum(class_predictions != "no_selection") # Report failed selections separately from valid ones.
        correct_count <- sum(class_predictions == class_name) # A correct prediction must match the generating model.
        metric_position <- metric_position + 1L # Place the current method and class in the metrics export.
        metrics_rows[[metric_position]] <- data.frame(method = method_name, class = class_name, histories_attempted = attempted_count, histories_selected = selected_count, no_selection = attempted_count - selected_count, correctly_classified = correct_count, accuracy_all = correct_count / attempted_count, accuracy_selected_only = if (selected_count > 0L) correct_count / selected_count else NA_real_, stringsAsFactors = FALSE)
    }
}
comparison_metrics <- do.call(rbind, metrics_rows) # Compare classification accuracy on the same held-out sample.
write.csv(comparison_metrics, file.path(results_directory, "classification_metrics_by_class.csv"), row.names = FALSE, na = "")
png(file.path(results_directory, "classification_confusion_matrices.png"), width = 1200L * length(methods), height = 1050L, res = 150L)
par(mfrow = c(1L, length(methods)), mar = c(6, 5, 4, 2)) # Draw comparable matrices with the same zero-to-one scale.
for (method_name in names(matrix_tables)) { # Keep class order consistent between the figure and numeric CSV matrices.
    plot_matrix <- matrix_tables[[method_name]][nrow(matrix_tables[[method_name]]):1L, , drop = FALSE] # Put p1 at the top.
    image(seq_along(prediction_levels), seq_along(classes), t(plot_matrix), axes = FALSE, xlab = "Selected model", ylab = "Generating model", main = toupper(method_name), col = colorRampPalette(c("white", "steelblue4"))(100L), zlim = c(0, 1))
    axis(1L, at = seq_along(prediction_levels), labels = prediction_levels, las = 2L) # Include the no-selection column.
    axis(2L, at = seq_along(classes), labels = rev(classes), las = 1L) # Match the reversed plotting row order.
    for (column_index in seq_along(prediction_levels)) { # Display the proportion of histories in each matrix cell.
        text(column_index, seq_along(classes), labels = sprintf("%.2f", plot_matrix[, column_index]), col = ifelse(plot_matrix[, column_index] > 0.55, "white", "black"))
    }
    box() # Keep the boundary of each matrix panel visible without changing the numeric classifications.
}
dev.off() # Complete the comparison figure after all method panels and labels have been drawn.
observed_results <- completed_summary[completed_summary$data_type == "observed", , drop = FALSE] # Separate empirical fits.
if (nrow(observed_results)) { # Observed data have model support values but no known true class or classification accuracy.
    write.csv(observed_results, file.path(results_directory, "observed_selection_summary.csv"), row.names = FALSE, na = "")
    write.csv(completed_models[completed_models$data_type == "observed", , drop = FALSE], file.path(results_directory, "observed_model_selection.csv"), row.names = FALSE, na = "")
}
cat(sprintf("\nCompleted RMark evaluation.\nResults: %s\n", results_directory)) # Report the newly created run directory.
print(comparison_metrics) # Review no-selection counts and full-denominator accuracy before comparing methods.
