# Generate Phyllotis capture-history images and paired individual histories entirely in native R.
# Install these packages once with install.packages(c("png", "jsonlite", "digest")) before running.
# Place PD_data_final.RData and PD.final.results.RData in the working directory, or edit their paths.
# Run with source("generate_phyllotis.R") after adjusting the settings in the first numbered step.
# Native R reproduces its own seeded experiments, while its random sequences differ from NumPy.

required_packages <- c("png", "jsonlite", "digest") # Three packages handle PNG, JSON and source-file checksums.
for (package_name in required_packages) { # Check dependencies before reading the input data or creating output.
    if (!requireNamespace(package_name, quietly = TRUE)) { # Identify missing packages with an actionable message.
        stop(sprintf("Install the missing R package with install.packages('%s').", package_name))
    }
}
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection") # Use fixed native R methods.

# 1. Edit the input paths, image dimensions and number of datasets to generate for each fitted model.
data_file <- "PD_data_final.RData" # Observed histories, monthly covariates and fixed-zero capture occasions.
results_file <- "PD.final.results.RData" # Saved RMark results with the fitted p1 to p5 POPAN models and model.table.
profiles_file <- "PD_top5_fitted_profiles.json" # Read a compatible existing cache, or create it from the fitted results.
output_directory <- "phyllotis_fitted_370_R" # Keep the native R experiment separate from datasets generated in Python.
images_per_class <- 3L # Set this to 5000 for a full run.
master_seed <- 20260921L # Derive a distinct integer seed for each demographic class and replicate in this experiment.
image_width <- 512L # Use at least 370 pixels so every monthly occasion retains a distinct horizontal position.
image_height <- 512L # Compress the ordered individuals into the full image height using event counts per pixel.
population_scale <- 1.0 # Multiply the fitted female and male superpopulation estimates before rounding to integers.
density_count_clip <- 8L # Saturate each color channel when eight or more same-sex events occupy one pixel.
density_min_intensity <- 105L # Set the lower positive intensity used by the original capture-density renderer.
report_every <- 25L # Print progress at this replicate interval and after completing each demographic class.
classes <- c("p1", "p2", "p3", "p4", "p5") # Preserve the five fitted demographic classes and their original labels.
simulation_version <- "sequential_phyllotis_fitted_paired_370_R_v1" # Identify the native R paired-export version.
profile_schema_version <- "phyllotis_rmark_fitted_profiles_v1" # Keep the profile cache compatible with the Python version.
class_rows <- data.frame(scenario = classes, rank = seq_along(classes), Phi = c("~yr*season_rain2", "~yr*season2", "~yr+season_rain2", "~yr+season2", "~yr+wet_dry2"), p = rep("~time", 5L), pent = rep("~yr*season_rain2", 5L), N = rep("~sex", 5L), stringsAsFactors = FALSE)
model_definitions <- setNames(vector("list", length(classes)), classes) # Store each class definition for the JSON records.
for (class_index in seq_along(classes)) { # Associate each named class with its fitted parameter formulas and rank.
    model_definitions[[classes[class_index]]] <- as.list(class_rows[class_index, c("rank", "Phi", "p", "pent", "N"), drop = FALSE])
}

# These two helpers perform only the repeated simulation and image-drawing tasks used by the generation loop.
simulate_scenario <- function(scenario, replicate, replicate_seed) { # Simulate from the fitted point-estimate profiles.
    set.seed(replicate_seed) # Reproduce this particular image/history pair independently of previous replicates.
    profile <- profiles[[scenario]] # Use the fitted model assigned to the current generating demographic class.
    n_females <- max(1L, as.integer(round(profile$Nstar_estimate$F * population_scale))) # Scale fitted female N.
    n_males <- max(1L, as.integer(round(profile$Nstar_estimate$M * population_scale))) # Scale fitted male N separately.
    nstar <- n_females + n_males # Add both fitted sex-specific superpopulation sizes after rounding and scaling.
    sex <- sample(c(rep("F", n_females), rep("M", n_males)), size = nstar, replace = FALSE) # Shuffle fixed sex counts.
    capture_probability <- profile$p # Use the fitted monthly capture probabilities, including the fixed-zero sessions.
    entry_probability <- profile$entry # Include initial entry and the subsequent multinomial pent probabilities.
    survival_probability <- profile$Phi # Use the fitted interval survival probabilities without sampling coefficients.
    entry_occasion <- sample.int(n_occasions, size = nstar, replace = TRUE, prob = entry_probability)
    cumulative_hazard <- c(0, cumsum(-log(pmin(pmax(survival_probability, 1e-12), 1)))) # Accumulate interval mortality.
    death_threshold <- cumulative_hazard[entry_occasion] + rexp(nstar, rate = 1) # Draw one survival threshold per individual.
    exit_boundary <- pmin(pmax(findInterval(death_threshold, cumulative_hazard), entry_occasion), n_occasions)
    captured <- matrix(0L, nrow = nstar, ncol = n_occasions) # Record captures only while individuals are available.
    for (occasion in seq_len(n_occasions)) { # Entry is one-based, while exit_boundary counts available monthly positions.
        present_indices <- which(entry_occasion <= occasion & exit_boundary >= occasion) # Select available individuals.
        if (length(present_indices) > 0L && capture_probability[occasion] > 0) { # Fixed-zero months remain uncaptured.
            captured[present_indices, occasion] <- as.integer(runif(length(present_indices)) < capture_probability[occasion])
        }
    }
    included <- rowSums(captured) > 0L # Retain only individuals with at least one capture in the simulated dataset.
    if (!any(included)) { # Stop rather than exporting an empty image or an unmatched individual capture-history table.
        stop(sprintf("Replicate %s/%d produced no captured individuals: review the fitted profiles.", scenario, replicate))
    }
    histories <- captured[included, , drop = FALSE] # Preserve all 370 occasions in the retained individual histories.
    observed_sex <- sex[included] # Keep each captured individual's sex aligned with its original capture-history row.
    observed_entry <- entry_occasion[included] # True entry may occur before an individual's first recorded capture.
    observed_exit <- exit_boundary[included] + 1L # Export the first unavailable session, with 371 meaning after the study.
    individual_ids <- sprintf("SIM%06d", seq_len(nrow(histories))) # Assign identifiers to the exported observed individuals.
    capture_counts <- rowSums(histories) # Count all capture events, including recaptures, for each retained individual.
    definition <- model_definitions[[scenario]] # Keep the model rank and all four generating formulas in the metadata.
    metadata <- list(simulation_version = simulation_version, scenario = scenario, model_rank = definition$rank, Phi_formula = definition$Phi, p_formula = definition$p, pent_formula = definition$pent, N_formula = definition$N)
    metadata$AICc_RMark <- if (is.null(profile$model_selection$AICc)) NA_real_ else profile$model_selection$AICc
    metadata$DeltaAICc_RMark <- if (is.null(profile$model_selection$DeltaAICc)) NA_real_ else profile$model_selection$DeltaAICc
    metadata$AICc_weight_RMark <- if (is.null(profile$model_selection$weight)) NA_real_ else profile$model_selection$weight
    metadata$replicate <- replicate # Record the replicate number used to name this matched image/history pair.
    metadata$seed <- replicate_seed # Record the native R seed so this same pair can be independently regenerated.
    metadata$n_sessions <- n_occasions # Keep the number of monthly sampling occasions explicit for downstream analysis.
    metadata$Nstar_total <- nstar # Record the whole simulated superpopulation before removing uncaptured individuals.
    metadata$Nstar_F <- n_females # Preserve the exact number of females used in the simulated superpopulation.
    metadata$Nstar_M <- n_males # Preserve the exact number of males used in the simulated superpopulation.
    metadata$fitted_Nstar_F <- profile$Nstar_estimate$F # Record the female estimate before population scaling and rounding.
    metadata$fitted_Nstar_M <- profile$Nstar_estimate$M # Record the male estimate before population scaling and rounding.
    metadata$population_scale <- population_scale # Preserve the multiplier applied to both fitted sex-specific estimates.
    metadata$n_obs <- nrow(histories) # Count captured individuals represented in the image and compressed history table.
    metadata$obs_F <- sum(observed_sex == "F") # Count captured females rather than the entire female superpopulation.
    metadata$obs_M <- sum(observed_sex == "M") # Count captured males rather than the entire male superpopulation.
    metadata$n_captures_total <- sum(capture_counts) # Include all initial captures and recaptures in the event total.
    metadata$n_recaptured_individuals <- sum(capture_counts > 1L) # Count individuals observed on multiple occasions.
    metadata$mean_captures_per_observed_individual <- mean(capture_counts) # Summarize capture frequency among observed rows.
    metadata$mean_residence_months_Nstar <- mean(exit_boundary - entry_occasion + 1L) # Include the entry occasion itself.
    metadata$p_min_fitted_nonzero <- min(capture_probability[capture_probability > 0]) # Exclude fixed-zero capture months.
    metadata$p_mean_fitted <- mean(capture_probability) # Include all calendar months in the mean capture probability.
    metadata$p_max_fitted <- max(capture_probability) # Preserve the largest fitted monthly capture probability.
    metadata$phi_min_fitted <- min(survival_probability) # Summarize the minimum fitted interval survival probability.
    metadata$phi_mean_fitted <- mean(survival_probability) # Summarize the mean fitted interval survival probability.
    metadata$phi_max_fitted <- max(survival_probability) # Summarize the maximum fitted interval survival probability.
    metadata$entry_probability_min_fitted <- min(entry_probability) # Summarize the smallest monthly entry probability.
    metadata$entry_probability_max_fitted <- max(entry_probability) # Summarize the largest monthly entry probability.
    metadata$entry_entropy_fitted <- -sum(entry_probability * log(pmax(entry_probability, 1e-300))) # Include initial entry.
    metadata$parameter_mode <- "fitted_point_estimate_parametric_bootstrap" # Sample individual histories at fixed parameters.
    metadata$fitted_results_sha256 <- fitted_results_sha256 # Link the generated pair to the source fitted RMark results.
    return(list(histories = histories, sex = observed_sex, individual_ids = individual_ids, entry = observed_entry, exit = observed_exit, metadata = metadata))
}

render_capture_history <- function(histories, sex, destination) { # Convert ordered individual histories into an RGB raster.
    if (!is.matrix(histories) || nrow(histories) < 1L || ncol(histories) != 370L || length(sex) != nrow(histories) || anyNA(histories) || any(!histories %in% c(0L, 1L)) || any(rowSums(histories) == 0L) || anyNA(sex) || any(!sex %in% c("F", "M"))) {
        stop("Rendering requires a nonempty binary N-by-370 capture matrix with one F/M sex value per row.")
    }
    first_capture <- max.col(histories, ties.method = "first") # Identify first capture before reordering the individuals.
    last_capture <- ncol(histories) + 1L - max.col(histories[, ncol(histories):1L, drop = FALSE], ties.method = "first")
    total_captures <- rowSums(histories) # Use capture frequency as the second sorting criterion for individual histories.
    row_order <- order(first_capture, -total_captures, last_capture, seq_len(nrow(histories))) # Resolve ties stably.
    histories <- histories[row_order, , drop = FALSE] # Sort a local copy without altering the paired history export.
    sex <- sex[row_order] # Apply the same ordering to the sex codes that determine each color-channel contribution.
    row_pixel <- image_height - 1L - pmin(floor((seq_len(nrow(histories)) - 1L) * image_height / nrow(histories)), image_height - 1L)
    column_pixel <- pmin(floor((seq_len(ncol(histories)) - 0.5) * image_width / ncol(histories)), image_width - 1L)
    capture_positions <- which(histories > 0L, arr.ind = TRUE) # Locate all capture events before counting pixel collisions.
    pixel_index <- row_pixel[capture_positions[, 1L]] * image_width + column_pixel[capture_positions[, 2L]] + 1L
    female_events <- sex[capture_positions[, 1L]] == "F" # Separate female and male contributions to each image pixel.
    female_counts <- matrix(tabulate(pixel_index[female_events], nbins = image_height * image_width), nrow = image_height, ncol = image_width, byrow = TRUE)
    male_counts <- matrix(tabulate(pixel_index[!female_events], nbins = image_height * image_width), nrow = image_height, ncol = image_width, byrow = TRUE)
    image <- array(0, dim = c(image_height, image_width, 3L)) # Initialize a black background with an unused green channel.
    for (channel in c(1L, 3L)) { # Female captures contribute to red, male captures to blue, and overlaps combine both.
        counts <- if (channel == 1L) female_counts else male_counts # Select the event counts for this color channel.
        intensity <- matrix(0, nrow = image_height, ncol = image_width) # Keep pixels without events fully black.
        positive <- counts > 0L # Apply intensity scaling only where this sex has one or more capture events.
        intensity[positive] <- density_min_intensity + (255 - density_min_intensity) * log1p(pmin(counts[positive], density_count_clip)) / log1p(density_count_clip)
        image[, , channel] <- pmin(pmax(round(intensity), 0), 255) / 255 # Match the original integer RGB channel values.
    }
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE) # Create the destination class directory.
    png::writePNG(image, target = destination) # Save the native R pixel array directly, without graphics antialiasing.
}

# 2. Read the RData directly in R, preserve text histories and reconstruct the 370-month observation calendar.
integer_settings <- c(images_per_class, master_seed, image_width, image_height, density_count_clip, density_min_intensity, report_every)
if (any(!is.finite(integer_settings)) || any(integer_settings != floor(integer_settings)) || images_per_class < 1L || images_per_class >= 1000000L || master_seed < 0L || master_seed + 5000000 + images_per_class > .Machine$integer.max || image_width < 370L || image_height < 128L || density_count_clip < 1L || density_min_intensity < 0L || density_min_intensity > 255L || report_every < 1L || !is.finite(population_scale) || population_scale <= 0) {
    stop("Use positive integer image/settings counts, fewer than 1,000,000 images per class, a valid seed, and positive population_scale.")
}
source_data <- normalizePath(path.expand(data_file), winslash = "/", mustWork = TRUE) # Resolve the observed-data input file.
data_environment <- new.env(parent = baseenv()) # Load the RData into its own environment without overwriting user objects.
loaded_names <- load(source_data, envir = data_environment) # Native R reads the original saved objects without conversion.
if (!all(c("PD.process", "PD.ddl", "PD.fixed", "cov2") %in% loaded_names)) { # Require the four original input objects.
    stop("PD_data_final.RData must contain PD.process, PD.ddl, PD.fixed and cov2 for this generator.")
}
process <- data_environment[["PD.process"]] # Use the original processed POPAN data and monthly interval information.
design_data <- data_environment[["PD.ddl"]] # Identify the parameter indexes associated with fixed-zero capture occasions.
fixed_parameters <- data_environment[["PD.fixed"]] # Use the capture probabilities fixed to zero in the original analysis.
covariates <- as.data.frame(data_environment[["cov2"]], stringsAsFactors = FALSE) # Retain the original covariate calendar.
n_occasions <- as.integer(process$nocc) # Require the same complete monthly observation design as the Python generator.
interval_months <- as.numeric(process$time.intervals) # Read the supplied survival intervals rather than inventing a calendar.
if (n_occasions != 370L || length(interval_months) != 369L || any(!is.finite(interval_months)) || any(abs(interval_months - 1) > 1e-8)) {
    stop("This generator requires 370 capture occasions and 369 consecutive intervals of one month each.")
}
observed <- as.data.frame(process$data, stringsAsFactors = FALSE) # Preserve all original individual histories and flags.
if (nrow(observed) < 1L || !all(c("id2", "ch", "sex", "freq", "trap_death") %in% names(observed))) {
    stop("PD.process$data requires id2, ch, sex, freq and trap_death, with nonempty individual capture histories.")
}
observed$id2 <- trimws(as.character(observed$id2)) # Keep original identifiers as text, including any leading zeros.
observed$ch <- trimws(as.character(observed$ch)) # Preserve binary history strings before converting them to a matrix.
observed$sex <- toupper(trimws(as.character(observed$sex))) # Standardize female and male codes to the exported F/M values.
observed$freq <- as.numeric(as.character(observed$freq)) # Preserve signed frequencies from the observed RMark input.
observed$trap_death <- as.numeric(as.character(observed$trap_death)) # Preserve the original trap-death indicator values.
if (anyNA(observed$id2) || any(observed$id2 == "") || anyDuplicated(observed$id2) || anyNA(observed$sex) || any(!observed$sex %in% c("F", "M")) || any(!is.finite(observed$freq)) || any(observed$freq != floor(observed$freq)) || anyNA(observed$trap_death) || any(!observed$trap_death %in% c(0, 1))) {
    stop("Observed data require unique IDs, F/M sex, finite integer frequencies and binary trap_death values.")
}
if (anyNA(observed$ch) || any(!grepl("^[01]{370}$", observed$ch, perl = TRUE)) || any(!grepl("1", observed$ch, fixed = TRUE))) {
    stop("Every observed capture history must be a nonzero binary text string of exactly 370 characters.")
}
observed_histories <- matrix(as.integer(unlist(strsplit(observed$ch, split = "", fixed = TRUE), use.names = FALSE)), nrow = nrow(observed), ncol = n_occasions, byrow = TRUE)
if (nrow(covariates) != n_occasions || !all(c("yr.mo", "time", "season2", "season_rain2", "rain_LMH2", "wet_dry2") %in% names(covariates))) {
    stop("cov2 must contain 370 rows and the columns yr.mo, time, season2, season_rain2, rain_LMH2 and wet_dry2.")
}
year_month <- as.numeric(as.character(covariates$yr.mo)) # Interpret original year.month values without reading decimal months as dates.
years <- floor(year_month + 1e-9) # Extract the integer year used by the fitted year-factor coefficients.
months <- round((year_month - years) * 100) # Recover the integer month from the two digits after the decimal point.
session_number <- as.numeric(as.character(covariates$time)) # Read factor labels rather than internal factor integer codes.
if (any(!is.finite(year_month)) || any(months < 1L | months > 12L) || anyNA(session_number) || !identical(session_number, as.numeric(seq_len(n_occasions)))) {
    stop("Use valid year.month covariates and consecutive cov2$time labels numbered 1 through 370.")
}
fixed_rows <- design_data$p[design_data$p$par.index %in% fixed_parameters$index, , drop = FALSE] # Match fixed p indexes.
fixed_zero_sessions <- sort(unique(as.integer(as.character(fixed_rows$time)))) # Keep each fixed-zero session once.
if (length(fixed_zero_sessions) == 0L || anyNA(fixed_zero_sessions) || any(fixed_zero_sessions < 1L | fixed_zero_sessions > n_occasions) || any(!is.finite(fixed_parameters$value)) || any(fixed_parameters$value != 0)) {
    stop("PD.fixed and PD.ddl$p must identify capture occasions whose probabilities are fixed to zero.")
}
sessions <- data.frame(session = session_number, date = sprintf("%04d-%02d-01", as.integer(years), as.integer(months)), year = years, month = months, year_index = match(years, sort(unique(years))) - 1L, season2 = as.character(covariates$season2), season_reproductive = ifelse(as.character(covariates$season2) == "B", "breeding", "nonbreeding"), season_rain2 = as.character(covariates$season_rain2), season_rain = ifelse(as.character(covariates$season_rain2) == "W", "wet", "dry"), rain_LMH2 = as.character(covariates$rain_LMH2), wet_dry2 = as.character(covariates$wet_dry2), annual_rain_class = ifelse(as.character(covariates$wet_dry2) == "W", "wet_year", "dry_year"), p_fixed_zero = session_number %in% fixed_zero_sessions, time_interval_months = c(NA_real_, interval_months), stringsAsFactors = FALSE)
if (anyNA(sessions) && anyNA(sessions[, setdiff(names(sessions), "time_interval_months"), drop = FALSE])) {
    stop("The reconstructed session calendar contains missing covariate values that prevent fitted-profile reconstruction.")
}
if (!setequal(sessions$season2, c("B", "N")) || !setequal(sessions$season_rain2, c("D", "W")) || !setequal(sessions$wet_dry2, c("D", "W"))) {
    stop("The fitted covariates require season2=B/N, season_rain2=D/W and wet_dry2=D/W factor labels.")
}
zero_capture_sessions <- which(colSums(observed_histories) == 0L) # Identify the zero-capture months in the observations.
if (!identical(as.integer(zero_capture_sessions), fixed_zero_sessions)) { # Preserve the original capture restrictions.
    stop("Zero-capture sessions in the observed histories do not match the fixed-zero occasions in PD.fixed.")
}
capture_counts_observed <- rowSums(observed_histories) # Keep original binary capture-event totals for the observed summary.
observed_summary <- list(n_occasions = n_occasions, first_month = sessions$date[1L], last_month = sessions$date[n_occasions], n_individuals = nrow(observed), n_females = sum(observed$sex == "F"), n_males = sum(observed$sex == "M"), male_fraction = mean(observed$sex == "M"), total_binary_capture_events = sum(capture_counts_observed), n_recaptured_individuals = sum(capture_counts_observed > 1L), mean_captures_per_individual = mean(capture_counts_observed), n_trap_death_rows = sum(observed$trap_death == 1L), n_negative_frequency_rows = sum(observed$freq < 0L), p_fixed_zero_sessions = fixed_zero_sessions)
source_data_sha256 <- digest::digest(file = source_data, algo = "sha256") # Identify the exact original input RData bytes.
rm(data_environment, process, design_data, fixed_parameters, covariates) # Release loaded design objects after extracting inputs.

# 3. Read a compatible profile cache, or reconstruct and validate probabilities from the saved RMark coefficients.
source_results <- path.expand(results_file) # The fitted RData is optional only when a validated profile cache already exists.
cache_file <- path.expand(profiles_file) # Use the same cache schema as the current Python generator for interoperability.
source_results_sha256 <- if (file.exists(source_results)) digest::digest(file = source_results, algo = "sha256") else NULL
payload <- NULL # Load cached point-estimate profiles only when they correspond to the current fitted-results file.
if (file.exists(cache_file)) { # Reuse the cache without invoking Python, RMark or the external MARK executable.
    candidate <- jsonlite::read_json(cache_file, simplifyVector = TRUE) # Recover probability vectors from the saved JSON.
    if (is.null(source_results_sha256) || identical(candidate$source_results_sha256, source_results_sha256)) {
        payload <- candidate # Use a matching cache, then validate its data identity, model definitions and probabilities.
        cat(sprintf("Fitted profiles loaded from cache: %s\n", cache_file)) # Report the selected source of model profiles.
    }
}
if (is.null(payload)) { # Extract point estimates directly from the saved native R objects when no matching cache is available.
    if (!file.exists(source_results)) { # Require either the fitted-results RData or a valid previously generated profile cache.
        stop("Provide PD.final.results.RData or a compatible PD_top5_fitted_profiles.json, or edit the corresponding paths.")
    }
    fitted_environment <- new.env(parent = baseenv()) # Isolate the fitted model objects from the user's current R workspace.
    fitted_names <- load(source_results, envir = fitted_environment) # Read saved results without refitting any model.
    if (!"PD.POPAN.final.results" %in% fitted_names) { # Require the saved fitted-model collection used by the original generator.
        stop("The fitted-results RData must contain the object PD.POPAN.final.results with p1 to p5 and model.table.")
    }
    fitted <- fitted_environment[["PD.POPAN.final.results"]] # Select the original fitted POPAN model collection.
    if (!all(c(classes, "model.table") %in% names(fitted))) { # Preserve the five fitted classes and their model-selection table.
        stop("PD.POPAN.final.results must contain fitted models p1 through p5 and the corresponding model.table.")
    }
    model_table <- as.data.frame(fitted$model.table, stringsAsFactors = FALSE) # Retain original AICc and model-weight values.
    payload <- list(schema_version = profile_schema_version, created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), source_results_file = basename(source_results), source_results_sha256 = source_results_sha256, source_data_file = basename(source_data), source_data_sha256 = source_data_sha256, n_occasions = n_occasions, parameter_mode = "RMark point estimates and conditional parametric bootstrap", uses_beta_vcv = FALSE, beta_vcv_reason = "p5 reported singular beta and numerical underflow in the variance-covariance matrix", models = list())
    for (scenario in classes) { # Reconstruct each fitted class once, before simulating any individual capture histories.
        definition <- model_definitions[[scenario]] # Retrieve the generating survival, capture, entry and N formulas.
        fitted_result <- fitted[[scenario]]$results # Read coefficients, real parameters and estimation diagnostics.
        beta <- fitted_result$beta # Coefficient row names include the parameter prefix followed by the factor term.
        real <- fitted_result$real # Use MARK's exported real estimates to validate every reconstructed probability.
        probability_vectors <- list() # Store the interval Phi and multinomial-entry linear predictors for this class.
        for (parameter_name in c("Phi", "pent")) { # Survival uses leading sessions, while entry uses sessions after the first.
            profile_rows <- if (parameter_name == "Phi") sessions[-n_occasions, , drop = FALSE] else sessions[-1L, , drop = FALSE]
            formula_text <- definition[[parameter_name]] # The class determines whether year and season interact or add.
            coefficient_rows <- startsWith(rownames(beta), paste0(parameter_name, ":")) # Select this parameter's coefficients.
            coefficients <- setNames(as.numeric(beta$estimate[coefficient_rows]), substring(rownames(beta)[coefficient_rows], nchar(parameter_name) + 2L))
            if (!"(Intercept)" %in% names(coefficients) || any(!is.finite(coefficients))) { # Require finite fitted predictors.
                stop(sprintf("The fitted %s coefficients for %s lack a finite intercept or contain nonfinite values.", parameter_name, scenario))
            }
            year_effect <- coefficients[paste0("yr", profile_rows$year)] # The reference year has no separate fitted coefficient.
            year_effect[is.na(year_effect)] <- 0 # Match the treatment-coded year effects used by the original reconstruction.
            predictor <- unname(coefficients[["(Intercept)"]] + year_effect) # Add the year effects to the common intercept.
            seasonal_factor <- if (grepl("season_rain2", formula_text, fixed = TRUE)) "season_rain2" else if (grepl("season2", formula_text, fixed = TRUE)) "season2" else "wet_dry2"
            nonreference <- if (seasonal_factor == "season2") "N" else "W" # Use the original treatment-coded factor contrasts.
            indicator <- profile_rows[[seasonal_factor]] == nonreference # Identify the nonreference seasonal or annual class.
            seasonal_coefficient <- coefficients[paste0(seasonal_factor, nonreference)] # Read the seasonal main-effect coefficient.
            seasonal_coefficient[is.na(seasonal_coefficient)] <- 0 # Reference or absent coefficients contribute zero.
            predictor <- predictor + unname(seasonal_coefficient) * indicator # Apply the seasonal main effect where appropriate.
            if (grepl("*", formula_text, fixed = TRUE)) { # Add year-by-season interactions only when the fitted formula includes them.
                interaction_effect <- coefficients[paste0("yr", profile_rows$year, ":", seasonal_factor, nonreference)]
                interaction_effect[is.na(interaction_effect)] <- 0 # Reference-year interaction terms are not fitted separately.
                predictor <- predictor + unname(interaction_effect) * indicator # Preserve each fitted year-specific seasonal effect.
            }
            probability_vectors[[parameter_name]] <- predictor # Retain the predictor before applying the parameter-specific link.
        }
        survival_probability <- plogis(pmin(pmax(probability_vectors$Phi, -35), 35)) # Match the original bounded-logit conversion.
        capture_rows <- startsWith(rownames(beta), "p:") # Extract monthly capture coefficients from the fitted time-factor model.
        capture_coefficients <- setNames(as.numeric(beta$estimate[capture_rows]), substring(rownames(beta)[capture_rows], 3L))
        if (!"(Intercept)" %in% names(capture_coefficients) || any(!is.finite(capture_coefficients))) {
            stop(sprintf("The fitted capture coefficients for %s require a finite intercept and finite estimates.", scenario))
        }
        capture_effect <- capture_coefficients[paste0("time", seq_len(n_occasions))] # Use treatment-coded session effects.
        capture_effect[is.na(capture_effect)] <- 0 # The first capture occasion is represented by the fitted intercept alone.
        capture_probability <- plogis(pmin(pmax(unname(capture_coefficients[["(Intercept)"]] + capture_effect), -35), 35))
        capture_probability[fixed_zero_sessions] <- 0 # Preserve the original occasions whose capture probabilities are fixed to zero.
        pent_predictor <- probability_vectors$pent # Entry uses the multinomial-logit link rather than independent logits.
        predictor_maximum <- max(0, pent_predictor) # Include the reference initial-entry predictor when stabilizing exponentiation.
        scaled_pent <- exp(pent_predictor - predictor_maximum) # Compute stable weights for the 369 subsequent entry occasions.
        scaled_first <- exp(-predictor_maximum) # Include the residual initial-entry probability as the multinomial reference.
        entry_probability <- c(scaled_first, scaled_pent) / (scaled_first + sum(scaled_pent)) # Normalize all 370 entry probabilities.
        pent_probability <- entry_probability[-1L] # MARK reports pent for entry after the first sampling occasion only.
        reconstructed <- list(Phi = survival_probability, p = capture_probability, pent = pent_probability)
        validation <- list() # Compare reconstructed probabilities with the original MARK real-parameter estimates.
        for (parameter_name in c("Phi", "p", "pent")) { # Verify all exported real probabilities, including both sex groups.
            real_rows <- startsWith(rownames(real), paste0(parameter_name, " ")) # Select this parameter's real estimates.
            labels <- rownames(real)[real_rows] # Retain the original session labels used by the MARK parameter table.
            if (length(labels) == 0L || any(!grepl(" t[0-9]+$", labels))) { # Reject unsupported real-parameter label formats.
                stop(sprintf("No valid session-labelled MARK real estimates were found for %s in %s.", parameter_name, scenario))
            }
            positions <- as.integer(sub(".* t([0-9]+)$", "\\1", labels)) - as.integer(parameter_name == "pent")
            if (any(positions < 1L | positions > length(reconstructed[[parameter_name]]))) { # Validate all session positions.
                stop(sprintf("MARK real estimates for %s contain out-of-range session labels in %s.", parameter_name, scenario))
            }
            maximum_difference <- max(abs(as.numeric(real$estimate[real_rows]) - reconstructed[[parameter_name]][positions]))
            if (!is.finite(maximum_difference) || maximum_difference > 1e-5) { # Use the original absolute-error tolerance.
                stop(sprintf("Reconstructed %s differs from MARK for %s, maximum absolute difference %.8g.", parameter_name, scenario, maximum_difference))
            }
            validation[[paste0(parameter_name, "_max_absolute_error_vs_MARK_real")]] <- maximum_difference
        }
        population_estimates <- list() # Read the separate female and male real N estimates instead of using new N ranges.
        for (sex_code in c("F", "M")) { # Preserve the sex-specific fitted superpopulation values from the original POPAN analysis.
            population_rows <- startsWith(rownames(real), paste0("N g", tolower(sex_code), " ")) # Match the female or male group.
            if (sum(population_rows) != 1L) { # Require one identifiable fitted superpopulation estimate for each sex.
                stop(sprintf("Cannot identify a unique fitted N estimate for sex %s in %s.", sex_code, scenario))
            }
            population_estimates[[sex_code]] <- as.numeric(real$estimate[population_rows]) # Retain the unrounded fitted N value.
        }
        table_position <- if ("model" %in% names(model_table)) match(scenario, as.character(model_table$model)) else match(scenario, classes)
        if (is.na(table_position)) { # Prevent a class from acquiring the model-selection statistics of another fitted model.
            stop(sprintf("The model-selection table has no matching row for class %s.", scenario))
        }
        selection <- as.list(model_table[table_position, setdiff(names(model_table), c("Phi", "p", "pent", "N", "model")), drop = FALSE])
        diagnostics <- c(list(singular_beta_numbers = as.integer(fitted_result$singular), npar_result = fitted_result$npar, npar_unadjusted = fitted_result$npar.unadjusted), validation)
        payload$models[[scenario]] <- list(definition = definition, model_selection = selection, diagnostics = diagnostics, Nstar_estimate = population_estimates, Phi = survival_probability, p = capture_probability, pent = pent_probability, entry = entry_probability)
        cat(sprintf("Validated fitted profiles for %s against MARK real estimates.\n", scenario)) # Report the completed reconstruction.
    }
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE) # Create the selected profile-cache directory.
    jsonlite::write_json(payload, path = cache_file, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
    rm(fitted_environment, fitted, fitted_result, beta, real, model_table) # Release the large fitted-results objects after extraction.
    invisible(gc()) # Recover memory before simulating many replicate histories and drawing their matched images.
    cat(sprintf("Fitted profiles saved in: %s\n", cache_file)) # The compatible cache can be reused in subsequent native R runs.
}
if (!identical(payload$schema_version, profile_schema_version) || !identical(payload$source_data_sha256, source_data_sha256) || payload$n_occasions != n_occasions || !setequal(names(payload$models), classes)) {
    stop("The profile cache must match this input RData, the expected schema, all five classes and 370 occasions.")
}
profiles <- payload$models # The generation helper uses the same validated fitted vectors stored in the compatible cache.
for (scenario in classes) { # Validate cached values before any observed or simulated image/history pair is written.
    profile <- profiles[[scenario]] # Check every model, regardless of whether it was reconstructed or loaded from JSON.
    if (length(profile$Phi) != 369L || length(profile$p) != 370L || length(profile$pent) != 369L || length(profile$entry) != 370L || any(!is.finite(profile$Phi)) || any(profile$Phi <= 0 | profile$Phi > 1) || any(!is.finite(profile$p)) || any(profile$p < 0 | profile$p > 1) || any(profile$p[fixed_zero_sessions] != 0) || any(!is.finite(profile$pent)) || any(!is.finite(profile$entry)) || any(profile$entry < 0) || abs(sum(profile$entry) - 1) > 1e-10 || any(abs(profile$pent - profile$entry[-1L]) > 1e-12)) {
        stop(sprintf("Invalid fitted probability vectors, fixed-zero captures or normalized entry probabilities for %s.", scenario))
    }
    if (!identical(profile$definition$Phi, model_definitions[[scenario]]$Phi) || !identical(profile$definition$p, "~time") || !identical(profile$definition$pent, "~yr*season_rain2") || !identical(profile$definition$N, "~sex")) {
        stop(sprintf("The profile cache contains different generating model formulas for %s.", scenario))
    }
    population_values <- as.numeric(unlist(profile$Nstar_estimate[c("F", "M")], use.names = FALSE)) # Validate both fitted sexes.
    if (length(population_values) != 2L || any(!is.finite(population_values)) || any(population_values <= 0) || sum(pmax(1, round(population_values * population_scale))) > .Machine$integer.max) {
        stop(sprintf("Fitted female and male N estimates for %s must be positive and fit the native R population range.", scenario))
    }
}
fitted_results_sha256 <- payload$source_results_sha256 # Retain the original fitted-results identity even when using only the cache.

# 4. Record the experiment, calendar and fitted profiles, then export the observed image and original histories.
configuration <- list(data_file = data_file, results_file = results_file, profiles_file = profiles_file, output_directory = output_directory, images_per_class = images_per_class, master_seed = master_seed, image_width = image_width, image_height = image_height, population_scale = population_scale, density_count_clip = density_count_clip, density_min_intensity = density_min_intensity, report_every = report_every)
expected <- list(simulation_version = simulation_version, source_data_sha256 = source_data_sha256, fitted_results_sha256 = fitted_results_sha256, count_per_class = images_per_class, config = configuration, models = model_definitions, rng_kind = RNGkind(), seed_rule = "master_seed + class_index * 1000000 + replicate", R_version = R.version.string)
configuration_signature <- digest::digest(jsonlite::toJSON(expected, auto_unbox = TRUE, digits = NA, null = "null"), algo = "sha256", serialize = FALSE)
root <- path.expand(output_directory) # Use the selected output directory for all image, history and metadata exports.
dir.create(root, recursive = TRUE, showWarnings = FALSE) # Create the experiment directory before recording its absolute path.
root <- normalizePath(root, winslash = "/", mustWork = TRUE) # Display an absolute destination in progress and completion messages.
metadata_file <- file.path(root, "experiment_info.json") # Prevent different input data or settings from sharing the same output.
if (file.exists(metadata_file)) { # Matching reruns regenerate every pair from its seed rather than skipping existing files.
    previous <- jsonlite::read_json(metadata_file, simplifyVector = TRUE) # Check the recorded experiment before replacing outputs.
    if (!identical(previous$configuration_signature, configuration_signature)) { # Preserve outputs from other experiments.
        stop("This output directory belongs to another experiment or changed configuration: choose another output_directory.")
    }
} else if (dir.exists(file.path(root, "images")) || dir.exists(file.path(root, "capture_histories"))) {
    stop("This output directory contains an untracked dataset: choose another output_directory before generating pairs.")
}
information <- c(expected, list(configuration_signature = configuration_signature, source_data_file = source_data, fitted_results_file = source_results, fitted_profiles_cache = cache_file, observed_summary = observed_summary, parameters_from_RMark = TRUE, parameter_statement = "Conditional parametric bootstrap using fixed fitted Phi, p, pent and sex-specific N estimates without sampling coefficient uncertainty.", entry_statement = "The 369 MARK pent estimates and residual initial-entry probability define a multinomial distribution over 370 occasions.", image_statement = "Ordered individual capture histories compressed into an RGB raster with female events in red and male events in blue."))
jsonlite::write_json(information, path = metadata_file, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
jsonlite::write_json(configuration, path = file.path(root, "config.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
jsonlite::write_json(observed_summary, path = file.path(root, "observed_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
write.csv(sessions, file = file.path(root, "sessions.csv"), row.names = FALSE, na = "") # Preserve all calendar covariates.
definition_rows <- vector("list", length(classes)) # Combine the model definitions with the original model-selection statistics.
fitted_rows <- vector("list", length(classes)) # Export one full 370-session fitted-probability table per class.
for (class_index in seq_along(classes)) { # Keep the class order unchanged in the definition and probability-profile exports.
    scenario <- classes[class_index] # Select the fitted class for these experiment-level metadata tables.
    profile <- profiles[[scenario]] # Use the same point estimates that will generate the individual histories below.
    definition_rows[[class_index]] <- as.data.frame(c(list(scenario = scenario), model_definitions[[scenario]], profile$model_selection), stringsAsFactors = FALSE)
    fitted_rows[[class_index]] <- data.frame(scenario = scenario, session = seq_len(n_occasions), date = sessions$date, Phi = c(profile$Phi, NA_real_), p = profile$p, pent = c(NA_real_, profile$pent), entry_probability = profile$entry, Nstar_F = profile$Nstar_estimate$F, Nstar_M = profile$Nstar_estimate$M, stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, definition_rows), file = file.path(root, "class_definitions.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, fitted_rows), file = file.path(root, "fitted_probability_profiles.csv"), row.names = FALSE, na = "")
render_capture_history(observed_histories, observed$sex, file.path(root, "observed.png")) # Draw the original observed histories.
observed_export <- observed[, intersect(c("id2", "freq", "sex", "ch", "grid", "trap_death"), names(observed)), drop = FALSE]
names(observed_export)[names(observed_export) == "id2"] <- "ID_individual" # Match the identifier column of simulated tables.
observed_connection <- gzfile(file.path(root, "observed_capture_histories.csv.gz"), open = "wt", encoding = "UTF-8", compression = 6L)
tryCatch(write.csv(observed_export, file = observed_connection, row.names = FALSE, na = ""), finally = close(observed_connection))
parameter_rows <- vector("list", length(classes) * images_per_class) # Reserve one metadata record per simulated image/history pair.
index_rows <- vector("list", length(classes) * images_per_class) # Reserve one file-mapping record per simulated pair.
pair_position <- 0L # Maintain a single export position across all generating classes and their replicate histories.
cat(sprintf("Observed data: %d individuals, %d sessions, fixed-zero captures: %s\n", nrow(observed), n_occasions, paste(fixed_zero_sessions, collapse = ", ")))

# 5. Simulate each class and replicate sequentially, then save the image and its unaltered individual histories.
# Each compressed history is kept as text and includes true entry and exit sessions for downstream RMark evaluation.
for (class_index in seq_along(classes)) { # Repeat the same transparent generation sequence for each of the five fitted classes.
    scenario <- classes[class_index] # Select the generating class before advancing through its individual replicates.
    cat(sprintf("[%d/5] %s: %d image/history replicates\n", class_index, scenario, images_per_class)) # Report the current class.
    for (replicate in seq_len(images_per_class)) { # One replicate corresponds to a complete simulated study and one image.
        replicate_seed <- as.integer(master_seed + class_index * 1000000 + replicate) # Avoid overflowing native R integer seeds.
        simulation_result <- simulate_scenario(scenario, replicate, replicate_seed) # Call the repeated simulation helper.
        histories <- simulation_result$histories # Use the same individual capture matrix for drawing and history export.
        sex <- simulation_result$sex # Preserve sex alignment in the unsorted capture-history table exported below.
        metadata <- simulation_result$metadata # Retain the fitted-model parameters and realized capture summaries.
        basename <- sprintf("%s_%05d", scenario, replicate) # Give matching image and history files a common basename.
        image_relative <- paste0(scenario, "/", basename, ".png") # Store portable paths relative to the images directory.
        history_relative <- paste0(scenario, "/", basename, ".csv.gz") # Store portable paths relative to capture_histories.
        render_capture_history(histories, sex, file.path(root, "images", image_relative)) # Draw this simulated capture study.
        first_capture <- max.col(histories, ties.method = "first") # Export the first recorded capture of every individual.
        last_capture <- n_occasions + 1L - max.col(histories[, n_occasions:1L, drop = FALSE], ties.method = "first")
        history_table <- data.frame(scenario = scenario, replicate = replicate, seed = replicate_seed, image = image_relative, ID_individual = simulation_result$individual_ids, freq = 1L, sex = sex, ch = apply(histories, MARGIN = 1L, FUN = paste0, collapse = ""), n_captures = rowSums(histories), first_capture = first_capture, last_capture = last_capture, entry_session_true = simulation_result$entry, exit_session_exclusive_true = simulation_result$exit, stringsAsFactors = FALSE)
        history_file <- file.path(root, "capture_histories", history_relative) # Locate the matched compressed history table.
        dir.create(dirname(history_file), recursive = TRUE, showWarnings = FALSE) # Create the history folder for this class.
        temporary_file <- paste0(history_file, ".part") # Complete the gzip table before assigning its final output filename.
        history_connection <- gzfile(temporary_file, open = "wt", encoding = "UTF-8", compression = 6L) # Use gzip level six.
        tryCatch(write.csv(history_table, file = history_connection, row.names = FALSE), finally = close(history_connection))
        if (file.exists(history_file) && !file.remove(history_file)) { # Replace only outputs from this recorded matching experiment.
            stop(sprintf("Could not replace the existing capture-history table: %s", history_file))
        }
        if (!file.rename(temporary_file, history_file)) { # Detect a failed final file move instead of recording an incomplete pair.
            stop(sprintf("Could not finalize the compressed capture-history table: %s", history_file))
        }
        metadata$image_file <- image_relative # Link the simulation parameters to this pair's relative image filename.
        metadata$history_file <- history_relative # Link the same metadata row to the compressed individual histories.
        pair_position <- pair_position + 1L # Advance only after the image and compressed history file have both been written.
        parameter_rows[[pair_position]] <- as.data.frame(metadata, stringsAsFactors = FALSE) # Retain the realized summary row.
        index_rows[[pair_position]] <- data.frame(scenario = scenario, replicate = replicate, seed = replicate_seed, image_file = image_relative, history_file = history_relative, n_obs = metadata$n_obs, n_captures_total = metadata$n_captures_total, n_recaptured_individuals = metadata$n_recaptured_individuals, stringsAsFactors = FALSE)
        if (replicate %% report_every == 0L || replicate == images_per_class) { # Report completed pairs at the requested interval.
            cat(sprintf("%s %d/%d, N*=%d, observed=%d, captures=%d\n", scenario, replicate, images_per_class, metadata$Nstar_total, metadata$n_obs, metadata$n_captures_total))
        }
    }
}

# 6. Export the paired-file index, per-image parameters and per-class summaries after completing generation.
parameters_table <- do.call(rbind, parameter_rows) # Combine the simulation parameters and realized capture summaries.
index_table <- do.call(rbind, index_rows) # Combine the portable image/history paths used by downstream analysis scripts.
write.csv(parameters_table, file = file.path(root, "parameters_per_image.csv"), row.names = FALSE, na = "")
write.csv(index_table, file = file.path(root, "capture_history_index.csv"), row.names = FALSE) # Save the complete file mapping.
summary_rows <- vector("list", length(classes)) # Reserve one summary row for each generating fitted demographic class.
for (class_index in seq_along(classes)) { # Summarize all replicates within each generating class without altering individual histories.
    scenario <- classes[class_index] # Keep the same class order as the fitted definitions and the generation loop.
    selected <- parameters_table[parameters_table$scenario == scenario, , drop = FALSE] # Select this class's replicate metadata.
    summary_rows[[class_index]] <- data.frame(scenario = scenario, n_images = nrow(selected), mean_Nstar = mean(selected$Nstar_total), mean_n_obs = mean(selected$n_obs), mean_captures = mean(selected$n_captures_total), mean_recaptured = mean(selected$n_recaptured_individuals), fitted_mean_phi = mean(selected$phi_mean_fitted), fitted_mean_p = mean(selected$p_mean_fitted), stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, summary_rows), file = file.path(root, "summary_by_class.csv"), row.names = FALSE) # Export class means.
expected_pairs <- length(classes) * images_per_class # Require one matched image and one history table for every replicate.
generated_images <- list.files(file.path(root, "images"), pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
generated_histories <- list.files(file.path(root, "capture_histories"), pattern = "\\.csv\\.gz$", recursive = TRUE, full.names = TRUE)
if (length(generated_images) != expected_pairs || length(generated_histories) != expected_pairs) { # Reject incomplete exports.
    stop(sprintf("Incomplete paired export: images=%d, histories=%d, expected=%d.", length(generated_images), length(generated_histories), expected_pairs))
}
cat(sprintf("Completed %d PNG images and %d history files, directory: %s\n", length(generated_images), length(generated_histories), root))
