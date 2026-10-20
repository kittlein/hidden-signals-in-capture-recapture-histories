# Generate paired rodent capture histories and images entirely in R from the selected workbook sheet.
# Install these packages once in the R console before running this script with source() or Rscript.
# install.packages(c("readxl", "png", "jsonlite", "digest"))
# Each executable statement occupies its own line, and only simulation and drawing use helper functions.

required_packages <- c("readxl", "png", "jsonlite", "digest") # These are the four external R dependencies.
for (package_name in required_packages) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
        stop(sprintf("Install the missing R package with install.packages('%s').", package_name))
    }
}

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection") # Fix the native R random-number methods.

# 1. Select the species and the sampling design, then choose the number of image/history pairs per class.
species <- "Akodon azarae" # Change to "Oligoryzomys flavescens" to select that species and its worksheet.
workbook_file <- "Rodents_Necochea.xlsx" # Place the workbook in the working directory or use its full path.
history_sheet <- species # The selected species name must match its worksheet name in the shared workbook.
calendar_sheet <- "sessions" # Both species use the same ten-session calendar from this workbook worksheet.
if (species == "Akodon azarae") {
    species_prefix <- "akodon" # Use a species-specific prefix to distinguish the generated output directories.
    base_population_min <- 280L # Minimum Akodon superpopulation size before applying the population multiplier.
    base_population_max <- 440L # Maximum Akodon superpopulation size before applying the population multiplier.
    p_base_min <- 0.18 # Lower bound of the baseline Akodon capture probability drawn for each simulated dataset.
    p_base_max <- 0.34 # Upper bound of the baseline Akodon capture probability drawn for each simulated dataset.
    master_seed <- 123L # Reproduce this native R experiment when its settings and random-number methods match.
} else if (species == "Oligoryzomys flavescens") {
    species_prefix <- "oligoryzomys" # Use a separate prefix for the outputs generated from this species worksheet.
    base_population_min <- 180L # Minimum Oligoryzomys superpopulation size before applying the multiplier.
    base_population_max <- 280L # Maximum Oligoryzomys superpopulation size before applying the multiplier.
    p_base_min <- 0.065 # Lower bound of the baseline Oligoryzomys capture probability drawn for each dataset.
    p_base_max <- 0.13 # Upper bound of the baseline Oligoryzomys capture probability drawn for each dataset.
    master_seed <- 20260920L # Use a reproducible master seed for the native R Oligoryzomys simulation experiment.
} else {
    stop("Set species to 'Akodon azarae' or 'Oligoryzomys flavescens', matching the workbook sheet names.")
}
n_occasions <- 10L # For Akodon with 30 sessions, change this value to 30 and population_multiplier to 3.
population_multiplier <- 1L # Change to 3 to triple each sampled population exactly in the Akodon 30-session design.
images_per_class <- 3L # Change to 5000 for a full run.
report_every <- 100L # Print progress after this many replicates and at the end of each demographic class.
image_width <- 512L # Horizontal image size in pixels, matching the capture-history images used in the paper.
image_height <- 512L # Vertical image size in pixels, with individuals distributed across the full image height.
dashed_season <- "warm" # Warm-season marks are dashed, and cold-season marks are continuous horizontal lines.
dash_length_px <- 6L # Length of each visible horizontal segment before rounding its inclusive pixel endpoints.
gap_length_px <- 4L # Horizontal gap between successive segments, using the same spacing as the paper renderer.
stroke_width_px <- 1L # Use the same one-pixel thickness for both the continuous and dashed capture marks.
background_color <- c(0L, 0L, 0L) # Use a black RGB background without adding seasonal bands behind the captures.
sex_colors <- list(M = c(30L, 91L, 255L), F = c(255L, 42L, 42L)) # Males are blue and females are red.
month_length_days <- 30.4375 # Convert the workbook's intervals in days into intervals expressed in months.
output_directory <- sprintf("%s_R_%d_sessions_Nx%d", species_prefix, n_occasions, population_multiplier) # Separate native R output.
classes <- c("p1", "p2", "p3", "p4", "p5") # Preserve the five demographic classes and their original meanings.
simulation_version <- "rodents_paired_R_v1" # Identify the native R simulation and paired export version.
class_rows <- data.frame(scenario = classes, Phi = c("~1", "~season", "~season", "~1", "~season+sex"), p = c("~1", "~season", "~1", "~season", "~1"), pent = c("~1", "~1", "~season", "~season", "~1"), N = rep("~1", 5L), stringsAsFactors = FALSE)
class_descriptions <- setNames(character(length(classes)), classes) # Store readable model descriptions by class.
for (class_index in seq_along(classes)) {
    class_descriptions[[classes[class_index]]] <- paste(sprintf("%s(%s)", c("Phi", "p", "pent", "N"), substring(as.character(unlist(class_rows[class_index, c("Phi", "p", "pent", "N")], use.names = FALSE)), 2L)), collapse = " ")
}

# These two helpers handle the repeated simulation and raster drawing tasks called by the generation loop.
simulate_capture_histories <- function(scenario, replicate_seed, seasons, interval_months) {
    set.seed(replicate_seed) # Each image/history pair has an independent, reproducible native R random sequence.
    base_population <- base_population_min - 1L + sample.int(base_population_max - base_population_min + 1L, size = 1L)
    parameters <- list() # Keep the sampled demographic parameters together for calculation and metadata export.
    parameters$phi_month_base <- runif(1L, min = 0.72, max = 0.90) # Draw the baseline monthly survival probability.
    parameters$p_base <- runif(1L, min = p_base_min, max = p_base_max) # Use the selected species capture range.
    parameters$beta_phi_season <- sample(c(-1, 1), size = 1L) * runif(1L, min = 0.30, max = 0.75)
    parameters$beta_p_season <- sample(c(-1, 1), size = 1L) * runif(1L, min = 0.30, max = 0.80)
    parameters$beta_pent_season <- sample(c(-1, 1), size = 1L) * runif(1L, min = 0.50, max = 1.25)
    parameters$beta_phi_sex <- sample(c(-1, 1), size = 1L) * runif(1L, min = 0.25, max = 0.70)
    n_total <- base_population * population_multiplier # Multiply the sampled population rather than its limits alone.
    sex <- ifelse(runif(n_total) < 0.50, "M", "F") # Assign the sex of each simulated individual independently.
    individual_ids <- sprintf("SIM%05d", seq_len(n_total)) # Keep stable identifiers for the individual history rows.
    male_indicator <- as.numeric(sex == "M") # The sex coefficient applies only to males in the fifth model.
    warm_survival <- as.numeric(seasons[-length(seasons)] == "warm") # Use the season at the start of each interval.
    warm_capture <- as.numeric(seasons == "warm") # Use the season of the corresponding capture or entry occasion.

    # Calculate monthly survival and apply only the seasonal and sex effects included in the current model.
    survival_predictor <- matrix(qlogis(parameters$phi_month_base), nrow = n_total, ncol = n_occasions - 1L)
    if (scenario %in% c("p2", "p3", "p5")) {
        survival_predictor <- sweep(survival_predictor, MARGIN = 2L, STATS = parameters$beta_phi_season * warm_survival, FUN = "+")
    }
    if (scenario == "p5") {
        survival_predictor <- sweep(survival_predictor, MARGIN = 1L, STATS = parameters$beta_phi_sex * male_indicator, FUN = "+")
    }
    monthly_survival <- plogis(survival_predictor) # Convert the survival linear predictor into monthly probabilities.
    interval_survival <- sweep(monthly_survival, MARGIN = 2L, STATS = interval_months, FUN = "^")

    # Calculate capture probabilities and normalize entry weights across all sampling occasions in this design.
    capture_predictor <- matrix(qlogis(parameters$p_base), nrow = n_total, ncol = n_occasions)
    if (scenario %in% c("p2", "p4")) {
        capture_predictor <- sweep(capture_predictor, MARGIN = 2L, STATS = parameters$beta_p_season * warm_capture, FUN = "+")
    }
    capture_probability <- plogis(capture_predictor) # Convert capture predictors into probabilities for all individuals.
    entry_predictor <- rep(0, n_occasions) # Constant entry weights are used unless the model includes season.
    if (scenario %in% c("p3", "p4")) {
        entry_predictor <- parameters$beta_pent_season * warm_capture # Include the seasonal entry coefficient here.
    }
    entry_weights <- exp(entry_predictor - max(entry_predictor)) # Stabilize exponentiation before normalization.
    entry_probability <- entry_weights / sum(entry_weights) # The probabilities across entry occasions sum to one.
    entry_occasion <- sample.int(n_occasions, size = n_total, replace = TRUE, prob = entry_probability)

    # Follow survival after entry, then record captures only while each individual is alive and available.
    alive <- matrix(FALSE, nrow = n_total, ncol = n_occasions) # Availability starts at the sampled entry occasion.
    captured <- matrix(0L, nrow = n_total, ncol = n_occasions) # Store each capture as a binary matrix value.
    alive[, 1L] <- entry_occasion == 1L # Individuals entering in the first session are available immediately.
    for (occasion in seq.int(2L, n_occasions)) {
        survivors <- alive[, occasion - 1L] & (runif(n_total) < interval_survival[, occasion - 1L])
        alive[, occasion] <- survivors | (entry_occasion == occasion) # Preserve new arrivals at their entry occasion.
    }
    if (any(alive & outer(entry_occasion, seq_len(n_occasions), FUN = ">")) || !all(alive[cbind(seq_len(n_total), entry_occasion)])) {
        stop("Invalid survival/entry recurrence: an individual is available before entry or absent at entry.")
    }
    for (occasion in seq_len(n_occasions)) {
        captured[, occasion] <- as.integer(alive[, occasion] & (runif(n_total) < capture_probability[, occasion]))
    }
    included <- rowSums(captured) > 0L # Retain only individuals with at least one recorded capture in this dataset.
    histories <- captured[included, , drop = FALSE] # Preserve every session, including sessions with no captures.
    observed_sex <- sex[included] # Keep sex values aligned with the retained individual history rows.
    observed_ids <- individual_ids[included] # Preserve individual identifiers for the paired history export.
    if (nrow(histories) == 0L) {
        stop(sprintf("No captured individuals in %s, seed=%d: review the parameter ranges.", scenario, replicate_seed))
    }
    observed_entry <- entry_occasion[included] # Entry among captured individuals may precede their first capture.
    first_capture <- max.col(histories, ties.method = "first") # Identify the first recorded capture for each row.
    metadata <- list(species = species, simulation_version = simulation_version, scenario = scenario, scenario_definition = unname(class_descriptions[[scenario]]), seed = replicate_seed, n_sessions = n_occasions)
    metadata$n_obs <- nrow(histories) # Number of distinct individuals retained in the capture-history table.
    metadata$n_captures_total <- sum(histories) # Total capture events, including recaptures, across all sessions.
    metadata$Nstar_total <- n_total # Keep the simulated superpopulation size before removing uncaptured individuals.
    metadata$Nstar_F <- sum(sex == "F") # Number of females in the entire simulated superpopulation for this dataset.
    metadata$Nstar_M <- sum(sex == "M") # Number of males in the entire simulated superpopulation for this dataset.
    metadata$obs_F <- sum(observed_sex == "F") # Number of captured females represented in the exported histories.
    metadata$obs_M <- sum(observed_sex == "M") # Number of captured males represented in the exported histories.
    metadata$base_Nstar_total <- base_population # Save the integer population sampled before applying the multiplier.
    metadata$population_multiplier <- population_multiplier # Record the exact population scaling used in this design.
    metadata$phi_month_base <- parameters$phi_month_base # Save the baseline monthly survival parameter for this pair.
    metadata$p_base <- parameters$p_base # Save the baseline capture probability sampled from the species-specific range.
    metadata$beta_phi_season_active <- if (scenario %in% c("p2", "p3", "p5")) parameters$beta_phi_season else 0
    metadata$beta_p_season_active <- if (scenario %in% c("p2", "p4")) parameters$beta_p_season else 0
    metadata$beta_pent_season_active <- if (scenario %in% c("p3", "p4")) parameters$beta_pent_season else 0
    metadata$beta_phi_sex_active <- if (scenario == "p5") parameters$beta_phi_sex else 0 # Other models omit sex effects.
    for (occasion in seq_len(n_occasions)) {
        metadata[[paste0("true_entry_probability_session_", occasion)]] <- entry_probability[occasion]
        metadata[[paste0("captures_session_", occasion)]] <- sum(histories[, occasion]) # Count all capture events here.
        metadata[[paste0("simulated_entries_session_", occasion)]] <- sum(entry_occasion == occasion)
        metadata[[paste0("observed_entries_session_", occasion)]] <- sum(observed_entry == occasion)
        metadata[[paste0("first_captures_session_", occasion)]] <- sum(first_capture == occasion)
        metadata[[paste0("true_p_F_session_", occasion)]] <- mean(capture_probability[sex == "F", occasion])
        metadata[[paste0("true_p_M_session_", occasion)]] <- mean(capture_probability[sex == "M", occasion])
    }
    for (occasion in seq_len(n_occasions - 1L)) {
        metadata[[paste0("true_phi_F_interval_", occasion)]] <- mean(interval_survival[sex == "F", occasion])
        metadata[[paste0("true_phi_M_interval_", occasion)]] <- mean(interval_survival[sex == "M", occasion])
    }
    return(list(histories = histories, sex = observed_sex, individual_ids = observed_ids, metadata = metadata))
}

render_capture_history <- function(histories, sex, destination, seasons) { # Draw one paired history dataset.
    if (!is.matrix(histories) || ncol(histories) != length(seasons) || nrow(histories) == 0L || length(sex) != nrow(histories) || anyNA(histories) || any(!histories %in% c(0L, 1L)) || any(rowSums(histories) == 0L) || anyNA(sex) || any(!sex %in% names(sex_colors))) {
        stop("Capture histories must be a nonempty binary matrix with matching F/M sex values and calendar.")
    }
    first_capture <- max.col(histories, ties.method = "first") # Find the first capture before sorting individual rows.
    row_order <- order(first_capture, -rowSums(histories), seq_len(nrow(histories))) # Preserve original order for ties.
    histories <- histories[row_order, , drop = FALSE] # Place early first captures near the bottom of the image.
    sex <- sex[row_order] # Apply the same row ordering to the sex values used to select the drawing color.
    image <- array(0, dim = c(image_height, image_width, 3L)) # Build a native RGB raster without antialiasing.
    for (channel in seq_len(3L)) {
        image[, , channel] <- background_color[channel] / 255 # Initialize each channel of the image background.
    }
    cell_width <- image_width / ncol(histories) # Give every capture occasion the same horizontal cell width.
    cell_height <- image_height / nrow(histories) # Use the full image height for the retained individual histories.
    half_mark <- max(1, 0.36 * cell_width) # Match the horizontal mark extent used by the original paper renderer.
    for (individual in seq_len(nrow(histories))) {
        center_y <- image_height - (individual - 0.5) * cell_height # Calculate the zero-based vertical pixel position.
        pixel_rows <- seq.int(round(center_y) - floor((stroke_width_px - 1L) / 2L), round(center_y) + floor(stroke_width_px / 2L))
        pixel_rows <- pixel_rows[pixel_rows >= 0L & pixel_rows < image_height] + 1L # Clip before converting to R indexing.
        ink <- sex_colors[[sex[individual]]] / 255 # Convert the selected RGB color into the range required by png.
        for (occasion in which(histories[individual, ] > 0L)) {
            center_x <- (occasion - 0.5) * cell_width # Calculate the zero-based horizontal center of the capture mark.
            left_edge <- center_x - half_mark # Use the same left endpoint for continuous marks and dashed segments.
            right_edge <- center_x + half_mark # Set the inclusive right endpoint of the complete horizontal mark.
            segment_starts <- left_edge # A continuous mark consists of a single segment extending to the right edge.
            if (seasons[occasion] == dashed_season) {
                segment_starts <- seq(from = left_edge, to = right_edge, by = dash_length_px + gap_length_px)
                segment_starts <- segment_starts[segment_starts < right_edge] # Exclude a segment starting at the last edge.
            }
            for (segment_start in segment_starts) {
                segment_end <- right_edge # Continuous captures keep the full horizontal extent of the mark.
                if (seasons[occasion] == dashed_season) {
                    segment_end <- min(segment_start + dash_length_px, right_edge) # Clip the final visible dash.
                }
                pixel_columns <- seq.int(round(segment_start), round(segment_end)) # Round both inclusive endpoints.
                pixel_columns <- pixel_columns[pixel_columns >= 0L & pixel_columns < image_width] + 1L
                if (length(pixel_rows) > 0L && length(pixel_columns) > 0L) {
                    for (channel in seq_len(3L)) {
                        image[pixel_rows, pixel_columns, channel] <- ink[channel] # Later marks overwrite earlier pixels.
                    }
                }
            }
        }
    }
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE) # Create the class folder if necessary.
    png::writePNG(image, target = destination) # Save exact RGB pixel values without using a graphics device.
}

# 2. Validate settings, then read the observed histories as text so their leading zeros are preserved.
integer_settings <- c(n_occasions, population_multiplier, images_per_class, report_every, image_width, image_height, dash_length_px, gap_length_px, stroke_width_px, base_population_min, base_population_max)
if (any(!is.finite(integer_settings)) || any(integer_settings < 1) || any(integer_settings != floor(integer_settings)) || !n_occasions %in% c(10L, 30L) || images_per_class >= 1000000L) {
    stop("Use 10 or 30 occasions, positive integer settings, and fewer than 1,000,000 images per class.")
}
if (!is.finite(master_seed) || master_seed < 0 || master_seed != floor(master_seed) || master_seed + length(classes) * 1000000 + images_per_class > .Machine$integer.max) {
    stop("Choose a nonnegative integer master_seed whose derived replicate seeds fit R's integer range.")
}
if (base_population_min > base_population_max || base_population_max * population_multiplier > .Machine$integer.max || !is.finite(p_base_min) || !is.finite(p_base_max) || !(0 < p_base_min && p_base_min <= p_base_max && p_base_max < 1)) {
    stop("Population limits must be positive and ordered, and capture probability limits must lie between zero and one.")
}
if (!dashed_season %in% c("cold", "warm") || !is.finite(month_length_days) || month_length_days <= 0) {
    stop("Use a cold/warm dashed season and a positive finite month_length_days conversion value.")
}
source_workbook <- normalizePath(path.expand(workbook_file), winslash = "/", mustWork = TRUE)
observed <- as.data.frame(readxl::read_excel(source_workbook, sheet = history_sheet, col_types = "text", .name_repair = "minimal"), stringsAsFactors = FALSE)
base_sessions <- as.data.frame(readxl::read_excel(source_workbook, sheet = calendar_sheet, col_types = "text", .name_repair = "minimal"), stringsAsFactors = FALSE)
if (nrow(observed) == 0L || !all(c("ID_individual", "freq", "sex", "ch") %in% names(observed))) {
    stop("The selected history sheet requires ID_individual, freq, sex and ch, with one row per individual.")
}
if (!all(c("session", "time_interval", "season") %in% names(base_sessions)) || nrow(base_sessions) != 10L) {
    stop("The sessions worksheet requires session, time_interval and season, with exactly ten rows.")
}
base_sessions$session <- as.numeric(base_sessions$session) # Convert session numbers after reading all columns as text.
base_sessions <- base_sessions[order(base_sessions$session), , drop = FALSE] # Use the chronological session order.
if (anyNA(base_sessions$session) || !identical(base_sessions$session, as.numeric(seq_len(10L)))) {
    stop("The workbook calendar must contain consecutive sessions numbered 1 through 10, without duplicates.")
}
observed$ch <- trimws(observed$ch) # Preserve each original ten-session binary string, including its leading zeros.
observed$sex <- toupper(trimws(observed$sex)) # Standardize the observed sex values before choosing image colors.
base_sessions$season <- tolower(trimws(base_sessions$season)) # Use the cold/warm labels expected by the models.
if (anyNA(observed$ch) || any(!grepl("^[01]{10}$", observed$ch)) || any(!grepl("1", observed$ch, fixed = TRUE))) {
    stop("Observed histories must be nonzero binary strings of length 10, stored as text in the workbook.")
}
if (anyNA(observed$ID_individual) || any(trimws(observed$ID_individual) == "") || anyDuplicated(observed$ID_individual) || anyNA(observed$sex) || any(!observed$sex %in% c("F", "M"))) {
    stop("Observed individuals require unique nonempty identifiers and sex values coded as F or M.")
}
observed$freq <- as.numeric(observed$freq) # Require one observed history row per individual rather than grouped counts.
if (anyNA(observed$freq) || any(observed$freq != 1)) {
    stop("One row per observed individual and freq=1 are required for capture-history image rendering.")
}
if (!is.na(base_sessions$time_interval[1L])) {
    stop("The first session has no preceding sampling interval: leave its time_interval cell empty.")
}
base_sessions$time_interval <- as.numeric(base_sessions$time_interval) # Keep the first interval missing as in the input.
base_intervals <- base_sessions$time_interval[-1L] # Extract the nine intervals between the ten observed sessions.
if (anyNA(base_sessions$season) || !setequal(base_sessions$season, c("cold", "warm")) || !is.na(base_sessions$time_interval[1L]) || any(!is.finite(base_intervals)) || any(base_intervals <= 0)) {
    stop("Use cold/warm seasons, a missing first interval, and positive finite day intervals after the first session.")
}

# 3. For Akodon with 30 sessions and N times three, change the two design settings in step 1 to 30 and 3.
# Repeat the ten-session calendar three times, using its median interval to connect successive calendar blocks.
# This expands the simulation design, while the empirical workbook histories remain ten-session observations.
sessions <- base_sessions # The ten-session design uses the input calendar without extending the observation period.
boundary_interval_days <- NULL # No connecting interval is needed when using the original ten-session calendar.
if (n_occasions == 30L) {
    boundary_interval_days <- median(base_intervals) # Use an explicit, reproducible interval between calendar blocks.
    expanded_intervals <- c(base_intervals, boundary_interval_days, base_intervals, boundary_interval_days, base_intervals)
    sessions <- data.frame(session = seq_len(30L), time_interval = c(NA_real_, expanded_intervals), season = rep(base_sessions$season, times = 3L), stringsAsFactors = FALSE)
}
sessions$label <- sprintf("S%02d", sessions$session) # Use the same session labels in the calendar and diagnostics.
sessions$interval_months <- sessions$time_interval / month_length_days # Convert day intervals for survival calculations.
seasons <- sessions$season # Keep the full design calendar available to both simulation and image rendering.
interval_months <- sessions$interval_months[-1L] # Survival uses the interval between successive sampling occasions.

# 4. Record the experiment configuration, class formulas and calendar, then export the observed ten-session data.
configuration <- list(species = species, workbook_file = workbook_file, history_sheet = history_sheet, calendar_sheet = calendar_sheet, n_occasions = n_occasions, population_multiplier = population_multiplier, images_per_class = images_per_class, master_seed = master_seed, image_width = image_width, image_height = image_height, base_population_min = base_population_min, base_population_max = base_population_max, p_base_min = p_base_min, p_base_max = p_base_max, dashed_season = dashed_season, dash_length_px = dash_length_px, gap_length_px = gap_length_px, stroke_width_px = stroke_width_px, background_color = background_color, sex_colors = sex_colors, month_length_days = month_length_days)
information <- list(species = species, simulation_version = simulation_version, config = configuration, source_sha256 = digest::digest(file = source_workbook, algo = "sha256"), n_sessions = n_occasions, Nstar_range = c(base_population_min, base_population_max) * population_multiplier, calendar_expansion = if (n_occasions == 30L) "three repetitions of ten sessions" else "workbook calendar", boundary_interval_days = boundary_interval_days, observed_export = n_occasions == 10L, models = class_rows, rng_kind = RNGkind(), seed_rule = "master_seed + class_index * 1000000 + replicate", R_version = R.version.string)
information$configuration_signature <- digest::digest(jsonlite::toJSON(information, auto_unbox = TRUE, digits = NA, null = "null", na = "null"), algo = "sha256", serialize = FALSE)
root <- path.expand(output_directory) # Resolve the output relative to the working directory selected for this run.
dir.create(root, recursive = TRUE, showWarnings = FALSE) # Create the root directory before recording its full path.
root <- normalizePath(root, winslash = "/", mustWork = TRUE) # Use an absolute output path in the progress messages.
metadata_file <- file.path(root, "experiment_info.json") # Use this record to prevent mixing different experiments.
if (file.exists(metadata_file)) {
    previous_information <- jsonlite::read_json(metadata_file, simplifyVector = TRUE) # Read the existing experiment signature.
    if (!identical(previous_information$configuration_signature, information$configuration_signature)) {
        stop("This output directory has a different configuration: choose another output_directory.")
    }
} else if (dir.exists(file.path(root, "images")) || dir.exists(file.path(root, "capture_histories"))) {
    stop("This output directory contains an untracked dataset: choose another output_directory.")
}
jsonlite::write_json(information, path = metadata_file, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
jsonlite::write_json(configuration, path = file.path(root, "config.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
write.csv(sessions, file = file.path(root, "sessions.csv"), row.names = FALSE, na = "") # Export the actual design calendar.
write.csv(class_rows, file = file.path(root, "class_definitions.csv"), row.names = FALSE) # Preserve model meanings for analysis.
if (n_occasions == 10L) {
    observed_histories <- matrix(as.integer(unlist(strsplit(observed$ch, split = "", fixed = TRUE), use.names = FALSE)), nrow = nrow(observed), ncol = 10L, byrow = TRUE)
    render_capture_history(observed_histories, observed$sex, file.path(root, "observed.png"), seasons)
    observed_connection <- gzfile(file.path(root, "observed_capture_histories.csv.gz"), open = "wt", encoding = "UTF-8")
    tryCatch(write.csv(observed, file = observed_connection, row.names = FALSE, na = ""), finally = close(observed_connection))
}
parameter_rows <- vector("list", length(classes) * images_per_class) # Reserve one metadata record per generated pair.
index_rows <- vector("list", length(classes) * images_per_class) # Reserve the image/history mapping in generation order.
pair_position <- 0L # Advance once per matched image and history table written by the sequential generation loop.
cat(sprintf("Generating %d-session %s datasets, N multiplier=%d, output=%s\n", n_occasions, species, population_multiplier, root))

# 5. Generate each class and replicate, then save its image and the same individual capture histories together.
# Native R uses different random sequences from NumPy, while preserving the model equations and distributions.
for (class_index in seq_along(classes)) {
    scenario <- classes[class_index] # Select the demographic formula corresponding to this generated image class.
    for (replicate in seq_len(images_per_class)) {
        replicate_seed <- as.integer(master_seed + class_index * 1000000 + replicate) # Unique integer seed for this pair.
        simulation_result <- simulate_capture_histories(scenario, replicate_seed, seasons, interval_months)
        histories <- simulation_result$histories # The renderer and compressed history export use the same binary matrix.
        sex <- simulation_result$sex # Keep the sex values aligned with individual rows before rendering sorts its copy.
        individual_ids <- simulation_result$individual_ids # Preserve identifiers in the individual history export.
        metadata <- simulation_result$metadata # Retain true parameters and counts for downstream model assessment.
        basename <- sprintf("%s_%05d", scenario, replicate) # Use matching file basenames for the image and history pair.
        image_relative <- paste0(scenario, "/", basename, ".png") # Store a portable image path relative to images/.
        history_relative <- paste0(scenario, "/", basename, ".csv.gz") # Store a portable path relative to capture_histories/.
        render_capture_history(histories, sex, file.path(root, "images", image_relative), seasons)
        history_file <- file.path(root, "capture_histories", history_relative) # Select the compressed individual history file.
        dir.create(dirname(history_file), recursive = TRUE, showWarnings = FALSE) # Create the history folder for this class.
        history_table <- data.frame(scenario = scenario, replicate = replicate, seed = replicate_seed, image = image_relative, ID_individual = individual_ids, freq = 1L, sex = sex, ch = apply(histories, MARGIN = 1L, FUN = paste0, collapse = ""), stringsAsFactors = FALSE)
        temporary_file <- paste0(history_file, ".part") # Finish writing the compressed table before naming the final file.
        history_connection <- gzfile(temporary_file, open = "wt", encoding = "UTF-8") # Save a gzip-compressed UTF-8 CSV.
        tryCatch(write.csv(history_table, file = history_connection, row.names = FALSE), finally = close(history_connection))
        if (file.exists(history_file) && !file.remove(history_file)) {
            stop(sprintf("Could not replace the existing capture-history table: %s", history_file))
        }
        if (!file.rename(temporary_file, history_file)) {
            stop(sprintf("Could not finalize the compressed capture-history table: %s", history_file))
        }
        metadata$replicate <- replicate # Record the replicate number together with the true simulation parameters.
        metadata$image_file <- image_relative # Link the metadata row to its image within the images directory.
        metadata$history_file <- history_relative # Link the same metadata row to the individual capture histories.
        pair_position <- pair_position + 1L # Maintain a single position across the five classes and their replicates.
        parameter_rows[[pair_position]] <- as.data.frame(metadata, stringsAsFactors = FALSE) # Store the parameter row.
        index_rows[[pair_position]] <- data.frame(scenario = scenario, replicate = replicate, seed = replicate_seed, image_file = image_relative, history_file = history_relative, n_obs = nrow(histories), n_captures_total = sum(histories), Nstar_total = metadata$Nstar_total, stringsAsFactors = FALSE)
        if (replicate %% report_every == 0L || replicate == images_per_class) {
            cat(sprintf("%s: %d/%d image/history pairs\n", scenario, replicate, images_per_class))
        }
    }
}

# 6. Export file mappings, true simulation parameters, per-class summaries and temporal entry diagnostics.
index_table <- do.call(rbind, index_rows) # Combine the individual mapping records without changing generation order.
write.csv(index_table, file = file.path(root, "capture_history_index.csv"), row.names = FALSE)
parameters_table <- do.call(rbind, parameter_rows) # Combine the true parameters and counts for every generated pair.
write.csv(parameters_table, file = file.path(root, "parameters_per_image.csv"), row.names = FALSE, na = "")
summary_rows <- vector("list", length(classes)) # Store one summary row for each demographic class in the experiment.
diagnostic_rows <- vector("list", length(classes) * n_occasions) # Store one temporal diagnostic row per class and session.
diagnostic_position <- 0L # Advance once for each class/session combination entered into the diagnostic table.
for (class_index in seq_along(classes)) {
    scenario <- classes[class_index] # Keep the same class order used for simulation and the exported file mapping.
    class_parameters <- parameters_table[parameters_table$scenario == scenario, , drop = FALSE]
    summary_rows[[class_index]] <- data.frame(scenario = scenario, n_images = nrow(class_parameters), mean_n_obs = mean(class_parameters$n_obs), mean_n_captures_total = mean(class_parameters$n_captures_total), mean_Nstar_total = mean(class_parameters$Nstar_total), stringsAsFactors = FALSE)
    for (occasion in seq_len(n_occasions)) {
        diagnostic_position <- diagnostic_position + 1L # Count this class/session combination in the diagnostic export.
        diagnostic_rows[[diagnostic_position]] <- data.frame(scenario = scenario, session = occasion, label = sprintf("S%02d", occasion), season = seasons[occasion], simulated_entries = sum(class_parameters[[paste0("simulated_entries_session_", occasion)]]), observed_entries = sum(class_parameters[[paste0("observed_entries_session_", occasion)]]), first_captures = sum(class_parameters[[paste0("first_captures_session_", occasion)]]), stringsAsFactors = FALSE)
    }
}
write.csv(do.call(rbind, summary_rows), file = file.path(root, "summary_by_class.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostic_rows), file = file.path(root, "temporal_entry_diagnostic.csv"), row.names = FALSE)
cat(sprintf("Completed %d matched images and history tables: %s\n", nrow(index_table), root))
