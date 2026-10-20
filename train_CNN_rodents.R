# ==================================================================================================
# Generic CNN training for generated rodent images in images/p1, images/p2, images/p3, images/p4, images/p5
# Run the whole script after restarting R, using one generated image dataset for each training run.
# ==================================================================================================

# 1. Select the generated dataset and keep or change the original CNN training settings.
# Set the folder containing images/p1 to images/p5, or use the images folder itself.
dataset_directory <- "akodon_R_10_sessions_Nx1/images"
python_executable <- "/usr/bin/python3.12" # Leave empty to use reticulate's Python, or enter the path to your TensorFlow Python.
observed_image_file <- "" # Optional full path, otherwise the script searches the dataset for an observed PNG.
cfg <- list(cnn_width = 128L, cnn_height = 32L, resize_method = "area", seed = 123L, batch_size = 32L, epochs = 50L, patience = 8L, learning_rate = 3e-4, validation_fraction = 0.20)
classes <- c("p1", "p2", "p3", "p4", "p5") # Class meanings are specific to the selected generated dataset.

if (!dir.exists(dataset_directory)) stop("Dataset directory does not exist: ", dataset_directory)
dataset_directory <- normalizePath(dataset_directory, winslash = "/", mustWork = TRUE) # Resolve the selected folder.
if (basename(dataset_directory) == "images" && all(dir.exists(file.path(dataset_directory, classes)))) {
    dataset_directory <- dirname(dataset_directory) # Accept a path pointing directly to the images folder.
}
images_directory <- file.path(dataset_directory, "images") # All generated training PNGs are read from this folder.
run_name <- format(Sys.time(), "run_%Y%m%d_%H%M%S") # Give each training run its own generic result folder.
results_directory <- file.path(dataset_directory, "cnn_results", run_name) # Save one selected dataset per run.
if (dir.exists(results_directory)) results_directory <- tempfile(pattern = paste0(run_name, "_"), tmpdir = file.path(dataset_directory, "cnn_results"))
class_directories <- file.path(images_directory, classes) # Each directory supplies its own integer class label.
if (any(!dir.exists(class_directories))) stop("Missing class folders: ", paste(class_directories[!dir.exists(class_directories)], collapse = ", "))
if (!is.finite(cfg$validation_fraction) || cfg$validation_fraction <= 0 || cfg$validation_fraction >= 1) stop("validation_fraction must be between zero and one.")
if (any(!is.finite(c(cfg$cnn_height, cfg$cnn_width, cfg$batch_size, cfg$epochs, cfg$patience, cfg$seed, cfg$learning_rate)))) stop("CNN settings must contain finite numeric values.")
if (any(c(cfg$cnn_height, cfg$cnn_width, cfg$batch_size, cfg$epochs, cfg$patience) < 1) || cfg$learning_rate <= 0) stop("CNN dimensions, batch size, epochs, patience and learning rate must be positive.")
if (!cfg$resize_method %in% c("area", "bilinear", "nearest", "bicubic", "lanczos3", "lanczos5", "gaussian", "mitchellcubic")) stop("Unsupported TensorFlow resize method.")

# 2. Read the model definitions when the generator exported them, otherwise retain the class labels.
model_definitions <- setNames(classes, classes) # No species-specific model formula is assumed by the CNN.
definitions_file <- file.path(dataset_directory, "class_definitions.csv") # This file is optional for training.
if (file.exists(definitions_file)) { # Read the definitions supplied by the selected image generator.
    definitions_table <- read.csv(definitions_file, stringsAsFactors = FALSE, check.names = FALSE)
    class_column <- intersect(c("scenario", "class", "class_name"), names(definitions_table)) # Accept generator labels.
    if (!length(class_column)) stop("class_definitions.csv must contain a scenario, class or class_name column.")
    definition_positions <- match(classes, definitions_table[[class_column[1L]]]) # Preserve the p1 to p5 order.
    if (anyNA(definition_positions) || anyDuplicated(definitions_table[[class_column[1L]]])) stop("class_definitions.csv must contain one unique definition for each class p1 to p5.")
    definitions_table <- definitions_table[definition_positions, , drop = FALSE] # Match definitions to CNN outputs.
    if (all(c("Phi", "p", "pent", "N") %in% names(definitions_table))) { # Use the supplied POPAN formulas.
        model_definitions <- setNames(paste0("Phi ", definitions_table$Phi, " | p ", definitions_table$p, " | pent ", definitions_table$pent, " | N ", definitions_table$N), classes)
    } else {
        definition_column <- intersect(c("definition", "model", "description"), names(definitions_table))
        if (length(definition_column)) model_definitions <- setNames(as.character(definitions_table[[definition_column[1L]]]), classes)
    }
}
model_table <- data.frame(class = classes, model = unname(model_definitions[classes]), stringsAsFactors = FALSE)
cat("\nClass definitions for the selected dataset:\n") # Print the model meaning associated with each CNN output.
print(model_table) # The labels p1 to p5 refer to the definitions supplied by this particular dataset.

# 3. Find the optional observed image before training, allowing datasets without an observed image.
if (nzchar(observed_image_file)) { # Check an explicitly selected observed image before starting training.
    if (!file.exists(observed_image_file)) stop("The configured observed image does not exist: ", observed_image_file)
    observed_image_file <- normalizePath(observed_image_file, winslash = "/", mustWork = TRUE)
} else {
    default_observed_file <- file.path(dataset_directory, "observed.png") # Preferred name used by the generators.
    observed_candidates <- list.files(dataset_directory, pattern = "(^observed|_observed)\\.png$", full.names = TRUE, ignore.case = TRUE)
    if (file.exists(default_observed_file)) observed_image_file <- default_observed_file # Use the standard file first.
    if (!nzchar(observed_image_file) && length(observed_candidates) == 1L) observed_image_file <- observed_candidates[1L]
    if (!nzchar(observed_image_file) && length(observed_candidates) > 1L) message("Several observed PNGs were found, set observed_image_file to select one. Training will continue.")
}

# 4. Split each class separately, keeping all available images and recording the exact file assignment.
set.seed(as.integer(cfg$seed)) # Use a reproducible R seed for the training and validation file partition.
manifest_rows <- vector("list", length(classes)) # Build one simple file table for each generating class.
for (class_index in seq_along(classes)) { # Make the training and validation partition within each class.
    class_files <- sort(list.files(class_directories[class_index], pattern = "\\.png$", full.names = TRUE, ignore.case = TRUE))
    class_files <- class_files[!file.info(class_files)$isdir] # Exclude directories whose names happen to end in PNG.
    image_count <- length(class_files) # Use the actual number of generated files rather than a configured count.
    if (image_count < 2L) stop("Class ", classes[class_index], " needs at least two PNGs for training and validation.")
    class_files <- normalizePath(class_files, winslash = "/", mustWork = TRUE) # Keep file paths portable on Windows.
    validation_count <- max(1L, min(image_count - 1L, as.integer(floor(image_count * cfg$validation_fraction))))
    validation_positions <- sort(sample.int(image_count, size = validation_count)) # Sample within this class only.
    class_manifest <- data.frame(image_file = class_files, class = classes[class_index], label = class_index - 1L, subset = "training", stringsAsFactors = FALSE)
    class_manifest$subset[validation_positions] <- "validation" # Reserve at least one image from every class.
    manifest_rows[[class_index]] <- class_manifest # Keep the selected paths in their original sorted order.
}
split_manifest <- do.call(rbind, manifest_rows) # Combine the five class tables without mixing their labels.
training_manifest <- split_manifest[split_manifest$subset == "training", , drop = FALSE] # Select training files.
validation_manifest <- split_manifest[split_manifest$subset == "validation", , drop = FALSE] # Keep validation order.
if (anyDuplicated(split_manifest$image_file)) stop("Repeated image paths were found in the selected class folders.")
if (length(intersect(training_manifest$image_file, validation_manifest$image_file))) stop("Training and validation files overlap.")
training_counts <- tabulate(training_manifest$label + 1L, nbins = length(classes)) # Count each training class.
validation_counts <- tabulate(validation_manifest$label + 1L, nbins = length(classes)) # Count each validation class.
distribution_table <- data.frame(class = classes, total_images = training_counts + validation_counts, training_images = training_counts, validation_images = validation_counts, validation_fraction = validation_counts / (training_counts + validation_counts))
if (any(training_counts == 0L) || any(validation_counts == 0L)) stop("Every class must occur in both data subsets.")
if (length(unique(distribution_table$total_images)) != 1L) warning("Classes contain different image counts, training will use all available images.")
cat("\nTraining and validation image counts:\n") # Check the class distribution before starting the CNN.
print(distribution_table) # Very small datasets use at least one validation image per class.
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE) # Create the selected dataset's result folder.
write.csv(split_manifest, file.path(results_directory, "split_manifest.csv"), row.names = FALSE)
write.csv(distribution_table, file.path(results_directory, "dataset_distribution.csv"), row.names = FALSE)
write.csv(model_table, file.path(results_directory, "class_definitions_used.csv"), row.names = FALSE)
saveRDS(list(dataset_directory = dataset_directory, observed_image_file = observed_image_file, config = cfg, classes = model_table), file.path(results_directory, "training_configuration.rds"))

# 5. Load the R interfaces and initialize the Python installation that provides TensorFlow and Keras.
required_packages <- c("reticulate", "keras3", "tensorflow", "tfdatasets", "ggplot2") # Required R packages.
for (package_name in required_packages) { # Report any missing R dependency before initializing Keras.
    if (!requireNamespace(package_name, quietly = TRUE)) stop("Install the missing R package before running this script: ", package_name)
}
if (nzchar(python_executable)) reticulate::use_python(python_executable, required = TRUE) # Set Python before loading Keras.
Sys.setenv(KERAS_BACKEND = "tensorflow") # Use the TensorFlow backend for the existing CNN architecture.
library(keras3) # Provides the original CNN layers, optimizer, callbacks and model training methods.
library(tensorflow) # Provides image decoding, resizing, random seeds and the TensorFlow dataset backend.
library(tfdatasets) # Converts the explicit file partition into streamed TensorFlow training datasets.
library(ggplot2) # Draws learning curves and the row-normalized validation confusion matrix.
keras3::set_random_seed(as.integer(cfg$seed)) # Initialize the R, Python and TensorFlow seeds used by Keras.
tf$random$set_seed(as.integer(cfg$seed)) # Keep the TensorFlow random seed explicit for the selected run.
cat("\nR version:", R.version.string, "\n") # Record the software used for this training run.
print(reticulate::py_config()) # Confirm that reticulate selected the intended TensorFlow Python installation.
cat("\nAvailable TensorFlow GPUs:\n") # The same script also works when TensorFlow uses the CPU.
print(tf$config$list_physical_devices("GPU")) # Report available devices without requiring a specific GPU.
runtime_information <- c(paste("R:", R.version.string), paste("keras3:", packageVersion("keras3")), paste("tensorflow R:", packageVersion("tensorflow")), capture.output(reticulate::py_config()))
writeLines(runtime_information, file.path(results_directory, "runtime_information.txt")) # Save runtime details.

# 6. Decode and resize images in one shared function, also used later for the optional observed image.
load_image <- function(image_file, class_label) { # Reuse the same image preprocessing for every prediction.
    image_bytes <- tf$io$read_file(image_file) # Read one generated PNG without loading the whole dataset into RAM.
    image_tensor <- tf$io$decode_png(image_bytes, channels = 3L) # Decode female, male and background RGB values.
    image_tensor <- tf$cast(image_tensor, dtype = tf$float32) # Keep values on the zero to 255 scale for rescaling.
    image_tensor <- tf$image$resize(image_tensor, size = as.integer(c(cfg$cnn_height, cfg$cnn_width)), method = cfg$resize_method)
    reticulate::tuple(image_tensor, tf$cast(class_label, dtype = tf$int32)) # Return the image and its zero-based label.
}
autotune <- tf$data$AUTOTUNE # Let TensorFlow choose the parallel image loading and prefetch settings.
training_dataset <- tensor_slices_dataset(reticulate::tuple(tf$constant(training_manifest$image_file, dtype = tf$string), tf$constant(as.integer(training_manifest$label), dtype = tf$int32)))
validation_dataset <- tensor_slices_dataset(reticulate::tuple(tf$constant(validation_manifest$image_file, dtype = tf$string), tf$constant(as.integer(validation_manifest$label), dtype = tf$int32)))
training_dataset <- dataset_shuffle(training_dataset, buffer_size = as.integer(nrow(training_manifest)), seed = as.integer(cfg$seed), reshuffle_each_iteration = TRUE)
training_dataset <- dataset_map(training_dataset, map_func = load_image, num_parallel_calls = autotune)
validation_dataset <- dataset_map(validation_dataset, map_func = load_image, num_parallel_calls = autotune)
dataset_options <- tf$data$Options() # Keep parallel mapping deterministic so validation paths remain aligned.
dataset_options$experimental_deterministic <- TRUE # Preserve the order recorded in validation_manifest.
training_dataset <- training_dataset$with_options(dataset_options) # Shuffle only the training input files.
validation_dataset <- validation_dataset$with_options(dataset_options) # Validation is never shuffled.
training_dataset <- dataset_batch(training_dataset, batch_size = as.integer(cfg$batch_size)) # Keep the final batch.
validation_dataset <- dataset_batch(validation_dataset, batch_size = as.integer(cfg$batch_size)) # Use every held-out file.
training_dataset <- dataset_prefetch(training_dataset, buffer_size = autotune) # Prepare the next training batch.
validation_dataset <- dataset_prefetch(validation_dataset, buffer_size = autotune) # Prepare the next validation batch.

# 7. Keep the original CNN architecture, building it one layer at a time so the script remains sequential.
model <- keras_model_sequential(input_shape = c(as.integer(cfg$cnn_height), as.integer(cfg$cnn_width), 3L))
model <- layer_rescaling(model, scale = 1 / 255) # Apply the same zero to one rescaling to every input image.
model <- layer_conv_2d(model, filters = 32L, kernel_size = c(3L, 3L), padding = "same", use_bias = FALSE)
model <- layer_batch_normalization(model) # Normalize the first convolution's feature activations.
model <- layer_activation(model, activation = "relu") # Use the original nonlinear activation function.
model <- layer_conv_2d(model, filters = 32L, kernel_size = c(3L, 3L), padding = "same", use_bias = FALSE)
model <- layer_batch_normalization(model) # Normalize the second convolution's feature activations.
model <- layer_activation(model, activation = "relu") # Complete the first convolutional feature block.
model <- layer_max_pooling_2d(model, pool_size = c(2L, 2L), padding = "same") # Reduce spatial dimensions.
model <- layer_conv_2d(model, filters = 64L, kernel_size = c(3L, 3L), padding = "same", use_bias = FALSE)
model <- layer_batch_normalization(model) # Normalize the next convolutional feature block.
model <- layer_activation(model, activation = "relu") # Retain positive nonlinear feature activations.
model <- layer_conv_2d(model, filters = 64L, kernel_size = c(3L, 3L), padding = "same", use_bias = FALSE)
model <- layer_batch_normalization(model) # Normalize the second convolution within the 64-filter block.
model <- layer_activation(model, activation = "relu") # Complete the second convolutional feature block.
model <- layer_max_pooling_2d(model, pool_size = c(2L, 2L), padding = "same") # Reduce spatial dimensions.
model <- layer_spatial_dropout_2d(model, rate = 0.10) # Use the original regularization for this block.
model <- layer_conv_2d(model, filters = 128L, kernel_size = c(3L, 3L), padding = "same", use_bias = FALSE)
model <- layer_batch_normalization(model) # Normalize the first convolution within the 128-filter block.
model <- layer_activation(model, activation = "relu") # Retain positive nonlinear feature activations.
model <- layer_conv_2d(model, filters = 128L, kernel_size = c(3L, 3L), padding = "same", use_bias = FALSE)
model <- layer_batch_normalization(model) # Normalize the second convolution within the 128-filter block.
model <- layer_activation(model, activation = "relu") # Complete the third convolutional feature block.
model <- layer_max_pooling_2d(model, pool_size = c(2L, 2L), padding = "same") # Reduce spatial dimensions.
model <- layer_spatial_dropout_2d(model, rate = 0.15) # Use the original regularization for this block.
model <- layer_conv_2d(model, filters = 256L, kernel_size = c(3L, 3L), padding = "same", use_bias = FALSE)
model <- layer_batch_normalization(model) # Normalize the first convolution within the 256-filter block.
model <- layer_activation(model, activation = "relu") # Retain positive nonlinear feature activations.
model <- layer_conv_2d(model, filters = 256L, kernel_size = c(3L, 3L), padding = "same", use_bias = FALSE)
model <- layer_batch_normalization(model) # Normalize the second convolution within the 256-filter block.
model <- layer_activation(model, activation = "relu") # Complete the fourth convolutional feature block.
model <- layer_max_pooling_2d(model, pool_size = c(2L, 2L), padding = "same") # Reduce spatial dimensions.
model <- layer_spatial_dropout_2d(model, rate = 0.20) # Use the original regularization for this block.
model <- layer_global_average_pooling_2d(model) # Summarize each convolutional feature map over its pixels.
model <- layer_dense(model, units = 256L, use_bias = FALSE) # Preserve the original hidden dense layer.
model <- layer_batch_normalization(model) # Normalize the hidden dense layer before its activation.
model <- layer_activation(model, activation = "relu") # Retain positive nonlinear hidden activations.
model <- layer_dropout(model, rate = 0.30) # Preserve the original dense-layer dropout regularization.
model <- layer_dense(model, units = as.integer(length(classes)), activation = "softmax") # Return five class scores.
compile(object = model, optimizer = optimizer_adam(learning_rate = as.numeric(cfg$learning_rate)), loss = "sparse_categorical_crossentropy", metrics = "accuracy")
summary(model) # Print the original CNN architecture for the selected input dimensions and five classes.

# 8. Train the CNN and save the model with the lowest validation loss in this dataset's result folder.
best_model_file <- file.path(results_directory, "best_model.keras") # Use a generic name for any rodent dataset.
callbacks <- list(callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5, patience = 3L, min_delta = 1e-4, min_lr = 1e-7, verbose = 1L), callback_model_checkpoint(filepath = best_model_file, monitor = "val_loss", mode = "min", save_best_only = TRUE, save_weights_only = FALSE, verbose = 1L), callback_early_stopping(monitor = "val_loss", mode = "min", patience = as.integer(cfg$patience), min_delta = 1e-4, restore_best_weights = TRUE, verbose = 1L))
cat("\nStarting CNN training for:", basename(dataset_directory), "\n") # Train one selected image design.
history <- fit(object = model, x = training_dataset, validation_data = validation_dataset, epochs = as.integer(cfg$epochs), callbacks = callbacks, verbose = 1L, view_metrics = "auto")
history_table <- data.frame(epoch = seq_along(history$metrics$loss), as.data.frame(history$metrics), check.names = FALSE)
write.csv(history_table, file.path(results_directory, "training_history.csv"), row.names = FALSE)
history_plot <- plot(history) + labs(x = "Epoch") + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 16), legend.text = element_text(size = 15), strip.text = element_text(size = 16))
print(history_plot) # Display the learning curves for the training and validation image subsets.
ggsave(filename = file.path(results_directory, "training_history.png"), plot = history_plot, width = 10, height = 6, dpi = 300)

# 9. Load the best saved model and predict validation images in the same order as their recorded paths.
if (!file.exists(best_model_file)) stop("Training did not create the expected best_model.keras file.")
best_model <- load_model(best_model_file, compile = TRUE) # Evaluate the actual checkpoint selected by val_loss.
validation_results <- evaluate(object = best_model, x = validation_dataset, verbose = 1L) # Use the full holdout subset.
print(validation_results) # Report Keras validation loss and accuracy for the selected best checkpoint.
evaluation_table <- data.frame(metric = names(validation_results), value = as.numeric(unlist(validation_results)), stringsAsFactors = FALSE)
write.csv(evaluation_table, file.path(results_directory, "validation_evaluation.csv"), row.names = FALSE)
validation_probabilities <- as.matrix(predict(object = best_model, x = validation_dataset, verbose = 1L))
if (nrow(validation_probabilities) != nrow(validation_manifest) || ncol(validation_probabilities) != length(classes)) stop("Validation prediction dimensions do not match the recorded image partition.")
predicted_classes <- classes[max.col(validation_probabilities, ties.method = "first")] # Select the largest softmax score.
colnames(validation_probabilities) <- paste0("probability_", classes) # Record the class order explicitly.
validation_predictions <- data.frame(image_file = validation_manifest$image_file, true_class = validation_manifest$class, predicted_class = predicted_classes, correct = validation_manifest$class == predicted_classes, validation_probabilities, stringsAsFactors = FALSE, check.names = FALSE)
write.csv(validation_predictions, file.path(results_directory, "validation_predictions.csv"), row.names = FALSE)

# 10. Calculate and save the count matrix, row-normalized matrix and classification performance by class.
true_factor <- factor(validation_manifest$class, levels = classes) # Rows represent the generating classes.
predicted_factor <- factor(predicted_classes, levels = classes) # Columns represent the CNN's predicted classes.
confusion_counts <- table(True = true_factor, Predicted = predicted_factor) # Include all five classes in the matrix.
confusion_proportions <- sweep(confusion_counts, 1L, rowSums(confusion_counts), FUN = "/") # Normalize by true class.
accuracy_score <- sum(diag(confusion_counts)) / sum(confusion_counts) # Overall validation classification accuracy.
class_metrics <- data.frame(class = classes, validation_images = as.integer(rowSums(confusion_counts)), correctly_classified = as.integer(diag(confusion_counts)), recall = as.numeric(diag(confusion_proportions)), precision = ifelse(colSums(confusion_counts) > 0, diag(confusion_counts) / colSums(confusion_counts), NA_real_), stringsAsFactors = FALSE)
cat("\nValidation confusion counts:\n") # Print the number of predictions for each true and predicted class.
print(confusion_counts) # Correct classifications occur on the main diagonal of the count matrix.
cat("\nRow-normalized confusion matrix:\n") # Print the validation proportions within each true generating class.
print(round(confusion_proportions, 3)) # Each row sums to one because every class contributes validation images.
cat("\nOverall validation accuracy:", sprintf("%.4f", accuracy_score), "\n") # Report the fraction classified correctly.
print(class_metrics) # Display validation recall and precision separately for each generating class.
write.csv(as.data.frame.matrix(confusion_counts), file.path(results_directory, "confusion_counts.csv"), row.names = TRUE)
write.csv(as.data.frame.matrix(confusion_proportions), file.path(results_directory, "confusion_proportions.csv"), row.names = TRUE)
write.csv(class_metrics, file.path(results_directory, "validation_metrics_by_class.csv"), row.names = FALSE)
confusion_data <- as.data.frame(as.table(confusion_proportions), responseName = "Proportion") # Prepare the plot table.
confusion_data$True <- factor(confusion_data$True, levels = rev(classes)) # Show p1 at the top of the plot.
confusion_data$Predicted <- factor(confusion_data$Predicted, levels = classes) # Show predicted classes from p1 to p5.
confusion_plot <- ggplot(confusion_data, aes(x = Predicted, y = True, fill = Proportion)) + geom_tile(color = "white", linewidth = 1.0) + geom_text(aes(label = sprintf("%.3f", Proportion)), size = 4.8, color = "black") + scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1), name = "Proportion") + coord_equal() + labs(title = "Row-normalized confusion matrix", subtitle = sprintf("%s, validation accuracy = %.3f", basename(dataset_directory), accuracy_score), x = "Predicted class", y = "True class") + theme_minimal(base_size = 15) + theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15), axis.title = element_text(size = 17), axis.text = element_text(size = 14, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
print(confusion_plot) # Display the same confusion-matrix style for any of the generated rodent datasets.
ggsave(filename = file.path(results_directory, "confusion_matrix.png"), plot = confusion_plot, width = 8.5, height = 7.5, dpi = 300)

# 11. Classify the optional observed image using exactly the same decoding, resizing and model rescaling.
if (nzchar(observed_image_file)) { # The observed prediction is optional for all supported image datasets.
    observed_data <- reticulate::py_to_r(load_image(tf$constant(observed_image_file, dtype = tf$string), tf$constant(0L, dtype = tf$int32)))
    observed_batch <- tf$expand_dims(observed_data[[1L]], axis = 0L) # Add the single-image batch dimension.
    observed_probabilities <- as.numeric(predict(object = best_model, x = observed_batch, verbose = 0L))
    observed_predictions <- data.frame(class = classes, model = unname(model_definitions[classes]), probability = observed_probabilities, percentage = 100 * observed_probabilities, stringsAsFactors = FALSE)
    observed_predictions <- observed_predictions[order(observed_predictions$probability, decreasing = TRUE), , drop = FALSE]
    rownames(observed_predictions) <- NULL # Print the class scores in decreasing softmax order.
    cat("\nObserved image classification:", observed_image_file, "\n") # Identify the image used for this prediction.
    print(observed_predictions) # These are CNN softmax scores, not Akaike weights from POPAN model selection.
    write.csv(observed_predictions, file.path(results_directory, "observed_predictions.csv"), row.names = FALSE)
} else {
    message("No observed image was selected, training and validation results were saved without this optional prediction.")
}
cat("\nTraining completed. Results directory:", normalizePath(results_directory, winslash = "/"), "\n")
