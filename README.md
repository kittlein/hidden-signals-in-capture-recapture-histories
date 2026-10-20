# Deep learning reveals hidden signals in capture-recapture histories: addressing data scarcity and computational complexity

Run the scripts from the directory containing the input files, or set their full paths at the beginning of each script. The examples use `generate_akodon.R`, `generate_phyllotis.R`, `train_CNN_rodents.R` and `evaluate_RMark_rodents.R` as the filenames. If the shared generator is named `generate_rodents.R`, use that filename in its execution command.

## Generate Akodon with ten sessions

Use `generate_akodon.R` with `Rodents_Necochea.xlsx`. Set the species and workbook at the beginning of the script:

```r
species <- "Akodon azarae" # Read the Akodon worksheet and select its population range, capture range and master seed.
workbook_file <- "Rodents_Necochea.xlsx" # Use the actual workbook filename, or replace this value with its full path.
n_occasions <- 10L # Use the ten sampling sessions recorded in the original workbook calendar for the simulations.
population_multiplier <- 1L # Keep each sampled base superpopulation size unchanged for the ten-session design.
images_per_class <- 5000L # Set the required number of matching image and capture-history datasets per demographic model.
```

The script reads observed histories from `Akodon azarae` and the calendar from `sessions`. The species sheets require `ID_individual`, `freq`, `sex` and `ch`. Histories must be ten-digit binary strings stored as text, with one row per individual and `freq = 1`. The calendar requires `session`, `time_interval` and `season`, with sessions numbered 1 to 10, a missing first interval and positive subsequent intervals in days.

For Akodon, the base superpopulation size is sampled from 280 to 440 and the baseline capture probability from 0.18 to 0.34. Both Akodon and Oligoryzomys simulations assign sex independently with probability 0.5 for each sex.

Captures use continuous lines for `cold` and dashed lines for `warm`, both one pixel thick. Females are red and males blue. The required R packages are `readxl`, `png`, `jsonlite` and `digest`.

Run in CMD, PowerShell or a terminal:

```powershell
Rscript generate_akodon.R
```

Images and matching histories are saved in `akodon_R_10_sessions_Nx1`.

## Settings for the other Akodon and Oligoryzomys designs

**Akodon: thirty sessions and N times three.** Keep the Akodon species setting and change:

```r
n_occasions <- 30L # Repeat the original ten-session calendar three times to obtain thirty simulated sampling sessions.
population_multiplier <- 3L # Multiply each sampled base superpopulation size by three before assigning individual sex.
```

The same workbook supplies the base calendar. The median original interval connects successive ten-session blocks. This expands the simulation design while the observed workbook histories remain ten-session observations. Output is saved in `akodon_R_30_sessions_Nx3`.

**Oligoryzomys: ten sessions.** Change the species setting and retain the original design:

```r
species <- "Oligoryzomys flavescens" # Read the Oligoryzomys worksheet and select its population range, capture range and seed.
n_occasions <- 10L # Use the common ten-session workbook calendar for the original Oligoryzomys simulation design.
population_multiplier <- 1L # Keep the sampled Oligoryzomys superpopulation size unchanged in the ten-session design.
```

The script automatically uses the `Oligoryzomys flavescens` sheet, N from 180 to 280, baseline capture probability from 0.065 to 0.13 and master seed 20260920. Output is saved in `oligoryzomys_R_10_sessions_Nx1`.

## Generate Phyllotis with 370 sessions

Use `generate_phyllotis.R` with `PD_data_final.RData` and a compatible `PD_top5_fitted_profiles.json`. The input data correspond to [Kelt et al. (2025)](https://doi.org/10.1002/ece3.72248). Download the data and original R code from the [Figshare repository cited in the paper](https://figshare.com/s/cf634e2c9e10d8855da3).

The JSON is a derived file containing fitted parameter profiles and model-selection information for the five leading POPAN models for *Phyllotis* reported in Table 3 of that paper, mapped to classes `p1` to `p5`. It must correspond to the same input RData. Set their filenames or full paths at the beginning of the script:

```r
data_file <- "PD_data_final.RData" # Read the observed histories, monthly covariates and fixed-zero capture occasions.
profiles_file <- "PD_top5_fitted_profiles.json" # Read the compatible fitted profiles associated with this input RData.
output_directory <- "phyllotis_fitted_370_R" # Save the Phyllotis images, histories and metadata in this experiment directory.
images_per_class <- 5000L # Generate this many matching images and capture-history datasets for each fitted demographic model.
population_scale <- 1.0 # Preserve the fitted female and male superpopulation sizes before rounding to integer individuals.
```

If the compatible JSON is available, `PD.final.results.RData` is not required. Otherwise, provide that RData and set `results_file` to its filename or full path. The script reconstructs the profiles from the saved RMark results, checks them against the MARK real-parameter estimates and creates the JSON cache. It does not refit the models.

The input RData must contain `PD.process`, `PD.ddl`, `PD.fixed` and `cov2`. The generator uses 370 monthly occasions and the fitted point estimates of survival, capture, entry and sex-specific superpopulation size for each model. Coefficient uncertainty is not sampled. Capture probabilities remain zero in the occasions specified by `PD.fixed`.

Each image represents ordered individual capture histories on a black RGB raster. Female captures contribute to the red channel and male captures to the blue channel. Capture events sharing a pixel are accumulated by sex and converted to color intensity. The required R packages are `png`, `jsonlite` and `digest`.

Run:

```powershell
Rscript generate_phyllotis.R
```

Output is saved in `phyllotis_fitted_370_R`. Each simulated image has a matching `.csv.gz` table containing individual histories, sex, capture counts, first and last captures, and true entry and exit sessions. The observed image and its original histories are also exported.

## Train the convolutional neural network

Use `train_CNN_rodents.R` to train a separate CNN for each species and sampling design. Set `dataset_directory` to the selected dataset root. The script also accepts its `images` subdirectory. It reads PNG files from `images/p1` to `images/p5` and reads the class meanings from `class_definitions.csv` when that file is available.

| Dataset | `dataset_directory` |
| --- | --- |
| Akodon, ten sessions | `akodon_R_10_sessions_Nx1` |
| Akodon, thirty sessions and N times three | `akodon_R_30_sessions_Nx3` |
| Oligoryzomys, ten sessions | `oligoryzomys_R_10_sessions_Nx1` |
| Phyllotis, 370 sessions | `phyllotis_fitted_370_R` |

For an Akodon training run, set the following at the beginning of the script:

```r
dataset_directory <- "akodon_R_10_sessions_Nx1" # Select the generated dataset containing images/p1 to images/p5.
python_executable <- "" # Leave empty for reticulate's Python, or enter the full path to the TensorFlow Python executable.
observed_image_file <- "" # Leave empty to find the observed PNG automatically, or provide the required full filename.
cfg <- list(cnn_width = 128L, cnn_height = 32L, resize_method = "area", seed = 123L, batch_size = 32L, epochs = 50L, patience = 8L, learning_rate = 3e-4, validation_fraction = 0.20)
```

The required R packages are `reticulate`, `keras3`, `tensorflow`, `tfdatasets` and `ggplot2`. The selected Python installation must provide TensorFlow and compatible Keras dependencies. Set `python_executable` before importing Python modules when a specific installation is needed.

Each class needs at least two simulated PNGs. The script reserves twenty percent within each class for validation, rounded down and limited so that at least one image remains in both subsets. Images are loaded as RGB, resized to 128 pixels wide by 32 pixels high using area interpolation, and scaled within the network. Early stopping monitors validation loss. The exported `split_manifest.csv` records the image path, class, integer label and training/validation membership.

Capture-history tables are not network inputs and are not required for training. Retain the generator's history folders and `capture_history_index.csv` for the paired RMark evaluation described below.

Run:

```powershell
Rscript train_CNN_rodents.R
```

Each run creates a new folder under `cnn_results`, for example `akodon_R_10_sessions_Nx1/cnn_results/run_20261007_120000`. Results include:

| Output | Contents |
| --- | --- |
| `best_model.keras` | Model selected by the lowest validation loss |
| `split_manifest.csv` and `dataset_distribution.csv` | Image assignment and class counts in training and validation |
| `class_definitions_used.csv` | Class labels and model meanings used in this run |
| `training_configuration.rds` and `runtime_information.txt` | Training settings and software information |
| `training_history.csv` and `training_history.png` | Training and validation loss and accuracy by epoch |
| `validation_evaluation.csv` | Overall validation loss and accuracy |
| `validation_predictions.csv` | Image paths, true classes, predicted classes and softmax scores |
| `confusion_counts.csv` and `confusion_proportions.csv` | Confusion counts and proportions normalized within each generating class |
| `confusion_matrix.png` and `validation_metrics_by_class.csv` | Confusion-matrix figure, recall and precision by class |
| `observed_predictions.csv` | Classification scores for the observed PNG, when available |

The ten-session Akodon and Oligoryzomys datasets and the Phyllotis dataset include an observed image. The expanded thirty-session Akodon design does not create an observed thirty-session image.

## Evaluate Akodon with RMark and compare it with the CNN

Use `evaluate_RMark_rodents.R` with the same generated experiment directory. The script fits POPAN models to the capture-history tables paired with the selected images. It reads the candidate formulas from `class_definitions.csv`, the calendar from `sessions.csv`, and the exact image/history associations from `capture_history_index.csv`. The current generators store image paths relative to `images` and history paths relative to `capture_histories`; the evaluator resolves both against the selected dataset root.

Install the R package and install [Program MARK](https://sites.warnercnr.colostate.edu/gwhite/program-mark/) separately:

```r
install.packages("RMark")
```

On Windows, set `mark_directory` to the directory containing `mark.exe`, `mark64.exe` or `mark32.exe` if the executable is not found automatically. Use forward slashes in R paths, for example `"C:/Program Files/Mark"`. Installing the R package alone does not install the MARK executable.

For the Akodon ten-session example, edit the settings at the beginning of the evaluation script:

```r
dataset_directory <- "akodon_R_10_sessions_Nx1" # Evaluate the histories generated for this Akodon image dataset.
evaluation_source <- "cnn_validation" # Restrict the paired comparison to validation images from one completed CNN run.
cnn_run_directory <- "" # Select the latest completed prediction export, or enter the full path to a specific CNN run.
histories_per_class <- 100L # Sample at most one hundred eligible histories per generating class, or use zero for all.
evaluation_seed <- 123L # Reproduce the selected validation sample without modifying its original CNN split membership.
evaluate_observed <- TRUE # Also evaluate the exported observed histories when this sampling design includes them.
require_all_models <- TRUE # Require usable AICc values for all five candidates before assigning a selected model class.
adjust_parameter_count <- TRUE # Use RMark's design-based parameter count and export the original MARK count as well.
use_initial_model <- TRUE # Initialize the candidates from a simple model fitted to the current capture histories.
mark_directory <- "" # Enter the MARK executable directory if it is not in PATH or a detected standard installation.
```

The script automatically selects the most recently modified CNN validation prediction export with a matching `split_manifest.csv`. Set `cnn_run_directory` explicitly to compare against a particular trained model. Its validation image paths must still resolve to the same generated files.

With the default limit, the script evaluates up to 100 histories from each class. If each class has at least 100 validation images, this gives 500 paired datasets and 2,500 candidate-model fits, plus 500 simple initialization fits when enabled. Set a small limit, such as `histories_per_class <- 2L`, for an initial execution check. Phyllotis fits can take substantially longer because their models contain 370 occasions and many time/year parameters.

For Akodon and Oligoryzomys, the exported model definitions are:

| Class | Survival `Phi` | Capture `p` | Entry `pent` | Superpopulation `N` |
| --- | --- | --- | --- | --- |
| `p1` | `~1` | `~1` | `~1` | `~1` |
| `p2` | `~season` | `~season` | `~1` | `~1` |
| `p3` | `~season` | `~1` | `~season` | `~1` |
| `p4` | `~1` | `~season` | `~season` | `~1` |
| `p5` | `~season+sex` | `~1` | `~1` | `~1` |

The actual CSV definitions are authoritative. Both sexes are processed as groups. Survival intervals use `interval_months` from the exported calendar. The script maps survival to the starting occasion of each interval, capture to the current occasion and entry to the subsequent occasion. Thus the evaluation uses the same seasonal timing and month units as generation.

For each history, all five candidates are fitted to the same `ch`, `freq` and `sex` records. Binary histories are read as text to preserve leading zeros. The initialization model uses those same records. Generating class labels, simulated population sizes, true entry/exit sessions and fitted generating profiles are not supplied to the likelihood or used to initialize the fits.

The candidate with the lowest usable AICc is the RMark prediction. The exports include delta AICc, Akaike weights, parameter counts, warning/error messages and reported singular parameters. With `require_all_models <- TRUE`, an incomplete candidate set produces `no_selection`, which remains in the full classification denominator. Finite AICc is not a convergence certificate; inspect the diagnostics, particularly when singular or boundary parameters are reported. The design-based count assumes the design is estimable; both adjusted and original MARK criteria are retained for review. Exact numerical AICc ties are flagged and resolved in `p1` to `p5` order.

The CNN confusion matrix is recalculated on exactly the sampled validation pairs evaluated by RMark. It therefore uses the same images, history tables and class counts as the RMark matrix, including tasks with unavailable RMark selections. CNN softmax scores and Akaike weights represent different quantities and are exported under separate names.

Run:

```powershell
Rscript evaluate_RMark_rodents.R
```

Each evaluation creates a new directory under `rmark_results`, for example `akodon_R_10_sessions_Nx1/rmark_results/run_20261007_130000`. The script exports progress after each history so completed comparisons remain available if execution is interrupted. A new execution starts a fresh run.

| Output | Contents |
| --- | --- |
| `evaluation_manifest.csv` | Selected image/history pairs, generating classes, replicates and evaluation IDs |
| `class_definitions_used.csv` and `sessions_used.csv` | Exact candidate formulas and calendar used for fitting |
| `evaluation_configuration.rds` and `runtime_information.txt` | Selection settings, model-count settings, intervals and software information |
| `model_selection.csv` | Five candidate rows per history with AICc, delta AICc, Akaike weights, parameter counts and diagnostics |
| `selection_summary.csv` and `simulation_predictions.csv` | Selected model, fit completeness, tie flags and simulation classification outcomes |
| `rmark_confusion_counts.csv` and `rmark_confusion_proportions.csv` | Classification counts and proportions with a sixth `no_selection` column |
| `classification_metrics_by_class.csv` | Class accuracy using all attempted histories, selected-only accuracy and no-selection counts |
| `classification_confusion_matrices.png` | RMark confusion matrix and paired CNN matrix when available |
| `cnn_predictions_selected.csv` and `cnn_rmark_paired_predictions.csv` | Selected CNN predictions and their corresponding RMark outcomes |
| `cnn_confusion_counts.csv` and `cnn_confusion_proportions.csv` | CNN results recomputed on the same selected pairs, with the same class denominator |
| `observed_model_selection.csv` and `observed_selection_summary.csv` | Observed-data model support and selected model, when observed histories exist |
| `fit_logs` | Captured fitting messages and a saved summary for each evaluated history |

CNN-specific files are produced only for `cnn_validation`. Observed histories are evaluated separately and do not contribute to a confusion matrix because their generating model is unknown.

## Use the same RMark script for the other generated designs

Change `dataset_directory` using the dataset table above. Keep the same evaluation script for Akodon with thirty sessions and N times three, Oligoryzomys with ten sessions, and Phyllotis with 370 sessions. No workbook, source RData or generating JSON is needed during this evaluation; the exported history tables, model definitions, index and calendar supply the inputs.

**Akodon, thirty sessions and N times three.** The script uses all thirty exported sessions, including the intervals connecting repeated ten-session blocks. It fits population size from the histories without using the simulated true N. There is no observed thirty-session history, so observed evaluation is skipped.

**Oligoryzomys, ten sessions.** Set `dataset_directory <- "oligoryzomys_R_10_sessions_Nx1"`. The exported formulas have the same structure as the Akodon candidates, and the script reads the Oligoryzomys histories and calendar from that experiment directory.

**Phyllotis, 370 sessions.** Set `dataset_directory <- "phyllotis_fitted_370_R"`. The script reads these formulas from the generator output:

| Class | Survival `Phi` | Capture `p` | Entry `pent` | Superpopulation `N` |
| --- | --- | --- | --- | --- |
| `p1` | `~yr*season_rain2` | `~time` | `~yr*season_rain2` | `~sex` |
| `p2` | `~yr*season2` | `~time` | `~yr*season_rain2` | `~sex` |
| `p3` | `~yr+season_rain2` | `~time` | `~yr*season_rain2` | `~sex` |
| `p4` | `~yr+season2` | `~time` | `~yr*season_rain2` | `~sex` |
| `p5` | `~yr+wet_dry2` | `~time` | `~yr*season_rain2` | `~sex` |

The exported year is treated as the `yr` factor. The script retains the `season2`, `season_rain2` and `wet_dry2` factors and the monthly intervals from `time_interval_months`. Capture probability is fixed at zero in both sexes for every occasion marked by `p_fixed_zero`, in the initialization model and all five candidates. Negative frequencies in original observed histories are preserved for MARK's trap-death handling. Simulation-only entry and exit columns remain outside the fit.

To evaluate generated histories before training a CNN, change only the selection mode:

```r
evaluation_source <- "generated_sample" # Sample the generator's indexed histories directly without requiring CNN outputs.
```

This mode still creates the RMark ranking, diagnostics, classifications and observed-data evaluation. It does not create a paired CNN comparison. To compare with the CNN afterwards, run a separate evaluation with `evaluation_source <- "cnn_validation"` and the required completed CNN run.

For the fitting API and parameter-count options, see the [RMark documentation](https://sjbonner.r-universe.dev/RMark/doc/manual.html).
