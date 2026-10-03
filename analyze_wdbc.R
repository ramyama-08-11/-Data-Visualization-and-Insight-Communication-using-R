options(stringsAsFactors = FALSE, width = 110)

data_dir <- file.path("data")
output_dir <- file.path("outputs")
evidence_dir <- file.path("evidence")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(evidence_dir, showWarnings = FALSE, recursive = TRUE)

data_path <- file.path(data_dir, "wdbc.data")
if (!file.exists(data_path)) {
  archive_path <- tempfile(fileext = ".zip")
  download.file(
    "https://archive.ics.uci.edu/static/public/17/breast+cancer+wisconsin+diagnostic.zip",
    archive_path,
    mode = "wb",
    quiet = TRUE
  )
  unzip(archive_path, files = "wdbc.data", exdir = data_dir)
  unlink(archive_path)
}

measurements <- c(
  "radius", "texture", "perimeter", "area", "smoothness", "compactness",
  "concavity", "concave_points", "symmetry", "fractal_dimension"
)
feature_names <- c(
  paste0(measurements, "_mean"),
  paste0(measurements, "_se"),
  paste0(measurements, "_worst")
)
column_names <- c("id", "diagnosis", feature_names)
data <- read.csv(data_path, header = FALSE, col.names = column_names)
data$diagnosis <- factor(data$diagnosis, levels = c("B", "M"))
predictors <- c("radius_mean", "texture_mean", "smoothness_mean", "concave_points_mean")

if (anyNA(data) || nrow(data) != 569L || ncol(data) != 32L) {
  stop("The downloaded data failed its expected completeness or dimension checks.")
}

set.seed(20261003)
train_rows <- unlist(lapply(levels(data$diagnosis), function(class_label) {
  class_rows <- which(data$diagnosis == class_label)
  sample(class_rows, floor(0.80 * length(class_rows)))
}))
train <- data[train_rows, ]
test <- data[-train_rows, ]

summary_rows <- lapply(predictors, function(variable) {
  values <- train[[variable]]
  data.frame(
    feature = variable,
    n = length(values),
    mean = mean(values),
    sd = sd(values),
    median = median(values),
    q1 = unname(quantile(values, 0.25)),
    q3 = unname(quantile(values, 0.75)),
    min = min(values),
    max = max(values)
  )
})
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, file.path(output_dir, "training_summary_statistics.csv"), row.names = FALSE)

primary_test <- t.test(concave_points_mean ~ diagnosis, data = train, var.equal = FALSE)
normality_tests <- do.call(rbind, lapply(levels(train$diagnosis), function(class_label) {
  result <- shapiro.test(train$concave_points_mean[train$diagnosis == class_label])
  data.frame(group = class_label, statistic = unname(result$statistic), p_value = result$p.value)
}))
hypothesis_table <- data.frame(
  test = "Welch two-sample t-test: mean concave points, M versus B",
  estimate_B = unname(primary_test$estimate[1]),
  estimate_M = unname(primary_test$estimate[2]),
  difference_M_minus_B = unname(primary_test$estimate[2] - primary_test$estimate[1]),
  confidence_low = -primary_test$conf.int[2],
  confidence_high = -primary_test$conf.int[1],
  p_value = primary_test$p.value
)
write.csv(hypothesis_table, file.path(output_dir, "hypothesis_test.csv"), row.names = FALSE)
write.csv(normality_tests, file.path(output_dir, "normality_tests.csv"), row.names = FALSE)

correlation_table <- do.call(rbind, lapply(combn(predictors, 2, simplify = FALSE), function(pair) {
  pearson <- cor.test(train[[pair[1]]], train[[pair[2]]], method = "pearson")
  spearman <- cor.test(train[[pair[1]]], train[[pair[2]]], method = "spearman", exact = FALSE)
  data.frame(
    feature_1 = pair[1], feature_2 = pair[2],
    pearson_r = unname(pearson$estimate), pearson_p = pearson$p.value,
    spearman_rho = unname(spearman$estimate), spearman_p = spearman$p.value
  )
}))
write.csv(correlation_table, file.path(output_dir, "correlation_tests.csv"), row.names = FALSE)

vif_values <- sapply(predictors, function(target) {
  other_predictors <- setdiff(predictors, target)
  auxiliary_model <- lm(reformulate(other_predictors, response = target), data = train)
  1 / (1 - summary(auxiliary_model)$r.squared)
})
write.csv(data.frame(feature = names(vif_values), VIF = unname(vif_values)),
          file.path(output_dir, "variance_inflation_factors.csv"), row.names = FALSE)

model_formula <- reformulate(predictors, response = "diagnosis")
fit_model <- function(training_data) {
  glm(model_formula, data = training_data, family = binomial())
}
auc_rank <- function(actual, probability) {
  positive <- actual == "M"
  n_positive <- sum(positive)
  n_negative <- sum(!positive)
  if (n_positive == 0L || n_negative == 0L) return(NA_real_)
  positive_ranks <- sum(rank(probability, ties.method = "average")[positive])
  (positive_ranks - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)
}
classification_metrics <- function(actual, probability, threshold = 0.5) {
  predicted <- factor(ifelse(probability >= threshold, "M", "B"), levels = c("B", "M"))
  actual <- factor(actual, levels = c("B", "M"))
  confusion <- table(actual = actual, predicted = predicted)
  tn <- unname(confusion["B", "B"])
  fp <- unname(confusion["B", "M"])
  fn <- unname(confusion["M", "B"])
  tp <- unname(confusion["M", "M"])
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  precision <- if (tp + fp == 0) NA_real_ else tp / (tp + fp)
  c(
    accuracy = (tp + tn) / sum(confusion), sensitivity = sensitivity,
    specificity = specificity, precision = precision,
    F1 = 2 * precision * sensitivity / (precision + sensitivity),
    balanced_accuracy = (sensitivity + specificity) / 2,
    AUC = auc_rank(actual, probability), Brier = mean((as.integer(actual == "M") - probability)^2)
  )
}

cv_folds <- 5L
cv_repeats <- 5L
cv_results <- vector("list", cv_folds * cv_repeats)
result_index <- 1L
for (repeat_index in seq_len(cv_repeats)) {
  fold_by_row <- integer(nrow(train))
  for (class_label in levels(train$diagnosis)) {
    class_rows <- which(train$diagnosis == class_label)
    fold_by_row[class_rows] <- sample(rep(seq_len(cv_folds), length.out = length(class_rows)))
  }
  for (fold_index in seq_len(cv_folds)) {
    assessment <- train[fold_by_row == fold_index, ]
    analysis <- train[fold_by_row != fold_index, ]
    fold_model <- fit_model(analysis)
    fold_probability <- predict(fold_model, newdata = assessment, type = "response")
    cv_results[[result_index]] <- data.frame(
      repeat_id = repeat_index,
      fold = fold_index,
      as.list(classification_metrics(assessment$diagnosis, fold_probability))
    )
    result_index <- result_index + 1L
  }
}
cv_table <- do.call(rbind, cv_results)
write.csv(cv_table, file.path(output_dir, "cross_validation_metrics.csv"), row.names = FALSE)
metric_names <- names(classification_metrics(train$diagnosis, rep(0.5, nrow(train))))
cv_summary <- do.call(rbind, lapply(metric_names, function(metric) {
  data.frame(metric = metric, mean = mean(cv_table[[metric]], na.rm = TRUE),
             sd = sd(cv_table[[metric]], na.rm = TRUE))
}))
write.csv(cv_summary, file.path(output_dir, "cross_validation_summary.csv"), row.names = FALSE)

final_model <- fit_model(train)
test_probability <- predict(final_model, newdata = test, type = "response")
test_metrics <- classification_metrics(test$diagnosis, test_probability)
test_prediction <- factor(ifelse(test_probability >= 0.5, "M", "B"), levels = c("B", "M"))
test_confusion <- table(actual = factor(test$diagnosis, levels = c("B", "M")), predicted = test_prediction)
test_predictions <- data.frame(
  id = test$id, actual = test$diagnosis, probability_M = test_probability,
  predicted = test_prediction
)
write.csv(test_predictions, file.path(output_dir, "holdout_predictions.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(test_confusion), file.path(output_dir, "holdout_confusion_matrix.csv"))
write.csv(data.frame(metric = names(test_metrics), value = unname(test_metrics)),
          file.path(output_dir, "holdout_metrics.csv"), row.names = FALSE)
write.csv(data.frame(term = names(coef(final_model)), estimate = unname(coef(final_model)),
                     standard_error = sqrt(diag(vcov(final_model))),
                     p_value = summary(final_model)$coefficients[, "Pr(>|z|)"]),
          file.path(output_dir, "logistic_coefficients.csv"), row.names = FALSE)

capture.output({
  cat("UCI WDBC analysis; R", as.character(getRversion()), "\n")
  cat("Rows:", nrow(data), "| training:", nrow(train), "| holdout:", nrow(test), "\n")
  cat("Class counts (B, M):\n")
  print(table(data$diagnosis))
  cat("\nTraining summary statistics:\n")
  print(summary_table, row.names = FALSE, digits = 4)
  cat("\nWelch hypothesis test:\n")
  print(primary_test)
  cat("\nShapiro-Wilk tests (exploratory):\n")
  print(normality_tests, row.names = FALSE, digits = 4)
  cat("\nCandidate-predictor correlations:\n")
  print(correlation_table, row.names = FALSE, digits = 4)
  cat("\nVariance inflation factors:\n")
  print(vif_values, digits = 4)
  cat("\nRepeated stratified 5-fold CV (5 repeats), mean and SD:\n")
  print(cv_summary, row.names = FALSE, digits = 4)
  cat("\nFinal model:\n")
  print(summary(final_model))
  cat("\nIndependent holdout metrics (threshold 0.50):\n")
  print(test_metrics, digits = 4)
  cat("\nHoldout confusion matrix (rows=actual, columns=predicted):\n")
  print(test_confusion)
  cat("\nSession information:\n")
  print(sessionInfo())
}, file = file.path(output_dir, "analysis_console_output.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

png(file.path(evidence_dir, "feature_distributions.png"), width = 1400, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (variable in predictors) {
  benign <- train[[variable]][train$diagnosis == "B"]
  malignant <- train[[variable]][train$diagnosis == "M"]
  bins <- pretty(range(c(benign, malignant)), n = 25)
  hist(benign, breaks = bins, col = adjustcolor("#2878B5", 0.6), border = "white",
       main = variable, xlab = variable, ylab = "Count")
  hist(malignant, breaks = bins, col = adjustcolor("#D9534F", 0.55), border = "white", add = TRUE)
  legend("topright", legend = c("Benign", "Malignant"), fill = c("#2878B5", "#D9534F"), bty = "n")
}
dev.off()

png(file.path(evidence_dir, "correlation_matrix.png"), width = 1100, height = 900, res = 150)
correlation_matrix <- cor(train[predictors])
par(mar = c(10, 12, 3, 2))
image(seq_along(predictors), seq_along(predictors), correlation_matrix,
      col = hcl.colors(50, "Blue-Red 3", rev = TRUE), zlim = c(-1, 1), axes = FALSE,
      xlab = "", ylab = "", main = "Training-predictor Pearson correlations")
axis(1, at = seq_along(predictors), labels = predictors, las = 2, cex.axis = 0.85)
axis(2, at = seq_along(predictors), labels = predictors, las = 2, cex.axis = 0.85)
text(rep(seq_along(predictors), each = length(predictors)),
     rep(seq_along(predictors), times = length(predictors)),
     labels = sprintf("%.2f", as.vector(correlation_matrix)), cex = 0.85)
dev.off()

roc_thresholds <- c(Inf, sort(unique(test_probability), decreasing = TRUE), -Inf)
roc_points <- do.call(rbind, lapply(roc_thresholds, function(threshold) {
  predicted_positive <- test_probability >= threshold
  actual_positive <- test$diagnosis == "M"
  c(FPR = sum(predicted_positive & !actual_positive) / sum(!actual_positive),
    TPR = sum(predicted_positive & actual_positive) / sum(actual_positive))
}))
roc_points <- roc_points[order(roc_points[, "FPR"], roc_points[, "TPR"]), , drop = FALSE]
png(file.path(evidence_dir, "holdout_roc_curve.png"), width = 1000, height = 850, res = 150)
plot(roc_points[, "FPR"], roc_points[, "TPR"], type = "l", lwd = 3, col = "#D9534F",
    xlim = c(0, 1), ylim = c(0, 1), xlab = "False-positive rate",
     ylab = "True-positive rate", main = sprintf("Independent holdout ROC (AUC = %.3f)", test_metrics["AUC"]))
abline(0, 1, lty = 2, col = "#666666")
dev.off()

png(file.path(evidence_dir, "holdout_confusion_matrix.png"), width = 950, height = 800, res = 150)
par(mar = c(5, 7, 4, 2))
image(1:2, 1:2, test_confusion, col = hcl.colors(30, "BluYl"), axes = FALSE,
      xlab = "Predicted diagnosis", ylab = "Actual diagnosis", main = "Independent holdout, threshold 0.50")
axis(1, at = 1:2, labels = c("Benign", "Malignant"))
axis(2, at = 1:2, labels = c("Benign", "Malignant"), las = 1)
text(rep(1:2, each = 2), rep(1:2, times = 2), labels = as.vector(test_confusion), cex = 1.7)
dev.off()

png(file.path(evidence_dir, "model_diagnostics.png"), width = 1400, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(fitted(final_model), residuals(final_model, type = "deviance"),
     xlab = "Fitted probability", ylab = "Deviance residual", main = "Residuals vs fitted")
abline(h = 0, lty = 2, col = "#777777")
qqnorm(rstandard(final_model, type = "deviance"), main = "Standardized deviance residuals")
qqline(rstandard(final_model, type = "deviance"), col = "#D9534F")
plot(cooks.distance(final_model), type = "h", xlab = "Training row", ylab = "Cook's distance",
     main = "Influence screening")
abline(h = 4 / nrow(train), lty = 2, col = "#D9534F")
calibration_group <- cut(rank(test_probability, ties.method = "first"),
                         breaks = seq(0, length(test_probability), length.out = 11), include.lowest = TRUE)
calibration <- aggregate(data.frame(predicted = test_probability,
                                    observed = as.integer(test$diagnosis == "M")),
                         by = list(group = calibration_group), FUN = mean)
plot(calibration$predicted, calibration$observed, pch = 19, cex = 1.2,
     xlim = c(0, 1), ylim = c(0, 1), asp = 1, xlab = "Mean predicted risk",
     ylab = "Observed malignant fraction", main = "Holdout calibration (10 groups)")
abline(0, 1, lty = 2, col = "#777777")
dev.off()

cat("Analysis complete. Outputs: outputs/; figures: evidence/\n")

render_text_panel <- function(lines, path, title) {
  png(path, width = 2200, height = max(850, 27 * (length(lines) + 3)), res = 160)
  par(mar = c(0.2, 0.4, 1.8, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.015, 0.99, title, adj = c(0, 1), family = "sans", cex = 1.15)
  y_positions <- seq(0.94, 0.02, length.out = length(lines))
  text(0.015, y_positions, labels = lines, adj = c(0, 1), family = "mono", cex = 0.76)
  dev.off()
}

source_lines <- readLines("analyze_wdbc.R", warn = FALSE)
code_anchors <- grep(
  "set.seed\\(|primary_test <- t.test|pearson <- cor.test|fit_model <- function|auc_rank <- function|classification_metrics <- function|cv_folds <-|for \\(repeat_index|final_model <- fit_model|test_metrics <-",
  source_lines
)
source_end <- match("render_text_panel <- function(lines, path, title) {", source_lines) - 1L
code_anchors <- code_anchors[code_anchors < source_end]
code_indices <- sort(unique(unlist(lapply(code_anchors, function(index) {
  seq(max(1L, index - 1L), min(length(source_lines), index + 1L))
}))))
code_excerpt <- sprintf("%3d | %s", code_indices, source_lines[code_indices])
render_text_panel(code_excerpt, file.path(evidence_dir, "r_code_excerpt.png"),
                  "Selected source lines from analyze_wdbc.R")

console_lines <- readLines(file.path(output_dir, "analysis_console_output.txt"), warn = FALSE)
extract_console_block <- function(start, end) {
  first <- grep(start, console_lines, fixed = TRUE)[1]
  last <- grep(end, console_lines, fixed = TRUE)[1]
  if (is.na(first) || is.na(last) || last <= first) stop("A console evidence section was not found.")
  console_lines[first:(last - 1L)]
}
console_excerpt <- c(
  extract_console_block("Rows:", "Shapiro-Wilk tests"),
  extract_console_block("Repeated stratified", "Final model:"),
  extract_console_block("Independent holdout metrics", "Session information:")
)
console_excerpt <- console_excerpt[nzchar(console_excerpt)]
render_text_panel(console_excerpt, file.path(evidence_dir, "r_console_results.png"),
                  "Selected results from the R 4.6.1 run log")