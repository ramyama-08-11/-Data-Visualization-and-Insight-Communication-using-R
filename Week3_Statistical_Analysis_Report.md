---
title: "Week 3: Statistical Analysis and Predictive Modeling"
subtitle: "Logistic regression using the Wisconsin Diagnostic Breast Cancer dataset"
author: "R Statistical Analysis Project"
date: "3 October 2026"
---

# Executive Summary

This project uses the public Wisconsin Breast Cancer (Diagnostic) dataset to examine whether four mean nuclear-image measurements distinguish malignant from benign diagnoses. The target is malignant (M); benign (B) is the reference class. A logistic regression model was evaluated with a stratified 80/20 train/holdout split. Five repeats of stratified five-fold cross-validation were run on training data only; the holdout was left untouched until final evaluation.

On the independent 115-case holdout, the model achieved 93.9% accuracy, 93.0% sensitivity, 94.4% specificity, and an ROC AUC of 0.994 at a fixed 0.50 probability threshold. It made seven errors: three malignant cases were classified as benign, and four benign cases were classified as malignant. Repeated cross-validation estimated mean accuracy of 93.5% (SD 2.2 percentage points) and mean AUC of 0.983 (SD 0.010). These results indicate strong discrimination on this benchmark, not clinical readiness. The dataset is small, historic, and derived from digitized images; external validation, decision-specific threshold selection, and clinical review would be required before any real-world use.

# Dataset and Rationale

The Breast Cancer Wisconsin (Diagnostic) dataset is hosted by the UCI Machine Learning Repository. UCI reports 569 observations, 30 real-valued predictors, no missing values, and a categorical diagnosis target. The predictors describe characteristics of cell nuclei calculated from digitized fine-needle-aspirate images. Ten base measurements (including radius, texture, smoothness, and concave points) are represented as mean, standard error, and worst values. The response has 357 benign and 212 malignant cases.

This dataset suits a supervised classification exercise because it is openly accessible, has a clearly defined binary outcome, continuous measurements suitable for distribution and correlation analysis, and enough cases for a stratified holdout plus repeated cross-validation. It is also a known limitation: it is a compact benchmark rather than a contemporary, representative clinical sample. The task is educational pattern classification, not diagnosis or treatment guidance.

The model uses four pre-specified mean measurements: radius, texture, smoothness, and concave points. Limiting the model to four fields keeps the analysis interpretable and avoids fitting every correlated feature variant. This is a teaching choice, not a claim that these are the uniquely optimal predictors. The dataset's ID column is retained for record tracking and excluded from predictors.

**Source and attribution:** Wolberg, W., Mangasarian, O., Street, N., & Street, W. (1993). *Breast Cancer Wisconsin (Diagnostic)* [Dataset]. UCI Machine Learning Repository. [https://doi.org/10.24432/C5DW2B](https://doi.org/10.24432/C5DW2B). UCI lists the data under CC BY 4.0. The R script downloads the official archive when `data/wdbc.data` is absent.

# Questions and Hypotheses

**Primary inferential question:** Does the mean concave-points measurement differ between malignant and benign cases in the training sample?

- Null hypothesis, $H_0$: the population mean is equal for malignant and benign groups, $\mu_M - \mu_B = 0$.
- Alternative hypothesis, $H_A$: the population means differ, $\mu_M - \mu_B \ne 0$.
- Significance level: $\alpha = 0.05$.

**Predictive question:** How well can a four-predictor logistic regression discriminate malignant from benign cases on observations not used for fitting? Malignant is treated as the positive class. The probability threshold is fixed at 0.50 for the reported confusion matrix; the ROC AUC evaluates ranking across thresholds.

# Reproducible Workflow

The analysis uses base R only (`stats`, `graphics`, and other standard R packages). It sets seed `20261003`, checks the expected 569-by-32 source data shape and absence of missing values, then performs class-wise random sampling. The resulting split contains 454 training observations (285 benign, 169 malignant) and 115 holdout observations (72 benign, 43 malignant). Cross-validation folds are stratified within the training set. Each fold's model is fitted without its assessment fold, preventing fold-level information leakage. The independent holdout is used only after cross-validation and final model fitting.

The overall workflow is:

1. Download/read and validate the UCI data; encode diagnosis as a two-level factor.
2. Create reproducible class-stratified training and holdout sets.
3. Summarize the four selected predictors and compare one pre-specified feature between diagnosis groups.
4. Examine within-class normality, Pearson and Spearman associations, and predictor VIFs.
5. Estimate a binomial logistic regression using repeated stratified cross-validation on training data.
6. Refit on all training data, then report discrimination, probability error, and classification errors on the holdout.
7. Export run logs, CSV tables, diagnostic plots, and code/output evidence images.

The analysis reports unadjusted p-values for the six exploratory pairwise correlation tests. These tests describe the selected predictors; they are not six independent confirmatory hypotheses. The Welch comparison is the single pre-specified inferential test.

# Exploratory Analysis and Hypothesis Test

Training-set summaries show that the selected fields have different scales and shapes. For example, mean radius was 14.08 (SD 3.57; median 13.28), while mean concave points was 0.0493 (SD 0.0392; median 0.0346). The selected predictors' distributions overlap to differing degrees: radius and concave points show visible separation by diagnosis, whereas texture and smoothness overlap more.

![Training-set distributions of the four model predictors, grouped by diagnosis.](evidence/feature_distributions.png)

The Shapiro-Wilk test rejected normality for mean concave points in both training groups: benign $W=0.947$, $p=1.18\times10^{-8}$; malignant $W=0.961$, $p=1.07\times10^{-4}$. The histograms also show skew and possible extreme values. The primary comparison therefore uses Welch's two-sample t-test, which does not assume equal group variances. Its mean-based inference still relies on independent observations and can be sensitive to strong outliers; the sizable group samples support approximate sampling-distribution reasoning, but a rank-based or bootstrap sensitivity analysis would strengthen the conclusion.

The training means were 0.0261 for benign and 0.0884 for malignant cases. The estimated malignant-minus-benign difference was 0.0623 measurement units (95% CI 0.0566 to 0.0680), Welch $t=-21.67$ for the test's benign-minus-malignant ordering, $df=212.32$, $p<2.2\times10^{-16}$. At $\alpha=0.05$, reject $H_0$: this sample provides strong evidence of a difference in mean concave-points measurements. This association alone does not establish causation or clinical utility.

![R console evidence: descriptive summaries, hypothesis test, validation metrics, and holdout confusion matrix.](evidence/r_console_results.png)

# Correlation and Predictor Checks

Pearson and Spearman tests were used because the predictors are continuous and some distributions are skewed. Radius and concave points were strongly positively associated (Pearson $r=0.825$, Spearman $\rho=0.762$). Smoothness and concave points had a moderate positive association (Pearson $r=0.548$), while texture and smoothness were nearly uncorrelated (Pearson $r=-0.050$). The rank-based results broadly support the direction of the linear correlations.

![Training-set Pearson correlation matrix for the four selected predictors.](evidence/correlation_matrix.png)

Variance inflation factors were 4.77 for radius, 1.16 for texture, 2.28 for smoothness, and 6.82 for concave points. The radius/concave-points association merits attention: coefficient estimates can be less stable even when predictions remain useful. This is why interpretation focuses on holdout predictive performance rather than treating individual coefficients as causal effects. Correlation and VIF are screening diagnostics, not automatic pass/fail tests.

# Model Specification and Validation

The model is a binomial generalized linear model with the logit link:

$$
\log\left(\frac{P(M)}{1-P(M)}\right) = \beta_0 + \beta_1(\text{radius mean}) + \beta_2(\text{texture mean}) + \beta_3(\text{smoothness mean}) + \beta_4(\text{concave points mean}).
$$

The model is fitted with base R's binomial `glm`:

```r
model_formula <- reformulate(predictors, response = "diagnosis")
fit_model <- function(training_data) {
	glm(model_formula, data = training_data, family = binomial())
}
```

The positive class is malignant. A 0.50 cutoff is a transparent default for the confusion matrix, not an optimized or clinically justified decision threshold. Accuracy is accompanied by sensitivity, specificity, precision, F1, balanced accuracy, AUC, and Brier score so class imbalance and probability error are visible. AUC is calculated from positive-class ranks with average ranks for ties.

Repeated stratified five-fold cross-validation (five repeats, 25 fold assessments) was conducted using only the 454 training cases. The reported mean and standard deviation across folds describe variation under these splits; repeated folds are not independent samples, so these standard deviations are descriptive and should not be read as formal confidence intervals.

Each class is assigned to folds separately, maintaining the diagnosis proportions in every assessment set:

```r
for (class_label in levels(train$diagnosis)) {
	class_rows <- which(train$diagnosis == class_label)
	fold_by_row[class_rows] <- sample(
		rep(seq_len(cv_folds), length.out = length(class_rows))
	)
}
```

After cross-validation, the final model is refitted on all training data and evaluated once on the holdout:

```r
final_model <- fit_model(train)
test_probability <- predict(final_model, newdata = test, type = "response")
test_metrics <- classification_metrics(test$diagnosis, test_probability)
```

| Metric | Cross-validation mean | Fold SD | Holdout |
|:--|--:|--:|--:|
| Accuracy | 0.935 | 0.022 | 0.939 |
| Sensitivity (malignant recall) | 0.897 | 0.060 | 0.930 |
| Specificity | 0.957 | 0.031 | 0.944 |
| Precision | 0.929 | 0.047 | 0.909 |
| F1 | 0.911 | 0.032 | 0.920 |
| Balanced accuracy | 0.927 | 0.027 | 0.937 |
| ROC AUC | 0.983 | 0.010 | 0.994 |
| Brier score | 0.050 | 0.014 | 0.032 |

On the holdout, the confusion matrix had 68 true negatives, 4 false positives, 3 false negatives, and 40 true positives. Thus the model detected 40 of 43 malignant cases and correctly classified 68 of 72 benign cases at the selected threshold. The high AUC indicates excellent ranking in this split, but only 43 positive cases contribute to holdout sensitivity, so its apparent precision should not obscure sampling uncertainty.

![ROC curve for the independent holdout; the dashed diagonal is chance-level ranking.](evidence/holdout_roc_curve.png)

![Independent holdout confusion matrix at a probability threshold of 0.50.](evidence/holdout_confusion_matrix.png)

# Diagnostics and Assumptions

Logistic regression does not require normally distributed predictors or normally distributed residuals. More relevant conditions include independent observations, a sufficiently appropriate log-odds form for continuous predictors, absence of complete separation, and manageable collinearity. The residual-versus-fitted plot displays the two expected outcome bands for binary data; the deviance-residual Q-Q plot shows tail departures, which are a prompt to inspect fit rather than a failed normal-residual assumption. Cook's-distance screening flags some observations above the reference line $4/n$; they should be checked for influence, not removed automatically. The ten-group holdout calibration plot gives a coarse probability check, but small-bin counts make it noisy.

![Training-model residual, influence, and holdout calibration diagnostics.](evidence/model_diagnostics.png)

The analysis does not formally test logit linearity, model separation, or calibration with a large independent sample. The modest predictor count and successful model fit support a useful baseline, not proof that all assumptions hold. The holdout is a single random split from the same source dataset, not external validation.

# Interpretation, Strengths, and Improvements

**Strengths:** The response, positive class, hypothesis, split, and threshold are explicit; the holdout remains separate from fitting and cross-validation; the fold construction is stratified and reproducible; performance is reported using threshold-dependent and threshold-independent metrics; and the report includes distribution, association, residual, influence, calibration, ROC, and confusion-matrix evidence.

**Limitations and next steps:**

- The data are a historic, relatively small benchmark. Evaluate a genuinely external, more recent cohort before considering generalization.
- The fixed 0.50 cutoff treats error types symmetrically. Choose a threshold from an explicit cost or sensitivity requirement using training-only predictions, then evaluate once on untouched data.
- Investigate the three false negatives and four false positives; inspect records and data provenance before drawing case-level conclusions.
- Check continuous-predictor logit linearity (for example, restricted splines compared within nested cross-validation), separation, and calibration more formally.
- Standardize predictors and compare a penalized logistic model (ridge or elastic net) to address correlated predictors. Keep feature selection and tuning inside the resampling loop to avoid optimistic estimates.
- Add bootstrap intervals for holdout metrics and a Mann-Whitney or bootstrap sensitivity analysis for the skewed concave-points comparison. Consider multiplicity adjustment if correlation tests become confirmatory.
- Compare a parsimonious model with the full 30-feature set only under the same nested resampling design. Report threshold selection, calibration, and interpretability trade-offs rather than choosing by accuracy alone.
- This dataset is not a clinical decision system. Predictions must not replace pathology, clinician judgment, or validated clinical procedures.

# Reproducibility and Evidence

Run the analysis from the project root in a current R installation:

```r
source("analyze_wdbc.R")
```

Or from a terminal:

```powershell
& "$env:LOCALAPPDATA\Programs\R\R-4.6.1\bin\Rscript.exe" analyze_wdbc.R
```

The script uses base R and writes machine-readable tables and the full captured output to `outputs/`, and plots plus source/output evidence panels to `evidence/`. The included PNG evidence panels are generated directly by R from the executed source and captured console log; they are code/output image exports, not photographs of the RStudio interface. `outputs/session_info.txt` records the runtime; this run used R 4.6.1 on Windows 11 x64.

![Selected lines from the executed R source.](evidence/r_code_excerpt.png)

The complete reproducible analysis is provided separately in `analyze_wdbc.R`; the full run log is in `outputs/analysis_console_output.txt`.

# Suggested 32-Hour Learning Plan

This is a suggested allocation for completing and extending the assignment; it is not a claim that 32 hours were spent producing this report.

| Work block | Hours | Learning activity |
|:--|--:|:--|
| Dataset and research question | 3 | Review the UCI data dictionary, provenance, target definition, and limitations; define the primary hypothesis and prediction goal. |
| R setup and data audit | 4 | Reproduce the download/import, inspect dimensions, labels, missingness, duplicates, and feature definitions. |
| Exploratory statistics | 6 | Build summaries and plots; compare distributions; assess outliers, normality, and effect-size interpretation. |
| Hypothesis testing and assumptions | 5 | Justify Welch testing; inspect correlation, rank association, variance, and multiplicity; run a nonparametric sensitivity analysis. |
| Modeling and resampling | 6 | Fit logistic regression, implement stratified cross-validation, define positive-class metrics, and compare threshold choices without holdout leakage. |
| Diagnostics and improvements | 4 | Inspect residuals, influence, calibration, logit form, and penalized-model alternatives. |
| Documentation and reproducibility | 4 | Rerun from a clean session, verify outputs, explain findings and caveats, and prepare the final report. |
| **Total** | **32** | |

# References

1. Wolberg, W., Mangasarian, O., Street, N., & Street, W. (1993). *Breast Cancer Wisconsin (Diagnostic)* [Dataset]. UCI Machine Learning Repository. [https://doi.org/10.24432/C5DW2B](https://doi.org/10.24432/C5DW2B).
2. UCI Machine Learning Repository. *Breast Cancer Wisconsin (Diagnostic): dataset page and variable information*. [https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic](https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic).
3. Street, W. N., Wolberg, W. H., & Mangasarian, O. L. (1993). Nuclear feature extraction for breast tumor diagnosis. In *Biomedical Image Processing and Biomedical Visualization* (pp. 861-870). SPIE. [https://doi.org/10.1117/12.148698](https://doi.org/10.1117/12.148698).