# ===========================================================================
# US-007: BRCA Molecular Subtype Classification
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# 方法：Random Forest + LASSO + XGBoost for PAM50 subtype prediction
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(glmnet)
})

cat("\n========== US-007: Molecular Subtype Classification ==========\n\n")

INPUT_DIR  <- "D:/Users/Desktop/R_Work/data/processed"
FIG_DIR    <- "D:/Users/Desktop/R_Work/results/figures"
TBL_DIR    <- "D:/Users/Desktop/R_Work/results/tables"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load and prepare data ----
cat("Step 1: Loading and preparing classification data...\n")

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
tpm_mat <- brca$tpm
clinical <- brca$clinical

cat(sprintf("  TPM: %d genes x %d samples\n", nrow(tpm_mat), ncol(tpm_mat)))

# Map expression samples to clinical
expr_patients <- substr(colnames(tpm_mat), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_patients, clinical$pid_short)

subtypes <- clinical$molecular_subtype[match_idx]
names(subtypes) <- colnames(tpm_mat)

# Remove NA subtypes
valid_idx <- !is.na(subtypes)
cat(sprintf("  Samples with known subtype: %d / %d\n", sum(valid_idx), length(valid_idx)))

tpm_sub <- tpm_mat[, valid_idx]
subtypes_sub <- subtypes[valid_idx]

cat(sprintf("  Subtype distribution:\n"))
for (lv in names(table(subtypes_sub))) {
  cat(sprintf("    %s: %d\n", lv, table(subtypes_sub)[lv]))
}

# Remove subtypes with too few samples
min_samples <- 15
subtype_counts <- table(subtypes_sub)
keep_subtypes <- names(subtype_counts[subtype_counts >= min_samples])
keep_idx <- subtypes_sub %in% keep_subtypes

tpm_sub <- tpm_sub[, keep_idx]
labels <- factor(subtypes_sub[keep_idx])

cat(sprintf("  After filtering (>%d samples): %d samples, %d subtypes\n",
            min_samples, ncol(tpm_sub), length(levels(labels))))
cat(sprintf("  Subtypes: %s\n", paste(levels(labels), collapse = ", ")))

# ---- 2. Feature selection: top variable genes ----
cat("\nStep 2: Feature selection (top 500 variable genes)...\n")

log_tpm <- log2(tpm_sub + 1)
gene_var <- apply(log_tpm, 1, var)
top500 <- names(sort(gene_var, decreasing = TRUE))[1:500]

X <- t(log_tpm[top500, ])
y <- labels

cat(sprintf("  Feature matrix: %d samples x %d genes\n", nrow(X), ncol(X)))

# ---- 3. Train/Test split ----
cat("\nStep 3: Creating train/test split...\n")

set.seed(42)
train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
y_train <- y[train_idx]
y_test  <- y[-train_idx]

cat(sprintf("  Train: %d | Test: %d\n", length(y_train), length(y_test)))

# ---- 4. Random Forest ----
cat("\nStep 4: Training Random Forest classifier...\n")

set.seed(42)
rf_model <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 500,
  importance = TRUE
)

rf_pred <- predict(rf_model, X_test)
rf_cm <- confusionMatrix(rf_pred, y_test)
cat(sprintf("  RF Accuracy: %.4f | Kappa: %.4f\n",
            rf_cm$overall["Accuracy"], rf_cm$overall["Kappa"]))

# RF feature importance
rf_imp <- importance(rf_model)
rf_imp_df <- data.frame(
  Gene = rownames(rf_imp),
  MeanDecreaseGini = rf_imp[, "MeanDecreaseGini"],
  stringsAsFactors = FALSE
) %>% arrange(desc(MeanDecreaseGini))
write.csv(rf_imp_df, file.path(TBL_DIR, "classification_rf_features.csv"), row.names = FALSE)

# ---- 5. LASSO (multinomial) ----
cat("\nStep 5: Training LASSO classifier...\n")

set.seed(42)
cv_lasso <- cv.glmnet(
  x = X_train,
  y = y_train,
  family = "multinomial",
  alpha = 1,
  nfolds = 5
)

lasso_pred <- predict(cv_lasso, X_test, s = "lambda.min", type = "class")
lasso_pred <- factor(as.vector(lasso_pred), levels = levels(y_train))
lasso_cm <- confusionMatrix(lasso_pred, y_test)
cat(sprintf("  LASSO Accuracy: %.4f | Kappa: %.4f\n",
            lasso_cm$overall["Accuracy"], lasso_cm$overall["Kappa"]))

# LASSO selected features
lasso_coef <- coef(cv_lasso, s = "lambda.min")
lasso_genes <- unique(unlist(lapply(lasso_coef, function(x) {
  nz <- which(x[-1] != 0)
  if (length(nz) > 0) colnames(X_train)[nz]
})))
cat(sprintf("  LASSO selected %d genes\n", length(lasso_genes)))
write.csv(data.frame(Gene = lasso_genes), file.path(TBL_DIR, "classification_lasso_features.csv"), row.names = FALSE)

# ---- 6. XGBoost ----
cat("\nStep 6: Training XGBoost classifier...\n")

has_xgboost <- requireNamespace("xgboost", quietly = TRUE)

xgboost_acc <- NA
if (has_xgboost) {
  library(xgboost)

  # Convert to numeric labels
  y_train_num <- as.numeric(y_train) - 1
  y_test_num <- as.numeric(y_test) - 1

  dtrain <- xgb.DMatrix(data = X_train, label = y_train_num)
  dtest  <- xgb.DMatrix(data = X_test, label = y_test_num)

  params <- list(
    objective = "multi:softprob",
    num_class = length(levels(y_train)),
    max_depth = 6,
    eta = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8
  )

  set.seed(42)
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 100,
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0
  )

  xgb_pred_prob <- predict(xgb_model, dtest)
  xgb_pred_prob <- matrix(xgb_pred_prob, ncol = length(levels(y_train)), byrow = TRUE)
  xgb_pred_class <- levels(y_train)[max.col(xgb_pred_prob)]
  xgb_pred_class <- factor(xgb_pred_class, levels = levels(y_train))

  xgb_cm <- confusionMatrix(xgb_pred_class, y_test)
  xgboost_acc <- xgb_cm$overall["Accuracy"]
  cat(sprintf("  XGBoost Accuracy: %.4f | Kappa: %.4f\n",
              xgb_cm$overall["Accuracy"], xgb_cm$overall["Kappa"]))

  # XGBoost importance
  xgb_imp <- xgb.importance(model = xgb_model, feature_names = colnames(X_train))
  write.csv(xgb_imp, file.path(TBL_DIR, "classification_xgb_features.csv"), row.names = FALSE)
} else {
  cat("  xgboost not installed, skipping.\n")
}

# ---- 7. Model comparison ----
cat("\nStep 7: Generating model comparison...\n")

model_metrics <- data.frame(
  Model = c("Random Forest", "LASSO"),
  Accuracy = c(rf_cm$overall["Accuracy"], lasso_cm$overall["Accuracy"]),
  Kappa = c(rf_cm$overall["Kappa"], lasso_cm$overall["Kappa"]),
  F1_Macro = c(mean(rf_cm$byClass[, "F1"], na.rm = TRUE),
               mean(lasso_cm$byClass[, "F1"], na.rm = TRUE)),
  stringsAsFactors = FALSE
)

if (has_xgboost) {
  model_metrics <- rbind(model_metrics, data.frame(
    Model = "XGBoost",
    Accuracy = xgb_cm$overall["Accuracy"],
    Kappa = xgb_cm$overall["Kappa"],
    F1_Macro = mean(xgb_cm$byClass[, "F1"], na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
}

write.csv(model_metrics, file.path(TBL_DIR, "classification_metrics.csv"), row.names = FALSE)
cat("  Model metrics:\n")
print(model_metrics)

# ---- 8. ROC curves for RF (multi-class) ----
cat("\nStep 8: Generating multi-class ROC curves...\n")

rf_prob <- predict(rf_model, X_test, type = "prob")

pdf(file.path(FIG_DIR, "classification_roc.pdf"), width = 10, height = 8)

# One-vs-all ROC for each class
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
par(mfrow = c(2, 2))

for (i in seq_along(levels(y_train))) {
  class_name <- levels(y_train)[i]
  binary_truth <- (y_test == class_name)

  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(binary_truth, rf_prob[, class_name], quiet = TRUE)
    pROC::plot.roc(roc_obj, main = paste("ROC:", class_name),
                   col = colors[i], lwd = 2, legacy.axes = TRUE)
    legend("bottomright",
           sprintf("AUC = %.3f", pROC::auc(roc_obj)),
           bty = "n", cex = 1.2)
  }
}

par(mfrow = c(1, 1))
dev.off()
cat("  ROC curves saved.\n")

# ---- 9. Confusion matrix heatmap ----
cat("\nStep 9: Generating confusion matrix heatmaps...\n")

pdf(file.path(FIG_DIR, "classification_cm.pdf"), width = 14, height = 5)

# RF confusion matrix
rf_cm_table <- as.data.frame.matrix(rf_cm$table)
rf_cm_prop <- sweep(rf_cm_table, 1, rowSums(rf_cm_table), "/")

library(pheatmap)
par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))

# Using pheatmap for better visualization
pheatmap(as.matrix(rf_cm_prop),
         main = "Random Forest: Confusion Matrix (Proportion)",
         display_numbers = TRUE,
         number_format = "%.2f",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#377EB8"))(50))

dev.off()
cat("  Confusion matrix heatmap saved.\n")

cat("\n========== US-007 Complete: Classification ==========\n")
cat(sprintf("  Best Model: %s (Accuracy=%.3f)\n",
            model_metrics$Model[which.max(model_metrics$Accuracy)],
            max(model_metrics$Accuracy)))
