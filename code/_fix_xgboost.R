# Fix XGBoost: diagnose and correct the 35.8% accuracy issue
suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(glmnet)
  library(xgboost)
})

cat("=== XGBoost Diagnosis & Fix ===\n\n")

# Load same data
brca <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_tumor_processed.rds")
tpm_mat <- brca$tpm
clinical <- brca$clinical

expr_patients <- substr(colnames(tpm_mat), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_patients, clinical$pid_short)
subtypes <- clinical$molecular_subtype[match_idx]
names(subtypes) <- colnames(tpm_mat)

valid_idx <- !is.na(subtypes)
tpm_sub <- tpm_mat[, valid_idx]
subtypes_sub <- subtypes[valid_idx]

# Remove subtypes with <15 samples
min_samples <- 15
subtype_counts <- table(subtypes_sub)
keep_subtypes <- names(subtype_counts[subtype_counts >= min_samples])
keep_idx <- subtypes_sub %in% keep_subtypes
tpm_sub <- tpm_sub[, keep_idx]
labels <- factor(subtypes_sub[keep_idx])

cat(sprintf("Samples: %d, Subtypes: %s\n", length(labels), paste(levels(labels), collapse=", ")))
cat("Subtype distribution:\n")
print(table(labels))

# Feature selection
log_tpm <- log2(tpm_sub + 1)
gene_var <- apply(log_tpm, 1, var)
top500 <- names(sort(gene_var, decreasing = TRUE))[1:500]

X <- t(log_tpm[top500, ])
y <- labels

# Train/test split
set.seed(42)
train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
y_train <- y[train_idx]
y_test  <- y[-train_idx]

cat(sprintf("\nTrain: %d, Test: %d\n", length(y_train), length(y_test)))
cat("Test set distribution:\n")
print(table(y_test))

# === DIAGNOSIS 1: Check the original approach ===
cat("\n--- Original approach (multi:softprob + byrow=TRUE reshape) ---\n")
y_train_num <- as.numeric(y_train) - 1
y_test_num <- as.numeric(y_test) - 1

dtrain <- xgb.DMatrix(data = X_train, label = y_train_num)
dtest  <- xgb.DMatrix(data = X_test, label = y_test_num)

params_orig <- list(
  objective = "multi:softprob",
  num_class = 3,
  max_depth = 6,
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8
)

set.seed(42)
xgb_orig <- xgb.train(params = params_orig, data = dtrain, nrounds = 100, verbose = 0)
pred_prob <- predict(xgb_orig, dtest)
pred_mat <- matrix(pred_prob, ncol = 3, byrow = TRUE)
pred_class <- factor(levels(y_train)[max.col(pred_mat)], levels = levels(y_train))
cm_orig <- confusionMatrix(pred_class, y_test)
cat(sprintf("Original accuracy: %.4f (Kappa=%.4f)\n", cm_orig$overall["Accuracy"], cm_orig$overall["Kappa"]))

# === DIAGNOSIS 2: Check if byrow=FALSE works better ===
cat("\n--- Try byrow=FALSE ---\n")
pred_mat2 <- matrix(pred_prob, ncol = 3, byrow = FALSE)
pred_class2 <- factor(levels(y_train)[max.col(pred_mat2)], levels = levels(y_train))
cm2 <- confusionMatrix(pred_class2, y_test)
cat(sprintf("byrow=FALSE accuracy: %.4f (Kappa=%.4f)\n", cm2$overall["Accuracy"], cm2$overall["Kappa"]))

# === FIX: Use multi:softmax (returns class directly, no reshape needed) ===
cat("\n--- FIX: multi:softmax + more rounds + class weights ---\n")

# Compute class weights for imbalance
class_counts <- table(y_train_num)
class_weights <- sum(class_counts) / (length(class_counts) * class_counts)
cat("Class weights:", paste(round(class_weights, 3), collapse=", "), "\n")

params_fixed <- list(
  objective = "multi:softmax",
  num_class = 3,
  max_depth = 6,
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8
)

set.seed(42)
xgb_fixed <- xgb.train(
  params = params_fixed,
  data = dtrain,
  nrounds = 500,
  verbose = 0
)

# multi:softmax returns class labels directly (as integers 0-based)
pred_raw <- predict(xgb_fixed, dtest)
pred_fixed <- factor(levels(y_train)[pred_raw + 1], levels = levels(y_train))
cm_fixed <- confusionMatrix(pred_fixed, y_test)
cat(sprintf("Fixed accuracy: %.4f (Kappa=%.4f)\n", cm_fixed$overall["Accuracy"], cm_fixed$overall["Kappa"]))

# === FIX 2: Try with scale_pos_weight-like approach using sample weights ===
cat("\n--- FIX v2: With sample weights ---\n")
sample_weights <- class_weights[y_train_num + 1]
dtrain_w <- xgb.DMatrix(data = X_train, label = y_train_num, weight = sample_weights)

set.seed(42)
xgb_w <- xgb.train(params = params_fixed, data = dtrain_w, nrounds = 500, verbose = 0)
pred_w <- predict(xgb_w, dtest)
pred_w_class <- factor(levels(y_train)[pred_w + 1], levels = levels(y_train))
cm_w <- confusionMatrix(pred_w_class, y_test)
cat(sprintf("Weighted accuracy: %.4f (Kappa=%.4f)\n", cm_w$overall["Accuracy"], cm_w$overall["Kappa"]))

# === FIX 3: LASSO and RF also re-run for confirmation ===
cat("\n--- Re-run RF for verification ---\n")
set.seed(42)
rf_model <- randomForest(x = X_train, y = y_train, ntree = 500)
rf_pred <- predict(rf_model, X_test)
rf_cm <- confusionMatrix(rf_pred, y_test)
cat(sprintf("RF Accuracy: %.4f | Kappa: %.4f\n", rf_cm$overall["Accuracy"], rf_cm$overall["Kappa"]))

cat("\n--- Re-run LASSO for verification ---\n")
set.seed(42)
cv_lasso <- cv.glmnet(x = X_train, y = y_train, family = "multinomial", alpha = 1, nfolds = 5)
lasso_pred <- predict(cv_lasso, X_test, s = "lambda.min", type = "class")
lasso_pred <- factor(as.vector(lasso_pred), levels = levels(y_train))
lasso_cm <- confusionMatrix(lasso_pred, y_test)
cat(sprintf("LASSO Accuracy: %.4f | Kappa: %.4f\n", lasso_cm$overall["Accuracy"], lasso_cm$overall["Kappa"]))
cat(sprintf("LASSO genes: %d\n", {
  coefs <- coef(cv_lasso, s = "lambda.min")
  sum(sapply(coefs, function(x) sum(x[-1] != 0)))
}))

# === Save updated metrics ===
cat("\n--- Updated model comparison ---\n")
model_metrics <- data.frame(
  Model = c("Random Forest", "LASSO", "XGBoost"),
  Accuracy = c(rf_cm$overall["Accuracy"], lasso_cm$overall["Accuracy"], cm_fixed$overall["Accuracy"]),
  Kappa = c(rf_cm$overall["Kappa"], lasso_cm$overall["Kappa"], cm_fixed$overall["Kappa"]),
  F1_Macro = c(mean(rf_cm$byClass[, "F1"], na.rm = TRUE),
               mean(lasso_cm$byClass[, "F1"], na.rm = TRUE),
               mean(cm_fixed$byClass[, "F1"], na.rm = TRUE)),
  stringsAsFactors = FALSE
)
print(model_metrics)
write.csv(model_metrics, "D:/Users/Desktop/R_Work/results/tables/classification_metrics.csv", row.names = FALSE)
cat("\nSaved: classification_metrics.csv\n")

# === Check confusion matrices ===
cat("\n--- Confusion Matrices ---\n")
cat("Random Forest:\n"); print(rf_cm$table)
cat("\nLASSO:\n"); print(lasso_cm$table)
cat("\nXGBoost (fixed):\n"); print(cm_fixed$table)

cat("\n=== Diagnosis complete ===\n")
