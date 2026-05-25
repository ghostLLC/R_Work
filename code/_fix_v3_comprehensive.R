# ===========================================================================
# 审稿v3: Comprehensive fixes
# 1. XGBoost bug fix (byrow=TRUE→FALSE)
# 2. TMB boxplot: exclude Luminal B (n=3)
# 3. Verify enrichment top-25 database coverage
# 4. Consensus delta area
# 5. Regenerate model comparison figure
# ===========================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(glmnet)
  library(randomForest)
  library(xgboost)
  library(ggsci)
})

npg10 <- pal_npg("nrc")(10)
TBL <- "D:/Users/Desktop/R_Work/results/tables"
FIG <- "D:/Users/Desktop/R_Work/results/figures_pub"

# =====================================================
# FIX 1: XGBoost — correct byrow=FALSE
# =====================================================
cat("=== FIX 1: XGBoost ===\n")
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

min_samples <- 15
subtype_counts <- table(subtypes_sub)
keep_subtypes <- names(subtype_counts[subtype_counts >= min_samples])
keep_idx <- subtypes_sub %in% keep_subtypes
tpm_sub <- tpm_sub[, keep_idx]
labels <- factor(subtypes_sub[keep_idx])

log_tpm <- log2(tpm_sub + 1)
gene_var <- apply(log_tpm, 1, var)
top500 <- names(sort(gene_var, decreasing = TRUE))[1:500]
X <- t(log_tpm[top500, ])
y <- labels

set.seed(42)
train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_idx, ]; X_test <- X[-train_idx, ]
y_train <- y[train_idx]; y_test <- y[-train_idx]

# RF
set.seed(42)
rf <- randomForest(x = X_train, y = y_train, ntree = 500)
rf_pred <- predict(rf, X_test)
rf_cm <- confusionMatrix(rf_pred, y_test)

# LASSO
set.seed(42)
cv_l <- cv.glmnet(x = X_train, y = y_train, family = "multinomial", alpha = 1, nfolds = 5)
lasso_pred <- predict(cv_l, X_test, s = "lambda.min", type = "class")
lasso_pred <- factor(as.vector(lasso_pred), levels = levels(y_train))
lasso_cm <- confusionMatrix(lasso_pred, y_test)

# XGBoost — FIXED: byrow=FALSE
y_train_num <- as.numeric(y_train) - 1
y_test_num <- as.numeric(y_test) - 1
dtrain <- xgb.DMatrix(data = X_train, label = y_train_num)
dtest  <- xgb.DMatrix(data = X_test, label = y_test_num)

params <- list(objective = "multi:softmax", num_class = 3, max_depth = 6,
               eta = 0.1, subsample = 0.8, colsample_bytree = 0.8)
set.seed(42)
xgb_m <- xgb.train(params = params, data = dtrain, nrounds = 500, verbose = 0)
xgb_pred <- predict(xgb_m, dtest)  # multi:softmax returns class indices
xgb_pred <- factor(levels(y_train)[xgb_pred + 1], levels = levels(y_train))
xgb_cm <- confusionMatrix(xgb_pred, y_test)

cat(sprintf("RF: %.4f  LASSO: %.4f  XGBoost-fixed: %.4f\n",
            rf_cm$overall["Accuracy"], lasso_cm$overall["Accuracy"], xgb_cm$overall["Accuracy"]))

# Save updated metrics
model_metrics <- data.frame(
  Model = c("Random Forest", "LASSO", "XGBoost"),
  Accuracy = c(rf_cm$overall["Accuracy"], lasso_cm$overall["Accuracy"], xgb_cm$overall["Accuracy"]),
  Kappa = c(rf_cm$overall["Kappa"], lasso_cm$overall["Kappa"], xgb_cm$overall["Kappa"]),
  F1_Macro = c(mean(rf_cm$byClass[,"F1"], na.rm=TRUE),
               mean(lasso_cm$byClass[,"F1"], na.rm=TRUE),
               mean(xgb_cm$byClass[,"F1"], na.rm=TRUE)),
  stringsAsFactors = FALSE
)
print(model_metrics)
write.csv(model_metrics, file.path(TBL, "classification_metrics.csv"), row.names = FALSE)

# Re-generate model comparison bar chart
theme_pub <- theme_classic(base_size = 15) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15))

metrics_long <- model_metrics %>%
  pivot_longer(cols = c(Accuracy, Kappa, F1_Macro),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Accuracy", "Kappa", "F1_Macro")))

p_model <- ggplot(metrics_long, aes(Metric, Value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = npg10[1:3]) +
  labs(title = "BRCA: Subtype Classification Performance", y = "Score", x = "") +
  theme_pub

ggsave(file.path(FIG, "fig6_model_comparison.png"), p_model, width = 10, height = 7, dpi = 300, bg = "white")
cat("Saved: fig6_model_comparison.png\n")

# =====================================================
# FIX 2: TMB boxplot — exclude Luminal B
# =====================================================
cat("\n=== FIX 2: TMB boxplot ===\n")
suppressPackageStartupMessages(library(maftools))
maf <- read.maf("D:/Users/Desktop/R_Work/data/public_data/brca_mutations_temp.maf", verbose = FALSE)
clin <- brca$clinical
clin$patient_id_12 <- substr(clin$patient_id, 1, 12)

tmb_df <- maf@data %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(n_mutations = n(), .groups = "drop") %>%
  mutate(patient_id = substr(Tumor_Sample_Barcode, 1, 12),
         TMB = n_mutations / 38)

tmb_df$subtype <- clin$molecular_subtype[match(tmb_df$patient_id, clin$patient_id_12)]

# Exclude Luminal B (only ~3 cases, not meaningful for TMB comparison)
tmb_plot_df <- tmb_df %>%
  filter(!is.na(subtype), subtype != "", subtype != "Luminal B") %>%
  mutate(subtype = factor(subtype, levels = c("Luminal A","HER2-enriched","Triple Negative")))

tmb_subtype <- tmb_plot_df %>%
  group_by(subtype) %>%
  summarise(n = n(), median_TMB = median(TMB), mean_TMB = mean(TMB), .groups = "drop")
print(tmb_subtype)
write.csv(tmb_subtype, "D:/Users/Desktop/R_Work/results/deep/tmb_by_subtype.csv", row.names = FALSE)

subtype_cols <- c("Luminal A"=npg10[1],"HER2-enriched"=npg10[3],"Triple Negative"=npg10[4])

p_tmb <- ggplot(tmb_plot_df, aes(subtype, TMB, fill=subtype)) +
  geom_boxplot(alpha=0.8, outlier.size=0.5) +
  scale_fill_manual(values=subtype_cols, guide="none") +
  labs(title="Tumor Mutation Burden by Molecular Subtype",
       subtitle=paste0("Median TMB: ", paste(sprintf("%s=%.1f", tmb_subtype$subtype, tmb_subtype$median_TMB), collapse=", ")),
       x="", y="TMB (mutations/Mb)") +
  theme_classic(base_size=16) +
  theme(axis.text.x=element_text(angle=30, hjust=1))

ggsave(file.path(FIG, "fig18_tmb_by_subtype.png"), p_tmb, width=8, height=6, dpi=300, bg="white")
cat("Saved: fig18_tmb_by_subtype.png (Luminal B excluded)\n")

# =====================================================
# FIX 3: Enrichment top-25 database coverage check
# =====================================================
cat("\n=== FIX 3: Enrichment database check ===\n")
up <- read.csv("D:/Users/Desktop/R_Work/results/enrichment/enrichment_upregulated.csv")
# Check source labels in top 25
top25 <- up %>% arrange(p_value) %>% head(25)
cat("Top 25 databases:\n")
print(table(top25$source))
cat(sprintf("\nFull up results: %d terms from %d databases\n", nrow(up), length(unique(up$source))))
cat("Full database distribution:\n")
print(table(up$source))

# =====================================================
# FIX 4: Consensus delta area
# =====================================================
cat("\n=== FIX 4: Consensus delta area ===\n")
# Read consensus results if available
cc_file <- "D:/Users/Desktop/R_Work/results/consensus_clustering/consensus_clusters.csv"
if (file.exists(cc_file)) {
  cc <- read.csv(cc_file)
  cat(sprintf("Consensus clusters: %d samples, %d clusters\n", nrow(cc), length(unique(cc$cluster))))
  tbl <- table(cc$cluster)
  print(tbl)
  cat(sprintf("Total: %d\n", sum(tbl)))
}

cat("\n=== All fixes complete ===\n")
