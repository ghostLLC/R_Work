# ===========================================================================
# US-011: Survival Analysis — KM Curves, Cox Models, Prognostic Signature
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(tidyverse)
  library(glmnet)
})

cat("\n========== US-011: Survival Analysis ==========\n\n")

INPUT_DIR  <- "D:/Users/Desktop/R_Work/data/processed"
FIG_DIR    <- "D:/Users/Desktop/R_Work/results/figures"
TBL_DIR    <- "D:/Users/Desktop/R_Work/results/tables"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load and prepare survival data ----
cat("Step 1: Loading and preparing survival data...\n")

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
clinical <- brca$clinical
tpm_mat <- brca$tpm

# Clean survival data
clinical <- clinical %>%
  mutate(
    os_time_years = os_time / 365.25,
    os_status = vital_status
  )

cat(sprintf("  Patients: %d | Events (death): %d\n",
            nrow(clinical), sum(clinical$os_status == 1, na.rm = TRUE)))
cat(sprintf("  OS time range: %.1f - %.1f years\n",
            min(clinical$os_time_years, na.rm = TRUE),
            max(clinical$os_time_years, na.rm = TRUE)))

# ---- 2. KM curves by molecular subtype ----
cat("\nStep 2: Kaplan-Meier curves by molecular subtype...\n")

clinical_subtype <- clinical %>%
  filter(!is.na(molecular_subtype), molecular_subtype != "",
         !is.na(os_time_years), os_time_years > 0)

fit_subtype <- survfit(
  Surv(os_time_years, os_status) ~ molecular_subtype,
  data = clinical_subtype
)

subtype_colors <- c("Luminal A" = "#1B9E77", "Luminal B" = "#D95F02",
                    "HER2-enriched" = "#7570B3", "Triple Negative" = "#E7298A")

pdf(file.path(FIG_DIR, "km_curves_subtype.pdf"), width = 9, height = 7)
p <- ggsurvplot(
  fit_subtype,
  data = clinical_subtype,
  pval = TRUE,
  palette = subtype_colors,
  legend.title = "Molecular Subtype",
  legend.labs = names(subtype_colors),
  xlab = "Time (Years)",
  ylab = "Overall Survival Probability",
  risk.table = TRUE,
  risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 14)
)
print(p)
dev.off()
cat("  KM curves by subtype saved.\n")

# Log-rank test
lr_subtype <- survdiff(Surv(os_time_years, os_status) ~ molecular_subtype,
                       data = clinical_subtype)
cat(sprintf("  Log-rank p-value: %.4e\n", lr_subtype$pvalue))

# ---- 3. KM curves by stage ----
cat("\nStep 3: KM curves by pathologic stage...\n")

clinical_stage <- clinical %>%
  filter(!is.na(stage_simple), stage_simple != "",
         !is.na(os_time_years), os_time_years > 0)

fit_stage <- survfit(
  Surv(os_time_years, os_status) ~ stage_simple,
  data = clinical_stage
)

stage_colors <- c("Stage I" = "#1B9E77", "Stage II" = "#D95F02",
                  "Stage III" = "#7570B3", "Stage IV" = "#E7298A")

pdf(file.path(FIG_DIR, "km_curves_stage.pdf"), width = 9, height = 7)
p <- ggsurvplot(
  fit_stage,
  data = clinical_stage,
  pval = TRUE,
  palette = stage_colors,
  legend.title = "Pathologic Stage",
  xlab = "Time (Years)",
  ylab = "Overall Survival Probability",
  risk.table = TRUE,
  risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 14)
)
print(p)
dev.off()
cat("  KM curves by stage saved.\n")

# ---- 4. Multivariate Cox regression ----
cat("\nStep 4: Multivariate Cox proportional hazards model...\n")

cox_data <- clinical %>%
  transmute(
    os_time = os_time_years,
    os_status = os_status,
    age = as.numeric(age_at_diagnosis),
    lymph_nodes = as.numeric(positive_lymph_nodes),
    stage_II   = as.numeric(stage_simple == "Stage II"),
    stage_III  = as.numeric(stage_simple == "Stage III"),
    stage_IV   = as.numeric(stage_simple == "Stage IV"),
    subtype_LumB = as.numeric(molecular_subtype == "Luminal B"),
    subtype_HER2 = as.numeric(molecular_subtype == "HER2-enriched"),
    subtype_TN   = as.numeric(molecular_subtype == "Triple Negative")
  ) %>%
  filter(os_time > 0, !is.na(os_status))

# Replace NAs with 0
cox_data[is.na(cox_data)] <- 0

cox_model <- coxph(
  Surv(os_time, os_status) ~ age + lymph_nodes + stage_II + stage_III + stage_IV +
    subtype_LumB + subtype_HER2 + subtype_TN,
  data = cox_data
)

cox_summary <- summary(cox_model)
cat(sprintf("  Concordance: %.3f | Likelihood ratio p=%.2e\n",
            cox_summary$concordance[1], cox_summary$logtest["pvalue"]))

# Filter out variables with extreme/infinite coefficients
cox_coefs <- cox_summary$coefficients
finite_coefs <- cox_coefs[, "coef"] > -20 & cox_coefs[, "coef"] < 20
cox_coefs_finite <- cox_coefs[finite_coefs, , drop = FALSE]

if (nrow(cox_coefs_finite) > 0) {
  # Build simplified model with only finite variables
  finite_vars <- rownames(cox_coefs_finite)
  formula_str <- paste("Surv(os_time, os_status) ~",
                       paste(finite_vars, collapse = " + "))
  cox_simple <- coxph(as.formula(formula_str), data = cox_data)

  # Forest plot with simplified model
  pdf(file.path(FIG_DIR, "cox_forest_plot.pdf"), width = 10, height = 6)
  tryCatch({
    ggforest(cox_simple, data = cox_data,
             main = "BRCA: Multivariate Cox Regression Forest Plot")
  }, error = function(e) {
    # Fallback: manual forest plot
    coef_table <- summary(cox_simple)$coefficients
    conf_int <- summary(cox_simple)$conf.int

    plot_data <- data.frame(
      Variable = rownames(coef_table),
      HR = conf_int[, "exp(coef)"],
      Lower = conf_int[, "lower .95"],
      Upper = conf_int[, "upper .95"],
      pvalue = coef_table[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )

    ggplot(plot_data, aes(x = HR, y = rev(Variable))) +
      geom_point(size = 3, color = "#E41A1C") +
      geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
      labs(title = "BRCA: Multivariate Cox Regression",
           x = "Hazard Ratio (95% CI)", y = "") +
      theme_bw(base_size = 14) +
      scale_x_log10()
  })
  dev.off()
  cat("  Cox forest plot saved.\n")
}

# Save Cox results
cox_table <- data.frame(
  Variable = rownames(cox_summary$coefficients),
  HR = cox_summary$coefficients[, "exp(coef)"],
  CI_lower = cox_summary$conf.int[, "lower .95"],
  CI_upper = cox_summary$conf.int[, "upper .95"],
  pvalue = cox_summary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)
write.csv(cox_table, file.path(TBL_DIR, "cox_regression.csv"), row.names = FALSE)

# ---- 5. Prognostic gene signature (LASSO-Cox) ----
cat("\nStep 5: Building prognostic gene signature with LASSO-Cox...\n")

# Match expression to survival data
expr_patients <- substr(colnames(tpm_mat), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_patients, clinical$pid_short)

# Filter to patients with complete survival data
surv_ok <- !is.na(clinical$os_time[match_idx]) & clinical$os_time[match_idx] > 0 &
           !is.na(clinical$os_status[match_idx])

cat(sprintf("  Patients with complete survival data: %d\n", sum(surv_ok)))

if (sum(surv_ok) > 500) {
  # Use top variable genes + DEGs
  deg_table <- read.csv(file.path(TBL_DIR, "brca_degs_deseq2.csv"))
  top_degs <- deg_table %>%
    filter(!is.na(padj), padj < 0.01, abs(log2FC) > 2) %>%
    pull(gene)

  degs_in_expr <- intersect(top_degs, rownames(tpm_mat))
  cat(sprintf("  Candidate prognostic genes: %d\n", length(degs_in_expr)))

  if (length(degs_in_expr) > 50) {
    # Subset expression
    X <- t(log2(tpm_mat[degs_in_expr, surv_ok] + 1))
    y <- Surv(
      clinical$os_time[match_idx][surv_ok] / 365.25,
      clinical$os_status[match_idx][surv_ok]
    )

    set.seed(42)
    cv_fit <- cv.glmnet(
      x = X, y = y,
      family = "cox",
      alpha = 1,
      nfolds = 5
    )

    # Extract selected genes
    coefs <- coef(cv_fit, s = "lambda.min")
    selected_idx <- which(as.matrix(coefs) != 0)
    selected_genes <- rownames(coefs)[selected_idx]
    selected_coefs <- coefs[selected_idx]

    cat(sprintf("  LASSO-Cox selected %d prognostic genes\n", length(selected_genes)))

    if (length(selected_genes) >= 3) {
      # Calculate risk score
      risk_score <- X[, selected_genes, drop = FALSE] %*% selected_coefs
      risk_score <- as.vector(risk_score)

      # Stratify by median risk
      risk_group <- ifelse(risk_score > median(risk_score), "High Risk", "Low Risk")
      risk_df <- data.frame(
        time = clinical$os_time[match_idx][surv_ok] / 365.25,
        status = clinical$os_status[match_idx][surv_ok],
        risk_group = risk_group,
        risk_score = risk_score,
        stringsAsFactors = FALSE
      )

      # KM by risk group
      fit_risk <- survfit(Surv(time, status) ~ risk_group, data = risk_df)

      pdf(file.path(FIG_DIR, "km_curves_signature.pdf"), width = 9, height = 7)
      p <- ggsurvplot(
        fit_risk,
        data = risk_df,
        pval = TRUE,
        palette = c("High Risk" = "#E41A1C", "Low Risk" = "#377EB8"),
        legend.title = "Prognostic Signature",
        xlab = "Time (Years)",
        ylab = "Overall Survival Probability",
        risk.table = TRUE,
        risk.table.height = 0.25,
        ggtheme = theme_bw(base_size = 14)
      )
      print(p)
      dev.off()
      cat("  KM curves by gene signature saved.\n")

      # Save prognostic genes
      prog_genes <- data.frame(
        Gene = selected_genes,
        Coefficient = selected_coefs,
        stringsAsFactors = FALSE
      ) %>% arrange(desc(abs(Coefficient)))
      write.csv(prog_genes, file.path(TBL_DIR, "prognostic_genes.csv"), row.names = FALSE)

      # ---- 6. Time-dependent ROC ----
      cat("\nStep 6: Time-dependent ROC curves...\n")

      if (requireNamespace("timeROC", quietly = TRUE)) {
        library(timeROC)

        roc_3yr <- timeROC(
          T = risk_df$time,
          delta = risk_df$status,
          marker = risk_df$risk_score,
          cause = 1,
          times = c(1, 3, 5),
          iid = TRUE
        )

        pdf(file.path(FIG_DIR, "time_roc.pdf"), width = 8, height = 7)
        plot(roc_3yr, time = 1, col = "#1B9E77", lwd = 2, title = "BRCA: Time-dependent ROC")
        plot(roc_3yr, time = 3, col = "#D95F02", lwd = 2, add = TRUE)
        plot(roc_3yr, time = 5, col = "#7570B3", lwd = 2, add = TRUE)
        legend("bottomright",
               legend = c(
                 sprintf("1-Year AUC = %.3f", roc_3yr$AUC[1]),
                 sprintf("3-Year AUC = %.3f", roc_3yr$AUC[2]),
                 sprintf("5-Year AUC = %.3f", roc_3yr$AUC[3])
               ),
               col = c("#1B9E77", "#D95F02", "#7570B3"), lwd = 2, bty = "n")
        dev.off()

        # Save AUC table
        roc_auc_df <- data.frame(
          Time = c("1-Year", "3-Year", "5-Year"),
          AUC = roc_3yr$AUC,
          stringsAsFactors = FALSE
        )
        write.csv(roc_auc_df, file.path(TBL_DIR, "time_roc_auc.csv"), row.names = FALSE)
        cat(sprintf("  Time-ROC: 1yr AUC=%.3f, 3yr AUC=%.3f, 5yr AUC=%.3f\n",
                    roc_3yr$AUC[1], roc_3yr$AUC[2], roc_3yr$AUC[3]))
      }
    }
  }
} else {
  cat("  Insufficient samples for prognostic signature (<500 with complete data)\n")
}

cat("\n========== US-011 Complete: Survival Analysis ==========\n")
