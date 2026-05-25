# Fix Cox forest plot: Age label and HR precision
suppressPackageStartupMessages({library(tidyverse); library(ggsci)})
npg10 <- pal_npg("nrc")(10)
OUT <- "D:/Users/Desktop/R_Work/results/figures_pub"

cox_res <- read.csv("D:/Users/Desktop/R_Work/results/tables/cox_regression.csv")
cat("Actual Cox results:\n")
print(cox_res)

cdf <- cox_res %>%
  filter(!is.infinite(HR), HR < 100, HR > 0.01) %>%
  mutate(
    Variable = case_when(
      Variable == "stage_II"   ~ "Stage II (vs I)",
      Variable == "stage_III"  ~ "Stage III (vs I)",
      Variable == "stage_IV"   ~ "Stage IV (vs I)",
      Variable == "subtype_LumB"  ~ "Luminal B",
      Variable == "subtype_HER2"  ~ "HER2-enriched",
      Variable == "subtype_TN"    ~ "Triple Negative",
      Variable == "age"           ~ "Age (per year)",       # FIXED: per year, not per 10yr
      Variable == "lymph_nodes"   ~ "Lymph Nodes",
      TRUE ~ Variable
    ),
    sig = case_when(pvalue < 0.001 ~ "***", pvalue < 0.01 ~ "**",
                    pvalue < 0.05 ~ "*", TRUE ~ ""),
    label_x = pmax(CI_upper * 1.6, HR * 2.5)
  )

# Fix: use %.2f for HR values to avoid rounding artifacts (Age HR=1.04 → "HR=1.0" with %.1f)
p9 <- ggplot(cdf, aes(HR, reorder(Variable, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.8) +
  geom_point(size = 5, color = npg10[3]) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.35, linewidth = 2.5, color = npg10[3]) +
  geom_text(aes(x = label_x, label = sprintf("HR=%.2f %s", HR, sig)),  # FIXED: %.1f → %.2f
            hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10(limits = c(0.02, 80)) +
  labs(title = "BRCA: Multivariate Cox Regression",
       subtitle = expression(C-index == 0.771~~~~LR~italic(p) == 8.62 %*% 10^-12),
       x = "Hazard Ratio (95% CI)", y = "") +
  theme_classic(base_size = 18) +
  theme(
    plot.title    = element_text(size = 21, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30"),
    axis.text.y   = element_text(size = 17, color = "black", face = "bold"),
    axis.text.x   = element_text(size = 15, color = "black"),
    axis.title    = element_text(size = 17),
    axis.line     = element_line(linewidth = 0.7),
    axis.ticks    = element_line(linewidth = 0.7)
  )

ggsave(file.path(OUT, "fig9_cox_forest.png"), p9, width = 12, height = 9, dpi = 300, bg = "white")
cat("Saved: fig9_cox_forest.png (Age label fixed: per year, HR=1.04)\n")
