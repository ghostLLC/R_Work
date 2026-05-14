# Fig3: Pure ggplot2 KM curve — NO survminer dependency for rendering
library(survival)
library(tidyverse)
library(ggsci)
library(showtext)
showtext_auto()
font_add("songti", regular = "C:/Windows/Fonts/simsun.ttc")

npg10 <- pal_npg("nrc")(10)
OUT <- "D:/Users/Desktop/R_Work/results/figures_pub"

brca <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_tumor_processed.rds")
clinical <- brca$clinical %>%
  filter(!is.na(stage_simple), stage_simple != "", os_time > 0) %>%
  mutate(os_years = os_time / 365.25)

fit <- survfit(Surv(os_years, vital_status) ~ stage_simple, data = clinical)

# Build KM data frame for ggplot
km_df <- data.frame(
  time     = fit$time,
  surv     = fit$surv,
  upper    = fit$upper,
  lower    = fit$lower,
  n_risk   = fit$n.risk,
  n_event  = fit$n.event,
  strata   = rep(names(fit$strata), fit$strata)
) %>%
  mutate(strata = gsub("stage_simple=", "", strata))

# Add t=0 rows
stages <- unique(km_df$strata)
start_df <- data.frame(
  time = 0, surv = 1, upper = 1, lower = 1, n_risk = NA, n_event = 0, strata = stages
)
km_df <- bind_rows(start_df, km_df)

stage_cols <- c("Stage I" = npg10[1], "Stage II" = npg10[2],
                "Stage III" = npg10[3], "Stage IV" = npg10[4])

# Log-rank p-value
lr <- survdiff(Surv(os_years, vital_status) ~ stage_simple, data = clinical)
pval <- round(1 - pchisq(lr$chisq, df = length(lr$n) - 1), 6)

# Number at risk table (build manually)
time_cuts <- c(0, 5, 10, 15, 20)
nrisk_df <- data.frame()
for (s in stages) {
  sub <- clinical %>% filter(stage_simple == s) %>%
    mutate(os_years = os_time / 365.25)
  for (t in time_cuts) {
    n_at_risk <- sum(sub$os_years >= t)
    nrisk_df <- rbind(nrisk_df, data.frame(Time = t, Stage = s, N = n_at_risk))
  }
}
nrisk_df$Time <- factor(paste0(nrisk_df$Time, "yr"), levels = paste0(time_cuts, "yr"))

# ---- PLOT ----
p <- ggplot(km_df, aes(x = time, y = surv, color = strata, fill = strata)) +
  geom_step(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.08, show.legend = FALSE) +
  scale_color_manual(values = stage_cols, name = "Stage") +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_x_continuous(breaks = seq(0, 25, 5), limits = c(0, 25), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.02), expand = c(0, 0),
                     labels = scales::percent) +
  labs(title = "BRCA: Overall Survival by Pathologic Stage",
       subtitle = sprintf("Log-rank p = %.4f  |  Events: 87 / 1,105 patients", pval),
       x = "Time (Years)", y = "Overall Survival") +
  annotate("text", x = 21, y = 0.95, label = sprintf("p = %.4f", pval),
           size = 5, family = "songti", hjust = 0, color = "grey30") +
  theme_classic(base_size = 16, base_family = "songti") +
  theme(
    plot.title    = element_text(size = 19, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey40"),
    axis.text     = element_text(size = 14, color = "black"),
    axis.title    = element_text(size = 15),
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 13),
    legend.position = c(0.75, 0.25),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    axis.line     = element_line(linewidth = 0.6),
    axis.ticks    = element_line(linewidth = 0.6)
  )

# Number at risk sub-table
p_nrisk <- ggplot(nrisk_df, aes(x = Time, y = Stage, label = N)) +
  geom_text(size = 4.5, family = "songti") +
  labs(title = "Number at Risk", x = "", y = "") +
  theme_minimal(base_size = 13, base_family = "songti") +
  theme(
    plot.title  = element_text(size = 12, hjust = 0, face = "bold"),
    panel.grid  = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 12, face = "bold")
  )

# Combine: main plot (85%) + nrisk table (15%) using grid
library(gridExtra)
library(grid)

png(file.path(OUT, "fig3_km_stage.png"), width = 10, height = 9,
    units = "in", res = 300, bg = "white")
grid.arrange(
  p,
  p_nrisk + theme(plot.margin = margin(t = -5, l = 40)),
  ncol = 1,
  heights = c(4, 0.9)
)
dev.off()

cat(sprintf("Saved: fig3_km_stage.png (10x9, %d strata, log-rank p=%.4f)\n",
            length(stages), pval))
