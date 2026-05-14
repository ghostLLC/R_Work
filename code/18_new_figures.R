# ===========================================================================
# 18: Figures from R_Analysis results — pub-quality versions
# Enrichment dotplots + Association rules + Consensus clustering
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE, bitmapType = "cairo")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
})

npg10 <- pal_npg("nrc")(10)
OUT   <- "D:/Users/Desktop/R_Work/results/figures_pub"
ENR   <- "D:/Users/Desktop/R_Work/results/enrichment"
ASC   <- "D:/Users/Desktop/R_Work/results/association"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# Pub theme (same as fig standards)
# =========================================================
theme_pub <- theme_classic(base_size = 15) +
  theme(
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    axis.title    = element_text(size = 15),
    axis.text     = element_text(size = 14, color = "black"),
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 12),
    axis.line     = element_line(linewidth = 0.7),
    axis.ticks    = element_line(linewidth = 0.7),
    panel.grid    = element_blank()
  )

# =========================================================
# Fig12: Enrichment dotplot — Top 25 upregulated pathways
# =========================================================
cat("Fig12: Enrichment dotplot\n")

up <- read.csv(file.path(ENR, "enrichment_upregulated.csv"))
down <- read.csv(file.path(ENR, "enrichment_downregulated.csv"))
all_degs <- read.csv(file.path(ENR, "enrichment_all_differentially_expressed.csv"))

# Clean source labels for display
clean_source <- function(df) {
  df %>% mutate(source_label = case_when(
    source == "GO:MF" ~ "GO: Molecular Function",
    source == "GO:BP" ~ "GO: Biological Process",
    source == "GO:CC" ~ "GO: Cellular Component",
    source == "KEGG"  ~ "KEGG Pathway",
    source == "REAC"  ~ "Reactome",
    source == "WP"    ~ "WikiPathway",
    TRUE ~ source
  ))
}

up <- clean_source(up)
down <- clean_source(down)
all_degs <- clean_source(all_degs)

# Top 25 upregulated
plot_up <- up %>%
  arrange(p_value) %>%
  head(25) %>%
  mutate(term_name = factor(term_name, levels = rev(term_name)),
         log10p = -log10(p_value))

p12 <- ggplot(plot_up, aes(x = log10p, y = term_name,
                            size = intersection_size, color = source_label)) +
  geom_point(alpha = 0.85) +
  scale_color_brewer(palette = "Set1", name = "Source") +
  scale_size_continuous(range = c(3, 10), name = "Gene Count") +
  labs(title = "GO/KEGG Enrichment: Upregulated Genes",
       subtitle = sprintf("Top 25 of %d significant terms (FDR<0.05)", nrow(up)),
       x = expression(-log[10](italic(p))), y = "") +
  theme_pub +
  theme(axis.text.y = element_text(size = 11),
        legend.position = "right")

ggsave(file.path(OUT, "fig12_enrichment_up.png"), p12,
       width = 14, height = 8, dpi = 300, bg = "white")
cat("  Saved: fig12_enrichment_up.png\n")

# =========================================================
# Fig13: Enrichment dotplot — Top 25 downregulated
# =========================================================
cat("Fig13: Enrichment downregulated\n")

plot_down <- down %>%
  arrange(p_value) %>%
  head(25) %>%
  mutate(term_name = factor(term_name, levels = rev(term_name)),
         log10p = -log10(p_value))

p13 <- ggplot(plot_down, aes(x = log10p, y = term_name,
                              size = intersection_size, color = source_label)) +
  geom_point(alpha = 0.85) +
  scale_color_brewer(palette = "Set1", name = "Source") +
  scale_size_continuous(range = c(3, 10), name = "Gene Count") +
  labs(title = "GO/KEGG Enrichment: Downregulated Genes",
       subtitle = sprintf("Top 25 of %d significant terms (FDR<0.05)", nrow(down)),
       x = expression(-log[10](italic(p))), y = "") +
  theme_pub +
  theme(axis.text.y = element_text(size = 11),
        legend.position = "right")

ggsave(file.path(OUT, "fig13_enrichment_down.png"), p13,
       width = 14, height = 8, dpi = 300, bg = "white")
cat("  Saved: fig13_enrichment_down.png\n")

# =========================================================
# Fig14: Enrichment — Top 5 sources (up vs down comparison)
# =========================================================
cat("Fig14: Enrichment source comparison\n")

source_comp <- bind_rows(
  up %>% mutate(Direction = "Upregulated"),
  down %>% mutate(Direction = "Downregulated")
) %>%
  count(Direction, source_label) %>%
  group_by(Direction) %>%
  mutate(pct = n / sum(n) * 100)

p14 <- ggplot(source_comp, aes(x = reorder(source_label, n), y = n, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.88) +
  scale_fill_manual(values = c("Upregulated" = npg10[3], "Downregulated" = npg10[5])) +
  labs(title = "Enrichment Results by Source Database",
       subtitle = sprintf("Up: %d terms | Down: %d terms (FDR<0.05)", nrow(up), nrow(down)),
       x = "", y = "Number of Significant Terms") +
  coord_flip() +
  theme_pub

ggsave(file.path(OUT, "fig14_enrichment_comparison.png"), p14,
       width = 10, height = 6, dpi = 300, bg = "white")
cat("  Saved: fig14_enrichment_comparison.png\n")

# =========================================================
# Fig15: Association rules scatter
# =========================================================
cat("Fig15: Association rules\n")

rules_csv <- file.path(ASC, "association_rules.csv")
if (file.exists(rules_csv)) {
  rules <- read.csv(rules_csv)
  if (nrow(rules) > 0) {
    # Take top 20 rules by lift
    top_rules <- rules %>%
      arrange(desc(lift)) %>%
      head(20) %>%
      mutate(LHS_short = substr(LHS, 1, 40))

    p15 <- ggplot(top_rules, aes(x = support, y = confidence,
                                  size = lift, color = lift)) +
      geom_point(alpha = 0.8) +
      scale_color_gradient(low = npg10[5], high = npg10[3], name = "Lift") +
      scale_size_continuous(range = c(3, 12), name = "Lift") +
      labs(title = "Association Rules: Gene Expression → Clinical Features",
           subtitle = sprintf("Top 20 of %d rules (lift-ranked)", nrow(rules)),
           x = "Support", y = "Confidence") +
      theme_pub

    ggsave(file.path(OUT, "fig15_association_rules.png"), p15,
           width = 10, height = 7, dpi = 300, bg = "white")
    cat("  Saved: fig15_association_rules.png\n")
  }
}

# =========================================================
# Fig16: Top enrichment terms summary (bar chart of most significant)
# =========================================================
cat("Fig16: Enrichment summary bars\n")

# Top 10 most significant terms from ALL DEGs
top_all <- all_degs %>%
  arrange(p_value) %>%
  head(10) %>%
  mutate(term_short = ifelse(nchar(term_name) > 60,
                              paste0(substr(term_name, 1, 57), "..."),
                              term_name),
         term_short = factor(term_short, levels = rev(term_short)),
         log10p = -log10(p_value))

p16 <- ggplot(top_all, aes(x = log10p, y = term_short, fill = source_label)) +
  geom_bar(stat = "identity", alpha = 0.88, width = 0.7) +
  scale_fill_brewer(palette = "Set1", name = "Source") +
  geom_text(aes(label = sprintf("p=%.1e", p_value)),
            hjust = -0.1, size = 4.5) +
  labs(title = "Top 10 Most Significant Enriched Terms",
       subtitle = "All DEGs (6,768 genes), FDR < 0.05",
       x = expression(-log[10](italic(p))), y = "") +
  theme_pub +
  theme(axis.text.y = element_text(size = 12))

ggsave(file.path(OUT, "fig16_enrichment_top10.png"), p16,
       width = 13, height = 6.5, dpi = 300, bg = "white")
cat("  Saved: fig16_enrichment_top10.png\n")

cat("\n========== 5 new figures generated ==========\n")
