suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(maftools)
})
npg10 <- pal_npg("nrc")(10)

maf <- read.maf("D:/Users/Desktop/R_Work/data/public_data/brca_mutations_temp.maf", verbose = FALSE)
clin <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_tumor_processed.rds")$clinical

tmb <- maf@data %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(n_mut = n(), .groups = "drop") %>%
  mutate(patient_id = substr(Tumor_Sample_Barcode, 1, 12), TMB = n_mut / 38)

clin$pid12 <- substr(clin$patient_id, 1, 12)
tmb$subtype <- clin$molecular_subtype[match(tmb$patient_id, clin$pid12)]

tmb_sub <- tmb %>%
  filter(!is.na(subtype), subtype != "") %>%
  group_by(subtype) %>%
  summarise(n = n(), median = round(median(TMB), 2),
            mean = round(mean(TMB), 2), sd = round(sd(TMB), 2),
            .groups = "drop")
print(tmb_sub)
write.csv(tmb_sub, "D:/Users/Desktop/R_Work/results/deep/tmb_by_subtype.csv", row.names = FALSE)

tmb_plot <- tmb %>%
  filter(!is.na(subtype), subtype != "") %>%
  mutate(subtype = factor(subtype, levels = c("Luminal A","Luminal B","HER2-enriched","Triple Negative")))
cols <- c("Luminal A"=npg10[1],"Luminal B"=npg10[2],"HER2-enriched"=npg10[3],"Triple Negative"=npg10[4])

p <- ggplot(tmb_plot, aes(subtype, TMB, fill = subtype)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.5) +
  scale_fill_manual(values = cols, guide = "none") +
  labs(title = "Tumor Mutation Burden by Molecular Subtype", x = "", y = "TMB (mutations/Mb)") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("D:/Users/Desktop/R_Work/results/figures_pub/fig18_tmb_by_subtype.png",
       p, width = 8, height = 6, dpi = 300, bg = "white")
cat("Done\n")
