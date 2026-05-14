# BRCA Paper Figure Standards (locked in 2026-05-14)

## Global defaults
- **Format**: PNG, 300 DPI, white background
- **Language**: English only (no Chinese)
- **Color palette**: NPG (Nature Publishing Group) via `ggsci::pal_npg("nrc")(10)`
- **R packages**: ggplot2, ggrepel, ggsci, pheatmap, survival, survminer, gridExtra, maftools

## Font sizes (all in pt)
| Element | Size |
|---------|------|
| `theme_classic(base_size)` | 15-17 |
| `plot.title` | 18-21 bold, centered |
| `plot.subtitle` | 12-14, grey30, centered |
| `axis.title` | 15-17 |
| `axis.text` | 14-15, black |
| `legend.title` | 14-15 |
| `legend.text` | 12-14 |
| `geom_text` (bar labels) | 5-6 bold |
| `geom_text` (forest plot HR labels) | 5-6 bold |
| `geom_text` (nrisk table) | 5-6 bold |
| `geom_text` (miRNA-mRNA tiles) | 5-7 bold |
| `pheatmap fontsize_row` | 10 |
| `pheatmap fontsize` | 14 |
| `annotate` text | 5-6 bold |

## Figure dimensions (inches)
| Figure | Width | Height |
|--------|-------|--------|
| Volcano | 10 | 8 |
| PCA | 10 | 8 |
| KM curve | 11 | 10 |
| Oncoplot | 12 | 6 |
| Heatmap | 10 | 8 |
| Model comparison | 10 | 8 |
| WGCNA SFT | 16 | 6.5 |
| WGCNA modules | 10 | 8 |
| Cox forest | 12 | 9 |
| miRNA-mRNA | 13 | 10 |
| Summary | 22 | 14 |

## Theme settings
```
theme_classic(base_size) +
  theme(
    plot.title       = element_text(size = 18-21, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 12-14, hjust = 0.5, color = "grey30"),
    axis.text        = element_text(size = 14-15, color = "black"),
    axis.title       = element_text(size = 15-17),
    legend.title     = element_text(size = 14-15),
    legend.text      = element_text(size = 12-14),
    panel.grid       = element_blank(),
    axis.line        = element_line(linewidth = 0.7),
    axis.ticks       = element_line(linewidth = 0.7)
  )
```

## Key decisions
- Survminer's ggsurvplot was abandoned — risk table doesn't render Chinese well, legend placement is buggy
- KM curves use pure ggplot2 with manual nrisk table via grid.arrange
- Oncoplot: top 12 genes, no colBar, big legendFontSize (1.5)
- Heatmap: top 25 genes, subsampled 80 tumors, thin dendrograms, no borders
- Volcano: ggrepel labels for top 10 up + top 10 down DEGs
- Model comparison: grouped bar chart (3 metrics × 3 models), not single-bar
