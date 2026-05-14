# ===========================================================================
# Generate Final Paper .docx
# 项目：BRCA多组学数据挖掘课程论文
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(officer)
  library(dplyr)
})

cat("\n========== Generate Course Paper .docx ==========\n\n")

OUTPUT  <- "D:/Users/Desktop/R_Work/results/BRCA_Paper_final.docx"
DRAFT   <- "D:/Users/Desktop/R_Work/docs/paper_final_20260514.md"

# ---- 1. Create document and read draft ----
cat("Step 1: Loading merged draft...\n")
doc <- read_docx()
lines <- readLines(DRAFT, encoding = "UTF-8", warn = FALSE)
cat(sprintf("  %d lines loaded\n", length(lines)))

# ---- 2. Add content ----
cat("Step 2: Adding content to document...\n")

for (line in lines) {
  line <- trimws(line)

  # Skip empty markdown formatting lines
  if (nchar(line) == 0) {
    doc <- body_add_par(doc, "", style = "Normal")
    next
  }

  # Remove markdown formatting for display
  clean <- line
  clean <- gsub("\\*\\*([^*]+)\\*\\*", "\\1", clean)  # bold
  clean <- gsub("\\*([^*]+)\\*", "\\1", clean)         # italic
  clean <- gsub("^#+\\s+", "", clean)                   # headers
  clean <- gsub("^[-*]\\s+", "", clean)                 # list bullets
  clean <- gsub("`([^`]+)`", "\\1", clean)              # inline code

  if (nchar(clean) < 2) {
    doc <- body_add_par(doc, "", style = "Normal")
    next
  }

  # Detect section headers
  is_h1 <- grepl("^摘\\s*要$|^关键词|^1\\s+绪\\s*论|^2\\s|^3\\s|^4\\s|^结\\s*论$|^参考文献$", clean)
  is_h2 <- grepl("^[1-5]\\.[0-9]+\\s+[^0-9]", clean) && nchar(clean) < 80
  is_h3 <- grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\s+[^0-9]", clean) && nchar(clean) < 80

  if (is_h1) {
    doc <- body_add_par(doc, clean, style = "heading 1")
  } else if (is_h2) {
    doc <- body_add_par(doc, clean, style = "heading 2")
  } else if (is_h3) {
    doc <- body_add_par(doc, clean, style = "heading 3")
  } else {
    doc <- body_add_par(doc, clean, style = "Normal")
  }
}

# ---- 3. Save ----
cat("\nStep 3: Saving document...\n")
dir.create(dirname(OUTPUT), showWarnings = FALSE, recursive = TRUE)
print(doc, target = OUTPUT)
cat(sprintf("  Saved: %s (%.1f KB)\n", OUTPUT, file.size(OUTPUT) / 1024))
cat("\n========== Done ==========\n")
