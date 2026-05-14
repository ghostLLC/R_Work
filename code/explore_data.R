Sys.setlocale("LC_ALL", "English")

data_dir <- "D:/Users/Desktop/R_Work/data/private_data"

cat("===== Exploring all 6 cancer expression datasets =====\n\n")

for (c in c("brca", "chol", "coad", "kirc", "luad", "stad")) {
  fname <- file.path(data_dir, paste0(c, "_exp.Rdata"))
  cat(sprintf("Loading: %s\n", fname))
  e <- new.env()
  load(fname, envir = e)
  nm <- ls(e)[1]
  x <- get(nm, envir = e)
  cat(sprintf("  %s: %d genes x %d samples\n", toupper(c), nrow(x), ncol(x)))
  cat(sprintf("  Class: %s, Type: %s\n", class(x), typeof(x)))
  cat(sprintf("  Size: %.1f MB\n", as.numeric(object.size(x)) / 1024^2))
  cat(sprintf("  First 2 rownames: %s\n", paste(head(rownames(x), 2), collapse=", ")))
  cat(sprintf("  First 2 colnames: %s\n", paste(head(colnames(x), 2), collapse=", ")))
  rm(e, x); gc()
  cat("\n")
}

# Also load clinical data
cat("\n===== Clinical datasets =====\n\n")
for (c in c("brca", "chol", "coad", "kirc", "luad", "stad")) {
  fname <- file.path(data_dir, paste0(c, "_clinical.Rdata"))
  cat(sprintf("Loading: %s\n", fname))
  e <- new.env()
  load(fname, envir = e)
  nm <- ls(e)[1]
  x <- get(nm, envir = e)
  cat(sprintf("  %s: %d patients x %d variables\n", toupper(c), nrow(x), ncol(x)))
  cat(sprintf("  Class: %s\n", class(x)))
  cat(sprintf("  Size: %.1f KB\n", as.numeric(object.size(x)) / 1024))
  cat(sprintf("  Colnames[1:5]: %s\n", paste(head(colnames(x), 5), collapse=", ")))
  rm(e, x); gc()
  cat("\n")
}
