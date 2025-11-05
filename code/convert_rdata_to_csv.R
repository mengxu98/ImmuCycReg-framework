# Convert RData files to CSV for comparison
# Run this script to convert all RData files to CSV format

base_dir <- getwd()
results_dir <- file.path(base_dir, "results")

# List of RData files to convert
rdata_files <- c(
  "TCGA-LUAD.Rdata",
  "GTEx-LUAD.Rdata",
  "CNV-LUAD.Rdata",
  "Peak-LUAD.Rdata",
  "ATAC-LUAD.Rdata",
  "Geneinfo_df.Rdata"
)

cat("============================================================\n")
cat("Converting RData files to CSV\n")
cat("============================================================\n")

for (file in rdata_files) {
  rdata_path <- file.path(results_dir, file)
  csv_path <- file.path(results_dir, gsub("\\.Rdata$", ".csv", file, ignore.case = TRUE))
  
  if (file.exists(rdata_path)) {
    cat(sprintf("\nConverting: %s\n", file))
    
    # Load the RData file
    load(rdata_path)
    
    # Get the variable name (assuming it's the same as filename without extension)
    var_name <- gsub("\\.Rdata$", "", file, ignore.case = TRUE)
    var_name <- gsub("-", "_", var_name)
    
    # Try to find the variable
    if (exists(var_name)) {
      data <- get(var_name)
      
      # Save as CSV
      if (is.data.frame(data) || is.matrix(data)) {
        write.csv(data, csv_path, row.names = TRUE)
        cat(sprintf("  ✓ Saved to: %s\n", csv_path))
      } else {
        cat(sprintf("  ⚠ %s is not a data frame or matrix\n", var_name))
      }
    } else {
      # Try common variable names
      common_names <- c("tcga_luad", "gtex_luad", "tcga_cnv", "peak_luad", 
                       "genes_adj_peak", "geneinfo_df")
      
      found <- FALSE
      for (name in common_names) {
        if (exists(name)) {
          data <- get(name)
          if (is.data.frame(data) || is.matrix(data)) {
            write.csv(data, csv_path, row.names = TRUE)
            cat(sprintf("  ✓ Saved to: %s (as %s)\n", csv_path, name))
            found <- TRUE
            break
          }
        }
      }
      
      if (!found) {
        cat(sprintf("  ✗ Could not find variable for %s\n", file))
        cat(sprintf("  Available variables: %s\n", paste(ls(), collapse = ", ")))
      }
    }
  } else {
    cat(sprintf("  ⚠ File not found: %s\n", rdata_path))
  }
}

cat("\n============================================================\n")
cat("Conversion complete\n")
cat("============================================================\n")

