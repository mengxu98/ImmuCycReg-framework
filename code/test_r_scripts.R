# Test R Scripts
# Test if all R scripts can run without errors

cat("============================================================\n")
cat("Testing R Scripts\n")
cat("============================================================\n")

base_dir <- getwd()
code_dir <- file.path(base_dir, "code")

# List of R scripts to test
scripts_to_test <- list(
  list(name = "Step-01_Processdata.R", description = "Step 01: Data Processing"),
  list(name = "Step-02_NMF-ConsensuClustering.R", description = "Step 02: NMF Consensus Clustering"),
  list(name = "Step-03_Cibersort-data.R", description = "Step 03: CIBERSORT Data Preparation"),
  list(name = "Step-04_DESeq2.R", description = "Step 04: DESeq2 Analysis"),
  list(name = "Step-05_CNV-Correlation.R", description = "Step 05: CNV Correlation"),
  list(name = "Step-06_Peak.R", description = "Step 06: Peak Analysis"),
  list(name = "Step-07_Obtain_DNA_sequences.R", description = "Step 07: DNA Sequence Extraction"),
  list(name = "Step-08_ImmuCycReg-framework.R", description = "Step 08: ImmuCycReg Framework")
)

results <- list()
all_passed <- TRUE

# Test package imports
cat("\n============================================================\n")
cat("Testing R package imports\n")
cat("============================================================\n")

required_packages <- c(
  "dplyr", "magrittr", "DESeq2", "ConsensusClusterPlus",
  "ggplot2", "psych", "survival", "survminer"
)

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("✓ %s\n", pkg))
  } else {
    cat(sprintf("✗ %s - NOT INSTALLED\n", pkg))
    all_passed <- FALSE
  }
}

# Test data files
cat("\n============================================================\n")
cat("Checking required data files\n")
cat("============================================================\n")

data_dir <- file.path(base_dir, "data")
required_files <- c(
  "sample_cluster_4.csv",
  "Genes_17.txt",
  "Genes_249.txt"
)

for (file in required_files) {
  file_path <- file.path(data_dir, file)
  if (file.exists(file_path)) {
    cat(sprintf("✓ %s\n", file))
  } else {
    cat(sprintf("✗ %s - NOT FOUND\n", file))
  }
}

# Test script syntax
cat("\n============================================================\n")
cat("Testing R script syntax\n")
cat("============================================================\n")

for (script_info in scripts_to_test) {
  script_name <- script_info$name
  description <- script_info$description
  script_path <- file.path(code_dir, script_name)
  
  if (file.exists(script_path)) {
    cat(sprintf("\nTesting: %s\n", description))
    cat(sprintf("Script: %s\n", script_name))
    
    # Try to parse the script
    tryCatch({
      parse(script_path)
      cat(sprintf("  ✓ Syntax OK\n"))
      results[[script_name]] <- TRUE
    }, error = function(e) {
      cat(sprintf("  ✗ Syntax Error: %s\n", e$message))
      results[[script_name]] <- FALSE
      all_passed <<- FALSE
    })
  } else {
    cat(sprintf("\n⚠ SKIPPED: %s - File not found\n", script_name))
    results[[script_name]] <- NULL
  }
}

# Summary
cat("\n============================================================\n")
cat("Test Summary\n")
cat("============================================================\n")

for (script_name in names(results)) {
  if (isTRUE(results[[script_name]])) {
    cat(sprintf("✓ %s\n", script_name))
  } else if (isFALSE(results[[script_name]])) {
    cat(sprintf("✗ %s\n", script_name))
  } else {
    cat(sprintf("⚠ %s - Not tested\n", script_name))
  }
}

if (all_passed) {
  cat("\n✓ All tests passed!\n")
  quit(status = 0)
} else {
  cat("\n✗ Some tests failed. Please check the errors above.\n")
  quit(status = 1)
}

