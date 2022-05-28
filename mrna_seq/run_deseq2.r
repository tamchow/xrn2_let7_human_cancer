init <- function() {
  if (!suppressMessages(require("DESeq2", quietly = TRUE))) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("DESeq2", lib = .libPaths()[1])
  }
  suppressMessages(library(tools))
}

init()
library("DESeq2")

do_deseq2 <- function(count_file,
                      count_meta_file,
                      groups = NA,
                      normalizer = NA,
                      # normalize "normalizer" RPM to this RPM
                      normalize_to = 1e4) {
  count_data <- read.csv(count_file, header = TRUE, row.names = 1)
  count_data_meta <- read.csv(
    count_meta_file,
    row.names = 1,
    colClasses = c("character", "factor", "factor", "factor")
  )
  count_data[] <- lapply(count_data, function(x) {
    if (is.numeric(x)) {
      as.integer(round(x))
    } else {
      x
    }
  })
  count_data <- count_data[, rownames(count_data_meta)]
  if (!is.na(normalizer)) {
    cat(sprintf("\nDetected normalizer row label: %s\n", normalizer))
    normalizer_value <- count_data[normalizer, ]
    normalize_factor <- as.integer(round(normalize_to / normalizer_value))
    cat("
    Normalization factor:
    ")
    print(normalize_factor)
    check_transitive <- all(normalize_factor == 1)
    if (check_transitive) {
      cat("
      Checking transitivity of normalization
      as normalization factor is unity
      ")
    }
    new_count_data <- data.frame(count_data)
    for (irow in seq_len(nrow(new_count_data))) {
      row <- new_count_data[irow, ]
      new_row <- row * normalize_factor
      if (check_transitive & !identical(row, new_row)) {
        R.utils::printf(
          "Transitivity check failed: row %d count data unequal",
          irow,
          file = stderr()
        )
      }
      new_count_data[irow, ] <- row * normalize_factor
    }
    if (check_transitive & !identical(count_data, new_count_data)) {
      write.csv(count_data, file = "./count_data.csv")
      write.csv(new_count_data, file = "./new_count_data.csv")
      print(count_data - new_count_data)
      print(
        "Transitivity check failed: table count data unequal",
        file = stderr()
      )
    }
    count_data <- new_count_data
  }
  sample_groups <- unique(count_data_meta$group)
  all_results <- c()
  for (sample_group in sample_groups) {
    sample_group_meta <- count_data_meta[
      count_data_meta$group == sample_group,
    ]
    sample_group_meta$set <- droplevels(sample_group_meta$set)
    count_data_group <- count_data[
      , names(count_data) %in% rownames(sample_group_meta)
    ]
    count_data_dds <- DESeqDataSetFromMatrix(
      countData = count_data_group,
      colData = sample_group_meta[, names(sample_group_meta) != "group"],
      design = ~ set + treatment
    )

    count_data_analysis <- tryCatch(
      DESeq(count_data_dds, quiet = TRUE),
      error = function(error) {
        tryCatch(
          DESeq(count_data_dds, sfType = "iterate", quiet = TRUE),
          error = function(error) {
            DESeq(count_data_dds, sfType = "poscounts", quiet = TRUE)
          }
        )
      }
    )
    count_data_results_names <- resultsNames(count_data_analysis)
    if (is.na(groups) || (length(groups) != 2)) {
      groups <- sort(unique(count_data_meta$treatment), decreasing = TRUE)
      if (length(groups) != 2) {
        stop(paste0(
          "Contrast groups could not be ",
          "unambiguously estimated from the metadata."
        ))
      }
    }
    contrast <- c("treatment", as.character(groups))
    cat(
      sprintf(
        "
        Detected experimental group: %s
        Detected control group: %s
        If this is incorrect, please specify groups manually after other args.
        ", groups[1], groups[2]
      )
    )
    if (!all(
      gettextf("treatment_%s_vs_%s", groups[1], groups[2])
      %in% count_data_results_names
    )) {
      results_names <- paste(count_data_results_names, collapse = ",")
      stop(
        sprintf(
          "contrast groups are not consistent with resultsNames: %s",
          results_names
        )
      )
    }
    results_file <- paste0(
      tools::file_path_sans_ext(count_file),
      "_", sample_group, "_deseq2", ".csv"
    )
    count_data_results <- results(count_data_analysis, contrast = contrast)
    all_results[[results_file]] <- count_data_results
  }
  all_results
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "
    No arguments provided.
    Usage:
    RScript run_deseq2.r --args <count file path> <count metadata file path> [<experimental treatment> <control treatment>]"
  )
}
count_file <- args[2]
count_meta_file <- args[3]
groups <- NA
normalizer <- NA
if (!(is.na(args[4]) || is.na(args[5]))) {
  groups <- c(args[4], args[5])
}
if (length(args) > 5) {
  normalizer <- args[6]
}
all_results <- do_deseq2(
  count_file = count_file,
  count_meta_file = count_meta_file,
  groups = groups,
  normalizer = normalizer
)
for (result_file in names(all_results)) {
  write.csv(all_results[[result_file]], file = result_file)
  cat(sprintf("\nResults saved to %s\n", result_file))
}