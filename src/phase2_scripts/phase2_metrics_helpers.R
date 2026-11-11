make_sample_ids <- function(metadata, grouping_columns) {
  missing_columns <- setdiff(grouping_columns, colnames(metadata))
  if (length(missing_columns) > 0) {
    stop("Sample ID columns are missing: ", paste(missing_columns, collapse = ", "))
  }
  values <- lapply(metadata[grouping_columns], function(column) {
    value <- as.character(column)
    vapply(value, function(component) {
      if (is.na(component)) {
        return("N")
      }
      bytes <- as.integer(charToRaw(enc2utf8(component)))
      paste0("V", paste(sprintf("%02X", bytes), collapse = ""))
    }, character(1))
  })
  do.call(paste, c(values, sep = "-"))
}

spreadsheet_safe <- function(value, field) {
  value <- as.character(value)
  value[is.na(value)] <- ""
  stripped <- trimws(value, whitespace = "[ \\t\\r\\n]")
  numeric_value <- field == "value" & grepl(
    "^-?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?%?$",
    stripped, perl = TRUE
  )
  formula_value <- grepl("^[=+@-]", stripped)
  value[formula_value & !numeric_value] <- paste0("'", value[formula_value & !numeric_value])
  value
}

acquire_metrics_lock <- function(metrics_path, timeout_seconds = 300,
                                 poll_seconds = 0.1) {
  lock_path <- paste0(metrics_path, ".lock")
  deadline <- Sys.time() + timeout_seconds
  repeat {
    if (dir.create(lock_path, showWarnings = FALSE)) {
      return(lock_path)
    }
    if (Sys.time() >= deadline) {
      stop("Timed out waiting for metrics lock: ", lock_path,
           ". Confirm no writer is running, then remove a stale lock directory if needed.")
    }
    Sys.sleep(poll_seconds)
  }
}

append_qc_metrics <- function(metrics_path, qc_summary, modality,
                              project_id, source_file, units = character(),
                              project_root = NULL) {
  required_columns <- c(
    "project_id", "run_id", "modality", "stage", "sample_id",
    "sample_scope", "metric", "value", "unit", "source_file"
  )
  metrics_path <- path.expand(metrics_path)
  metrics_parent <- normalizePath(dirname(metrics_path), winslash = "/",
                                  mustWork = TRUE)
  metrics_path <- file.path(metrics_parent, basename(metrics_path))
  if (nzchar(Sys.readlink(metrics_path))) {
    stop("Sample metrics destination must not be a symlink: ", metrics_path)
  }
  if (is.null(project_root)) project_root <- dirname(metrics_parent)
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  parent_with_separator <- paste0(metrics_parent, "/")
  root_with_separator <- paste0(project_root, "/")
  if (!identical(metrics_parent, project_root) &&
      !startsWith(parent_with_separator, root_with_separator)) {
    stop("Sample metrics destination must be inside project root: ", metrics_path)
  }
  lock_path <- acquire_metrics_lock(metrics_path)
  on.exit(unlink(lock_path, recursive = TRUE), add = TRUE)
  if (!file.exists(metrics_path)) {
    stop("Sample metrics CSV not found: ", metrics_path)
  }

  existing <- read.csv(
    metrics_path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = "character", na.strings = character()
  )
  missing_columns <- setdiff(required_columns, colnames(existing))
  if (length(missing_columns) > 0) {
    stop("Sample metrics CSV has an unsupported schema; missing columns: ",
         paste(missing_columns, collapse = ", "))
  }
  if (!all(c("sample_id", "run_id") %in% colnames(qc_summary))) {
    stop("QC summary must include sample_id and run_id columns")
  }

  timestamp <- gsub("[^0-9T]", "", format(Sys.time(), "%Y%m%dT%H%M%OS6"))
  backup_path <- sub("\\.csv$", paste0(".before_", tolower(modality), "_qc_",
                                         timestamp, ".csv"), metrics_path)
  backup_index <- 1L
  while (file.exists(backup_path)) {
    backup_path <- sub("\\.csv$", paste0(".before_", tolower(modality),
                                           "_qc_", timestamp, "_", backup_index,
                                           ".csv"), metrics_path)
    backup_index <- backup_index + 1L
  }
  if (!file.copy(metrics_path, backup_path)) {
    stop("Could not back up sample metrics CSV to: ", backup_path)
  }

  metric_columns <- setdiff(colnames(qc_summary), c("sample_id", "run_id"))
  new_rows <- do.call(rbind, lapply(seq_len(nrow(qc_summary)), function(row_index) {
    metrics <- metric_columns[vapply(metric_columns, function(metric) {
      !is.na(qc_summary[[metric]][row_index])
    }, logical(1))]
    if (length(metrics) == 0) {
      return(NULL)
    }
    data.frame(
      project_id = project_id,
      run_id = as.character(qc_summary$run_id[row_index]),
      modality = toupper(modality),
      stage = "qc",
      sample_id = as.character(qc_summary$sample_id[row_index]),
      sample_scope = "qc_sample",
      metric = metrics,
      value = vapply(metrics, function(metric) {
        as.character(qc_summary[[metric]][row_index])
      }, character(1)),
      unit = vapply(metrics, function(metric) {
        if (metric %in% names(units)) unname(units[[metric]]) else ""
      }, character(1)),
      source_file = source_file,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(new_rows)) {
    stop("QC summary contains no metrics to append")
  }
  new_rows <- new_rows[, required_columns, drop = FALSE]
  existing <- existing[, required_columns, drop = FALSE]
  is_previous_run <- existing$stage == "qc" &
    toupper(existing$modality) == toupper(modality) &
    existing$source_file == source_file
  existing <- existing[!is_previous_run, , drop = FALSE]
  updated <- rbind(existing, new_rows)
  for (column in required_columns) {
    updated[[column]] <- spreadsheet_safe(updated[[column]], column)
  }
  temp_path <- tempfile(pattern = paste0(".", basename(metrics_path), "."),
                        tmpdir = dirname(metrics_path))
  on.exit(unlink(temp_path), add = TRUE)
  write.csv(updated, temp_path, row.names = FALSE, na = "")
  if (!file.rename(temp_path, metrics_path)) {
    stop("Could not atomically update sample metrics CSV: ", metrics_path)
  }
  message("Backed up original metrics to ", backup_path)
  message("Appended ", nrow(new_rows), " ", toupper(modality), " QC metrics")
  invisible(backup_path)
}

plot_depth_cell_correlation <- function(qc_summary, modality, output_dir) {
  plot_data <- qc_summary[, c("filtered_cells", "mean_depth_after")]
  names(plot_data) <- c("filtered_cells", "mean_depth")
  plot_data <- plot_data[is.finite(plot_data$filtered_cells) &
                           is.finite(plot_data$mean_depth), , drop = FALSE]
  n_samples <- nrow(plot_data)
  estimate <- p_value <- NA_real_
  method <- "Spearman"

  if (n_samples >= 3 && length(unique(plot_data$filtered_cells)) > 1 &&
      length(unique(plot_data$mean_depth)) > 1) {
    test <- cor.test(plot_data$filtered_cells, plot_data$mean_depth,
                     method = "spearman", exact = FALSE)
    estimate <- unname(test$estimate)
    p_value <- test$p.value
  }

  stats <- data.frame(
    modality = toupper(modality),
    metric_x = "filtered_cells",
    metric_y = "mean_depth_after",
    method = method,
    n_samples = n_samples,
    correlation = estimate,
    p_value = p_value
  )
  write.csv(stats,
            file.path(output_dir, paste0("correlation_", tolower(modality),
                                         "_filtered_cells_vs_depth.csv")),
            row.names = FALSE, na = "")

  annotation <- if (is.finite(estimate)) {
    sprintf("Spearman rho = %.2f; p = %.3g; n = %d", estimate, p_value, n_samples)
  } else {
    sprintf("Spearman test needs >=3 non-constant samples; n = %d", n_samples)
  }
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = !!rlang::sym("filtered_cells"),
                 y = !!rlang::sym("mean_depth"))
  ) +
    ggplot2::geom_point() +
    ggplot2::annotate("text", x = Inf, y = Inf, label = annotation,
                      hjust = 1.05, vjust = 1.5, size = 3.5) +
    ggplot2::labs(
      title = paste(toupper(modality), "filtered cell count vs sequencing depth"),
      x = "Filtered cells per sample",
      y = if (toupper(modality) == "RNA") {
        "Mean RNA counts per cell"
      } else {
        "Mean peak fragments per cell"
      }
    ) +
    ggplot2::theme_classic()
  if (n_samples >= 2 && length(unique(plot_data$filtered_cells)) > 1) {
    plot <- plot + ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                                        se = TRUE, na.rm = TRUE)
  }
  ggplot2::ggsave(
    file.path(output_dir, paste0("correlation_", tolower(modality),
                                 "_filtered_cells_vs_depth.pdf")),
    plot, width = 7, height = 5
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0("correlation_", tolower(modality),
                                 "_filtered_cells_vs_depth.png")),
    plot, width = 7, height = 5, dpi = 300
  )
  invisible(stats)
}