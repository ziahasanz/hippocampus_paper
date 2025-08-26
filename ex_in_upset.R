library(UpSetR)
library(dplyr)
library(tools)

# ---- Desired set order (top-level groups) ----
set_order <- c("lERM", "lERF", "lSTRM", "lSTRF",
               "eERM", "eERF", "eSTRM", "eSTRF",
               "cE14M", "cE14F")

# ---- Colors: Female = red, Male = black ----
female_col <- "#ff0000"
male_col   <- "#000000"
order_colors <- ifelse(grepl("F$", set_order), female_col, male_col)
names(order_colors) <- set_order

# ---- Path with Cluster_* folders ----
main_path <- "/your_file_path/ex_DEGs/" #for inhibiotory just change the folder name to in_DEGs

subfolders <- list.dirs(path = main_path, full.names = TRUE, recursive = FALSE)

for (folder in subfolders) {
  cluster_name <- basename(folder)                   # e.g., "Cluster_0"
  cluster_id   <- sub("^Cluster_", "C", cluster_name)
  
  deg_files <- list.files(path = folder, pattern = "\\.csv$", full.names = TRUE)
  combined_deg_list <- list()
  
  for (file in deg_files) {
    condition <- tools::file_path_sans_ext(basename(file))  # e.g., "lERM_C0"
    deg_data  <- read.csv(file, check.names = FALSE)
    filtered  <- deg_data %>% filter(p_val_adj <= 0.05, abs(logFC) >= 0.5)
    combined_deg_list[[condition]] <- filtered$gene_names
  }
  
  if (!length(combined_deg_list)) {
    cat("Skipping", cluster_name, "- No data after filtering.\n")
    next
  }
  
  ordered_sets <- paste0(set_order, "_", cluster_id)
  present_sets <- ordered_sets[ordered_sets %in% names(combined_deg_list)]
  if (!length(present_sets)) {
    cat("Skipping", cluster_name, "- No matching files for expected sets.\n")
    next
  }
  
  present_base_groups <- sub(paste0("_", cluster_id, "$"), "", present_sets)
  present_colors <- unname(order_colors[present_base_groups])
  
  # ---- draw & store in Global env ----
  assign(
    cluster_name,
    upset(
      fromList(combined_deg_list),
      sets = present_sets,
      keep.order = TRUE,
      sets.bar.color = present_colors,
      sets.x.label = "Total DEGs",
      mainbar.y.label = paste("DEGs (Distinct and Overlap) -", cluster_name)
    )
  )
  
  # ---- show on screen ----
  cat("Displaying UpSet plot for", cluster_name, "\n")
  print(get(cluster_name))
  
  # ---- record the drawn plot, then replay into a PDF (7x6 inches) ----
  rec <- grDevices::recordPlot()  # capture what was just drawn
  pdf_file <- file.path(folder, paste0(cluster_name, "_UpSet.pdf"))
  grDevices::cairo_pdf(pdf_file, width = 7, height = 6, onefile = FALSE)
  grDevices::replayPlot(rec)
  grDevices::dev.off()
  cat("Saved:", pdf_file, "\n")
  
  # Pause if you want to inspect interactively
  cat("Press [Enter] to continue to the next cluster...\n")
  readline()
}