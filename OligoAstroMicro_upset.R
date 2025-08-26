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

# ---- Root with glia subfolders ----
main_path <- "your_file_path/OligoAstroMicro_DEGs/"

# Map subfolder -> cell-type suffix used in filenames
suffix_map <- c(
  astro = "Astrocytes",
  micro = "Microglia",
  oligo = "Oligodendrocytes"
)

# List immediate subfolders (astro, micro, oligo)
subfolders <- list.dirs(path = main_path, full.names = TRUE, recursive = FALSE)

for (folder in subfolders) {
  cell_key <- basename(folder)                 # "astro" | "micro" | "oligo"
  if (!cell_key %in% names(suffix_map)) {
    cat("Skipping unknown subfolder:", cell_key, "\n")
    next
  }
  cell_suffix <- suffix_map[[cell_key]]        # "Astrocytes" | "Microglia" | "Oligodendrocytes"
  
  # Read CSVs in this cell-type folder
  deg_files <- list.files(path = folder, pattern = "\\.csv$", full.names = TRUE)
  combined_deg_list <- list()
  
  for (file in deg_files) {
    condition <- tools::file_path_sans_ext(basename(file))  # e.g., "lERM_Astrocytes"
    deg_data  <- read.csv(file, check.names = FALSE)
    filtered  <- deg_data %>% filter(p_val_adj <= 0.05, abs(logFC) >= 0.5)
    combined_deg_list[[condition]] <- filtered$gene_names
  }
  
  # Skip if nothing after filtering
  if (!length(combined_deg_list)) {
    cat("Skipping", cell_suffix, "- No data after filtering.\n")
    next
  }
  
  # Expected sets for this cell type in fixed order (e.g., "lERM_Astrocytes", ..., "cE14F_Astrocytes")
  ordered_sets <- paste0(set_order, "_", cell_suffix)
  
  # Keep only sets that actually exist (preserves order)
  present_sets <- ordered_sets[ordered_sets %in% names(combined_deg_list)]
  if (!length(present_sets)) {
    cat("Skipping", cell_suffix, "- No matching files for expected sets.\n")
    next
  }
  
  # Colors by sex for the present sets
  present_base_groups <- sub(paste0("_", cell_suffix, "$"), "", present_sets)
  present_colors <- unname(order_colors[present_base_groups])
  
  # ---- draw & store in Global env (object named by subfolder: astro/micro/oligo) ----
  obj_name <- cell_key
  assign(
    obj_name,
    upset(
      fromList(combined_deg_list),
      sets = present_sets,
      keep.order = TRUE,
      sets.bar.color = present_colors,
      sets.x.label = "Total DEGs",
      mainbar.y.label = paste("DEGs (Distinct and Overlap) -", cell_suffix)
    )
  )
  
  # ---- show on screen ----
  cat("Displaying UpSet plot for", cell_suffix, "\n")
  print(get(obj_name))
  
  # ---- record the drawn plot, then replay into a PDF (7x6 inches) ----
  rec <- grDevices::recordPlot()  # capture what was just drawn
  pdf_file <- file.path(folder, paste0(cell_suffix, "_UpSet.pdf"))
  grDevices::cairo_pdf(pdf_file, width = 7, height = 6, onefile = FALSE)
  grDevices::replayPlot(rec)
  grDevices::dev.off()
  cat("Saved:", pdf_file, "\n")
  
  # Pause to inspect interactively
  cat("Press [Enter] to continue to the next cell type...\n")
  readline()
}
