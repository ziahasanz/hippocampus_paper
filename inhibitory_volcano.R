library(EnhancedVolcano)
library(tidyverse)
library(patchwork)

# Path to directory with CSV files
path <- "~/you_file_path/in_DEGs/Cluster_#" #load the data from the in_DEGs, for examples Cluster_1

# Desired 10-panel order: 5 female (row 1), 5 male (row 2)
order_female <- c("cE14F","eSTRF","eERF","lSTRF","lERF")
order_male   <- c("cE14M","eSTRM","eERM","lSTRM","lERM")
desired_base <- c(order_female, order_male)

# List all CSV files and parse the SampleBase from "<SampleBase>_Something.csv"
files_tbl <- tibble(path = list.files(path, pattern = "\\.csv$", full.names = TRUE)) %>%
  mutate(
    fname      = basename(path),
    root       = tools::file_path_sans_ext(fname),
    SampleBase = sub("_[^_]+$", "", root)  # take everything before last underscore
  )

# Reorder to your exact desired sequence
ord_tbl <- tibble(SampleBase = desired_base, ord = seq_along(desired_base)) %>%
  left_join(files_tbl, by = "SampleBase") %>%
  arrange(ord)

# (Optional) warn if any expected files are missing
if (any(is.na(ord_tbl$path))) {
  warning("Missing files for: ",
          paste(ord_tbl$SampleBase[is.na(ord_tbl$path)], collapse = ", "))
}

# Volcano factory
generate_volcano_plot <- function(data, title) {
  EnhancedVolcano(
    data,
    lab = NA,
    x = "logFC",
    y = "p_val_adj",
    pCutoff = 0.05,
    FCcutoff = 0.5,
    xlim = c(-4, 4),
    ylim = c(0, 10),
    labSize = 0,
    maxoverlapsConnectors = 0,
    drawConnectors = FALSE,
    legendPosition = "none",
    pointSize = 0.8,
    colAlpha = 4/5,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    cutoffLineType = "blank",
    border = "full",
    borderWidth = 0.5,
    borderColour = "black",
    xlab = expression(Log[2]~FC)
  ) +
    ggtitle(title) +
    theme(
      axis.text.x  = element_text(size = 10),
      axis.text.y  = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      plot.title   = element_text(size = 8, hjust = 0.5)
    )
}

# Build plots in the specified order (skip missing gracefully)
plots <- lapply(seq_len(nrow(ord_tbl)), function(i) {
  if (is.na(ord_tbl$path[i])) {
    ggplot() + theme_void() + ggtitle(paste0(ord_tbl$SampleBase[i], " (missing)")) +
      theme(plot.title = element_text(size = 8, hjust = 0.5),
            panel.border = element_rect(colour = "grey80", fill = NA))
  } else {
    df <- read.csv(ord_tbl$path[i], check.names = FALSE)
    generate_volcano_plot(df, ord_tbl$root[i])
  }
})

# Arrange: 5 columns × 2 rows (row1 = female, row2 = male)
grid_plot <- wrap_plots(plots, ncol = 5, nrow = 2, byrow = TRUE)
print(grid_plot)

#ggsave("in_cluster_#.pdf",  plot = grid_plot, width = 14, height = 6, dpi = 400)
#ggsave("in_cluster_#.tiff", plot = grid_plot, width = 14, height = 6, dpi = 600, device = "tiff")
