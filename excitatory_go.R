library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Rn.eg.db)
library(ggnewscale)
library(RColorBrewer)

# ------------------------------ Parameters ------------------------------------
deg_dir   <- "~/your_file_path/ex_G0_DEGs_c012"
pval_thr  <- 0.05
logfc_thr <- 0.5
topN      <- 10

# Faceting (clusters) and x-axis ordering
cluster_order <- c("C0","C1","C2")  # facets
base_order    <- c("cE14F","cE14M","eSTRF","eSTRM","eERF","eERM","lSTRF","lSTRM","lERF","lERM")
desired_order <- as.vector(outer(base_order, cluster_order, paste, sep = "_")) # optional

# ---------------------------- Helper functions --------------------------------
strip_ver <- function(x) gsub("\\.\\d+$", "", x)  # remove Ensembl version suffix
norm <- function(x) {
  x |>
    stringr::str_replace_all("[\u2010-\u2015\u2212]", "-") |>  # normalize dashes
    stringr::str_squish() |>
    stringr::str_to_lower()
}

# ---------------------- 1) Load DEGs and split UP/DOWN ------------------------
deg_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE)
combined_up   <- list()
combined_down <- list()

for (file in deg_files) {
  condition <- tools::file_path_sans_ext(basename(file))  # e.g., "cE14F_C0"
  df <- read.csv(file, check.names = FALSE)
  
  gene_col <- dplyr::case_when(
    "Gene"       %in% names(df) ~ "Gene",
    "gene_names" %in% names(df) ~ "gene_names",
    TRUE ~ NA_character_
  )
  if (is.na(gene_col)) next
  
  up_genes <- df %>%
    filter(p_val_adj < pval_thr, logFC >  logfc_thr) %>%
    pull(!!sym(gene_col)) %>% unique()
  
  dn_genes <- df %>%
    filter(p_val_adj < pval_thr, logFC < -logfc_thr) %>%
    pull(!!sym(gene_col)) %>% unique()
  
  combined_up[[condition]]   <- up_genes
  combined_down[[condition]] <- dn_genes
}

# Strip Ensembl version suffixes (if present)
combined_up   <- lapply(combined_up,   strip_ver)
combined_down <- lapply(combined_down, strip_ver)

# ------------------- 2) Build compareCluster input & run ----------------------
gene_sets <- c(
  setNames(combined_up,   paste0(names(combined_up),   "_UP")),
  setNames(combined_down, paste0(names(combined_down), "_DOWN"))
)

go_all <- compareCluster(
  geneCluster = gene_sets,
  fun     = "enrichGO",
  OrgDb   = org.Rn.eg.db,
  keyType = "ENSEMBL",
  ont     = "BP",
  readable = TRUE
)

# ---------------- 3) Tidy results: parse sample & cluster ---------------------
# Cluster like "cE14F_C0_UP" or "..._DOWN"
go_df <- as_tibble(go_all@compareClusterResult) %>%
  mutate(
    Regulation = if_else(grepl("_UP$", Cluster), "Up", "Down"),
    Sample     = sub("_(UP|DOWN)$", "", Cluster),           # drop _UP/_DOWN
    Cluster3   = sub("^.*_(C[0-2])$", "\\1", Sample),       # take suffix C0/C1/C2
    SampleBase = sub("_(C[0-2])$", "", Sample),             # prefix before _C#
    negLog10Q  = -log10(qvalue + 1e-300)
  ) %>%
  filter(!is.na(Cluster3), pvalue < 0.05, qvalue < 0.05, Count >= 3) %>%
  mutate(
    Cluster3   = factor(Cluster3,   levels = cluster_order),
    SampleBase = factor(SampleBase, levels = base_order),
    Regulation = factor(Regulation, levels = c("Up","Down"))
  ) %>%
  droplevels()

#for topN GO terms run 4A # for your selected terms run 4B # then continue from step 5

# ---------------- 4A) Option: keep topN per panel (default) -------------------
sel <- go_df %>%
  group_by(Regulation, Cluster3, Sample) %>%
  slice_min(order_by = p.adjust, n = topN, with_ties = FALSE) %>%
  ungroup()

global_levels <- sel %>%
  group_by(Description) %>%
  summarise(median_padj = median(p.adjust, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_padj) %>%
  pull(Description)

plot_df <- go_df %>%
  filter(Description %in% global_levels) %>%
  mutate(Description = factor(Description, levels = rev(global_levels))) %>%
  droplevels()

# ---------------- 4B) Option: selcted GO term from top list for main figure-> GO IDs -----------
keep_terms <- c(
  "oxidative phosphorylation",
  "aerobic respiration",
  "cellular respiration",
  "cytoplasmic translation",
  "ribosome biogenesis",
  "proton motive force−driven ATP synthesis",
  "ribosomal small subunit assembly",
  "mitochondrial respiratory chain complex assembly",
  "energy derivation by oxidation of organic compounds",
  "ribonucleoprotein complex biogenesis",
  "NADH dehydrogenase complex assembly",
  "mitochondrial respiratory chain complex I assembly",
  "ribosomal small subunit biogenesis",
  "ribosome assembly",
  "generation of precursor metabolites and energy",
  "proton transmembrane transport",
  "long-term memory",
  "memory",
  "neural nucleus development",
  "regulation of synaptic vesicle endocytosis",
  "protein folding",
  "vesicle-mediated transport in synapse",
  "autophagy",
  "process utilizing autophagic mechanism",
  "presynaptic endocytosis",
  "midbrain development",
  "synaptic vesicle endocytosis",
  "axon ensheathment",
  "ensheathment of neurons",
  "learning or memory",
  "cognition",
  "regeneration",
  "gliogenesis",
  "ephrin receptor signaling pathway",
  "positive regulation of cell projection organization",
  "positive regulation of neuron projection development",
  "peripheral nervous system axon regeneration",
  "neurotransmitter receptor internalization"
)

all_terms <- as_tibble(go_all@compareClusterResult) %>%
  distinct(ID, Description) %>% mutate(Desc_clean = norm(Description))

keep_terms_clean <- norm(keep_terms)
keep_ids <- all_terms %>%filter(Desc_clean %in% keep_terms_clean) %>% pull(ID)
plot_df <- go_df %>% filter(ID %in% keep_ids)
id2desc <- plot_df %>% distinct(ID, Description)
desc_levels <- id2desc$Description[match(keep_ids, id2desc$ID)]
desc_levels <- rev(unique(na.omit(desc_levels)))
plot_df <- plot_df %>%
  mutate(Description = factor(Description, levels = desc_levels)) %>%
  droplevels()

# ---------------- 5) Pad empty samples so ticks stay visible ------------------
pad_grid <- tidyr::expand_grid(
  Cluster3   = factor(cluster_order, levels = cluster_order),
  SampleBase = factor(base_order,    levels = base_order)
)

pad_df <- pad_grid %>%
  mutate(
    Description = factor(NA, levels = levels(plot_df$Description)),
    Count       = NA_integer_,
    Regulation  = factor("Up", levels = c("Up","Down")),  # won't draw without Description
    negLog10Q   = NA_real_
  )

plot_df_pad <- bind_rows(plot_df, pad_df) %>%
  mutate(SampleBase = factor(SampleBase, levels = base_order))

# ---------------------------- 6) Plot figure ----------------------------------
down_df <- filter(plot_df_pad, Regulation == "Down")
up_df   <- filter(plot_df_pad, Regulation == "Up")

p_go <- ggplot(plot_df_pad) +
  # DOWN (Blues)
  geom_point(
    data  = down_df,
    aes(x = SampleBase, y = Description, size = Count, fill = negLog10Q),
    shape = 21, color = "#08519C", stroke = 0.6, alpha = 0.95, na.rm = TRUE
  ) +
  scale_fill_gradientn(
    name    = "Down: -log10(q)",
    colours = brewer.pal(9, "Blues"),
    na.value = NA
  ) +
  ggnewscale::new_scale_fill() +
  # UP (Reds)
  geom_point(
    data  = up_df,
    aes(x = SampleBase, y = Description, size = Count, fill = negLog10Q),
    shape = 21, color = "#A50F15", stroke = 0.6, alpha = 0.95, na.rm = TRUE
  ) +
  scale_fill_gradientn(
    name    = "Up: -log10(q)",
    colours = brewer.pal(9, "Reds"),
    na.value = NA
  ) +
  # Place Gene Count legend first
  scale_size(name = "Gene Count", range = c(2, 8)) +
  guides(size = guide_legend(order = 1)) +
  facet_grid(. ~ Cluster3, scales = "free_x", space = "free_x") +
  scale_x_discrete(limits = base_order, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title = "GO BP enrichment per sample (Excitatory clusters C0 / C1 / C2)",
    x = "Samples", y = "GO terms"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 7),
    strip.text   = element_text(face = "bold"),
    panel.grid.major.y = element_line(size = 0.2, linetype = 3)
  )

p_go

# ---------------------------- 7) Save outputs ---------------------------------
saveRDS(go_all, file = "excitatory_GO_up_down_compareCluster.rds")

raw_df <- as.data.frame(go_all@compareClusterResult)
write.csv(raw_df, file = "excitatory_GO_BP_results_raw.csv", row.names = FALSE)

write.csv(go_df,   file = "excitatory_GO_BP_results_tidy.csv",  row.names = FALSE)
write.csv(plot_df, file = "excitatory_GO_BP_results_plotdf.csv", row.names = FALSE)

ggsave("excitatory_GO_BP_plot.pdf",  plot = p_go, width = 11, height = 8.5, units = "in")
ggsave("excitatory_GO_BP_plot.svg",  plot = p_go, width = 11, height = 8.5, units = "in", bg = "transparent")
ggsave("excitatory_GO_BP_plot.jpg",  plot = p_go, width = 11, height = 8.5, units = "in", dpi = 300)