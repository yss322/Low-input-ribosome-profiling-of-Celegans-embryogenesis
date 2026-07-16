library(data.table)
library(ribor)
library(edgeR)
library(sva)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyr)
library(dplyr)

# --- Helper functions ---

robust_clr_transform <- function(x, pseudo_count = 1) {
  x_adjusted <- x + pseudo_count
  geo_mean <- exp(mean(log(x_adjusted)))
  log(x_adjusted / geo_mean)
}

apply_robust_clr <- function(df, pseudo_count = 1) {
  genes <- df[, 1]
  clr_data <- apply(df[, -1], 2, robust_clr_transform, pseudo_count = pseudo_count)
  clr_data <- cbind(genes, as.data.frame(clr_data))
  return(clr_data)
}

filter_genes_by_cpm_with_transcripts <- function(count_matrix, min_cols, desired_cpm) {
  transcript_col_name <- colnames(count_matrix)[1]
  transcript_names <- count_matrix[, 1]
  count_matrix_without_transcripts <- count_matrix[, -1]
  cpm_matrix <- cpm(count_matrix_without_transcripts)
  genes_meet_criteria <- rowSums(cpm_matrix >= desired_cpm) >= min_cols
  filtered_count_matrix <- count_matrix_without_transcripts[genes_meet_criteria, ]
  filtered_transcript_names <- transcript_names[genes_meet_criteria]
  final_matrix <- cbind(TranscriptName = filtered_transcript_names, filtered_count_matrix)
  colnames(final_matrix)[1] <- transcript_col_name
  return(final_matrix)
}

my_custom_theme <- function() {
  theme_cowplot() +
    theme(
      plot.subtitle = element_text(size = unit(8, "pt")),
      text = element_text(family = "Helvetica"),
      legend.title = element_text(size = unit(8, "pt")),
      legend.text = element_text(size = unit(8, "pt")),
      axis.title = element_text(size = unit(8, "pt")),
      axis.text = element_text(size = unit(8, "pt")),
      axis.line = element_line(size = 1, color = "black"),
      axis.ticks = element_line(size = 1, color = "black"),
      axis.ticks.length = unit(5, "pt"),
      axis.title.x = element_text(margin = margin(t = 5, unit = "pt")),
      axis.title.y = element_text(margin = margin(r = 5, unit = "pt"))
    )
}

# --- Colors ---
ribo_orange <- rgb(228, 88, 10, maxColorValue = 255)
rna_blue <- rgb(55, 135, 192, maxColorValue = 255)

# ===========================================================================
# SECTION 1: Data Pipeline
# ===========================================================================

cat("Loading .ribo files...\n")
input_dir <- normalizePath("../../input_data", mustWork = TRUE)
ITP.ribo <- Ribo(file.path(input_dir, "ITP.ribo"))
RNASEQ.ribo <- Ribo(file.path(input_dir, "RNASEQ.ribo"))

# --- Extract WT ribo CDS counts ---
cat("Extracting WT ribo CDS counts...\n")

one_cell_B7_1 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 26, range.upper = 40,
  experiment = "WT_1-cell_1",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

one_cell_B10_2 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 25, range.upper = 39,
  experiment = "WT_1-cell_2",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

two_cell_B1_1 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 26, range.upper = 40,
  experiment = "WT_2-cell_1",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

two_cell_B2_2 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 25, range.upper = 40,
  experiment = "WT_2-cell_2",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

two_cell_B3_3 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 25, range.upper = 39,
  experiment = "WT_2-cell_3",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

two_cell_B7_4 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 22, range.upper = 40,
  experiment = "WT_2-cell_4",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

four_cell_B2_1 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 23, range.upper = 40,
  experiment = "WT_4-cell_1",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

four_cell_B3_2 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 24, range.upper = 38,
  experiment = "WT_4-cell_2",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

four_cell_B3_3 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 22, range.upper = 38,
  experiment = "WT_4-cell_3",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

four_cell_B7_4 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 23, range.upper = 40,
  experiment = "WT_4-cell_4",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

eight_cell_B3_1 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 22, range.upper = 40,
  experiment = "WT_8-cell_1",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

eight_cell_B7_2 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 24, range.upper = 40,
  experiment = "WT_8-cell_2",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

eight_cell_B8_3 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 22, range.upper = 40,
  experiment = "WT_8-cell_3",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

eight_cell_B9_4 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 23, range.upper = 40,
  experiment = "WT_8-cell_4",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

eight_cell_B10_5 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 23, range.upper = 40,
  experiment = "WT_8-cell_5",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

# --- Extract OMA ribo CDS counts (needed for batch correction) ---
cat("Extracting OMA ribo CDS counts...\n")

oma_1_1cell_B6_1 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 23, range.upper = 40,
  experiment = "OMA-1_1-cell_1",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_1cell_B10_2 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 25, range.upper = 40,
  experiment = "OMA-1_1-cell_2",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_1cell_B10_3 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 24, range.upper = 40,
  experiment = "OMA-1_1-cell_3",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_2cell_B4_1 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 26, range.upper = 40,
  experiment = "OMA-1_2-cell_1",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_2cell_B5_2 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 24, range.upper = 40,
  experiment = "OMA-1_2-cell_2",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_2cell_B6_3 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 29, range.upper = 40,
  experiment = "OMA-1_2-cell_3",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_2cell_B10_4 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 25, range.upper = 40,
  experiment = "OMA-1_2-cell_4",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_2cell_B8_5 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 23, range.upper = 40,
  experiment = "OMA-1_2-cell_5",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_2cell_B9_6 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 24, range.upper = 40,
  experiment = "OMA-1_2-cell_6",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_4cell_B4_1 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 27, range.upper = 40,
  experiment = "OMA-1_4-cell_1",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_4cell_B5_2 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 26, range.upper = 39,
  experiment = "OMA-1_4-cell_2",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_4cell_B6_3 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 25, range.upper = 40,
  experiment = "OMA-1_4-cell_3",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

oma_1_4cell_B9_4 <- as.data.table(get_region_counts(ITP.ribo,
  range.lower = 25, range.upper = 40,
  experiment = "OMA-1_4-cell_4",
  length = TRUE, transcript = FALSE, tidy = F, region = "CDS"))

# --- Combine ribo counts and reshape ---
cat("Combining ribo counts...\n")

list_of_datasets <- list(
  one_cell_B7_1, one_cell_B10_2,
  two_cell_B1_1, two_cell_B2_2, two_cell_B3_3, two_cell_B7_4,
  four_cell_B2_1, four_cell_B3_2, four_cell_B3_3, four_cell_B7_4,
  eight_cell_B3_1, eight_cell_B7_2, eight_cell_B8_3, eight_cell_B9_4, eight_cell_B10_5,
  oma_1_1cell_B6_1, oma_1_1cell_B10_2, oma_1_1cell_B10_3,
  oma_1_2cell_B4_1, oma_1_2cell_B5_2, oma_1_2cell_B6_3,
  oma_1_2cell_B10_4, oma_1_2cell_B8_5, oma_1_2cell_B9_6,
  oma_1_4cell_B4_1, oma_1_4cell_B5_2, oma_1_4cell_B6_3, oma_1_4cell_B9_4
)

combined <- rbindlist(list_of_datasets)
rcw_ribo <- dcast(combined, transcript ~ experiment)

colnames(rcw_ribo) <- c("transcript",
  "one_cell_B7_1.RIBO", "one_cell_B10_2.RIBO",
  "two_cell_B1_1.RIBO", "two_cell_B2_2.RIBO", "two_cell_B3_3.RIBO", "two_cell_B7_4.RIBO",
  "four_cell_B2_1.RIBO", "four_cell_B3_2.RIBO", "four_cell_B3_3.RIBO", "four_cell_B7_4.RIBO",
  "eight_cell_B3_1.RIBO", "eight_cell_B7_2.RIBO", "eight_cell_B8_3.RIBO", "eight_cell_B9_4.RIBO", "eight_cell_B10_5.RIBO",
  "oma_1_1cell_B6_1.RIBO", "oma_1_1cell_B10_2.RIBO", "oma_1_1cell_B10_3.RIBO",
  "oma_1_2cell_B4_1.RIBO", "oma_1_2cell_B5_2.RIBO", "oma_1_2cell_B6_3.RIBO",
  "oma_1_2cell_B10_4.RIBO", "oma_1_2cell_B8_5.RIBO", "oma_1_2cell_B9_6.RIBO",
  "oma_1_4cell_B4_1.RIBO", "oma_1_4cell_B5_2.RIBO", "oma_1_4cell_B6_3.RIBO", "oma_1_4cell_B9_4.RIBO"
)

# --- Extract RNA-seq counts ---
cat("Extracting RNA-seq counts...\n")

rnaseq_neb <- get_rnaseq(ribo.object = RNASEQ.ribo, tidy = F, compact = F, region = "CDS")
rnaseq_neb <- as.data.table(rnaseq_neb)
rnaseq_w_diff <- dcast(rnaseq_neb, transcript ~ experiment)

colnames(rnaseq_w_diff) <- c("transcript",
  "OMA_1cell_B11_1.RNA", "OMA_1cell_B12_2.RNA", "OMA_1cell_B12_3.RNA",
  "OMA_2cell_B11_1.RNA", "OMA_2cell_B12_2.RNA", "OMA_2cell_B12_3.RNA",
  "OMA_4cell_B11_1.RNA", "OMA_4cell_B12_1.RNA", "OMA_4cell_B12_2.RNA",
  "WT_1cell_B11_1.RNA", "WT_1cell_B12_2.RNA", "WT_1cell_B12_3.RNA",
  "WT_2cell_B11_1.RNA", "WT_2cell_B12_2.RNA", "WT_2cell_B12_3.RNA",
  "WT_4cell_B11_1.RNA", "WT_4cell_B12_2.RNA", "WT_4cell_B12_3.RNA",
  "WT_8cell_B11_1.RNA", "WT_8cell_B12_2.RNA", "WT_8cell_B12_3.RNA"
)

# --- Clean up: extract gene names, remove histones ---
cat("Cleaning and filtering...\n")

ribo_counts <- as.data.table(rcw_ribo)
ribo_counts[, transcript := as.character(transcript)]
ribo_counts[, gene_name := sapply(strsplit(transcript, "gene_symbol:"), function(x) tail(x, n = 1))]
ribo_counts <- ribo_counts[!startsWith(gene_name, "his-")]
ribo_counts[, transcript := NULL]
setcolorder(ribo_counts, c("gene_name", setdiff(names(ribo_counts), "gene_name")))

rnaseq_neb_table <- as.data.table(rnaseq_w_diff)
rnaseq_neb_table[, transcript := as.character(transcript)]
rnaseq_neb_table[, gene_name := sapply(strsplit(transcript, "gene_symbol:"), function(x) tail(x, n = 1))]
rnaseq_neb_table <- rnaseq_neb_table[!startsWith(gene_name, "his-")]
rnaseq_neb_table[, transcript := NULL]
setcolorder(rnaseq_neb_table, c("gene_name", setdiff(names(rnaseq_neb_table), "gene_name")))

# --- CPM filter ---
rnaseq_neb_table_filter <- filter_genes_by_cpm_with_transcripts(rnaseq_neb_table, 18, 10)
ribo_counts_table_filter <- filter_genes_by_cpm_with_transcripts(ribo_counts, 10, 3)

# --- Batch correction ---
cat("Running batch correction...\n")

rnaseq_neb_counts <- as.matrix(rnaseq_neb_table_filter[, -1])
covar_1 <- c(rep(1, 3), rep(2, 3), rep(4, 3), rep(1, 3), rep(2, 3), rep(4, 3), rep(8, 3))
covar_2 <- c(rep(1, 9), rep(0, 12))
covar_matrix <- cbind(covar_1, covar_2)
batch_rna <- c(11, 12, 12, 11, 12, 12, 11, 12, 12, 11, 12, 12, 11, 12, 12, 11, 12, 12, 11, 12, 12)
adjusted_RNA <- ComBat_seq(rnaseq_neb_counts, batch = batch_rna, group = NULL, covar_mod = covar_matrix)
adjusted_table_RNA <- as.data.table(cbind(rnaseq_neb_table_filter[, 1], adjusted_RNA))

ribo_counts_table <- as.matrix(ribo_counts_table_filter[, -1])
batch_ribo <- c(7, 10, 2, 2, 3, 7, 2, 3, 3, 7, 3, 7, 8, 9, 10, 6, 10, 10, 4, 5, 6, 10, 8, 9, 4, 5, 6, 9)
covar_1_ribo <- c(1, 1, 2, 2, 2, 2, 4, 4, 4, 4, 8, 8, 8, 8, 8, 1, 1, 1, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4)
covar_2_ribo <- c(rep(0, 15), rep(1, 13))
cov_matrix_ribo <- cbind(covar_1_ribo, covar_2_ribo)
adjusted_ribo <- ComBat_seq(ribo_counts_table, batch = batch_ribo, group = NULL, covar_mod = cov_matrix_ribo)
adjusted_table_ribo <- cbind(ribo_counts_table_filter[, 1], adjusted_ribo)

# --- Merge and CLR transform ---
cat("CLR transforming...\n")

all_counts <- merge.data.table(
  as.data.table(adjusted_table_ribo),
  adjusted_table_RNA,
  by = "gene_name"
)

clr_transformed_data <- apply_robust_clr(all_counts)

# --- Average replicates ---
averaged_clr_data <- data.frame(
  gene_name = clr_transformed_data$gene_name,
  one_cell.ribo = rowMeans(clr_transformed_data[, c("one_cell_B7_1.RIBO", "one_cell_B10_2.RIBO")]),
  two_cell.ribo = rowMeans(clr_transformed_data[, c("two_cell_B1_1.RIBO", "two_cell_B2_2.RIBO", "two_cell_B3_3.RIBO", "two_cell_B7_4.RIBO")]),
  four_cell.ribo = rowMeans(clr_transformed_data[, c("four_cell_B2_1.RIBO", "four_cell_B3_2.RIBO", "four_cell_B3_3.RIBO", "four_cell_B7_4.RIBO")]),
  eight_cell.ribo = rowMeans(clr_transformed_data[, c("eight_cell_B3_1.RIBO", "eight_cell_B7_2.RIBO", "eight_cell_B8_3.RIBO", "eight_cell_B9_4.RIBO", "eight_cell_B10_5.RIBO")]),
  one_cell.rna = rowMeans(clr_transformed_data[, c("WT_1cell_B11_1.RNA", "WT_1cell_B12_2.RNA", "WT_1cell_B12_3.RNA")]),
  two_cell.rna = rowMeans(clr_transformed_data[, c("WT_2cell_B11_1.RNA", "WT_2cell_B12_2.RNA", "WT_2cell_B12_3.RNA")]),
  four_cell.rna = rowMeans(clr_transformed_data[, c("WT_4cell_B11_1.RNA", "WT_4cell_B12_2.RNA", "WT_4cell_B12_3.RNA")]),
  eight_cell.rna = rowMeans(clr_transformed_data[, c("WT_8cell_B11_1.RNA", "WT_8cell_B12_2.RNA", "WT_8cell_B12_3.RNA")])
)

averaged_clr_TE <- averaged_clr_data %>%
  mutate(
    WT_one_cell_TE = one_cell.ribo - one_cell.rna,
    WT_two_cell_TE = two_cell.ribo - two_cell.rna,
    WT_four_cell_TE = four_cell.ribo - four_cell.rna,
    WT_eight_cell_TE = eight_cell.ribo - eight_cell.rna
  ) %>%
  dplyr::select(gene_name, WT_one_cell_TE, WT_two_cell_TE, WT_four_cell_TE, WT_eight_cell_TE)

# ===========================================================================
# SECTION 2: Figure
# ===========================================================================

cat("Building figure...\n")

# --- Panel A: Translational Efficiency ---
target_gene <- "par-3"
reference_genes <- c("lem-2", "gpd-4", "nos-2")

te_long <- averaged_clr_TE %>%
  filter(gene_name %in% c(target_gene, reference_genes)) %>%
  pivot_longer(cols = -gene_name, names_to = "cell_stage", values_to = "TE") %>%
  mutate(cell_stage = factor(cell_stage,
    levels = c("WT_one_cell_TE", "WT_two_cell_TE", "WT_four_cell_TE", "WT_eight_cell_TE"),
    labels = c("1-cell", "2-cell", "4-cell", "8-cell")))

te_target <- te_long %>% filter(gene_name == target_gene)
te_refs <- te_long %>% filter(gene_name %in% reference_genes)

ref_labels <- te_refs %>%
  filter(cell_stage == "1-cell") %>%
  mutate(label = case_when(
    gene_name == "lem-2" ~ "lem-2\n(High TE)",
    gene_name == "gpd-4" ~ "gpd-4\n(House\nkeeping)",
    gene_name == "nos-2" ~ "nos-2\n(Low TE)"
  ))

target_label <- te_target %>% filter(cell_stage == "8-cell")

panel_A <- ggplot() +
  geom_line(data = te_refs, aes(x = cell_stage, y = TE, group = gene_name),
            linetype = "dashed", color = "grey55", linewidth = 1.2) +
  geom_point(data = te_refs, aes(x = cell_stage, y = TE, group = gene_name),
             shape = 1, color = "grey55", size = 3, stroke = 1) +
  geom_line(data = te_target, aes(x = cell_stage, y = TE, group = gene_name),
            color = "red", linewidth = 1.5) +
  geom_point(data = te_target, aes(x = cell_stage, y = TE, group = gene_name),
             color = "red", size = 3) +
  geom_text(data = ref_labels, aes(x = cell_stage, y = TE, label = label),
            hjust = 1, nudge_x = -0.15, size = 2.5, color = "grey30", lineheight = 0.85) +
  geom_text(data = target_label, aes(x = cell_stage, y = TE, label = target_gene),
            hjust = 0, nudge_x = 0.15, size = 3.5, color = "red", fontface = "bold") +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 2.5)) +
  labs(y = "Normalized Translational efficiency", x = "",
       title = paste0("WBGene00003918 (", target_gene, ")")) +
  my_custom_theme() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "none"
  )

# --- Panel B: Ribosome Occupancy and mRNA abundance ---
ribo_rna_par3 <- averaged_clr_data %>%
  filter(gene_name == target_gene) %>%
  dplyr::select(gene_name, one_cell.ribo, two_cell.ribo, four_cell.ribo, eight_cell.ribo,
                one_cell.rna, two_cell.rna, four_cell.rna, eight_cell.rna) %>%
  pivot_longer(cols = -gene_name, names_to = "variable", values_to = "value") %>%
  mutate(
    data_type = ifelse(grepl("\\.ribo$", variable), "Ribosome\noccupancy", "mRNA\nabundance"),
    cell_stage = case_when(
      grepl("one_cell", variable) ~ "1-cell",
      grepl("two_cell", variable) ~ "2-cell",
      grepl("four_cell", variable) ~ "4-cell",
      grepl("eight_cell", variable) ~ "8-cell"
    ),
    cell_stage = factor(cell_stage, levels = c("1-cell", "2-cell", "4-cell", "8-cell"))
  )

panel_B_labels <- ribo_rna_par3 %>%
  filter(cell_stage == "8-cell") %>%
  mutate(nudge_y = ifelse(data_type == "Ribosome\noccupancy", 0.4, -0.4))

panel_B <- ggplot(ribo_rna_par3, aes(x = cell_stage, y = value, color = data_type, group = data_type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_text(data = panel_B_labels, aes(label = data_type, y = value + nudge_y),
            hjust = 0, nudge_x = 0.12, size = 2.8, lineheight = 0.85, show.legend = FALSE) +
  scale_color_manual(values = c("Ribosome\noccupancy" = ribo_orange, "mRNA\nabundance" = rna_blue)) +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 2.5)) +
  labs(y = "Normalized Values", x = "") +
  my_custom_theme() +
  theme(legend.position = "none")

# --- Combine panels ---
final_figure <- panel_A / panel_B +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# --- Save ---
ggsave("figure_par3_WBGene00003918.pdf", final_figure,
       width = 140, height = 180, units = "mm", dpi = 600)
ggsave("figure_par3_WBGene00003918.png", final_figure,
       width = 140, height = 180, units = "mm", dpi = 300)

# --- Export full dataset for Observable dashboard ---
dashboard_data <- averaged_clr_data %>%
  mutate(
    WT_one_cell_TE = one_cell.ribo - one_cell.rna,
    WT_two_cell_TE = two_cell.ribo - two_cell.rna,
    WT_four_cell_TE = four_cell.ribo - four_cell.rna,
    WT_eight_cell_TE = eight_cell.ribo - eight_cell.rna
  )
dashboard_csv_path <- "../../dashboard/src/data/averaged_clr_all.csv"
dir.create(dirname(dashboard_csv_path), recursive = TRUE, showWarnings = FALSE)
write.csv(dashboard_data, dashboard_csv_path, row.names = FALSE)
cat("Exported", nrow(dashboard_data), "genes to", dashboard_csv_path, "\n")

cat("Done! Saved figure_par3_WBGene00003918.pdf and .png\n")
