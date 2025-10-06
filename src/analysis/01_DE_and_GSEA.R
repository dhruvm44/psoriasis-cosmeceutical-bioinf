#!/usr/bin/env Rscript

###############################################################################
# Project: Psoriasis Cosmeceutical Bioinformatics Validation
# Script: 01_DE_and_GSEA.R
# Purpose: Differential expression (DE) + GSEA enrichment analysis
# Dataset: GSE13355 (psoriasis skin biopsies)
#
# Author: Dhruv Mishra
# Email: dhruvmishra9977@gmail.com
# Program: Masters of Bioinformatics
###############################################################################

# In this script, I want to:
#   1. Download the GSE13355 dataset from GEO
#   2. Run differential expression (lesional vs non-lesional skin) using limma
#   3. Save the DE table and a volcano plot
#   4. Perform GSEA enrichment with Hallmark + Reactome
#   5. Save GSEA results and enrichment plots

###############################################################################
# Helper: install packages if missing
###############################################################################
install_if_missing <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else {
        install.packages(p, repos = "https://cloud.r-project.org")
      }
    }
  }
}

# CRAN packages
install_if_missing(c("ggplot2", "data.table", "dplyr"))

# Bioconductor packages
install_if_missing(c("GEOquery", "limma", "fgsea", "msigdbr"), bioc = TRUE)

###############################################################################
# Load libraries
###############################################################################
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(dplyr)
})

###############################################################################
# Set output directories
###############################################################################
out_dir <- "results"
tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# Step 1: Download GSE13355 from GEO
###############################################################################
message(">>> Downloading GSE13355...")

gse <- getGEO("GSE13355", GSEMatrix = TRUE)
expr_set <- gse[[1]]

###############################################################################
# Step 2: Extract phenotype and expression data
###############################################################################
pheno <- pData(expr_set)
exprs_data <- exprs(expr_set)

# In this dataset, the "title" column encodes sample type:
#   *_NN_* = non-lesional skin
#   *_PP_* = psoriatic lesional skin
# I want to use this to define my groups cleanly.

group <- ifelse(grepl("_PP_", pheno$title, ignore.case = TRUE), "Lesional",
                ifelse(grepl("_NN_", pheno$title, ignore.case = TRUE),
                       "Non_lesional", NA))

# Keep only valid samples (ignore any that don't fit NN/PP)
keep_samples <- !is.na(group)
exprs_data <- exprs_data[, keep_samples]
group <- factor(group[keep_samples], levels = c("Non_lesional", "Lesional"))

# Quick sanity check
message(">>> Sample counts by group:")
print(table(group))

message(">>> Expression matrix dimensions (probes x samples):")
print(dim(exprs_data))


###############################################################################
# Step 3: Differential expression with limma
###############################################################################
message(">>> Running differential expression...")

design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Lesional_vs_Non")

fit <- lmFit(exprs_data, design)
fit <- eBayes(fit)

de_results <- topTable(fit, coef = "Lesional_vs_Non", number = Inf, sort.by = "P")

# Save DE results
fwrite(de_results, file.path(tables_dir, "GSE13355_DE_results.csv"))

###############################################################################
# Step 4: Final Polished Volcano Plot
###############################################################################
message(">>> Creating volcano plot...")

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("hgu133plus2.db")

library(ggrepel)
library(hgu133plus2.db)

# Map probes → gene symbols
de_results$gene <- mapIds(hgu133plus2.db,
                          keys = rownames(de_results),
                          column = "SYMBOL",
                          keytype = "PROBEID",
                          multiVals = "first")

# Categorize genes by thresholds
de_results$threshold <- "Not Sig"
de_results$threshold[de_results$logFC > 1 & de_results$adj.P.Val < 0.05] <- "Up"
de_results$threshold[de_results$logFC < -1 & de_results$adj.P.Val < 0.05] <- "Down"

# Top 5 up and down genes by fold change
top_up <- de_results %>%
  filter(threshold == "Up") %>%
  top_n(5, wt = logFC)

top_down <- de_results %>%
  filter(threshold == "Down") %>%
  top_n(-5, wt = logFC)

# Psoriasis hallmark genes (always highlight for biological meaning)
highlight_genes <- c("IL17A", "TNF", "CXCL8", "KRT16", "DEFB4A", "S100A7A", "SERPINB4")

# Combine label set
label_genes <- unique(c(top_up$gene, top_down$gene, highlight_genes))

# Make label dataset (remove duplicates, only genes with symbols)
label_data <- de_results %>%
  filter(!is.na(gene)) %>%
  filter(gene %in% label_genes) %>%
  distinct(gene, .keep_all = TRUE)

# Volcano plot
volcano <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = threshold), alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Down" = "#1f78b4", "Up" = "#e31a1c", "Not Sig" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot: Psoriasis Lesional vs Non-lesional Skin",
    subtitle = "Significant genes (FDR < 0.05, |log2FC| > 1) highlighted",
    x = "log2 Fold Change (Lesional vs Non-lesional)",
    y = "-log10 Adjusted p-value",
    color = "Regulation"
  ) +
  theme(
    legend.position.inside = c(0.9, 0.9),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11, color = "grey30")
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.3,
    nudge_y = 3
  )

# Show in RStudio
print(volcano)

# Save high-resolution figure
ggsave(file.path(figures_dir, "GSE13355_volcano_FINAL.png"),
       volcano, width = 8, height = 6, dpi = 600)

# Interpretation of Volcano Plot:
# This volcano plot (GSE13355) compares lesional vs non-lesional psoriasis skin.
# - Upregulated in lesional skin (red): IL17A, TNF, CXCL8, KRT16, DEFB4A, S100A7A, SERPINB4
#   These are hallmark psoriasis drivers: IL-17/TNF pathway activation, keratinocyte hyperproliferation,
#   neutrophil chemotaxis, and induction of antimicrobial peptides.
# - Downregulated in lesional skin (blue): IL37, CCL27, WIF1, BTC, THRSP
#   These reflect loss of anti-inflammatory cytokines (IL-37), reduced skin immune surveillance (CCL27),
#   and impaired lipid/epidermal homeostasis.
#
# Literature context:
# - IL-23/IL-17 axis (Lowes et al., 2014, Annu Rev Immunol; Nestle et al., 2009, NEJM).
# - TNF as therapeutic target (Ettehadi et al., 1994, Br J Dermatol).
# - CXCL8-driven neutrophil infiltration (Homey et al., 2000, J Immunol).
# - Keratinocyte hyperproliferation and AMP induction (Haider et al., 2006, J Invest Derm; Harder & Schröder, 2005).
# - Loss of IL-37 and CCL27 in psoriasis lesions (Nold et al., 2010, Nat Immunol; van der Fits et al., 2004).

###############################################################################
# Step 4b: Supplementary Table - Top 20 DE genes
###############################################################################
# make sure we have a gene column (character)
de_results$gene <- as.character(de_results$gene)

# order by adjusted p-value and take top 20
ord <- order(de_results$adj.P.Val)
top20 <- de_results[ord, c("gene","logFC","adj.P.Val")]
top20 <- head(top20, 20)

# save
data.table::fwrite(top20, file.path(tables_dir, "GSE13355_top20_DE_genes.csv"))

message(">>> Volcano plot and top 20 DE genes table created.")

###############################################################################
# Step 5: GSEA Enrichment Analysis (gene-symbol ranks, msigdbr API)
###############################################################################
message(">>> Running GSEA enrichment...")

suppressPackageStartupMessages({
  library(hgu133plus2.db)
  library(msigdbr)
  library(fgsea)
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

# 5.1 Map probes -> symbols
symbols <- mapIds(
  hgu133plus2.db,
  keys     = rownames(de_results),
  column   = "SYMBOL",
  keytype  = "PROBEID",
  multiVals = "first"
)

# 5.2 Collapse probes to a single stat per gene (choose probe with largest |t|)
t_tbl <- data.frame(gene = symbols, t = de_results$t, stringsAsFactors = FALSE) %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  summarise(t = t[which.max(abs(t))], .groups = "drop")

ranks <- t_tbl$t
names(ranks) <- t_tbl$gene

# 5.3 Get MSigDB pathways (category / subcategory API)
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark <- split(x = hallmark_df$gene_symbol, f = hallmark_df$gs_name)

reactome_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome <- split(x = reactome_df$gene_symbol, f = reactome_df$gs_name)

# 5.4 Run fgsea
set.seed(42)
fgsea_h <- fgsea(pathways = hallmark, stats = ranks, nperm = 1000)
fgsea_r <- fgsea(pathways = reactome, stats = ranks, nperm = 1000)

# 5.5 Save results
fwrite(fgsea_h, file.path(tables_dir, "GSE13355_GSEA_hallmark.csv"))
fwrite(fgsea_r, file.path(tables_dir, "GSE13355_GSEA_reactome.csv"))
###############################################################################
# Step 5b (alternative): Classic multi-panel GSEA enrichment plots
###############################################################################
message(">>> Creating Broad Institute–style enrichment plots...")

library(patchwork)   # for arranging multiple plots
library(ggplot2)
library(stringr)

# Pick the top 3 pathways by FDR and |NES|
top3 <- fgsea_h %>%
  arrange(padj, desc(abs(NES))) %>%
  slice_head(n = 3)

# Function to make a styled enrichment plot
make_enrichment_plot <- function(pathway_name, ranks, pathways_list,
                                 up_label = "Lesional (up)",
                                 down_label = "Non-lesional (down)") {
  
  df <- plotEnrichment(pathways_list[[pathway_name]], ranks) %>%
    ggplot_build() %>% .$data[[1]]
  
  # Running enrichment score curve
  p <- ggplot(df, aes(x, y)) +
    geom_line(color = "forestgreen", linewidth = 1) +
    geom_hline(yintercept = 0, color = "grey40") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Enrichment Plot: ",
                     str_to_title(gsub("_", " ",
                                       gsub("^HALLMARK_", "", pathway_name)))),
      x = "Rank in Ordered Dataset",
      y = "Enrichment Score (ES)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11)
    )
  
  # Add labels for up/down phenotypes
  p <- p +
    annotate("text", x = Inf, y = max(df$y) * 0.9, label = up_label,
             hjust = 1.1, vjust = 1.5, size = 3.5, color = "red") +
    annotate("text", x = Inf, y = min(df$y) * 0.9, label = down_label,
             hjust = 1.1, vjust = -0.5, size = 3.5, color = "blue")
  
  return(p)
}

# Generate the three plots
plots <- lapply(top3$pathway, make_enrichment_plot,
                ranks = ranks, pathways_list = hallmark)

# Combine into a single multi-panel layout
combined_plot <- wrap_plots(plots, ncol = 2)

# Save as PDF and PNG
ggsave(file.path(figures_dir, "GSE13355_GSEA_top3_enrichment_classic.pdf"),
       combined_plot, width = 10, height = 6, dpi = 600)
ggsave(file.path(figures_dir, "GSE13355_GSEA_top3_enrichment_classic.png"),
       combined_plot, width = 10, height = 6, dpi = 600)

###############################################################################


# ---------------------------------------------------------------------------
# Enrichment plots for Top 3 pathways
# ---------------------------------------------------------------------------
message(">>> Creating enrichment plots for top 3 pathways...")

# Get top 3 pathways by padj + |NES|
top3 <- fgsea_h %>%
  arrange(padj, desc(abs(NES))) %>%
  slice_head(n = 3)

# Save combined into PDF
pdf(file.path(figures_dir, "GSE13355_GSEA_top3_enrichment.pdf"), width = 7, height = 5)
for (p in top3$pathway) {
  p_title <- gsub("^HALLMARK_", "", p) %>%
    gsub("_", " ", .) %>%
    str_to_title()
  print(
    plotEnrichment(hallmark[[p]], ranks) +
      labs(title = paste0("Enrichment Plot: ", p_title)) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11)
      )
  )
}
dev.off()
print(top3)

# Save individually as PNGs too
for (p in top3$pathway) {
  p_title <- gsub("^HALLMARK_", "", p) %>%
    gsub("_", " ", .) %>%
    str_to_title()
  
  p_plot <- plotEnrichment(hallmark[[p]], ranks) +
    labs(title = paste0("Enrichment Plot: ", p_title)) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11)
    )
  
  ggsave(file.path(figures_dir, paste0("GSE13355_GSEA_", p_title, ".png")),
         p_plot, width = 7, height = 5, dpi = 600)
}
###############################################################################
