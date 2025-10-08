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
#install_if_missing <- function(pkgs, bioc = FALSE) {
#  for (p in pkgs) {
#    if (!requireNamespace(p, quietly = TRUE)) {
#      if (bioc) {
#        if (!requireNamespace("BiocManager", quietly = TRUE)) {
#          install.packages("BiocManager")
#        }
#        BiocManager::install(p, ask = FALSE, update = FALSE)
#      } else {
#        install.packages(p, repos = "https://cloud.r-project.org")
#      }
#    }
#  }
#}

# CRAN packages
#install_if_missing(c("ggplot2", "data.table", "dplyr"))

# Bioconductor packages
#install_if_missing(c("GEOquery", "limma", "fgsea", "msigdbr"), bioc = TRUE)

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

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("hgu133plus2.db")

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
# Step 5b: GSEA Visualization (Bubble + Bar Plot Versions)
###############################################################################
message(">>> Creating GSEA bubble and bar plots...")

library(ggplot2)
library(dplyr)
library(stringr)

# Prepare top 10 enriched pathways
top_pathways <- fgsea_h %>%
  arrange(padj, desc(abs(NES))) %>%
  slice_head(n = 10) %>%
  mutate(
    GeneRatio = leadingEdge %>% sapply(length) / size,
    Count = sapply(leadingEdge, length),
    pathway_clean = gsub("^HALLMARK_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean),
    pathway_clean = str_to_title(pathway_clean),
    pathway_clean = str_wrap(pathway_clean, width = 25),
    pathway_clean = factor(pathway_clean, levels = rev(unique(pathway_clean)))
  )

# ---------------------------------------------------------------------------
# 1️⃣ Bubble (dot) plot — similar to enrichplot::dotplot style
# ---------------------------------------------------------------------------
bubble_plot <- ggplot(top_pathways, aes(x = GeneRatio, y = pathway_clean)) +
  geom_point(aes(size = Count, color = padj)) +
  scale_color_gradient(low = "red", high = "blue", name = "FDR q-value") +
  scale_size(range = c(3, 10), name = "Gene Count") +
  theme_bw(base_size = 14) +
  labs(
    title = "Top 10 Enriched Hallmark Pathways",
    subtitle = "Lesional vs Non-lesional (Psoriasis)\nColor = FDR; Size = Gene Count",
    x = "Gene Ratio (Leading Edge / Pathway Size)",
    y = "Pathway"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.position = "right"
  )

print(bubble_plot)
ggsave(file.path(figures_dir, "GSE13355_GSEA_bubble_plot.png"),
       bubble_plot, width = 8, height = 6, dpi = 600)

# ---------------------------------------------------------------------------
# 2️⃣ Horizontal bar plot — like KEGG/clusterProfiler style
# ---------------------------------------------------------------------------
bar_plot <- ggplot(top_pathways, aes(x = Count, y = pathway_clean, fill = padj)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue", name = "FDR q-value") +
  theme_bw(base_size = 14) +
  labs(
    title = "Top 10 Enriched Hallmark Pathways",
    subtitle = "Lesional vs Non-lesional (Psoriasis)\nColor = FDR q-value",
    x = "Gene Count (Leading Edge)",
    y = "Pathway"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.position = "right"
  )

print(bar_plot)
ggsave(file.path(figures_dir, "GSE13355_GSEA_bar_plot.png"),
       bar_plot, width = 8, height = 6, dpi = 600)


# ---------------------------------------------------------------------------
#  DOT (BUBBLE) PLOT — Compact Visualization of Enrichment Results
# ---------------------------------------------------------------------------
# What this plot shows:
#   - Summarizes the same top 10 enriched Hallmark pathways but adds more information per point.
#   - Each dot = one pathway.
#   - Dot size = number of genes contributing to enrichment (Leading Edge Count).
#   - Dot color = FDR-adjusted significance (q-value).
#   - X-axis = Gene Ratio = (Leading Edge genes) / (Total genes in pathway),
#               showing the *proportion* of a pathway involved in enrichment.
#
# Why we did this:
#   - The dot plot integrates *magnitude* (gene ratio), *significance* (FDR),
#     and *gene count* (dot size) in a single view.
#   - This allows easy comparison of enrichment intensity and significance across pathways.
#   - It’s especially useful for presentations and publications to convey both
#     scale and statistical weight at once.
#
# How to interpret:
#   - Bigger and darker dots = stronger, more significant enrichment.
#   - Pathways clustered near higher Gene Ratios have a large fraction of their genes
#     contributing to the enrichment signal.
#   - Consistent with the bar plot, "E2F Targets", "MYC Targets V1", and
#     "Interferon Alpha Response" dominate, showing that these biological processes
#     are both statistically significant and gene-rich.
#
# Psoriasis-specific biological meaning:
#   - High Gene Ratios in E2F and MYC pathways reinforce widespread activation of
#     proliferation programs in psoriatic lesions.
#   - Interferon signaling enrichment highlights chronic immune stimulation.
#   - These processes collectively capture the dual proliferative–immune nature
#     of psoriatic skin inflammation.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
#  BAR PLOT — Top 10 Enriched Hallmark Pathways
# ---------------------------------------------------------------------------
# What this plot shows:
#   - Displays the top 10 most significantly enriched Hallmark pathways (from fgsea_h).
#   - Each horizontal bar represents one pathway.
#   - Bar length corresponds to the "Leading Edge Gene Count" — the number of genes
#     driving that pathway’s enrichment.
#   - Color (magenta gradient) represents the adjusted FDR q-value; darker shades
#     indicate more statistically significant pathways.
#
# Why we did this:
#   - The bar plot provides a *straightforward summary* of enrichment results.
#   - It allows us to visually rank biological pathways by both significance and gene involvement.
#   - This gives a quick overview of which cellular processes dominate the differential signal
#     between Lesional and Non-lesional psoriasis skin samples.
#
# How to interpret:
#   - Pathways with longer bars and darker color are more significantly enriched.
#   - For example, "E2F Targets" and "MYC Targets V1" are highly enriched and contain
#     large sets of contributing genes.
#   - These pathways represent:
#       • Cell-cycle and proliferation genes (E2F Targets)
#       • Transcriptional and metabolic upregulation (MYC Targets)
#       • Immune signaling (Interferon Alpha/Gamma Response)
#       • Growth and metabolic control (mTORC1 Signaling, Oxidative Phosphorylation)
#
# Psoriasis-specific biological meaning:
#   - Lesional psoriasis tissue shows elevated keratinocyte proliferation and
#     strong activation of immune and interferon pathways.
#   - This plot visually confirms that proliferative and immune responses dominate
#     the transcriptomic signature of psoriatic lesions.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Enrichment plots for Top 3 pathways (Broad Institute–style multi-panel)
# ---------------------------------------------------------------------------
message(">>> Creating Broad Institute–style enrichment plots...")

library(patchwork)
library(stringr)
library(dplyr)
library(fgsea)
library(ggplot2)

# Select top 3 most significant Hallmark pathways
top3 <- fgsea_h %>%
  arrange(padj, desc(abs(NES))) %>%
  slice_head(n = 3)

# ---------------------------------------------------------------------------
# Improved Broad Institute–style enrichment plot function (final)
# ---------------------------------------------------------------------------
make_enrichment_plot <- function(pathway_name, ranks, pathways_list,
                                 up_label = "Lesional (up)",
                                 down_label = "Non-lesional (down)") {
  gene_ranks <- ranks
  pathway_genes <- pathways_list[[pathway_name]]
  
  # Identify pathway gene positions
  hits <- which(names(gene_ranks) %in% pathway_genes)
  hits <- sort(hits)
  
  # Compute running enrichment score
  N <- length(gene_ranks)
  Nh <- length(hits)
  Phit <- cumsum(ifelse(seq_len(N) %in% hits, abs(gene_ranks[seq_len(N)])^1, 0))
  Phit <- Phit / sum(abs(gene_ranks[hits])^1)
  Pmiss <- cumsum(ifelse(seq_len(N) %in% hits, 0, 1 / (N - Nh)))
  runningES <- Phit - Pmiss
  enrichment_df <- data.frame(Position = seq_along(runningES), ES = runningES)
  
  # Build publication-style plot
  p <- ggplot(enrichment_df, aes(x = Position, y = ES)) +
    geom_line(color = "#009E73", linewidth = 1) +  # green curve
    geom_segment(
      data = data.frame(x = hits),
      aes(x = x, xend = x, y = 0, yend = 0.05),
      color = "black", linewidth = 0.25, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
    labs(
      title = paste0(
        "Enrichment Plot: ",
        str_to_title(gsub("_", " ", gsub("^HALLMARK_", "", pathway_name)))
      ),
      x = "Rank in Ordered Dataset",
      y = "Enrichment Score (ES)"
    ) +
    annotate("text", x = N * 0.98, y = max(enrichment_df$ES) * 0.9,
             label = up_label, hjust = 1, vjust = 0, size = 3.8, color = "#D55E00") +
    annotate("text", x = N * 0.98, y = min(enrichment_df$ES) * 0.9,
             label = down_label, hjust = 1, vjust = 1, size = 3.8, color = "#0072B2") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85"),
      axis.line = element_line(color = "black", linewidth = 0.3)
    )
  
  return(p)
}

# ---------------------------------------------------------------------------
# Generate plots and save combined output
# ---------------------------------------------------------------------------
plots <- lapply(top3$pathway, make_enrichment_plot,
                ranks = ranks, pathways_list = hallmark)

# Match y-axis across all for comparability
common_ylim <- range(unlist(lapply(plots, function(p) ggplot_build(p)$data[[1]]$y)))
plots <- lapply(plots, function(p) p + ylim(common_ylim))

# Combine plots (2-column layout)
combined_plot <- wrap_plots(plots, ncol = 2)

# Ensure output directory exists
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Save combined figure
ggsave(file.path(figures_dir, "GSE13355_GSEA_top3_enrichment_classic_f.png"),
       combined_plot, width = 10, height = 6, dpi = 600)

message(">>> Classic enrichment plots saved successfully.")
# ---------------------------------------------------------------------------
#  ENRICHMENT CURVE PLOTS — Detailed View for Top 3 Pathways
# ---------------------------------------------------------------------------
# What these plots show:
#   - Each panel shows the running Enrichment Score (ES) for one top Hallmark pathway.
#   - The green line shows the running sum of enrichment as we move through all genes
#     ranked by differential expression (lesional vs non-lesional).
#   - Black tick marks represent genes from that pathway along the ranked list.
#   - The grey dashed line at ES = 0 separates up- and down-regulated regions.
#   - “Lesional (up)” and “Non-lesional (down)” annotations mark the direction
#     of enrichment for clarity.
#
# Why we did this:
#   - These enrichment plots are the *core diagnostic output* of GSEA.
#   - They show *where* in the ranked list the pathway genes occur — concentrated
#     at the top (upregulated) or bottom (downregulated).
#   - This visualization is essential to validate that the enrichment signal
#     is not random but driven by coordinated gene behavior.
#
# How to interpret:
#   - A curve peaking *above 0* → pathway enriched in Lesional samples (positive NES).
#   - A curve dipping *below 0* → pathway enriched in Non-lesional samples (negative NES).
#   - The more extreme the peak, the stronger the enrichment.
#
# Pathway-specific findings:
#   1. E2F Targets → Strong positive ES: over-activation of cell cycle and replication.
#   2. MYC Targets V1 → Positive ES: transcriptional upregulation via MYC signaling.
#   3. Interferon Alpha Response → Positive ES: immune activation and antiviral response.
#
# Psoriasis-specific biological meaning:
#   - Lesional skin is marked by rapid epidermal turnover (E2F/MYC) and immune activation
#     (interferon pathways), confirming known pathogenic signatures.
#   - Together, these plots visually confirm that psoriasis combines both
#     hyperproliferation and chronic inflammation at the molecular level.
# ---------------------------------------------------------------------------