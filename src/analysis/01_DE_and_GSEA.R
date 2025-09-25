#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(data.table)
  library(limma)
  library(fgsea)
})

option_list <- list(
  make_option(c("--data"), type="character", default="data/raw", help="raw data dir"),
  make_option(c("--out"), type="character", default="results", help="output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(file.path(opt$out, "tables"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out, "figures"), showWarnings=FALSE, recursive=TRUE)

# ---- PLACEHOLDER PIPELINE ----
# Replace this with actual GEO loading and preprocessing for your chosen accession(s).

message(">>> This is a scaffold. Insert GEO loading and normalization here.")

# Example dummy DE table
de <- tibble(gene=c("IL17A","TNF","KRT6A","CXCL8","PPARG","CNR2"),
             logFC=c(2.1,1.8,1.6,1.3,-0.8,-0.6),
             P.Value=c(1e-8,3e-7,2e-6,1e-5,0.01,0.03),
             adj.P.Val=p.adjust(c(1e-8,3e-7,2e-6,1e-5,0.01,0.03), method="BH"))
write_csv(de, file.path(opt$out, "tables", "de_top.csv"))

# Rank for GSEA
ranked <- setNames(de$logFC, de$gene)

# Example fgsea with hallmark (user to supply gmt)
# hallmark <- gmtPathways("data/raw/h.all.v7.5.1.symbols.gmt")
# fg <- fgsea(pathways = hallmark, stats = ranked, nperm = 1000)
# write_csv(fg, file.path(opt$out, "tables", "gsea_hallmark.csv"))

# Save simple plot placeholder
png(file.path(opt$out, "figures", "volcano_placeholder.png"), width=900, height=700)
with(de, plot(logFC, -log10(P.Value), pch=19, main="Volcano (placeholder)"))
text(de$logFC, -log10(de$P.Value), labels=de$gene, pos=3, cex=0.8)
dev.off()

message(">>> Wrote placeholder DE table and volcano plot. Replace with real analysis.")