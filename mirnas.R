################################################################################
# ANALYSIS OF MIRNA TARGETS AND SIGNALING PATHWAYS INVOLVED.
################################################################################

# Load required packages
suppressPackageStartupMessages({
library(multiMiR)
library(VennDiagram)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(BiocParallel)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(grid)
})

##[optional: define a working directory]
# setwd("~/Documents/master_genomica/genomica_transcriptomica/tarea_1")

# miRNAs to analyze
mirnas_file <- file.path("data", "mirnas.txt")
mirnas <- readLines(mirnas_file, warn = FALSE)
cat("\nmiRNAs:")
print(mirnas)

# output directory
results_dir <- "results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}


################################################################################
# SECTION 1: searching for miRNA targets
################################################################################

# Function to obtain targets for a specific miRNA
get_targets <- function(mirna) {
  cat("\nRetrieving targets for", mirna, "...\n")
  
  # multiMiR query with predictive and validated databases
  targets <- suppressWarnings(
    get_multimir(
      org = "hsa",
      mirna = mirna,
      table = "all",
      predicted.cutoff = 80,
      predicted.site = "conserved",
      summary = TRUE
    )
  )
  
  if (length(targets@data) > 0) {
    return(list(
      complete_data = targets@data,
      unique_targets = unique(targets@data$target_symbol)
    ))
  } else {
    return(list(
      complete_data = data.frame(),
      unique_targets = character(0)
    ))
  }
}

# Obtain and save targets for each miRNA
targets_list <- list() # Unique targets
complete_targets_list <- list() # All targets

for(mirna in mirnas) {
  results <- get_targets(mirna)
  
  # Complete list of targets
  complete_targets_list[[mirna]] <- results$complete_data
  
  # Save only the unique names of target genes
  targets_list[[mirna]] <- results$unique_targets
  targets_list[[mirna]] <- results$unique_targets[results$unique_targets != "NA" & results$unique_targets != ""]
  
  cat("miRNA", mirna, "has", length(targets_list[[mirna]]), "unique targets\n")
}

# Save all targets in a single file
all_targets_df <- data.frame()
for(mirna in mirnas) {
  if(nrow(complete_targets_list[[mirna]]) > 0) {
    # Add column to identify miRNA
    temp_df <- complete_targets_list[[mirna]]
    temp_df$source_mirna <- mirna
    all_targets_df <- rbind(all_targets_df, temp_df)
  }
}

if(nrow(all_targets_df) > 0) {
  all_targets_file <- file.path(results_dir, "all_mirnas_complete_targets.csv")
  write.csv(all_targets_df, all_targets_file, row.names = FALSE)
  cat("Complete data for all miRNAs saved at", all_targets_file, "\n\n")
}

# Comment which miRNA has the most and the least targets
target_counts <- sapply(targets_list, length)
str(target_counts)
max_targets_mirna <- names(target_counts)[which.max(target_counts)]
min_targets_mirna <- names(target_counts)[which.min(target_counts)]

cat("\nmiRNA with the most targets:", max_targets_mirna, "with", max(target_counts), "targets\n")
cat("miRNA with the fewest targets:", min_targets_mirna, "with", min(target_counts), "targets\n\n")

# Identify common targets across all miRNAs
common_targets <- Reduce(intersect, targets_list)
cat("Common targets for all miRNAs:", length(common_targets), "\n")

# Save common targets to file
common_targets_file <- file.path(results_dir, "common_targets.txt")
writeLines(common_targets, common_targets_file)
cat("Common targets saved at:", common_targets_file, "\n\n")


################################################################################
# SECTION 2: identify common targets of miRNAs
################################################################################

# Create Venn diagram to visualize common targets
venn_output <- file.path(results_dir, "miRNA_targets_venn.pdf")
pdf(venn_output, width = 10, height = 10) 
venn_obj <- venn.diagram(
  x = targets_list,
  category.names = mirnas,
  filename = NULL,
  output = FALSE, 
  lwd = 2,
  col = "black",
  fill = c("cornflowerblue", "green", "yellow", "purple", "red"),
  alpha = 0.50,
  main.cex = 2,
  margin = 0.1,
  print.params = FALSE
)

grid.draw(venn_obj)
invisible(dev.off())
cat("Venn diagram saved at:", venn_output, "\n\n")


################################################################################
# SECTION 3: functional enrichment analysis to identify overrepresented 
# cellular processes
################################################################################

# Convert gene symbols to Entrez IDs
entrez_ids <- suppressMessages(
  suppressWarnings(
    clusterProfiler::bitr(common_targets, 
                          fromType = "SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = "org.Hs.eg.db")
  )
)

cat("Performing functional enrichment analysis...\n")

# GO Biological Process
go_bp <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# KEGG Pathways
suppressMessages({
  kegg_pathway <- enrichKEGG(
    gene = entrez_ids$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
})

# Save enrichment results
go_file <- file.path(results_dir, "GO_enrichment_results.csv")
write.csv(as.data.frame(go_bp), go_file, row.names = FALSE)
cat("GO results saved at:", go_file, "\n")

# Context-relevant terms for atherosclerosis
context_file <- file.path("data", "context_keywords.txt")
context_keywords <- readLines(context_file, warn = FALSE) %>% 
  grep(pattern = "^\\s*#|^$", invert = TRUE, value = TRUE) %>% 
  trimws() %>% 
  .[nzchar(.)]

if (length(context_keywords) == 0) {
  stop("Error: No valid keywords were found in ", context_file)
}


# Filter GO results relevant to atherosclerosis context
go_bp_athero <- go_bp %>%
  filter(grepl(paste(context_keywords, collapse="|"), Description, ignore.case=TRUE))

# Results visualization
cat("\nSignificant GO BP terms found before filtering:", nrow(go_bp))
cat("\nSignificant GO BP terms found after filtering:", nrow(go_bp_athero), "\n")

# Bar plot of top 15 GO terms
go_output <- file.path(results_dir, "GO_athero_barplot.pdf")
pdf(go_output, width = 10, height = 8)
print(barplot(go_bp_athero, showCategory=15))
invisible(dev.off())
cat("GO plot saved at:", go_output, "\n")

# GO term network
go_network <- file.path(results_dir, "GO_athero_network.pdf")
pdf(go_network, width = 24, height = 20)
print(emapplot(pairwise_termsim(go_bp_athero), showCategory = 200))
invisible(dev.off())
cat("GO term network saved at:", go_network, "\n")

# Filter KEGG results relevant to atherosclerosis context
kegg_pathway_athero <- kegg_pathway %>%
  filter(grepl(paste(context_keywords, collapse="|"), Description, ignore.case=TRUE))

# Results visualization
cat("\nSignificant KEGG terms found before filtering:", nrow(kegg_pathway))
cat("\nSignificant KEGG terms found after filtering:", nrow(kegg_pathway_athero), "\n")

# Bar plot of top 15 KEGG terms
kegg_output <- file.path(results_dir, "KEGG_barplot.pdf")
pdf(kegg_output, width = 10, height = 8)
print(barplot(kegg_pathway_athero, showCategory=15))
invisible(dev.off())
cat("KEGG plot saved at:", kegg_output, "\n")

# KEGG term network
kegg_network <- file.path(results_dir, "KEGG_network.pdf")
pdf(kegg_network, width = 12, height = 10)
print(emapplot(pairwise_termsim(kegg_pathway_athero), showCategory = 200))
invisible(dev.off())
cat("KEGG term network saved at:", kegg_network, "\n")

# Save contextualized results
context_file <- file.path(results_dir, "athero_relevant_terms.txt")
con <- file(context_file, "w")
write_results <- function(data, title) {
  writeLines(title, con)
  writeLines(paste(rep("-", nchar(title)), collapse = ""), con)
  
  if (nrow(data) > 0) {
    for (i in 1:nrow(data)) {
      writeLines(sprintf("%s: %s (p-adj=%.5f, Count=%d)", 
                         data$ID[i], 
                         data$Description[i], 
                         data$p.adjust[i], 
                         data$Count[i]), con)
    }
  } else {
    writeLines("No relevant terms found.", con)
  }
  
  writeLines("\n", con)
}
write_results(as.data.frame(go_bp_athero), "RELEVANT GO TERMS FOR ATHEROSCLEROSIS")
write_results(as.data.frame(kegg_pathway_athero), "RELEVANT KEGG PATHWAYS FOR ATHEROSCLEROSIS")

close(con)
cat("Atherosclerosis-relevant terms saved at:", context_file, "\n\n")


################################################################################
# ANALYSIS SUMMARY
################################################################################

summary_file <- file.path(results_dir, "analysis_summary.txt")
con <- file(summary_file, "w")

writeLines("miRNA ANALYSIS SUMMARY", con)
writeLines("===================================================", con)
writeLines(paste("Analysis date:", format(Sys.Date(), "%B %d, %Y")), con)

writeLines("\n1. Analyzed miRNAs:", con)
writeLines(paste("   -", mirnas), con)

writeLines("\n2. Number of targets per miRNA:", con)
for(mirna in mirnas) {
  writeLines(sprintf("   - %-15s: %d targets", mirna, length(targets_list[[mirna]])), con)
}

writeLines(sprintf("\n3. miRNA with most targets: %s (%d targets)", max_targets_mirna, max(target_counts)), con)
writeLines(sprintf("4. miRNA with fewest targets: %s (%d targets)", min_targets_mirna, min(target_counts)), con)

writeLines("\n5. Common targets across all miRNAs:", con)
if(length(common_targets) > 0) {
  writeLines(paste("   Total:", length(common_targets), "common targets"), con)
  writeLines("   Top 10 common targets (alphabetically):", con)
  writeLines(paste("   -", sort(common_targets)[1:min(10, length(common_targets))]), con)
} else {
  writeLines("   No common targets found for all miRNAs", con)
}

writeLines("\n6. Enrichment analysis results:", con)
writeLines(sprintf("   - Significant GO (Biological Process) terms: %d", nrow(go_bp)), con)
writeLines(sprintf("   - Atherosclerosis-relevant GO terms: %d", nrow(go_bp_athero)), con)
writeLines(sprintf("   - Significant KEGG pathways: %d", nrow(kegg_pathway)), con)
writeLines(sprintf("   - Atherosclerosis-relevant KEGG pathways: %d", nrow(kegg_pathway_athero)), con)

writeLines("\n7. Top enriched biological processes (top 5):", con)
top_go <- head(go_bp_athero, 5)
for(i in 1:nrow(top_go)) {
  writeLines(sprintf("   %d. %s (p-adj=%.3e)", i, top_go$Description[i], top_go$p.adjust[i]), con)
}

writeLines("\n8. Top enriched KEGG pathways (top 5):", con)
top_kegg <- head(kegg_pathway_athero, 5)
for(i in 1:nrow(top_kegg)) {
  writeLines(sprintf("   %d. %s (p-adj=%.3e)", i, top_kegg$Description[i], top_kegg$p.adjust[i]), con)
}

writeLines("\n===================================================", con)
writeLines("End of summary", con)

close(con)
cat("Analysis summary saved at:", summary_file, "\n")