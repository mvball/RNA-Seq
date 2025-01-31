########## DESeq Analysis ########

# Function to check and install packages
install_if_missing <- function(package, source = "CRAN") {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (source == "CRAN") {
      install.packages(package)
    } else if (source == "Bioconductor") {
      BiocManager::install(package)
    }
  }
}

# Check if BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of CRAN and Bioconductor packages needed 
cran_packages <- c("ggplot2", "tidyr")
bioc_packages <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db")

# CRAN packages installed?
for (pkg in cran_packages) {
  install_if_missing(pkg, source = "CRAN")
}

# Bioconductor packages installed? 
for (pkg in bioc_packages) {
  install_if_missing(pkg, source = "Bioconductor")
}

# ggrepel package installed?
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

# Install and load the necessary package for heatmap
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# Load necesary libraries: 
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(tidyr)
library(pheatmap)

#### Load data and perform DEseq - 5 Steps ####

# Step 1: Read the data
Counts <- read.delim("~/Desktop/cleaned_featureCounts copy.txt", sep="\t", header = TRUE)

# Step 2: Set Geneid as row names and remove the Geneid column from the data frame
rownames(Counts) <- Counts$Geneid  
Counts <- Counts[, -1]

# Step 3: Set counts as numeric data type. Recreate sample names as columns. 
Counts_numeric <- Counts[, sapply(Counts, is.numeric)]  # Select only numeric columns (count data)
sample <- c("HER21", "HER22", "HER23", "NonTNBC1", "NonTNBC2", "NonTNBC3", "Normal1", "Normal2", "Normal3", "TNBC1", "TNBC2", "TNBC3")
colnames(Counts) <- sample

# Step 4: Assign groups to samples for coldata and create dataframe coldata for DEseq
group <- factor(c("HER2", "HER2", "HER2", "NONTNBC", "NONTNBC", "NONTNBC", "Normal", "Normal", "Normal", "TNBC", "TNBC", "TNBC"))
coldata <- data.frame(row.names = colnames(Counts), group)


# Step 5: Run DEseq analysis
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~group)
dds <- DESeq(dds)

#### Produce a PCA Plot ####
#Assign DESeq2 results to vsdata and do PCA plot
vsdata <- vst(dds, blind=TRUE)

plotPCA(vsdata, intgroup = "group")

##### Get statistics of DEGs, with filtering #####
results_total <- results(dds)
results_df <-as.data.frame(results_total)
results_df <- results_df[results_df$pvalue < 0.05, ]

# Assign gene symbols using org.Hs.eg.db package. Samples used "ENSEMBL" as keytype.
library("org.Hs.eg.db")
results_df$gene <- mapIds(org.Hs.eg.db, keys = rownames(results_df), keytype = "ENSEMBL", column = "SYMBOL")

# Calculate number of DE genes with padj < 0.05
num_DE_genes_total <- nrow(results_df)
num_DE_genes_total
# Remove rows where log2FoldChange is NA
results_df_no_na <- results_df[!is.na(results_df$log2FoldChange), ]
results_df_no_na_total <- nrow(results_df_no_na)
#Perform the up-regulated and down-regulated calculations
up_regulated <- sum(results_df_no_na$log2FoldChange > 1)
down_regulated <- sum(results_df_no_na$log2FoldChange < 1)
# Print results
up_regulated
down_regulated

###### Produce a heatmap of the top 30 genes #######

# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(org.Hs.eg.db)

# Select the top 30 genes with the smallest adjusted p-values
top_genes <- head(results_df_no_na[order(results_df_no_na$padj), ], 50)

# Debug: Check if ESR1 is in the top genes
if (!"ESR1" %in% top_genes$gene) {
  message("ESR1 not found in the top genes list. Forcing inclusion.")
  esr1_row <- significant_genes[significant_genes$gene == "ESR1", ]
  top_genes <- rbind(esr1_row, top_genes[1:(nrow(top_genes) - 1), ])
}

# Extract ENSEMBL IDs for the top genes
ensembl_ids <- rownames(top_genes)

# Ensure the ENSEMBL IDs exist in the DESeq2 object
valid_ids <- ensembl_ids[ensembl_ids %in% rownames(dds)]
if (length(valid_ids) < length(ensembl_ids)) {
  message("Warning: Some ENSEMBL IDs are missing in the DESeq2 object.")
  print(ensembl_ids[!ensembl_ids %in% rownames(dds)])
}

# Extract normalized counts for the valid ENSEMBL IDs
normalized_counts <- counts(dds, normalized = TRUE)[valid_ids, ]

# Replace rownames with gene symbols for better readability
rownames(normalized_counts) <- top_genes$gene[rownames(top_genes) %in% valid_ids]

# Scale data for visualization (z-score normalization)
scaled_counts <- t(scale(t(normalized_counts)))

# Create the heatmap
pheatmap(
  scaled_counts,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = coldata,
  main = "Heatmap of Top 50 DEGs",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)



####### Contrast "NONTNBC" vs "TNBC" ############
res <- results(dds, contrast = c("group", "NONTNBC", "TNBC"))

# convert the results to a data frame
res_df <- as.data.frame(res)

########## Filter results and Create Visualizations ##########

# Filter by p-value < 0.05
res_df <- res_df[res_df$pvalue < 0.05, ]

# Assign gene symbols using org.Hs.eg.db package. Samples used "ENSEMBL" as keytype.
library("org.Hs.eg.db")
res_df$gene <- mapIds(org.Hs.eg.db, keys = rownames(res_df), keytype = "ENSEMBL", column = "SYMBOL")

# Calculate number of DE genes with padj < 0.05
num_DE_genes <- nrow(res_df)
num_DE_genes


# Remove rows where log2FoldChange is NA
res_df_no_na <- res_df[!is.na(res_df$log2FoldChange), ]

# Now perform the up-regulated and down-regulated calculations
up_regulated <- sum(res_df_no_na$log2FoldChange > 1)
down_regulated <- sum(res_df_no_na$log2FoldChange < 1)
# Print results
up_regulated
down_regulated

num_DE_genes <- nrow(res_df)
num_DE_genes


# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Map ENSEMBL IDs to gene symbols and assign them as rownames
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = rownames(normalized_counts), 
                       keytype = "ENSEMBL", 
                       column = "SYMBOL")

# Assign gene symbols as rownames of the matrix
rownames(normalized_counts) <- gene_symbols

# Select a few genes of interest
genes_of_interest <- c("FOXA1", "AGR2", "AGR3", "TFF1", "ESR1", "ELF5", "FSIP1", "ACTG2", "ACADSB")  
normalized_counts_selected <- normalized_counts[genes_of_interest, ]

# View the normalized counts for the selected genes
normalized_counts_selected


#FOXA1 = It is downregulated or lost in tnbc, which is often ER-negative and HER2-negative.
#AGR2 = breast cancer oncogene, highly expressed in tnbc, leads to poor prognosis
#AGR3 = overexpression = contribute to tnbc aggressiveness. FOXA1 is coreg of AGR3
#TFF1 = TFF1 is a well-established marker of estrogen receptor (ER)-positive breast cancers. Low or absent TFF1 expression is typical in most TNBC cases. 
#ESR1 =  TNBC is defined by the absence of ESR1 (Estrogen Receptor 1)
#ELF5 expression is increased in basal-like or TNBC cell lines, contributing to their stem-like features.
#FSIP1 = association with key processes in cell growth, differentiation, and resistance, FSIP1 may be a potential biomarker for TNBC or a therapeutic target.
#ACTG2 = studies suggest that smooth muscle actin (like ACTG2) is involved in fibroblast activation and tumor stroma remodeling, which are key components of TNBC metastasis
#ACADSB = could be involved in the regulation of metabolic flexibility, helping TNBC cells adapt to nutrient deprivation and therapeutic stress


##### Volcano Plot of NonTNBC and TNBC Sig Diff Exp Genes #####

# Add library for volcano plot
library(ggplot2)
library(ggrepel)

# Prepare data for volcano plot
volcano_data <- as.data.frame(res)
volcano_data$gene <- mapIds(org.Hs.eg.db, 
                            keys = rownames(volcano_data), 
                            keytype = "ENSEMBL", 
                            column = "SYMBOL")

# Remove rows with NA gene symbols
volcano_data <- volcano_data[!is.na(volcano_data$gene), ]

# Add significance column based on padj and log2FoldChange thresholds
volcano_data$significance <- "Not Significant"
volcano_data$significance[volcano_data$padj < 0.05 & volcano_data$log2FoldChange > 1] <- "Upregulated"
volcano_data$significance[volcano_data$padj < 0.05 & volcano_data$log2FoldChange < -1] <- "Downregulated"

# Highlight the top 15 genes based on absolute log2FoldChange
top_genes <- volcano_data %>%
  dplyr::arrange(desc(abs(log2FoldChange))) %>%
  dplyr::filter(!is.na(padj)) %>%
  head(15)

# Volcano plot
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 15) +
  theme_minimal() +
  labs(
    title = "Significant DEGs in non-TNBC vs TNBC",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Significance"
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )


#### Heatmap of NonTNBC and TNBC, with genes of interest ####

# Define the columns (samples) you are interested in
samples_of_interest <- c("NonTNBC1", "NonTNBC2", "NonTNBC3", "TNBC1", "TNBC2", "TNBC3")

# Subset the normalized_counts_selected matrix to include only the desired columns
normalized_counts_subset <- normalized_counts_selected[, samples_of_interest]

# Remove rows with N/A gene symbols
normalized_counts_subset <- normalized_counts_subset[!is.na(rownames(normalized_counts_subset)), ]

# Check the dimensions of the subsetted data (it should only have the samples of interest)
dim(normalized_counts_subset)

# Load the pheatmap library
library(pheatmap)

# Create a heatmap with the selected genes (normalized counts)
pheatmap(
  normalized_counts_subset, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  show_rownames = TRUE, 
  show_colnames = TRUE, 
  scale = "row", # Scaling by rows (genes) so that the gene expression is comparable across samples
  main = "Heatmap of Top DE Genes (NonTNBC vs TNBC)"
)


#### Overrepresentation analysis ####

#BiocManager::install("clusterProfiler")
#library(clusterProfiler)
#library(org.Hs.eg.db) 

# Define DE genes and the universe of all genes (both in ENSEMBL IDs)
DE_genes <- rownames(res_df)  # Differentially expressed genes
all_genes <- rownames(Counts) # All genes measured

# Perform GO enrichment analysis (BP, MF, or CC, "ALL" for all subontologies)
go_results <- enrichGO(
  gene = DE_genes,
  universe = all_genes,
  OrgDb = org.Hs.eg.db,  # Human genome annotation
  ont = "ALL",            # Specify GO subontology 
  keyType = "ENSEMBL",   # Gene ID format
  pvalueCutoff = 0.05,   # P-value threshold
  qvalueCutoff = 0.05,# Adjusted p-value threshold
  readable = TRUE        # Convert ENSEMBL to gene symbols
)

# Visualize the results
if (nrow(go_results) > 0) {
  # Dot plot (Top 10 GO terms by default)
  dotplot(go_results, showCategory = 10)
} else {
  message("No significant GO terms found.")
}


##### Gene set enrichment analysis #####

# Load necessary libraries
library(fgsea)
library(ggplot2)

# Load gene sets (e.g., MSigDB collections or custom pathways)
File_path = "~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/rnaseq_course/h.all.v2024.1.Hs.symbols.gmt"

# Step 1: Load Gene Set Collection (MSigDB)
gmt_file <- File_path  
pathways <- gmtPathways(gmt_file)

# Step 2: Ensure Unique Gene Names in gene_list
# Check and remove duplicate names from gene_list
if (any(duplicated(names(gene_list)))) {
  gene_list <- gene_list[!duplicated(names(gene_list))]
}

# Optional: Add small random noise to handle tied ranks
if (sum(duplicated(gene_list)) > 0) {
  gene_list <- gene_list + runif(length(gene_list), min = -1e-6, max = 1e-6)
}

# Step 3: Run fgseaMultilevel
fgsea_results <- fgsea(
  pathways = pathways,
  stats = gene_list,
  minSize = 15,    # Minimum gene set size
  maxSize = 500    # Maximum gene set size
)

# Step 4: Order and Inspect Results
# Order by adjusted p-value
fgsea_results <- fgsea_results[order(fgsea_results$padj), ]

# Display top pathways
head(fgsea_results)

# Step 5: Visualize Results
# Example 1: Enrichment plot for a specific pathway
plotEnrichment(pathways[["HALLMARK_HYPOXIA"]], gene_list) +
  labs(title = "Enrichment plot for HALLMARK_HYPOXIA")

# Example 2: Summary plot for top pathways
top_pathways <- fgsea_results[1:10, ]  # Top 10 pathways
ggplot(top_pathways, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0)) +
  coord_flip() +
  labs(title = "Top Pathways", x = "Pathway", y = "Normalized Enrichment Score")

# Assume fgsea_results is your GSEA results data frame
# Select top pathways (e.g., based on adjusted p-value or NES)
topPathways <- fgsea_results[order(fgsea_results$padj), ]
topPathways <- topPathways[c(1:10, (nrow(topPathways)-9):nrow(topPathways)), ] # Top and bottom 10 pathways

# Factor pathways to control order in the plot
topPathways$pathway <- factor(topPathways$pathway, levels = rev(topPathways$pathway))

# Create the summary plot
ggplot(topPathways, aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
  theme_minimal() +
  labs(
    title = "GSEA Summary Plot",
    x = "Normalized Enrichment Score (NES)",
    y = "Pathways"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

#Positive NES (e.g., red): Indicates enrichment among upregulated genes.
#Negative NES (e.g., blue): Indicates enrichment among downregulated genes.
#The pathways with the highest absolute NES are the most significantly enriched in your data.
#A high positive NES suggests that the pathway is strongly activated (upregulated).
#A high negative NES suggests that the pathway is strongly suppressed (downregulated).

##### Dot plot of GSEA Results #####

ggplot(fgsea_results, aes(x = NES, y = reorder(pathway, NES), size = -log10(padj), color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    title = "Dot Plot of GSEA Results",
    x = "Normalized Enrichment Score (NES)",
    y = "Pathways",
    size = "-log10 Adjusted P-value"
  )