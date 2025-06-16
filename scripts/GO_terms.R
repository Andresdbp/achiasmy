# Andres Barboza
# Jun 11, 2024
# andresdbp@tamu.edu

library(biomaRt)
library(dplyr)
library(tidyr)
library(tibble)
library(topGO)
library(ggplot2)

# # Databases
# humandb <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# 
# # Get gene annotations including chromosome location
# human <- getBM(attributes = c('go_id', 'ensembl_gene_id', 'chromosome_name'),
#                           mart = ensembl)
# 


# Step 1: Connect to Ensembl BioMart
## Human
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
## Fruit Fly
ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")


# Step 2: Verify available attributes
attributes <- listAttributes(ensembl)
head(attributes)

# Verify that 'ensembl_gene_id' and 'chromosome_name' are present in the attributes list
if (!("ensembl_gene_id" %in% attributes$name) | !("chromosome_name" %in% attributes$name)) {
  stop("The required attributes 'ensembl_gene_id' or 'chromosome_name' are not available in the dataset.")
}

# Step 3: Fetch gene annotations
gene_annotations <- getBM(
  attributes = c('ensembl_gene_id', 'chromosome_name'),
  mart = ensembl
)

# Step 4: Filter to keep only numbered chromosomes and X, Y
valid_chromosomes <- c(as.character(1:99), "X", "Y", "W", "Z")
gene_annotations <- gene_annotations %>%
  filter(chromosome_name %in% valid_chromosomes)

# View the first few rows to ensure data was retrieved correctly and filtered
head(gene_annotations)

# Step 5: Fetch GO annotations
go_annotations <- getBM(
  attributes = c('ensembl_gene_id', 'go_id'),
  mart = ensembl
)

# View the first few rows to ensure data was retrieved correctly
head(go_annotations)

# Step 6: Map genes to GO terms
gene_go <- merge(gene_annotations, go_annotations, by = 'ensembl_gene_id')
head(gene_go)

# Step 7: Count GO term occurrences
go_counts <- gene_go %>%
  group_by(chromosome_name, go_id) %>%
  summarize(count = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = chromosome_name, values_from = count, values_fill = list(count = 0))
head(go_counts)

# Ensure numeric columns in go_counts and set go_id as row names
go_counts <- go_counts %>%
  column_to_rownames("go_id") %>%
  mutate(across(everything(), as.numeric))
go_counts <- go_counts[-1,]

# Function to perform Fisher's exact test for enrichment analysis
perform_fishers_test <- function(go_id, chromosome, go_counts, total_genes) {
  go_counts_chromosome <- go_counts[go_id, , drop = FALSE]
  
  a <- go_counts_chromosome[1, chromosome] # Count of genes with the GO term on the specific chromosome
  b <- sum(go_counts_chromosome[1, ]) - a # Count of genes with the GO term on other chromosomes
  c <- total_genes[chromosome] - a # Count of genes without the GO term on the specific chromosome
  d <- sum(total_genes) - total_genes[chromosome] - b # Count of genes without the GO term on other chromosomes
  
  # Ensure values are non-negative and finite
  a <- max(0, a)
  b <- max(0, b)
  c <- max(0, c)
  d <- max(0, d)
  if (!is.finite(a) || !is.finite(b) || !is.finite(c) || !is.finite(d)) {
    return(NA)
  }
  
  matrix <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  test <- fisher.test(matrix, alternative = "greater")
  return(test$p.value)
}

# Step 8: Prepare total genes count per chromosome
total_genes <- gene_annotations %>%
  group_by(chromosome_name) %>%
  summarize(total_genes = n(), .groups = 'drop') %>%
  as.data.frame()
total_genes <- setNames(total_genes$total_genes, total_genes$chromosome_name)

# Step 9: Enrichment analysis for a specific chromosome (e.g., Y chromosome)
chromosome_of_interest <- "X"

p_values <- sapply(rownames(go_counts), function(go_id) {
  perform_fishers_test(go_id, chromosome_of_interest, go_counts, total_genes)
})

adjusted_p_values <- p.adjust(p_values, method = "fdr")

enrichment_results <- data.frame(
  go_id = rownames(go_counts),
  p_value = p_values,
  adjusted_p_value = adjusted_p_values
) %>%
  arrange(adjusted_p_value)

# View the top GO terms most associated with the chromosome of interest
head(enrichment_results)

# Step 10: Visualize results
significant_terms <- enrichment_results %>%
  filter(adjusted_p_value < 0.05)

ggplot(significant_terms, aes(x = reorder(go_id, adjusted_p_value), y = -log10(adjusted_p_value))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = paste("Significant GO Terms Associated with Chromosome", chromosome_of_interest),
       x = "GO Term",
       y = "-log10(adjusted p-value)")

go_descriptions <- getBM(
  attributes = c("go_id", "name_1006"),
  filters = "go",
  values = significant_terms$go_id,
  mart = ensembl
)

sigf_go_desc <- merge(go_descriptions, significant_terms, by = 'go_id')

human_y <- sigf_go_desc
human_x <- sigf_go_desc
dmel_y <- sigf_go_desc
dmel_x <- sigf_go_desc

#TEST
ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

gene_annotations <- getBM(
  attributes = c('ensembl_gene_id', 'chromosome_name'),
  mart = ensembl
)
valid_chromosomes <- c(as.character(1:99), "X", "Y", "W", "Z")
gene_annotations <- gene_annotations %>%
  filter(chromosome_name %in% valid_chromosomes)
go_annotations <- getBM(
  attributes = c('ensembl_gene_id', 'go_id', 'name_1006'),
  mart = ensembl
)
gene_go_name <- merge(gene_annotations, go_annotations, by = 'ensembl_gene_id')
y_go_terms <- unique(gene_go_name[which(gene_go_name$chromosome_name == "Y"),c(3,4)])
