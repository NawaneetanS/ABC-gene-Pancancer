# Function to extract significance and regulation direction for multiple genes
library(TCGAplot)
paired_tpm <- get_all_paired_tpm()

extract_gene_significance <- function(genes, pcancers = c("BLCA","BRCA","COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")) {
  
  library(dplyr)
  library(ggpubr)
  
  results <- list()
  
  for(gene in genes) {
    for(cancer in pcancers) {
      # Filter data for specific cancer type
      data <- subset(paired_tpm, Cancer == cancer) %>%
        dplyr::select(Cancer, Group, ID, all_of(gene))
      
      # Skip if no data for this gene/cancer combination
      if(nrow(data) == 0) next
      
      # Statistical comparison
      fml <- as.formula(paste0(gene, "~Group"))
      cmp <- ggpubr::compare_means(fml, data = data, paired = TRUE)
      
      # Calculate regulation direction
      means <- aggregate(data[[gene]], by = list(data$Group), FUN = mean, na.rm = TRUE)
      colnames(means) <- c("Group", "Mean")
      
      direction <- "No_data"
      tumor_mean <- NA
      normal_mean <- NA
      
      if(nrow(means) == 2) {
        tumor_mean <- means$Mean[means$Group == "Tumor"]
        normal_mean <- means$Mean[means$Group == "Normal"]
        
        if(tumor_mean > normal_mean) {
          direction <- "Upregulated"
        } else {
          direction <- "Downregulated"
        }
      }
      
      # Store results
      results[[length(results) + 1]] <- data.frame(
        Gene = gene,
        Cancer = cancer,
        P_value = cmp$p,
        P_signif = cmp$p.signif,
        Direction = direction,
        Tumor_mean = tumor_mean,
        Normal_mean = normal_mean,
        Log2_fold_change = log2(tumor_mean / normal_mean),
        Sample_size = nrow(data)/2,  # Divided by 2 since it's paired data
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all results
  df_results <- do.call(rbind, results)
  return(df_results)
}

# Define genes you want to analyze
genes_to_test <- c(
  "ABCA1", "ABCA2", "ABCA3", "ABCA4", "ABCA5", "ABCA6", "ABCA7", "ABCA8", 
  "ABCA9", "ABCA10", "ABCA12", "ABCA13", "ABCB1", "ABCB4", "ABCB5", "ABCB6", 
  "ABCB7", "ABCB8", "ABCB9", "ABCB10", "ABCB11", "ABCC1", "ABCC2", "ABCC3", 
  "ABCC4", "ABCC5", "ABCC6", "ABCC8", "ABCC9", "ABCC10", "ABCC11", "ABCC12", 
  "ABCD1", "ABCD2", "ABCD3", "ABCD4", "ABCF1", "ABCF2", "ABCF3", "ABCG1", 
  "ABCG2", "ABCG4", "ABCG5", "ABCG8"
)

# Extract significance data
significance_df <- extract_gene_significance(genes_to_test)

# View the results
head(significance_df)

# Export to CSV for Excel
write.csv(significance_df, file = "gene_cancer_significance.csv", row.names = FALSE)

# You can also filter for significant results only
significant_only <- significance_df[significance_df$P_signif != "ns", ]
write.csv(significant_only, file = "significant_genes_only.csv", row.names = FALSE)

