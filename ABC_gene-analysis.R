library(janitor)
library(dplyr)
library(data.table)
library(biomaRt)
library(DESeq2)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(survival)
library(survminer)

setwd("/mnt/Linux_storage/ABC")

#==== 1. Read the file ====
pc_rawcounts <- fread(file = "merged_rawcounts.csv", nThread = 12)

#==== 2. Retrieve hugo symbols ====
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

##==== Get gene IDs ====
pc_rawcounts$Gene_ID <- gsub("\\..*$", "", pc_rawcounts$Gene)
geneID <- pc_rawcounts$Gene_ID
pc_rawcounts <- pc_rawcounts %>% 
  dplyr::select(Gene, Gene_ID, everything())
##==== Get Hugo symbols ====
geneHugo <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = geneID,
  mart = mart
)

#==== 3. Merge hugo symbols to the raw_counts ====
pc_rawcounts <- merge.data.frame(pc_rawcounts, geneHugo, 
                                 by.x = "Gene_ID", by.y = "ensembl_gene_id")
pc_rawcounts <- as.data.frame(pc_rawcounts)
pc_rawcounts <- pc_rawcounts %>% 
  dplyr::rename(Gene_Name = hgnc_symbol) %>% 
  dplyr::select(Gene_Name, everything(), -Gene, -Gene_ID)

#==== 4. Filter out for only ABC genes ====
abc_genes <- c(
  # ABCA
  "ABCA1", "ABCA2", "ABCA3", "ABCA4", "ABCA5", "ABCA6",
  "ABCA7", "ABCA8", "ABCA9", "ABCA10", "ABCA12", "ABCA13",
  
  # ABCB
  "ABCB1", "ABCB2", "ABCB3", "ABCB4", "ABCB5", "ABCB6",
  "ABCB7", "ABCB8", "ABCB9", "ABCB10", "ABCB11",
  
  # ABCC (including CFTR)
  "ABCC1", "ABCC2", "ABCC3", "ABCC4", "ABCC5", "ABCC6",
  "CFTR",   # functionally ABCC7
  "ABCC8", "ABCC9", "ABCC10", "ABCC11", "ABCC12",
  
  # ABCD
  "ABCD1", "ABCD2", "ABCD3", "ABCD4",
  
  # ABCE, ABCF
  "ABCE1",
  "ABCF1", "ABCF2", "ABCF3",
  
  # ABCG
  "ABCG1", "ABCG2", "ABCG4", "ABCG5", "ABCG8"
)

pc_rawcounts <- pc_rawcounts %>% 
  filter(Gene_Name %in% abc_genes) %>% 
  tibble::column_to_rownames("Gene_Name")

#==== 5. Clean the patient ids and make DESeq2 ready ====
##==== Retrieve the sample type code ====
samples <- colnames(pc_rawcounts)

sample_type_code <- str_sub(samples, 14, 15)

coldata <- data.frame(
  sample = samples,
  type = factor(
    case_when(
      sample_type_code %in% c("01", "02") ~ "Tumor",
      sample_type_code == "11" ~ "Normal",
      
      TRUE                     ~ "Other"
    ),
    levels = c("Normal", "Tumor")
  ),
  tcga_code = sample_type_code,
  row.names = samples
)

coldata <- coldata[!is.na(coldata$type), ]

count_mat <- as.matrix(pc_rawcounts) # Converting counts file to a matrix data type

count_mat <- count_mat[, colnames(count_mat) %in% rownames(coldata)] # Remove NA samples from count matrix

coldata <- coldata[colnames(count_mat), ] #Rearrange the sample ids in same order in both tables

#==== 6. Perform DESeq2 analysis ====
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ type
)

dds <- DESeq(dds)
res <- results(dds)

result <- as.data.frame(res)

#==== 7. Plot data ====
##==== Volcano plot ====
png("plots/ABC_volcano plot.png", width = 10, height = 10, res = 300, units = "in")
EnhancedVolcano(result, 
                lab = rownames(result),
                x = "log2FoldChange",
                y = "padj"
                )
dev.off()

#==================================
#       Further Exploration
#==================================

#==== 1. Split the raw counts file into different cancers ====
cancers <- c("BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

#==== 2. Read the Clinical Data ====
clinDat <- openxlsx::read.xlsx("TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR")

patient_lists <- list()
for (can in cancers) {
  patient_lists[[can]] <- clinDat$bcr_patient_barcode[clinDat$type == can]
}

#==== 3. Run DESeq2 for all cancers individually ====
all_results <- list()
can_vst <- list()
for (can in cancers) {
  
  message(can)
  
  samples <- colnames(pc_rawcounts)[str_sub(colnames(pc_rawcounts), 1, 12) %in% patient_lists[[can]]]
  
  sample_type_code <- str_sub(samples, 14, 15)

  coldata <- data.frame(
    sample = samples,
    type = factor(
      case_when(
        sample_type_code %in% c("01", "02") ~ "Tumor",
        sample_type_code == "11" ~ "Normal",

        TRUE                     ~ "Other"
      ),
      levels = c("Normal", "Tumor")
    ),
    tcga_code = sample_type_code,
    row.names = samples
  ) %>% 
    filter(type != "Other")
  
  count_mat <- as.matrix(pc_rawcounts[, rownames(coldata)])
  
  # Perform DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = ~ type
  )
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  result_df <- as.data.frame(res) %>% 
    tibble::rownames_to_column(var = "gene") %>% 
    dplyr::select(gene, log2FoldChange, pvalue, padj)
  ## This command is used in case vst doesnt  work due to very less number of genes with at least 5 counts
  vst_df <- tryCatch(vst(dds, blind = FALSE),
                     error = function(e) varianceStabilizingTransformation(dds, blind = FALSE))
  vst_df <- assay(vst_df) ## To get the normalised counts
  
  ## Make rownames into patient IDs (first 12 chars of rownames)
  vst_df <- vst_df %>%
    as.data.frame() %>%
    # transpose so patients are rows and genes are columns
    t() %>%
    as.data.frame() %>%
    # create patient_id from rownames
    tibble::rownames_to_column(var = "patient_id")
  
  ## If you really need to trim rownames to 12 characters:
  vst_df$patient_id <- stringr::str_sub(vst_df$patient_id, 1, 12)
  
  all_results[[can]] <- result_df
  can_vst[[can]] <- as.data.frame(vst_df)
}

 #==== 4. Merge and filter deseq data ====
pan_cancer_results <- bind_rows(all_results, .id = "cancer")
pan_cancer_results <- pan_cancer_results[pan_cancer_results$padj <= 0.05,]

#==== 5. Plot heatmap ====
heatmap_mat <- pan_cancer_results %>% 
  dplyr::select(gene, cancer, log2FoldChange) %>% 
  tidyr::pivot_wider(names_from = cancer, values_from = log2FoldChange) %>% 
  tibble::column_to_rownames("gene") %>% 
  as.matrix()

heatmap_mat[is.na(heatmap_mat)] <- 0 # Convert all NA to 0
heatmap_mat <- t(heatmap_mat)        # Transpose: cancers as rows, genes as columns

# Filter: keep genes where ANY cancer has |value| > 1.5
heatmap_mat <- heatmap_mat[, apply(abs(heatmap_mat) > 1.5, 2, any)]

heatmap_mat_display <- round(heatmap_mat, 1)

# Create directory if needed
if (!dir.exists("plots")) dir.create("plots")

png("plots/heatmap_pancan.png", width = 25, height = 13, res = 300, units = "in")
Heatmap(heatmap_mat_display,
        name = "log2FC",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", heatmap_mat_display[i, j]), 
                    x, y, gp = gpar(fontsize = 14))
        },
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 18),
          labels_gp = gpar(fontsize = 15),
          legend_height = unit(8, "cm"),
          legend_width = unit(2, "cm")
        ),
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 18),
        row_names_gp = gpar(fontsize = 18))
dev.off()

#=================================
#   Per cancer survival graph
#=================================

#==== 1. Read the Clinical Data ====
clinDat <- openxlsx::read.xlsx("TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR")

clinDat <- clinDat %>% 
  dplyr::select(bcr_patient_barcode, type, gender, vital_status, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time)

##==== 1.1 Gender specific filtering ====
male_clin <- clinDat %>% 
  filter(gender == "MALE")

female_clin <- clinDat %>% 
  filter(gender == "FEMALE")

##==== 1.2 Cancer specific filtering ====
cancers_clin <- c("BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

pancan_clin <- list()
for (can in cancers_clin) {
  
  message(can)
  
  pancan_clin[[can]] <- clinDat %>% 
    filter(type == can)
}


##==== 1.3 Chemotherapy resistant vs non resistant ====


#==== 2. Plot survival graphs for each cancer ====

##==== 2.1 Prepare the dataset for plotting ====
surv_mat <- list()

for (can in cancers_clin) {
  
  message(can)
  
  genes <- colnames(can_vst[[can]])[2:ncol(can_vst[[can]])]
  
  # initialize with the correct number of rows
  surv_mat[[can]] <- data.frame(
    patient_id = can_vst[[can]]$patient_id,
    stringsAsFactors = FALSE
  )
  
  for (gene in genes) {
    gene_median <- median(can_vst[[can]][[gene]], na.rm = TRUE)
    
    surv_mat[[can]][[gene]] <- ifelse(
      can_vst[[can]][[gene]] >= gene_median,
      "High",
      "Low"
    )
    surv_mat[[can]][[gene]] <- as.factor(surv_mat[[can]][[gene]])
  }
  
  # Merge clinical data
  surv_cols <- c("OS", "DSS", "DFI", "PFI")
  for (surv_type in surv_cols){
    col_time <- paste0(surv_type, ".time")
    pancan_clin[[can]][[col_time]] <- pancan_clin[[can]][[col_time]] / 30
  }
  
  surv_mat[[can]] <- merge.data.frame(surv_mat[[can]], pancan_clin[[can]], 
                                         by.x = "patient_id", by.y = "bcr_patient_barcode")
  
}

##==== 2.2 Create a fit for plotting and plot ====
survivalDat <- data.frame(
  Cancer = character(),
  Gene = character(),
  SurvType = character(),
  hazard_ratio = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (can in cancers_clin) {
  message(can)
  genes <- colnames(can_vst[[can]])[2:ncol(can_vst[[can]])]
  
  for (gene in genes) {
    for (surv_type in surv_cols) {
      
      timeCol  <- paste0(surv_type, ".time")
      eventCol <- surv_type
      
      # ensure gene grouping is a factor with >= 2 levels
      grp <- surv_mat[[can]][[gene]]
      grp <- as.factor(grp)
      
      # skip if only one group present
      if (nlevels(grp) < 2) {
        message("Skipping ", can, " - ", gene, " - ", surv_type,
                " (only ", nlevels(grp), " level)")
        next
      }
      
      surv_mat[[can]][[gene]] <- grp
      
      surv_formula <- as.formula(
        paste0("Surv(", timeCol, ", ", eventCol, ") ~ ", gene)
      )
      
      # KM fit
      surv_obj <- surv_fit(surv_formula, data = surv_mat[[can]])
      
      # Cox model
      cox_model   <- coxph(surv_formula, data = surv_mat[[can]])
      cox_summary <- summary(cox_model)
      
      hr   <- round(cox_summary$coefficients[1, "exp(coef)"], 3)
      pval <- signif(cox_summary$coefficients[1, "Pr(>|z|)"], 3)
      
      survivalDat <- rbind(survivalDat, data.frame(
        Cancer = can,
        Gene = gene,
        SurvType = surv_type,
        hazard_ratio = hr,
        pvalue = pval
      ))
    }
  }
}

write.csv(survivalDat, "Survival_Data.csv", row.names = FALSE)

# Filter survival data
survivalDat <- survivalDat %>% 
  filter(pvalue < 0.05) %>% 
  filter(hazard_ratio > 1.5 | hazard_ratio < 0.6)

survivalDat_sub <- survivalDat %>% 
  filter(Gene %in% colnames(heatmap_mat))

# Plot filtered results
for (i in 1:nrow(survivalDat)) {
  can <- survivalDat_sub$Cancer[i]
  gene <- survivalDat_sub$Gene[i]
  surv_type <- survivalDat_sub$SurvType[i]
  hr <- survivalDat_sub$hazard_ratio[i]
  pval <- survivalDat_sub$pvalue[i]
  
  message("Plotting: ", can, " - ", gene, " - ", surv_type)
  
  # Create directory: plots/Survival_plots/Cancer/SurvType/
  cancer_dir <- file.path("plots/Survival_plots", can)
  surv_dir <- file.path(cancer_dir, surv_type)
  
  if (!dir.exists(surv_dir)) {
    dir.create(surv_dir, recursive = TRUE)
  }
  
  # Survival formula
  timeCol <- paste0(surv_type, ".time")
  eventCol <- surv_type
  surv_formula <- as.formula(paste0("Surv(", timeCol, ", ", eventCol, ") ~ ", gene))
  
  # Fit survival model
  surv_obj <- surv_fit(surv_formula, data = surv_mat[[can]])
  
  # Plot title
  plot_title <- paste0(
    can, " | ", gene, "\n",
    "HR = ", hr, ", p = ", format(pval, scientific = TRUE)
  )
  
  # Save plot filename
  plot_name <- paste0(gene, ".png")
  png_file <- file.path(surv_dir, plot_name)
  
  png(png_file, width = 10, height = 8, units = "in", res = 300)
  
  print(
    ggsurvplot(
      surv_obj,
      data = surv_mat[[can]],
      risk.table = TRUE,
      pval = TRUE,
      conf.int = TRUE,
      title = plot_title,
      legend.title = "Expression",
      legend.labs = c("Low", "High"),
      xlab = paste(surv_type, "(months)"),
      ylab = "Survival probability",
      ggtheme = theme_minimal()
    )
  )
  dev.off()
}
