Create a rmd file with a name preferably ‘qc_batchcorrection’ and keep adding chunks for each sub-step.

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(R.utils)  
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(openxlsx)
library(dplyr)
library(scales)
library(HGNChelper)
library(dittoSeq)
library(harmony)
library(batchelor)
library(zellkonverter)
library(SingleCellExperiment)


# Firstly, decide a path where all outputs will be stored. Run in console.

setwd("C:/Users/smriti/Desktop/demo/WGSE279086/W2/")
getwd()

# how does raw data looks like ?? You may look for any sample.
```{r}
barcodes <- read.delim(gzfile("raw_geo_data/GSM8561110/barcodes.tsv.gz"), header = FALSE)
head(barcodes)

features <- read.delim(gzfile("raw_geo_data/GSM8561110/features.tsv.gz"), header = FALSE)
head(features)

matrix_lines <- readLines(gzfile("raw_geo_data/GSM8561110/matrix.mtx.gz"), n = 10)
matrix_lines
```

#Load all required libraries. Do R.version.string to ensure the latest version is activated.

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(R.utils)  
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(openxlsx)
library(scales)
library(HGNChelper)
library(dittoSeq)
library(harmony)
library(batchelor)
library(zellkonverter)
library(SingleCellExperiment)
```

#Manually create output and plots folder n define their paths below.
```{r Define the directories}
inputDir <- "input/"
outputDir <- "output/"
plotDir <- "plots/"
```

# This chunk looks inside a folder (input) and finds all files that contain single-cell data saved as Seurat objects (.rds files). Each of these files corresponds to one sample. It then loads all these Seurat objects into R and stores them together in a list, so that multiple samples can be handled at once. After loading, each Seurat object is given a meaningful name (usually based on sample) so that it is easy to identify which object belongs to which sample. Finally, it prints the names to confirm that all samples were loaded correctly.

```{r Read all the individual seurat objects}
rds_files <- list.files(inputDir, pattern = "_seurat\\.rds$", full.names = TRUE)
seurat_list <- lapply(rds_files, readRDS)
names(seurat_list) <- gsub("_seurat\\.rds$", "", basename(rds_files))
print(names(seurat_list))
```
# sanity check
```{r}
sapply(seurat_list, function(obj) {
  sum(duplicated(rownames(GetAssayData(obj, layer = "counts"))))
})

# Got 0s in all samples (meaning no duplicate genes are present in any of the samples) and thus no need to run and save clean_seurat_list. Proceed with seurat_list. But because in downstream, we have 'clean_seurat_list' everywhere, we'll do

clean_seurat_list <- seurat_list
print(names(clean_seurat_list))

# clean_seurat_list <- lapply(seurat_list, function(obj) {
#   # Get unique genes from RNA assay
#   counts_mat <- GetAssayData(obj, layer = "counts")
#   unique_genes <- !duplicated(rownames(counts_mat))
#   counts_mat_unique <- counts_mat[unique_genes, ]
#   obj[["RNA"]] <- CreateAssayObject(counts_mat_unique)
#   return(obj)
# })
```

# Merging all seurat objects into one and integrating metadata variables and QC 
```{r Merge all Seurat objects into one}
# Update orig.ident to GSM IDs inside each Seurat object
# for (s in names(clean_seurat_list)) {
#   clean_seurat_list[[s]]$orig.ident <- s  # s is GSM ID
# }

seurat_combined <- merge(clean_seurat_list[[1]], y = clean_seurat_list[-1], add.cell.ids = names(clean_seurat_list))
unique(seurat_combined$orig.ident)

table(seurat_combined$sample)
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)
metadata

# make a named vector: names = Sample, values = Condition
condition_map <- setNames(metadata$Condition, metadata$Sample)

seurat_combined@meta.data$Condition <- unname(condition_map[as.character(seurat_combined@meta.data$sample)])

table(seurat_combined$Condition)
sum(is.na(seurat_combined$Condition))

saveRDS(seurat_combined, file =  paste0(outputDir,"01_GSE279086_seurat_combined.rds"))
#seurat_combined <- readRDS(file = paste0(outputDir, "01_GSE279086_seurat_combined.rds"))

# Explore Combined Seurat Object}
class(seurat_combined)            # 'Seurat'
Assays(seurat_combined)           # (RNA) assays available in the rds
DefaultAssay(seurat_combined)     # default/active (RNA) assay in the rds
dim(seurat_combined)              # 28317 genes/features x 181249 cells
unique(seurat_combined$Condition)
```

# How is the data looking now??
```{r}
View(seurat_combined@meta.data)
Layers(seurat_combined[["RNA"]])
mat1 <- LayerData(seurat_combined, assay = "RNA", layer = "counts.GSM8561110")
dim(mat1)
mat1[1:10, 1:10]
```


#.................Before QC..............#
```{r}
study_id <- "GSE279086"

save_violin_plots_separate <- function(
  seurat_obj,
  plotDir,
  study_id,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
) {

  if (!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = TRUE)
  }

  meta <- seurat_obj@meta.data
  meta$plot_sample <- as.factor(meta$sample)

  for (feat in features) {

    if (!feat %in% colnames(meta)) {
      warning(paste("Skipping", feat, "- not found in metadata"))
      next
    }

    p <- ggplot(meta, aes(x = plot_sample, y = .data[[feat]])) +
      geom_violin(trim = TRUE, fill = "steelblue", alpha = 0.7) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
      labs(
        title = feat,
        x = "Sample",
        y = feat
      ) +
      theme_bw(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5)
      )

    ggsave(
      filename = file.path(
        plotDir,
        paste0(study_id, "_preQC_", feat, "_violin.png")
      ),
      plot  = p,
      width = 16,
      height = 9,
      dpi   = 400,
      bg    = "white"
    )
  }
}

save_violin_plots_separate(seurat_combined, plotDir, study_id)

# Scatter plot of nCount vs nFeature with marginal histograms
density_scatter_plot <- function(seurat_obj, filename){
  df <- data.frame(
    log1p_nCount_RNA = log1p(seurat_obj$nCount_RNA),
    log1p_nFeature_RNA = log1p(seurat_obj$nFeature_RNA),
    Condition = seurat_obj$Condition
  )
  
  p <- ggplot(df, aes(x = log1p_nCount_RNA, y = log1p_nFeature_RNA, colour = Condition)) +
    geom_point(alpha = 0.3, size = 0.5) +
    theme_minimal() +
    theme(legend.position.inside = c(0.05, 0.95), legend.justification = c("left", "top"), legend.key.size = unit(0.5, 'cm')) +
    guides(color = guide_legend(override.aes = list(size = 5))) + # 'size' here controls the symbol size
    labs(x = "log1p(nCount_RNA)", y = "log1p(nFeature_RNA)")
  
  p <- ggMarginal(p, type = "histogram", fill = "skyblue", bins = 40)
  ggsave(
  filename,
  plot = p,
  width = 8,
  height = 10,
  dpi = 400,
  bg = "white"
)

}

density_scatter_plot(seurat_combined, file.path(plotDir, paste0(study_id, "_preQC_density-scatter.png")))
```

#................. QC Functions ..............#
```{r Writing QC Functions}
fil_nCounts <- function(seurat_obj, threshold_mahalanobis){
  ncounts_data <- data.frame(nCount_RNA = seurat_obj@meta.data$nCount_RNA)
  
  mahalanobis_dist <- mahalanobis(
    x = ncounts_data,
    center = mean(ncounts_data$nCount_RNA),
    cov = var(ncounts_data$nCount_RNA)
  )
  
  seurat_obj$mahal_dist_nCount <- mahalanobis_dist
  mahal.fil <- qchisq(threshold_mahalanobis, df = 1)
  message(paste0("Mahalanobis threshold: ", round(mahal.fil, 3)))
  
  cells_to_remove <- rownames(seurat_obj@meta.data)[seurat_obj$mahal_dist_nCount > mahal.fil]
  message(paste0("Removing ", length(cells_to_remove), " cells (outliers by nCount_RNA)"))
  
  subset(seurat_obj, cells = setdiff(colnames(seurat_obj), cells_to_remove))
}


filter_nFeatures <- function(seurat_obj, min_feat, max_feat, ...) {
  return(subset(seurat_obj, subset = nFeature_RNA > min_feat & nFeature_RNA < max_feat))
}

filter_mt <- function(seurat_obj, max_mt) {
  return(subset(seurat_obj, subset = percent.mt < max_mt))
}

filter_rb <- function(seurat_obj, max_rb) {
  return(subset(seurat_obj, subset = percent.rb < max_rb))
}
```


# lets explore what are the most suitable QC parameters for our data
```{r}
qc_grid_search <- function(
  seurat_obj,
  min_features = 200,
  max_features = 7000,
  mt_values = c(20, 30),
  rb_values = c(10, 20),
  mahal_values = c(0.95, 0.98)
) {
  
  # before QC per sample
  df_pre_qc <- as.data.frame(table(seurat_obj$sample))
  colnames(df_pre_qc) <- c("Sample_ID", "Cells_Before_QC")
  
  results_list <- list()
  filtered_objects <- list()
  counter <- 1
  
  for (mt_cut in mt_values) {
    for (rb_cut in rb_values) {
      for (mahal_cut in mahal_values) {
        
        message(
          paste0(
            "Running QC with mt<", mt_cut,
            ", rb<", rb_cut,
            ", mahal=", mahal_cut
          )
        )
        
        # apply filters step by step
        tmp <- filter_nFeatures(seurat_obj, min_features, max_features)
        tmp <- filter_mt(tmp, mt_cut)
        tmp <- fil_nCounts(tmp, mahal_cut)
        tmp <- filter_rb(tmp, rb_cut)
        
        # after QC per sample
        df_post_qc <- as.data.frame(table(tmp$sample))
        colnames(df_post_qc) <- c("Sample_ID", "Cells_After_QC")
        
        # merge sample-wise counts
        qc_table <- merge(df_pre_qc, df_post_qc, by = "Sample_ID", all.x = TRUE)
        qc_table$Cells_After_QC[is.na(qc_table$Cells_After_QC)] <- 0
        
        # add parameter columns
        qc_table$min_features <- min_features
        qc_table$max_features <- max_features
        qc_table$max_percent_mt <- mt_cut
        qc_table$max_percent_rb <- rb_cut
        qc_table$threshold_mahalanobis <- mahal_cut
        qc_table$Cells_Removed <- qc_table$Cells_Before_QC - qc_table$Cells_After_QC
        qc_table$Percent_Removed <- round(
          100 * qc_table$Cells_Removed / qc_table$Cells_Before_QC, 2
        )
        
        # add total summary row separately
        total_row <- data.frame(
          Sample_ID = "TOTAL",
          Cells_Before_QC = sum(qc_table$Cells_Before_QC, na.rm = TRUE),
          Cells_After_QC = sum(qc_table$Cells_After_QC, na.rm = TRUE),
          min_features = min_features,
          max_features = max_features,
          max_percent_mt = mt_cut,
          max_percent_rb = rb_cut,
          threshold_mahalanobis = mahal_cut,
          Cells_Removed = sum(qc_table$Cells_Removed, na.rm = TRUE),
          Percent_Removed = round(
            100 * sum(qc_table$Cells_Removed, na.rm = TRUE) /
              sum(qc_table$Cells_Before_QC, na.rm = TRUE), 2
          )
        )
        
        qc_table_full <- rbind(qc_table, total_row)
        
        # save
        combo_name <- paste0(
          "mt", mt_cut, "_rb", rb_cut, "_mahal", mahal_cut
        )
        results_list[[combo_name]] <- qc_table_full
        filtered_objects[[combo_name]] <- tmp
        
        counter <- counter + 1
      }
    }
  }
  
  # combine all results into one big table
  final_results <- do.call(rbind, results_list)
  rownames(final_results) <- NULL
  
  return(list(
    qc_summary_table = final_results,
    filtered_objects = filtered_objects
  ))
}

qc_results <- qc_grid_search(
  seurat_obj = seurat_combined,
  min_features = 200,
  max_features = 7000,
  mt_values = c(20, 30),
  rb_values = c(10, 20),
  mahal_values = c(0.95, 0.98)
)

qc_summary <- qc_results$qc_summary_table
qc_totals <- qc_summary[qc_summary$Sample_ID == "TOTAL", ]

qc_totals[, c(
  "max_percent_mt",
  "max_percent_rb",
  "threshold_mahalanobis",
  "Cells_Before_QC",
  "Cells_After_QC",
  "Cells_Removed",
  "Percent_Removed"
)]
```


#.................Perform QC..............#
# These parameters define quality control rules that help distinguish real, healthy single cells from empty droplets, dying cells, or technical artifacts.

```{r QC thresholds}
study_id <- "GSE279086"
min_features    <- 200  # A cell must express at least 200 genes to be considered real. <200 --> empty droplets/cause most noise
max_features    <- 7000 # Cells with more than 7000 genes are removed. Extremely high gene counts often indicate doublets
max_percent_mt  <- 20 
max_percent_rb  <- 20 # Remove cells dominated by ribosomal expression. High rRNA --> low info content/technical bias
threshold_mahalanobis <- 0.95  # Remove the worst 5% of cells that look abnormal overall (catches subtle outliers)

seurat_filtered <- filter_nFeatures(seurat_combined, min_features, max_features)

seurat_filtered <- filter_mt(seurat_filtered, max_percent_mt)

seurat_filtered <- fil_nCounts(seurat_filtered, threshold_mahalanobis)

seurat_filtered <- filter_rb(seurat_filtered,max_percent_rb)

# Cell counts before & after QC
df_pre_qc <- as.data.frame(table(seurat_combined$sample))
colnames(df_pre_qc) <- c("Sample_ID", "Cells_Before_QC")

df_post_qc <- as.data.frame(table(seurat_filtered$sample))
colnames(df_post_qc) <- c("Sample_ID", "Cells_After_QC")

# Merge into final QC table (left join behavior)
final_qc_table <- merge(
  df_pre_qc,
  df_post_qc,
  by = "Sample_ID",
  all.x = TRUE
)

final_qc_table
total_cells_before_qc <- sum(final_qc_table$Cells_Before_QC, na.rm = TRUE)
total_cells_after_qc  <- sum(final_qc_table$Cells_After_QC, na.rm = TRUE)

total_cells_before_qc
total_cells_after_qc
```


```{r Save Summary Table}
qc_metrics_file <- file.path(outputDir, paste0(study_id, "_QC_metrics.csv"))
write.csv(final_qc_table, qc_metrics_file, row.names = FALSE)
cat("Saved QC summary table to:", qc_metrics_file, "\n")
```


#................. After QC..............#
```{r}
study_id <- "GSE279086"

save_violin_plots_separate <- function(
  seurat_obj,
  plotDir,
  study_id,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
) {

  if (!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = TRUE)
  }

  meta <- seurat_obj@meta.data
  meta$plot_sample <- as.factor(meta$sample)

  for (feat in features) {

    if (!feat %in% colnames(meta)) {
      warning(paste("Skipping", feat, "- not found in metadata"))
      next
    }

    p <- ggplot(meta, aes(x = plot_sample, y = .data[[feat]])) +
      geom_violin(trim = TRUE, fill = "steelblue", alpha = 0.7) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
      labs(
        title = feat,
        x = "Sample",
        y = feat
      ) +
      theme_bw(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5)
      )

    ggsave(
      filename = file.path(
        plotDir,
        paste0(study_id, "_postQC_", feat, "_violin.png")
      ),
      plot  = p,
      width = 16,
      height = 9,
      dpi   = 400,
      bg    = "white"
    )
  }
}

save_violin_plots_separate(seurat_filtered, plotDir, study_id)

# Scatter plot of nCount vs nFeature with marginal histograms
density_scatter_plot <- function(seurat_obj, filename){
  df <- data.frame(
    log1p_nCount_RNA = log1p(seurat_obj$nCount_RNA),
    log1p_nFeature_RNA = log1p(seurat_obj$nFeature_RNA),
    Condition = seurat_obj$Condition
  )
  
  p <- ggplot(df, aes(x = log1p_nCount_RNA, y = log1p_nFeature_RNA, colour = Condition)) +
    geom_point(alpha = 0.3, size = 0.5) +
    theme_minimal() +
    theme(legend.position.inside = c(0.05, 0.95), legend.justification = c("left", "top"), legend.key.size = unit(0.5, 'cm')) +
    guides(color = guide_legend(override.aes = list(size = 5))) + # 'size' here controls the symbol size
    labs(x = "log1p(nCount_RNA)", y = "log1p(nFeature_RNA)")
  
  p <- ggMarginal(p, type = "histogram", fill = "skyblue", bins = 40)
  ggsave(
  filename,
  plot = p,
  width = 8,
  height = 10,
  dpi = 400,
  bg = "white"
)

}

density_scatter_plot(seurat_filtered, file.path(plotDir, paste0(study_id, "_postQC_density-scatter.png")))
```
