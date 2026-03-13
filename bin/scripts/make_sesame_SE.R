library(sesame)
library(SummarizedExperiment)
library(sesameData)
library(ggrepel)
library(cowplot)
library(tidyverse)

# =========== Import betas and metadata, create SummarizedExperiment =====
betas <- openSesame(snakemake@input[['idat_dir']], 
  prep=snakemake@params[['sesame_prep']]) # BPPARAM = BiocParallel::MulticoreParam(snakemake@threads))

## Read in sample sheet:
samplesheet <- read_csv(snakemake@input[['samplesheet']], 
  skip = snakemake@params[['samplesheet_skip']]) 
# add checking for column names, make sure we have three essential columns
if(!all(c("Sample_Name", "Sentrix_ID", "Sentrix_Position") %in% colnames(samplesheet)))
stop("missing one or more of the essential columns: Sample_Name, Sentrix_ID, Sentrix_Position in Samplesheet")

samplesheet <- samplesheet %>%
    dplyr::mutate(sample = paste0(Sentrix_ID, "_", Sentrix_Position))

## check for duplicated samples in sample sheet:
if (anyDuplicated(samplesheet$sample)) {
  dup <- unique(samplesheet$sample[duplicated(samplesheet$sample)])
  stop("Duplicate constructed sample IDs in samplesheet: ", paste(dup, collapse = ", "))
}

## Read in meta sheet:
# check for duplicated samples in meta sheet: 
meta <- read_tsv(snakemake@input[['meta']]) 
if(!"Sample_Name" %in% colnames(meta))
stop("Error: missing column Sample_Name in the meta sheet")

if (anyDuplicated(meta$sample)) {
  dup <- unique(meta$sample[duplicated(meta$sample)])
  stop("Duplicate sample IDs in merged metadata: ", paste(dup, collapse = ", "))
}

meta <- meta %>%
    dplyr::left_join(samplesheet, by="Sample_Name") %>%
    tibble::column_to_rownames("sample")

## Create summarized experiment
if(setequal(rownames(meta), colnames(betas))) {         
  
  se <- SummarizedExperiment(assays = list(betas = betas),
                             colData = meta[colnames(betas), ])
} else if(length(intersect(rownames(meta), colnames(betas))) >= 1) {
  
  common_samples <- intersect(rownames(meta), colnames(betas))
  se <- SummarizedExperiment(assays = list(betas = betas[, common_samples]),
                             colData = meta[common_samples, ])
  
} else {
  stop ("Can not create SummarizedExperiment: no overlapping sample IDs between meta and beta matrix.")
}

## Get probe locations from manifest
rowdat <- sesameData_getManifestGRanges(snakemake@params[['manifest_id']])
stopifnot(sum(is.na(start(rowdat))) == 0)
stopifnot(sum(is.na(seqnames(rowdat))) == 0)

common_probes <- intersect(rownames(se), names(rowdat))
message(length(common_probes), " probes out of ", nrow(se), " have matches in the array manifest.")

if (length(common_probes) == 0) {
  stop("No overlapping probes between SummarizedExperiment and manifest.")
}

se <- se[common_probes, ]
rowRanges(se) <- rowdat[common_probes]

# save the SummarizedExperiment object:
saveRDS(se, snakemake@output[[1]])


# ================= Generate PCA plot =================
plot_pca_from_methylation_array <- function(se = se, group = "group") {
  noNA_rows <- rowSums(is.na(assay(se, "betas"))) == 0
  noNA_data <- se[noNA_rows, ]

  pca <- prcomp(t(assay(noNA_data, "betas")))
  pr_comps <- data.frame(pca$x)
  pr_comps$combo <- rownames(pr_comps)
  pr_comps <- pr_comps %>% left_join(colData(se) %>% as.data.frame() %>% rownames_to_column("combo"))

  ## now merge with sample names

  prop_var <- data.frame(t(summary(pca)$importance))
  names(prop_var) = c('sd', 'prop', 'cum')
  prop_var$num = 1:nrow(prop_var)

  pca_plot <- ggplot(pr_comps, aes(x=PC1, y=PC2, color = .data[[group]], label = Sample_Name)) +
    geom_point(size = 3, alpha = 0.8)  +
    xlab(paste0("PC1 (", prop_var$prop[1] * 100, "% variance)" )) +
    ylab(paste0("PC2 (", prop_var$prop[2] * 100 , "% variance)" )) +
    geom_text_repel(size = 3, max.overlaps = Inf, show.legend = F) +
    theme_cowplot()

  return(pca_plot)
}

## save the plot 
pca_p <- plot_pca_from_methylation_array(se = se, group = "group")
dir.create("analysis/plots")
save_plot("analysis/plots/pca_plot.pdf", pca_p, base_height = 6, base_width = 7.5)
