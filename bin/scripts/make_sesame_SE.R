library(sesame)
library(SummarizedExperiment)
library(sesameData)
library(ggrepel)
library(cowplot)
library(tidyverse)

betas <- openSesame(snakemake@params[['idat_dir']], prep=snakemake@params[['sesame_prep']]) # BPPARAM = BiocParallel::MulticoreParam(snakemake@threads))

samplesheet <- read_csv(snakemake@input[['samplesheet']], skip = snakemake@params[['samplesheet_skip']]) %>%
    dplyr::mutate(sample = paste0(Sentrix_ID, "_", Sentrix_Position))

meta <- read_tsv(snakemake@input[['meta']]) %>%
    dplyr::left_join(., samplesheet, by="Sample_Name") %>%
    tibble::column_to_rownames("sample")

se <- SummarizedExperiment(assays=list(betas=betas),
        colData = meta[colnames(betas), ])

# get probe locations from manifest
rowdat <- sesameData_getManifestGRanges(snakemake@params[['manifest_id']])
stopifnot(sum(is.na(start(rowdat))) == 0)
stopifnot(sum(is.na(seqnames(rowdat))) == 0)

# check that all the manifest probes are in experiment probes
stopifnot(all(names(rowdat) %in% rownames(se)))

# print out how many of the experiment probes have a match in the manifest
message(sum(rownames(se) %in% names(rowdat)), " probes out of ", nrow(se), " have matches in the array manifest.")

# subset for probes with manifest match
se <- se[names(rowdat), ]
rowRanges(se) <- rowdat

saveRDS(se, snakemake@output[[1]])


### generate PCA plot
plot_pca_from_methylation_array <- function(se = se, group = "group") {
  noNA_rows <- rowSums(is.na(assay(se))) == 0
  noNA_data <- se[noNA_rows, ]

  pca <- prcomp(t(assay(noNA_data)))
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
