library(sesame)
library(dplyr)
library(readr)
library(stringr)
library(SummarizedExperiment)

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
