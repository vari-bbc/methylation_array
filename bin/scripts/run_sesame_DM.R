library(sesame)
library(dplyr)
library(readr)
library(stringr)
library(SummarizedExperiment)

se <- readRDS(snakemake@input[[1]])
se <- se[, se$group %in% c(snakemake@wildcards[['group1']], snakemake@wildcards[['group2']])]
colData(se)$group <- relevel(factor(colData(se)$group), snakemake@wildcards[['group2']])
message(paste0("Keeping only samples from the groups: ", paste(unique(colData(se)$group), collapse=", ")))
se_ok <- (checkLevels(assay(se), colData(se)$group))
message(paste0(sum(se_ok), " out of ", nrow(se), " CpGs passed."))
se <- se[se_ok, ]

message("Running DML.")
smry <- DML(se, ~group)
test_result <- summaryExtractTest(smry)

message("Running DMR.")
merged <- DMR(se, smry, paste0("group", snakemake@wildcards[['group1']]), platform=snakemake@params[['array_type']])

saveRDS(list(se=se, dml_smry=smry, dml_res=test_result, dmr=merged), snakemake@output[[1]])
