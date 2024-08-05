library(sesame)
library(tidyverse)
library(SummarizedExperiment)
library(cowplot)

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


### plot volcano
plot_volcano_DMP_greater_than_2levels <- function(df = df, name = name, eff_thresh = eff_thresh, title_name=title_name) {
  ## below is a functional check
  stopifnot(all(str_detect(colnames(df), "Pval|Est")))
  colnames(df) <- str_split_fixed(colnames(df), "_", 2)[, 1]
  df <- df %>% mutate(pval_thresh = case_when(Pval < 0.05 ~ "pval_sig",
                                              TRUE  ~ "pval_not"))
  df <- df %>% mutate(eff_thresh = case_when(abs(Est) > eff_thresh ~ "eff_sig",
                                             TRUE  ~ "eff_not"))
  df <- df %>% mutate(color = case_when(pval_thresh == "pval_sig" & eff_thresh == "eff_sig" ~"purple",
                                        pval_thresh =="pval_sig" & eff_thresh =="eff_not" ~ "goldenrod",
                                        pval_thresh =="pval_not" & eff_thresh =="eff_sig" ~ "brown3",
                                        TRUE ~ "grey"))
  volcano <- ggplot(df) +
    geom_point(aes(x=Est, y = -log10(Pval), color = color), size = 0.75) +
    xlab(paste0("Beta value difference of: ", name)) +
    ylab("-log10(p-value)") +
    scale_color_manual(values = c("brown3", "goldenrod", "grey", "cornflowerblue"),
                       labels = c("beta diff", "pvalue", "NS", "pvalue & beta diff"))  +
    theme_classic(base_size = 13) +
    theme(legend.title = element_blank()) +
    ggtitle(paste0("Volcano plot for ", name))
  return(volcano)
}


volcano <- plot_volcano_DMP_greater_than_2levels(df = test_result %>% dplyr::select(matches("^Est_group|^Pval_group")),
                                                 title_name = paste0(snakemake@wildcards[['group1']],"-", snakemake@wildcards[['group2']]),
                                                 name = paste0(snakemake@wildcards[['group1']],"-", snakemake@wildcards[['group2']]),
                                                 eff_thresh = 0.05)

dir.create("analysis/plots/")
save_plot(paste0("analysis/plots/volcano_", snakemake@wildcards[['group1']], "-", snakemake@wildcards[['group2']],".jpeg"),
          volcano, base_height = 7, base_width = 7.5)
