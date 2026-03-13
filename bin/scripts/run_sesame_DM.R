library(sesame)
library(tidyverse)
library(SummarizedExperiment)
library(cowplot)
library(ggrastr)

start_time <- Sys.time()
message("Starting sesame DM analysis at ", format(start_time, "%Y-%m-%d %H:%M:%S"))

# ===== Configuration =====
PVAL_THRESHOLD <- snakemake@params[['pval_threshold']]
EFF_THRESHOLD <- snakemake@params[['eff_threshold']]

PLOT_HEIGHT <- 7
PLOT_WIDTH <- 7.5

# ===== Input Validation =====
group1 <- snakemake@wildcards[['group1']]
group2 <- snakemake@wildcards[['group2']]

if (group1 == group2) {
  stop("Error: group1 and group2 cannot be the same (", group1, ")")
}

# ===== Load and Prepare Data =====
se <- readRDS(snakemake@input[[1]])

if (!exists("se") || is.null(se)) {
  stop("Error: Failed to load SummarizedExperiment from ", snakemake@input[[1]])
}

available_groups <- unique(colData(se)$group)
if (!all(c(group1, group2) %in% available_groups)) {
  stop("Error: Groups '", group1, "' and/or '", group2, 
       "' not found in data. Available groups: ", paste(available_groups, collapse=", "))
}

se <- se[, se$group %in% c(group1, group2)]
message(paste0("Setting the reference to group ", group2))
colData(se)$group <- relevel(factor(colData(se)$group), group2)
message(paste0("Keeping only samples from the groups: ", paste(unique(colData(se)$group), collapse=", ")))
message(paste0("Sample count: ", ncol(se)))

# ===== Quality Control =====
se_ok <- (checkLevels(assay(se), colData(se)$group))
n_ok <- sum(se_ok)
n_total <- nrow(se)
message(paste0(n_ok, " out of ", n_total, " CpGs passed quality control in differential methylation (", 
               round(100 * n_ok / n_total, 1), "%)"))

if (n_ok < 100) {
  warning("Warning: Less than 100 CpGs passed checkLevels function from SeSAMe. Results may be unreliable.")
}

se <- se[se_ok, ]

# ===== Differential Methylation Analysis =====
message("Running DML (differential methylation locus analysis)...")
smry <- DML(se, ~group, BPPARAM = BiocParallel::MulticoreParam(snakemake@threads))
test_result <- summaryExtractTest(smry)

if (nrow(test_result) == 0) {
  stop("Error: No results from DML analysis. Check input data.")
}

message("Running DMR (differential methylation region analysis)...")
merged <- DMR(se, smry, paste0("group", group1), platform=snakemake@params[['array_type']])

# ===== Summary Statistics =====
## specify column names for the summary:
unadjusted_p_col <- paste0("Pval_group", group1)
eff_col <- paste0("Eff_group")
n_DMR_seg <- merged |>
  filter(!is.na(Seg_Pval_adj), Seg_Pval_adj < 0.05) |>
  summarise(n = dplyr::n_distinct(Seg_ID)) |>
  pull(n)

#--- Print some basic summary --- 
message("Analysis Summary for DMP:")
message("Note - p-values in DMP is not adjusted for mulitple testing.")

message(paste0("  - Total CpGs analyzed: ", nrow(se)))
message(paste0("  - Significant DMPs (p < ", PVAL_THRESHOLD, "): ", 
               sum(test_result[[unadjusted_p_col]] < PVAL_THRESHOLD, na.rm = T)))

message(paste0("  - DMPs with effect size > ", EFF_THRESHOLD, ": ", 
               sum(abs(test_result[[eff_col]]) > EFF_THRESHOLD, na.rm = T)))
message(paste0(
  "  - DMPs passing both p-value < ", PVAL_THRESHOLD, 
  " and |effect size| > ", EFF_THRESHOLD, ": ",
  sum(
    (test_result[[unadjusted_p_col]] < PVAL_THRESHOLD) &
      (abs(test_result[[eff_col]]) > EFF_THRESHOLD),
    na.rm = TRUE
  )
))

message("----------------------------------------------------")

message("Analysis Summary for DMR:")
message("Note - p-values in DMR is adjusted for mulitple testing.")
message(paste0("  - Significant DMRs identified: ", n_DMR_seg))

# ===== Save Results =====
message("Saving results to ", snakemake@output[[1]])
saveRDS(list(se = se, 
             dml_smry = smry, 
             dml_res = test_result, 
             dmr = merged), 
        snakemake@output[[1]])


# ===== Plotting Functions =====
plot_volcano_DMP_greater_than_2levels <- function(df, 
                                                  name, 
                                                  eff_thresh, 
                                                  pval_thresh, 
                                                  title_name) {
  ## Validate input
  if (nrow(df) == 0) {
    stop("Error: Empty dataframe provided to plot_volcano")
  }
  
  stopifnot(all(str_detect(colnames(df), "Pval|Est")))
  colnames(df) <- str_split_fixed(colnames(df), "_", 2)[, 1]
  
  ## Categorize by significance
  df <- df %>% 
    mutate(
      pval_sig = Pval < pval_thresh,
      eff_sig  = abs(Est) > eff_thresh,
      color = case_when(
        pval_sig & eff_sig  ~ "Both significant",
        pval_sig & !eff_sig ~ "P-value only",
        !pval_sig & eff_sig ~ "Effect size only",
        TRUE                ~ "Not significant"
      ),
      color = factor(color, levels = c(
        "Both significant", 
        "P-value only", 
        "Effect size only",
        "Not significant"
      ))
    )
  
  ## Count categories for legend
  color_counts <- df %>% 
    group_by(color) %>% 
    summarize(n = n(), .groups = 'drop')
  
  message("Volcano plot summary for ", name, ":")
  print(color_counts)
  
  ## Define consistent color mapping
  color_map <- c(
    "Both significant" = "#D53E4F",
    "P-value only" = "#F1A340",
    "Effect size only" = "#3288BD",
    "Not significant" = "grey"
  )
  # Only use colors for categories that exist in data
  colors_used <- color_map[intersect(names(color_map), unique(df$color))]
  
  volcano <- ggplot(df) +
    #geom_point(aes(x=Est, y = -log10(Pval), color = color), size = 0.75, alpha = 0.7) +
    geom_point_rast(
      aes(x=Est, y = -log10(Pval), color = color), 
      size = 0.75, 
      alpha = 0.7) +
    xlab(paste0("Beta value difference of: ", name)) +
    ylab("-log10(p-value)") +
    geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "grey50", size = 0.5) +
    geom_vline(xintercept = c(-eff_thresh, eff_thresh), linetype = "dashed", color = "grey50", size = 0.5) +
    scale_color_manual(values = colors_used, name = "") +
    theme_classic(base_size = 13) +
    theme(legend.position = "right") +
    ggtitle(paste0("Volcano plot for ", title_name))
  
  return(volcano)
}

# ===== Generate Volcano Plot =====
message("\nGenerating volcano plot...")

# Create output directory
dir.create("analysis/plots/", showWarnings = FALSE)

group_comparison <- paste0(group1, " vs ", group2)
plot_data <- test_result %>% 
  dplyr::select(matches("^Est_group|^Pval_group"))

if (nrow(plot_data) > 0) {
  volcano <- plot_volcano_DMP_greater_than_2levels(
    df = plot_data,
    title_name = group_comparison,
    name = group_comparison,
    eff_thresh = EFF_THRESHOLD,
    pval_thresh = PVAL_THRESHOLD
  )
  
  output_file <- paste0("analysis/plots/volcano_", group1, "-", group2, ".pdf")
  message("Saving volcano plot to ", output_file)
  save_plot(output_file, volcano, base_height = PLOT_HEIGHT, base_width = PLOT_WIDTH)
} else {
  warning("No data available for volcano plot")
}

end_time <- Sys.time()
message("Finished sesame DM analysis at ", format(end_time, "%Y-%m-%d %H:%M:%S"))
message("Total runtime: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes")