# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(stringr)
  library(rstatix) 
  library(ggpubr)  
  library(pheatmap)
})

#-----------------------------------------------------
option_list <- list(
  make_option(c("--deconv_method"), type = "character",
              help = "Deconvolution method (xcell or quantiseq)"),
  make_option(c("--clin_file"), type = "character",
              help = "histologies file (.tsv)"),
  make_option(c("--results_dir"), type = "character",
              help = "Directory containing deconvolution results"),
  make_option(c("--output_dir"), type = "character", 
              help = "Output directory")
)


opt <- parse_args(OptionParser(option_list = option_list))

# Validate input
if (!(opt$deconv_method %in% c("xcell", "quantiseq")))
  stop("Error: deconv_method must be one of xcell or quantiseq")

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

plots_dir <- file.path(opt$output_dir, "plots")
results_dir <- file.path(opt$output_dir)

# Create directories if they don't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


#-----------------------------------------------------

#expr_mat <- readRDS("data/gene-expression-rsem-tpm-collapsed.rds")
#histologies <- read_tsv("data/histologies.tsv")
#input
histologies <- read_tsv(opt$clin_file, guess_max = 10000)
deconv_results <- readRDS(file.path(opt$results_dir, paste0(opt$deconv_method, "_output.rds")))

#-----------------------------------------------------
#functions
prepare_mb_data <- function(deconv_results, method, keep_uncharacterized = TRUE) {
  mb_samples <- histologies %>%
    filter(cancer_group == "Medulloblastoma") %>%
    select(Kids_First_Biospecimen_ID, molecular_subtype, cancer_group)
  
  mb_data <- deconv_results %>%
    inner_join(mb_samples, 
               by = c("Kids_First_Biospecimen_ID", "cancer_group", "molecular_subtype")) %>%
    mutate(molecular_subtype = str_replace(molecular_subtype, "MB, ", ""),
           molecular_subtype = factor(molecular_subtype, 
                                      levels = c("WNT", "SHH", "Group3", "Group4")))
  
  if(method == "quantiseq" && !keep_uncharacterized) {
    mb_data <- mb_data %>% 
      filter(cell_type != "uncharacterized cell") %>%
      group_by(Kids_First_Biospecimen_ID) %>%
      mutate(fraction = fraction/sum(fraction)) %>%  # Rescale remaining fractions to sum to 1
      ungroup()
  }
  
  return(mb_data)
}

perform_statistical_tests <- function(mb_data) {
  # Kruskal-Wallis test
  kw_results <- mb_data %>%
    group_by(cell_type) %>%
    summarise(
      kw_pval = kruskal.test(fraction ~ molecular_subtype)$p.value
    ) %>%
    mutate(
      padj = p.adjust(kw_pval, method = "BH"),
      significant = padj < 0.05
    )
  
  # Dunn test for significant cell types
  dunn_results <- mb_data %>%
    group_by(cell_type) %>%
    filter(cell_type %in% (kw_results %>% filter(significant) %>% pull(cell_type))) %>%
    dunn_test(fraction ~ molecular_subtype, p.adjust.method = "BH") %>%
    mutate(
      p.adj = formatC(p.adj, format = "e", digits = 2),
      significance = case_when(
        as.numeric(p.adj) < 0.001 ~ "***",
        as.numeric(p.adj) < 0.01 ~ "**",
        as.numeric(p.adj) < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      comparison = paste(group1, "vs", group2)
    ) %>%
    select(
      `Cell Type` = cell_type,
      `Comparison` = comparison,
      `Statistic` = statistic,
      `Adjusted P-value` = p.adj,
      `Significance` = significance
    )
  
  return(list(kw = kw_results, dunn = dunn_results))
}

create_violin_plot <- function(mb_data, method, stat_results) {
  plot_title <- if (method == "xcell") {
    "xCell Enrichment Scores Across MB Subtypes"
  } else {
    "Immune Cell Composition Across MB Subtypes"
  }
  
  y_label <- if (method == "xcell") {
    "xCell Enrichment Score"
  } else {
    "Absolute Cell Fraction"
  }
  
  significance_labels <- stat_results$kw %>%
    filter(significant) %>%
    mutate(label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    )) %>%
    inner_join(
      mb_data %>%
        group_by(cell_type) %>%
        summarise(max_fraction = max(fraction, na.rm = TRUE)),
      by = "cell_type"
    ) %>%
    mutate(y_position = max_fraction * 1.05)
  
  plot <- ggplot(mb_data, aes(x = molecular_subtype, y = fraction, fill = molecular_subtype)) +
    geom_violin(alpha = 0.7, scale = "width", width = 0.8) +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.7, outlier.size = 0.5) +
    facet_wrap(~cell_type, scales = "free", ncol = 4) +
    geom_text(
      data = significance_labels,
      aes(x = 2.5, y = y_position, label = label),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      panel.spacing = unit(0.5, "lines"),
      strip.background = element_rect(fill = "gray95", color = "black")
    ) +
    labs(
      title = plot_title,
      x = "Molecular Subtype",
      y = y_label
    ) +
    scale_fill_manual(values = c(
      "WNT" = "#E69F00",
      "SHH" = "#56B4E9",
      "Group3" = "#009E73",
      "Group4" = "#CC79A7"
    ))
  
  return(plot)
}



#-----------------------------------------------------


#Execution
method <- opt$deconv_method
results_file <- file.path("analyses/immune-deconv/results", paste0(method, "_output.rds"))
deconv_results <- readRDS(results_file)

mb_data <- prepare_mb_data(deconv_results, method)

stat_results <- perform_statistical_tests(mb_data)

plot <- create_violin_plot(mb_data, method, stat_results = stat_results)
ggsave(
  file.path(plots_dir, paste0(method, "_mb_violin_plot.pdf")),
  plot,
  width = 15,
  height = 20,
  limitsize = FALSE
)


write_tsv(
  stat_results$kw,
  file.path(results_dir, paste0(method, "_mb_kw_results.tsv"))
)

write_tsv(
  stat_results$dunn,
  file.path(results_dir, paste0(method, "_mb_dunn_results.tsv"))
)
