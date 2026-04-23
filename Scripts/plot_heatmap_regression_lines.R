library(data.table)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(grid)
library(ggplot2)
library(reshape2)

path = "/work/users/h/e/hemilla/Nanopore/original_data/"
setwd(path)

lm_res <- fread("ont_m_matrix_original_nominal.tsv", sep = "\t")
elastic <- fread("ont_m_selected_cpgs_elastic_net.csv", sep = ",", header = TRUE)
drift_with_DNAm <- fread("drift_with_DNAm.csv", sep = ",", header = TRUE)
entropy <- fread("entropy_matrix_significant.tsv", sep = "\t")

plot_methylation_heatmap <- function(ont_m_matrix, input_name = "data") {

  # Sort by correlation (descending order)
  # sorted_data <- ont_m_matrix_filtered[order(ont_m_matrix_filtered$slope, decreasing = TRUE), ]
  # filter for NA values and by correlation desc as well
  
  # for drift: sorted_data <- ont_m_matrix[order(ont_m_matrix$drift_std, decreasing = TRUE), ]
  # for linear regression and entropy:
  sorted_data <- ont_m_matrix[order(ont_m_matrix$slope, decreasing = TRUE), ]
  sorted_data <- sorted_data[complete.cases(sorted_data[, c("5", "13", "24", "46", "69", "70", "91")]), ]
  
  # Extract methylation columns
  meth_matrix <- sorted_data[, c("5", "13", "24", "46", "69", "70", "91")]
  
  # Set row names to CpG identifiers (ont_cpg column)
  rownames(meth_matrix) <- sorted_data$ont_cpg
  meth_matrix <- as.matrix(meth_matrix)
  meth_matrix <- apply(meth_matrix, 2, as.numeric)
  
  # Define output filename
  file_name <- paste0("/work/users/h/e/hemilla/Nanopore/original_data/methylation_heatmap_", input_name, "_sorted_by_correlation.png")
  
  # Save to PNG with correct rendering
  png(file_name, width = 1000, height = 1200, res = 150)
  
  # Create and print the heatmap object to the PNG file
  ph <- pheatmap(meth_matrix,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 display_numbers = FALSE,
                 na_col = "white",
                 color = colorRampPalette(c("blue", "red"))(50),
                 main = paste("Methylation Heatmap -", input_name),
                 show_rownames = TRUE,
                 #      fontsize_row = 6,  # Reduce font size if there are many rows
                 #      fontsize_col = 10,  # Optional: Adjust column label size
                 #      cellwidth = 1,    # Optional: Adjust cell width for better fit
                 #      cellheight = 1,   # Optional: Adjust cell height for better fit
                 silent = TRUE)
  
  # Print the heatmap into the PNG device
  print(ph)
  
  # Turn off the device to save the image
  dev.off()
  
  message("Heatmap saved to ", file_name)
}

plot_methylation_heatmap(ont_m_matrix = lm_res, input_name = "linear_regression") # plots DNAm sorted by slope
plot_methylation_heatmap(ont_m_matrix = entropy, input_name = "entropy") # this plots entropy sorted by entropy association with age

drift_significant <- filter(drift_with_DNAm, whiteP < 0.05) # filter for significant drift sites
plot_methylation_heatmap(ont_m_matrix = drift_significant, input_name = "drift") # plot drift sites

plot_methylation_heatmap(ont_m_matrix = elastic, input_name = "elastic_net") # plot elastic net sites

# Regression line with points
plot_methylation_regression_lines_with_points <- function(ont_m_matrix_filtered, input_name = "data",
                                                          min_slope = 0, min_r_squared = -Inf,
                                                          min_p_value = Inf) {

  # Slope filter: ascending or descending only. We separate CpGs best direct correlated with ageing from the ones inverse correlated.
  if (min_slope >= 0) {
    slope_filter <- ont_m_matrix_filtered$slope >= min_slope
  } else {
    slope_filter <- ont_m_matrix_filtered$slope <= min_slope
  }
  
  # Apply all filters
  filtered_data <- subset(ont_m_matrix_filtered,
                          slope_filter &
                            r_squared >= min_r_squared &
                            p_value <= min_p_value)
  
  # Define age positions
  positions <- c(5, 13, 24, 46, 69, 70, 91)
  
  # Get top CpGs sorted by absolute slope
  top_cpgs_df <- head(filtered_data[order(-abs(filtered_data$slope)), c("ont_cpg", "slope")], 20)
  top_cpgs <- top_cpgs_df$ont_cpg
  top_cpgs_levels <- top_cpgs[order(-abs(top_cpgs_df$slope))]  # ensure legend is ordered
  
  # Build predicted values
  predicted_list <- lapply(1:nrow(filtered_data), function(i) {
    intercept <- filtered_data$intercept[i]
    slope <- filtered_data$slope[i]
    cpg <- filtered_data$ont_cpg[i]
    legend_label <- if (cpg %in% top_cpgs) cpg else NA
    
    data.frame(
      Position = positions,
      Methylation = intercept + slope * positions,
      ont_cpg = cpg,
      legend_label = legend_label
    )
  })
  predicted_df <- do.call(rbind, predicted_list)
  
  # Factor for legend ordering
  predicted_df$legend_label <- factor(predicted_df$legend_label, levels = top_cpgs_levels)
  
  # Get actual methylation data and melt it
  actual_meth <- filtered_data[, c("5", "13", "24", "46", "69", "70", "91")]
  actual_meth$ont_cpg <- filtered_data$ont_cpg
  actual_meth$legend_label <- ifelse(actual_meth$ont_cpg %in% top_cpgs, actual_meth$ont_cpg, NA)
  actual_meth$legend_label <- factor(actual_meth$legend_label, levels = top_cpgs_levels)
  
  actual_long <- reshape2::melt(actual_meth, id.vars = c("ont_cpg", "legend_label"), variable.name = "Position", value.name = "Methylation")
  actual_long$Position <- as.numeric(gsub("X", "", actual_long$Position))
  
  # Filter out non-top CpGs entirely
  predicted_df_top <- subset(predicted_df, !is.na(legend_label))
  actual_long_top <- subset(actual_long, !is.na(legend_label))
  
  # Output filename
  file_name <- paste0("/work/users/h/e/hemilla/Nanopore/original_data/methylation_regression_lines_with_points_", input_name, ".png")
  
  # Save plot
  png(file_name, width = 2000, height = 800, res = 150)
  
  p <- ggplot() +
    geom_line(data = predicted_df_top, aes(x = Position, y = Methylation, group = ont_cpg, color = legend_label), size = 0.3) +
    geom_point(data = actual_long_top, aes(x = Position, y = Methylation, group = ont_cpg, color = legend_label), size = 0.5) +
    theme_minimal() +
    ggtitle(paste("Methylation Regression Lines with Points - ", input_name, " (Top CpGs)")) +
    xlab("Age") +
    ylab("Methylation Level") +
    scale_x_continuous(breaks = positions) +
    scale_color_discrete(name = "Top CpGs (sorted by slope)") +
    guides(color = guide_legend(override.aes = list(linewidth = 3), ncol = 1)) +
    theme(legend.key.size = unit(0.6, "lines"))
  
  print(p)
  dev.off()
  
  message("Plot with regression lines and actual points saved to ", file_name)
}

lm_res_noNA <- lm_res[complete.cases(lm_res), ]
plot_methylation_regression_lines_with_points(lm_res_noNA, input_name = "Linear Regression",
                       min_slope = 0.005,
                       min_r_squared = 0.9,
                       min_p_value = 0.01)

entropy_noNA <- na.omit(entropy)
plot_methylation_regression_lines_with_points(entropy_noNA, input_name = "Entropy",
                                              min_slope = 0.005,
                                              min_r_squared = 0.9,
                                              min_p_value = 0.01)

# create Venn diagram to determine overlapping CpGs across different methods
library(eulerr)

cpgs <- list(Drift = c(drift_significant$ont_cpg), Entropy = c(entropy_noNA$ont_cpg),
             Elastic_net = c(elastic$ont_cpg), Linear = c(lm_res_noNA$ont_cpg))
plot(venn(cpgs))

# Compare with original results
og_en <- fread("ont_m_matrix_elasticnet.tsv", sep = "\t", header = TRUE) # elastic net results
og_lm <- fread("ont_m_matrix_filtered.tsv", sep = "\t", header = TRUE) # linear regression results

cpgs1 <- list(Găitănaru = c(og_lm$ont_cpg), Ours = c(lm_res_noNA$ont_cpg))
cpgs2 <- list(Găitănaru = c(og_en$ont_cpg), Ours = c(elastic$ont_cpg))

plot(venn(cpgs1))
plot(venn(cpgs2))

