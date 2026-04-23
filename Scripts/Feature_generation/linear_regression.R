################# Generate features for analysis ####################
#install.packages("data.table")
#install.packages("tidyverse")
#install.packages("dplyr")
#install.packages("readxl")
library(data.table)
library(tidyverse)
library(dplyr)

path = "/work/users/h/e/hemilla/Nanopore/original_data/"
setwd(path)

#1. Data prep ==========================
ont_m <- fread("ont_m.tsv", sep = "\t")
pheno <- fread("pheno_cell_type_original.csv", sep = ",")

pheno <- pheno %>% subset(select = c("sample", "sex")) # only adjusting for sex for now

# 1.Replicate linear regression and elastic net =========================
# Code adapted from Găitănaru et al. 2025

# Pivot and sort age columns numerically
pivot_sorted <- function(df) {
  df %>%
    select(age = 20, ont_cpg, fraction_modified) %>%
    pivot_wider(names_from = age, values_from = fraction_modified, values_fn = sum, values_fill = NA) %>%
    select(ont_cpg, all_of(as.character(sorted_ages)))  # Ensure column order
}

# Get sorted list of all unique ages (column 20)
sorted_ages <- sort(unique(ont_m[[20]]))
age_colnames <- as.character(sorted_ages)
ont_m_matrix  <- pivot_sorted(ont_m)
age_values <- as.numeric(colnames(ont_m_matrix)[-1])

# Check coverage
# How many CpGs are present at every age?
n_ages_total <- n_distinct(ont_m$age)
ont_m %>%
  group_by(ont_cpg) %>%
  summarise(n_ages = n_distinct(age)) %>%
  filter(n_ages == n_ages_total) %>%
  nrow()
# 2643694

identical(as.character(pheno$sample), colnames(ont_m_matrix)[-1]) # verify order

glimpse(ont_m_matrix)
colnames(ont_m_matrix)
head(ont_m_matrix$ont_cpg)

#2. Linear regression ==========================
# Function to compute intercept and slope per row 
fit_row_lm <- function(row, pheno_data) {
  y    <- as.numeric(row)
  mask <- !is.na(y)
  
  if (sum(mask) < 6) return(c(NA, NA, NA, NA, NA))
  
  y        <- y[mask]
  x        <- age_values[mask]
  res_x    <- res_x_full[mask]        # precomputed, just subset
  pheno_sub <- pheno_data[mask, ]
  
  model     <- lm(y ~ x + sex, data = pheno_sub)
  coefs     <- summary(model)$coefficients
  res_y     <- resid(lm(y ~ sex, data = pheno_sub))
  
  c(coef(model)[1],             # intercept
    coef(model)[2],             # slope
    summary(model)$r.squared,
    coefs["x", "Pr(>|t|)"],
    cor(res_y, res_x))          # partial r
}


res_x_full <- resid(lm(age_values ~ sex, data = pheno))

ont_m_matrix1 <- ont_m_matrix[-1] # exclude ont_cpg column

# Apply the function row-wise (excluding ont_cpg column)
regression_stats <- t(apply(
  ont_m_matrix1,
  1,
  function(x) {
    fit_row_lm(
      row = x,
      pheno_data = pheno
    )
  }
))

# add regression stats to matrix
ont_m_matrix$intercept <- regression_stats[, 1]
ont_m_matrix$slope <- regression_stats[, 2]
ont_m_matrix$r_squared <- regression_stats[, 3]
ont_m_matrix$p_value <- regression_stats[, 4]
ont_m_matrix$correlation <- regression_stats[, 5]

# write file
fwrite(ont_m_matrix, "ont_m_matrix_original.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#26875114

# filter for significant p-values and write file
ont_m_matrix_filtered <- ont_m_matrix %>%
  filter(p_value < 0.05, r_squared >= 0.8)
fwrite(ont_m_matrix_filtered, "ont_m_matrix_original_nominal.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# filter for FDR significant p-values and write file
ont_m_matrix$p_fdr <- p.adjust(ont_m_matrix$p_value, method = "BH")
ont_m_matrix_fdr <- ont_m_matrix %>%
  filter(p_fdr < 0.05, r_squared >= 0.8)
fwrite(ont_m_matrix_fdr, "ont_m_matrix_original_fdr.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# filter to remove sites with NA values and write file
ont_m_matrix_filtered_noNA <- ont_m_matrix_filtered[complete.cases(ont_m_matrix_filtered), ]
fwrite(ont_m_matrix_filtered_noNA, "ont_m_matrix_original_noNA.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

