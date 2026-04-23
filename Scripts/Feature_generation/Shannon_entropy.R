# Compute Shannon entropy

library(data.table)
library(tidyverse)
library(dplyr)

path = "/work/users/h/e/hemilla/Nanopore/original_data/"
setwd(path)

#1. Data prep ==========================
ont_m <- fread("ont_m.tsv", sep = "\t")
pheno <- fread("pheno_cell_type_original.csv")

pheno <- pheno %>% subset(select = c("sample", "sex")) # only use sex as covariate for now

# Pivot and sort age columns numerically (code adapted from Găitănaru et al. 2025)
pivot_sorted <- function(df) {
  df %>%
    select(age = 20, ont_cpg, fraction_modified) %>%
    pivot_wider(names_from = age, values_from = fraction_modified, values_fn = sum, values_fill = NA) %>%
    select(ont_cpg, all_of(as.character(sorted_ages)))  # Ensure column order
}

sorted_ages <- sort(unique(ont_m[[20]]))
age_colnames <- as.character(sorted_ages)
ont_m_matrix  <- pivot_sorted(ont_m)
age_values <- as.numeric(colnames(ont_m_matrix)[-1])

beta_matrix <- na.omit(ont_m_matrix) %>% as.data.frame()
rownames(beta_matrix) = beta_matrix$ont_cpg # name rows by cpg
beta_matrix$ont_cpg = NULL # remove cpg column
dim(beta_matrix) 

beta_matrix_t <- t(beta_matrix)
dim(beta_matrix_t) #[1]       7 2643694

# compute entropy using Shannon entropy equation - DNAm level as input
compute_entropy <- function(m) {
  eps <- 1e-6
  m_clamped <- pmax(pmin(m, 1 - eps), eps)   # clips to [1e-6, 0.999999]
  -m_clamped * log2(m_clamped) - (1 - m_clamped) * log2(1 - m_clamped)
} # the clamping is so the function can handle 0 and 1 values

# Apply across entire matrix (rows=ages, cols=CpGs)
entropy_matrix <- apply(beta_matrix, 2, compute_entropy) # apply across CpGs (cols)

# entropy_matrix: same dims as beta_matrix
dim(entropy_matrix) #7 2643694
identical(colnames(beta_matrix), as.character(age_values))

entropy_df <- entropy_matrix %>% as.data.frame()

# recreate ont_cpg column
entropy_df$ont_cpg <- rownames(entropy_df)
entropy_df <- relocate(entropy_df, ont_cpg) # move to first column

dim(entropy_df)        # 2643694       8
head(entropy_df)
glimpse(entropy_df)

fwrite(entropy_df, "Shannon_entropy_values.csv", sep = ",", row.names = FALSE) # save entropy values

#2. Run linear regression entropy ~ age ===============

library(data.table)
library(tidyverse)
library(dplyr)

path = "/work/users/h/e/hemilla/Nanopore/original_data/"
setwd(path)

#1. Data prep ==========================
pheno <- fread("pheno_cell_type_original.csv")
pheno <- pheno %>% subset(select = c("sample", "sex"))

entropy_matrix <- fread("Shannon_entropy_values.csv", sep = ",", header = TRUE)
dim(entropy_matrix) # [1] 2643694       8

age_values <- as.numeric(colnames(entropy_matrix)[-1])

# residualize age on sex to adjust for sex
res_x <- resid(lm(age_values ~ sex, data = pheno))

# Remove CpGs where entropy has near-zero variance across the 7 samples
entropy_var <- apply(entropy_matrix[-1,], 1, var)

# Inspect distribution
# summary(entropy_var)
# hist(log10(entropy_var + 1e-10), 
#      main = "Distribution of entropy variance across CpGs",
#      xlab = "log10(variance)")

# Separate the ID column before computation
ont_cpg_ids <- entropy_matrix[, "ont_cpg"]  # save IDs
numeric_matrix <- entropy_matrix[, -1, drop = FALSE]  # numeric only

# Compute variance on numeric columns only
entropy_var <- apply(numeric_matrix, 1, var)
min_var_threshold <- 1e-6 # threshold selected based on the histogram of entropy variance

# Filter for rows with high enough variance
keep <- entropy_var > min_var_threshold

entropy_matrix_filtered <- entropy_matrix[keep, ] # filter for high-variance CpGs

head(colnames(entropy_matrix_filtered)) # check colnames

cat("CpGs retained:", nrow(entropy_matrix_filtered), 
    "of", nrow(entropy_matrix), "\n")
# CpGs retained: 2023231 of 2643694 

# define linear regression function
get_lm <- function(row, cov_dat){
  y = as.numeric(row)
  x = age_values
  
  model <- lm(y ~ x + sex, data = cov_dat)
  
  coefs     <- summary(model)$coefficients
  res_y     <- resid(lm(y ~ sex, data = cov_dat))
  
  c(coef(model)[1],             # intercept
    coef(model)[2],             # slope
    summary(model)$r.squared,
    coefs["x", "Pr(>|t|)"],
    cor(res_y, res_x)) 
}

dim(entropy_matrix_filtered) # make sure rows are CpGs and cols are samples
# [1] 2023231       8

entropy_lm <- t(apply(
  entropy_matrix_filtered [, -1, drop = FALSE],
  1,
  function(x) {
    get_lm(
      row = x,
      cov_dat = pheno
    )
  }
))

# Add linear regression statistics to matrix
entropy_matrix_filtered$intercept <- entropy_lm[, 1]
entropy_matrix_filtered$slope <- entropy_lm[, 2]
entropy_matrix_filtered$r_squared <- entropy_lm[, 3]
entropy_matrix_filtered$p_value <- entropy_lm[, 4]
entropy_matrix_filtered$correlation <- entropy_lm[, 5]

# multiple test correction
entropy_matrix_filtered$p_fdr <- p.adjust(entropy_matrix_filtered$p_value, method = "BH")

# write file
fwrite(entropy_matrix_filtered, "entropy_matrix_stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# create matrix filtered only for significant values
# entropy_matrix_significant <- entropy_matrix_filtered %>%
#   filter(p_value < 0.05, r_squared >= 0.8)
# fwrite(entropy_matrix_significant, "entropy_matrix_significant.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
 

