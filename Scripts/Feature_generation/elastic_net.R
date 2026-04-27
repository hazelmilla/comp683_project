library(data.table)
library(tidyverse)
library(dplyr)
library(glmnet)

path = "/work/users/h/e/hemilla/Nanopore/original_data/"
setwd(path)

#1. Data prep ==========================
ont_m <- fread("ont_m_matrix_original.tsv", sep = "\t")
pheno <- fread("pheno_cell_type_original.csv", sep = ",")

# Format data for model
pheno <- pheno %>% subset(select = c("sample", "sex"))

identical(as.character(pheno$sample), colnames(ont_m)[2:8]) # verify order

ont_m_matrix <- ont_m %>% subset(select = c(colnames(ont_m)[1:8])) # select only CpG name and DNAm levels

pheno$sex <- ifelse(pheno$sex == "F", 0, 1) # code as numeric

#2. Run elastic net ==========================
# elastic net function
run_elastic_net <- function(input_beta, covariates, alpha = 0.5, n_columns = 7) {
  
  beta_matrix <- input_beta[, colSums(is.na(input_beta)) == 0]
  
  ages <- as.numeric(rownames(beta_matrix))
  if (any(is.na(ages))) stop("Non-numeric or NA values found in sample age row names.")
  stopifnot(length(ages) == nrow(beta_matrix))
  
  # Align covariates to sample order; exclude sample name column
  cov_cols <- setdiff(colnames(covariates), "sample")
  cov_matrix <- as.matrix(
    covariates[match(rownames(beta_matrix), covariates$sample), ..cov_cols]
  )
  
  beta_matrix <- matrix(as.numeric(beta_matrix), 
                        nrow = nrow(beta_matrix), 
                        dimnames = dimnames(beta_matrix))
  
  # Build combined matrix and penalty vector
  full_matrix <- cbind(cov_matrix, beta_matrix)
  penalty_vec <- c(rep(0, ncol(cov_matrix)), rep(1, ncol(beta_matrix)))
  
  cv_fit <- cv.glmnet(full_matrix, ages, alpha = alpha, penalty.factor = penalty_vec, nfolds = 7)
  plot(cv_fit)
  
  selected_coefs <- coef(cv_fit, s = "lambda.min")
  selected_cpgs <- rownames(selected_coefs)[which(selected_coefs != 0)]
  selected_cpgs <- setdiff(selected_cpgs, c("(Intercept)", cov_cols))
  
  return(selected_cpgs)
}

# Make sure samples are rows, CpGs are columns 
ont_m_matrix_noNA <- na.omit(ont_m_matrix)
ont_m_matrix_t <- t(ont_m_matrix_noNA[,-1]) # remove ont_cpg, transpose
dim(ont_m_matrix_t) #[1]        7 2643694

colnames(ont_m_matrix_t) <- ont_m_matrix_noNA$ont_cpg # name columns by cpg name

ont_m_selected_cpgs <- run_elastic_net(ont_m_matrix_t,
                                       covariates = pheno) # run elastic net function

# format and filter
ont_m_selected_cpgs <- data.frame(ont_m_selected_cpgs) 

ont_m_matrix_elasticnet <- ont_m_matrix %>%
  filter(ont_cpg %in% ont_m_selected_cpgs$ont_m_selected_cpgs)

# write output file
data.table::fwrite(ont_m_matrix_elasticnet,
                   file = "ont_m_selected_cpgs_elastic_net.csv", quote = FALSE, sep = ",")
