# Compute & assess epigenetic drift ===================================
# code adapted from Fan et al. 2025

library("lmtest") # code for drift calculation
library("sandwich") # code for drift calculation
library(dplyr)
library(tidyr)
library(data.table)

path = "/work/users/h/e/hemilla/Nanopore/original_data/"
setwd(path)

pheno <- fread("pheno_cell_type_original.csv", sep = ",")
ont_m <- fread("ont_m.tsv", sep = "\t")

pheno <- pheno %>% subset(select = c("sample", "sex")) # only controlling for sex for now

# Pivot code from Găitănaru et al. 2025
pivot_sorted <- function(df) {
  df %>%
    select(age = 20, ont_cpg, fraction_modified) %>%
    pivot_wider(names_from = age, values_from = fraction_modified, values_fn = sum, values_fill = NA) %>%
    select(ont_cpg, all_of(as.character(sorted_ages)))  # Ensure column order
}

sorted_ages <- sort(unique(ont_m[[20]]))
age_colnames <- as.character(sorted_ages)
ont_m_matrix  <- pivot_sorted(ont_m) # apply pivot function
age_values <- as.numeric(colnames(ont_m_matrix)[-1]) # get age values

identical(as.character(pheno$sample), colnames(ont_m_matrix)[-1]) # verify that names match up

ont_m_matrix[1:5,1:7] # check matrix

# if you are getting ont_m_matrix with linear regression data, run:
ont_m_matrix <- ont_m_matrix %>% subset(select = c("ont_cpg", as.character(age_values)))

cpg_dat <- ont_m_matrix[complete.cases(ont_m_matrix),] %>% as.data.frame()
dim(cpg_dat)
#[1] 2643694       8

rownames(cpg_dat) <- cpg_dat$ont_cpg # name rownames by CpG
cpg_dat$ont_cpg = NULL # remove ont_cpg column

# Prepare data for drift function
cpg_dat <- cpg_dat %>% as.matrix() %>% t() 
head(colnames(cpg_dat))
rownames(cpg_dat)
dim(cpg_dat)
#[1]       7 x 2643694

cgdat <- data.frame(cpg_dat)[,1:5]
dim(cgdat) # Run on cgdat to test subset

pheno$sex <- ifelse(pheno$sex == "F", 0, 1) # code as numeric

###### White method (robust) ########
# code adapted from Fan et al. 2025
getWhitePdriftCpG <- function(method, alpha, nSamples, nIterations, age, cpgdata, covdata) {
  
  whiteP <- betaX <- seX <- tX <- pX <- 
    betaX2 <- seX2 <- tX2 <- pX2 <- 
    betaInter <- seInter <- tInter <- pInter <-
    drift_abs <- drift_fold <- drift_std <- rep(NA, nIterations)
  
  age_young <- quantile(age, 0.1)
  age_old   <- quantile(age, 0.9)
  
  set.seed(111)
  for(i in 1:nIterations){
    x <- age
    y <- cpgdata[, i]
    fo <- as.formula(paste0("y~x+", paste0(colnames(covdata)[-1], collapse="+")))
    model1 <- lm(fo, data = covdata)
    
    if (method == "white") {
      model2 <- lm((model1$residuals)^2 ~ x + I(x^2))
      
      # White test statistic
      r2       <- summary(model2)$r.squared
      stat     <- nSamples * r2
      whiteP[i] <- pchisq(q = stat, df = 2, lower.tail = FALSE)
      names(whiteP)[i] <- names(cpgdata)[i]
      
      # Robust coefficients (HC0)
      robustcoef <- coeftest(model2, vcov = vcovHC(model2, type = "HC0"))
      
      betaInter[i] <- robustcoef['(Intercept)', 'Estimate']
      seInter[i]   <- robustcoef['(Intercept)', 'Std. Error']
      tInter[i]    <- robustcoef['(Intercept)', 't value']
      pInter[i]    <- robustcoef['(Intercept)', 'Pr(>|t|)']
      
      betaX[i] <- robustcoef['x', 'Estimate']
      seX[i]   <- robustcoef['x', 'Std. Error']
      tX[i]    <- robustcoef['x', 't value']
      pX[i]    <- robustcoef['x', 'Pr(>|t|)']
      
      betaX2[i] <- robustcoef['I(x^2)', 'Estimate']
      seX2[i]   <- robustcoef['I(x^2)', 'Std. Error']
      tX2[i]    <- robustcoef['I(x^2)', 't value']
      pX2[i]    <- robustcoef['I(x^2)', 'Pr(>|t|)']
      
      # --- Degree of drift at this CpG site ---
      # Predicted variance = intercept + β₁·age + β₂·age²
      var_young    <- betaInter[i] + betaX[i]*age_young + betaX2[i]*age_young^2
      var_old      <- betaInter[i] + betaX[i]*age_old   + betaX2[i]*age_old^2
      
      drift_abs[i]  <- var_old - var_young               # absolute variance gained
      drift_fold[i] <- var_old / var_young               # fold change in variance
      drift_std[i]  <- (var_old - var_young) / var(model1$residuals)  # fraction of total variance
    }
  }
  
  power <- length(whiteP[whiteP < alpha]) / length(whiteP)
  
  return(list(
    power    = power,
    whiteP   = whiteP,
    # Intercept (baseline variance)
    betaInter = betaInter, seInter = seInter, tInter = tInter, pInter = pInter,
    # Linear age effect on variance
    betaX    = betaX,  seX  = seX,  tX  = tX,  pX  = pX,
    # Quadratic age effect on variance
    betaX2   = betaX2, seX2 = seX2, tX2 = tX2, pX2 = pX2,
    # Drift magnitude
    drift_abs  = drift_abs,    # absolute variance change (10th → 90th percentile age)
    drift_fold = drift_fold,   # fold change in variance
    drift_std  = drift_std     # fraction of total CpG variance explained by age
  ))
}

driftCpG_coef <- getWhitePdriftCpG(method = "white", 
                                   alpha = 0.05, 
                                   nSamples = 7, 
                                   nIterations = ncol(cpg_dat), 
                                   age = as.numeric(rownames(cpg_dat)), 
                                   cpgdata = cpg_dat,
                                   covdata = pheno) # apply drift calculation

# Format for export
driftCpG_coef <- as.data.frame(bind_rows(driftCpG_coef)) 
driftCpG_coef$ont_cpg <- colnames(cpg_dat)
driftCpG_coef <- driftCpG_coef %>% relocate("ont_cpg")

# apply multiple test correction
driftCpG_coef$p_fdr <- p.adjust(driftCpG_coef$whiteP, method = "BH")
fwrite(driftCpG_coef, "driftCpG_robust_coef.csv", row.names = FALSE) # write file

driftCpG_significant <- filter(driftCpG_coef, whiteP < 0.05) # filter for significant drift sites
fwrite(driftCpG_significant, "driftCpG_robust_significant.csv", row.names = FALSE) # write file
