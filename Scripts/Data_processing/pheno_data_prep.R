library(data.table)
library(tidyverse)
library(dplyr)
library(readxl)

path = "~/Desktop/Nanopore_files"
setwd(path)

pheno_data <- read_excel("Supplementary_Material_S1_v4.xlsx", sheet = 9) # read in sex data

file_names <- list.files(path)
file_names_read <- file_names[grep("\\_imputed_deconv.tsv$", file_names, ignore.case = T)] # get file names

ctp_dfs <- lapply(file_names_read, function(x){fread(x, sep = "\t")}) # read in files
sample_names <- sapply(c(file_names_read), function(x){substr(x, 3, 4)}) # get sample names from files

sample_names[which(sample_names=="5_")] = "5" # rename sample to remove extra "_"

names(ctp_dfs) <- sample_names # name cell type proportion dataframes
ctp_df <- bind_rows(ctp_dfs, .id = "sample_names") # create dataframe with all sample cell type proportions

ctp_df <- dcast(
  ctp_df,
  sample_names ~ cell_type,
  value.var = "proportion",
  fill = 0
) # generate sample x cell type dataframe

pheno <- data.frame(sample_names = pheno_data$...2[-1], Sex = pheno_data$...4[-1]) # get sample name and sex data
pheno$sample_names = substr(pheno$sample_names, 2, 3) # edit sample name format

pheno <- left_join(pheno, as.data.frame(ctp_df), by = "sample_names") # merge with cell type proportions
fwrite(pheno, "pheno_cell_type_imputed.csv", row.names = FALSE) # write file with sex and cell type proportion
