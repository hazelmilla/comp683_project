###### Run Nanomix for 7 Nanopore samples after preprocessing ######
# 1. Data prep ======================
# Check data in tsv file
# https://github.com/Jonbroad15/nanomix/blob/main/README.md 

head -1 ~/ont_m.tsv | tr '\t' ‘\n’ # colnames
head -2 ~/ont_m.tsv | tail -1 | tr '\t' '\n' # first row
cut -f1 ~/ont_m.tsv | tail -n +2 # no rownames

# We need a tsv file with the following columns:
# {chr, start, end, total_calls, modified_calls}
# We have chrom, start, end, …. Nvalid_cov, Nmod
# This adjustment has been built into the deconvolution loop below

# Make adjustments to atlas so it only captures immune cell type proportions
awk 'BEGIN{FS=OFS="\t"}
NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="chr" || $i=="start" || $i=="end" || $i=="B-cell" || $i=="NK-cell" || $i=="T-cell" || $i=="monocyte" || $i=="granulocyte") {
      keep[i]=1
    }
  }
}
{
  out=""
  for (i=1; i<=NF; i++) {
    if (i in keep) {
      out = out (out=="" ? "" : OFS) $i
    }
  }
  print out
}' ~/39Bisulfite.tsv > ~/atlas_immune_cell.tsv

# 2. Install nanomix ======================
# Create nanomix environment

conda create -n nanomix python=3.10 -y
conda activate ~/miniconda3/envs/nanomix

which python                            # should show miniconda3/envs/nanomix/...
python3 --version                        # should be 3.10.x

# this needs to happen before maturin
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source $HOME/.cargo/env

# install nanomix from the source using maturin
pip install maturin
git clone https://github.com/Jonbroad15/nanomix.git
cd nanomix

rustup target add x86_64-apple-darwin
maturin develop --target x86_64-apple-darwin

#cd nanomix
#pip install -e .

# set PYTHONPATH
export PYTHONPATH="/Users/hazelmilla/nanomix/python:$PYTHONPATH"

# copy .so with expected name
cp ~/nanomix/python/nanomix/nanomix.cpython-310-darwin.so \
   ~/nanomix/python/nanomix/_nanomix.cpython-310-darwin.so

# verify Python can find nanomix
python -c "import sys; print('\n'.join(sys.path))"

# remove package resource, move to correct environment
rm /Users/hazelmilla/miniconda3/envs/nanomix/lib/python3.10/site-packages/pkg_resources

cp -r /Users/hazelmilla/miniconda3/envs/nanomix/lib/python3.10/site-packages/pip/_vendor/pkg_resources \
      /Users/hazelmilla/miniconda3/envs/nanomix/lib/python3.10/site-packages/pkg_resources
/Users/hazelmilla/miniconda3/envs/nanomix/bin/python -c "import pkg_resources; print('found at:', pkg_resources.__file__)"

# 3. Run nanomix deconvolution ======================
for f in ~/*_filtered.tsv; do
    sample=$(basename "$f" _filtered.tsv)
    echo "Processing $sample..."

    # filter and rename columns
    nanomix_input="${f/_filtered.tsv/_nanomix.tsv}"
    awk -F'\t' '
    NR==1{for(i=1;i<=NF;i++) col[$i]=i; print "chr\tstart\tend\ttotal_calls\tmodified_calls"; next}
    {print $col["chrom"]"\t"$col["start"]"\t"$col["end"]"\t"$col["Nvalid_cov"]"\t"$col["Nmod"]}
    ' "$f" > "$nanomix_input"
    echo "  Filtered: $nanomix_input"

    # deconvolute
    deconv_output="${f/_filtered.tsv/_deconv.tsv}"
    /Users/hazelmilla/miniconda3/envs/nanomix/bin/nanomix deconvolute -a atlas_immune_cell.tsv "$nanomix_input" > "$deconv_output"
    echo "  Deconvoluted: $deconv_output"
done

######## Run Nanomix on imputed data #########
#4. format tsv files in python ================
python
import pandas as pd
import os

os.chdir("/Users/hazelmilla/Desktop/Nanopore_files/")
df = pd.read_csv("ont_m_imputed_wide.csv", sep=",")

names_list = df.columns.tolist()
sample_names = names_list[1:] # get list of sample names w/o ont_cpg col
print(sample_names) 

dfs = {}

# run loop to generate 1 file per sample with correct formatting
# must include columns: {chr, start, end, total_calls, modified_calls}
for sample in sample_names: 
    # Load one filtered file at a time
    file_name = str(sample) + '_filtered.tsv'
    df_with_chr = pd.read_csv(file_name, sep="\t")
    
    df_subset = df[['ont_cpg', sample]].dropna().copy()
    DEPTH = 10
    df_subset["total_calls"] = DEPTH
    df_subset["modified_calls"] = (df_subset[sample] * DEPTH).round().astype(int)
    
    df_merged = pd.merge(df_with_chr, df_subset, on='ont_cpg')
    
    df_out = df_merged[['chrom', 'start', 'end', 'total_calls', 'modified_calls']]
    df_out = df_out.rename(columns={'chrom': 'chr'}) # rename 'chr' column
    
    df_out.to_csv(str(sample) + '_imputed_nanomix.tsv', sep="\t", index=False)
    
    # Explicitly free memory before next iteration
    del df_with_chr, df_subset, df_merged, df_out


exit()

# 5. Get everything set up for deconvolution again ================
# See step 2: install nanomix
# activate environment
conda activate ~/miniconda3/envs/nanomix

# verify that Python can find nanomix environment
python -c "import sys; print('\n'.join(sys.path))"

# 6. Run deconvolution ================
for f in "$HOME"/Desktop/Nanopore_files/*_imputed_nanomix.tsv; do
    sample=$(basename "$f" _imputed_nanomix.tsv)
    echo "Processing $sample..."

    deconv_output="${f/_imputed_nanomix.tsv/_imputed_deconv.tsv}"
    /Users/hazelmilla/miniconda3/envs/nanomix/bin/nanomix deconvolute \
        -a atlas_immune_cell.tsv "$f" > "$deconv_output"
    echo "  Deconvoluted: $deconv_output"
done

# Troubleshooting
cd ~/Desktop/Nanopore_files

bedtools intersect \
    -a xx5_imputed_nanomix.tsv \
    -b atlas_immune_cell.tsv \
    -wa -u | wc -l

# 66506 / 19282076 sites present in atlas

########## Running again on original data ##########

cd ~/Desktop/Nanopore_files
# convert bed files to tsv
# column count
head -1 xx5.bed | awk '{print NF}' # 18

# see each column on its own line with index
head -1 xx5.bed | tr '\t' '\n' | nl

for f in *.bed; do
    sample=$(basename "$f" .bed)
    echo "Processing $sample..."
    
    awk 'BEGIN{OFS="\t"; print "chr\tstart\tend\ttotal_calls\tmodified_calls"}
         {
             total = $5
             modified = int($11/100 * total + 0.5)
             print $1, $2, $3, total, modified
         }' "$f" > "${sample}_nanomix.tsv"
    
    echo "  Written: ${sample}_nanomix.tsv"
done

# activate environment and make sure Python can find it
conda activate ~/miniconda3/envs/nanomix
python -c "import sys; print('\n'.join(sys.path))"

cd ~/Desktop/Nanopore_files

# run nanomix
for f in "$HOME"/Desktop/Nanopore_files/*_aggregated_nanomix.tsv; do
    sample=$(basename "$f" _aggregated_nanomix.tsv)
    echo "Processing $sample..."

    deconv_output="${f/_aggregated_nanomix.tsv/_deconv.tsv}"
    /Users/hazelmilla/miniconda3/envs/nanomix/bin/nanomix deconvolute \
        -a atlas_immune_cell.tsv "$f" > "$deconv_output"
    echo "  Deconvoluted: $deconv_output"
done

cd ~/Desktop/Nanopore_files
nanomix simulate -a atlas_immune_cell.tsv xx5_deconv.tsv > nano_sim.tsv
nanomix evaluate -a atlas_immune_cell.tsv nano_sim.tsv > eval_output.tsv
