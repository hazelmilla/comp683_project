{\rtf1\ansi\ansicpg1252\cocoartf2867
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 ArialMT;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Run Nanomix to get cell type proportions\
# https://github.com/Jonbroad15/nanomix/blob/main/README.md\
\
# 1. Set up\
conda deactivate # deactivate any active environment\
\
# verify base is clean\
echo "Current env: $CONDA_DEFAULT_ENV"  # should be empty or base \
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 # Activate nanomix env\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 conda activate nanomix \
\
# Verify correct environment and Python version\
echo "Active env: $CONDA_DEFAULT_ENV"   # should be nanomix\
which python                             # should show miniconda3/envs/nanomix/...\
python --version                         # should be 3.10.x\
\
# Set PYTHONPATH\
export PYTHONPATH="/Users/hazelmilla/nanomix/python:$PYTHONPATH"\
\
# Copy .so with expected name\
cp ~/nanomix/python/nanomix/nanomix.cpython-310-darwin.so \\\
   ~/nanomix/python/nanomix/_nanomix.cpython-310-darwin.so\
\
# Verify Python can find nanomix\
python -c "import sys; print('\\n'.join(sys.path))"\
\
# needed to uninstall and reinstall pkg_resources in the correct location\
rm /Users/hazelmilla/miniconda3/envs/nanomix/lib/python3.10/site-packages/pkg_resources\
\
cp -r /Users/hazelmilla/miniconda3/envs/nanomix/lib/python3.10/site-packages/pip/_vendor/pkg_resources \\\
      /Users/hazelmilla/miniconda3/envs/nanomix/lib/python3.10/site-packages/pkg_resources\
/Users/hazelmilla/miniconda3/envs/nanomix/bin/python -c "import pkg_resources; print('found at:', pkg_resources.__file__)"\
\
# make sure nanomix is installed\
python -c "import nanomix; print('success:', nanomix.__file__)"\
\
#2. Make adjustments to the atlas - only need immune cell types\
\pard\pardeftab720\partightenfactor0
\cf0 \expnd0\expndtw0\kerning0
awk 'BEGIN\{FS=OFS="\\t"\}\
NR==1 \{\
  for (i=1; i<=NF; i++) \{\
    if ($i=="chr" || $i=="start" || $i=="end" || $i=="B-cell" || $i=="NK-cell" || $i=="T-cell" || $i=="monocyte" || $i=="granulocyte") \{\
      keep[i]=1\
    \}\
  \}\
\}\
\{\
  out=""\
  for (i=1; i<=NF; i++) \{\
    if (i in keep) \{\
      out = out (out=="" ? "" : OFS) $i\
    \}\
  \}\
  print out\
\}' ~/39Bisulfite.tsv > ~/atlas_immune_cell.tsv\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf0 \kerning1\expnd0\expndtw0 \
# 3. Check that deconvolution is working\
# check one of the tsv files\
\pard\pardeftab720\partightenfactor0
\cf0 \expnd0\expndtw0\kerning0
head -1 ~/xx5_filtered.tsv | tr '\\t' \'91\\n\'92 # colnames\kerning1\expnd0\expndtw0 \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 \
# run nanomix on 1 file to see if it\'92s working properly\
/Users/hazelmilla/miniconda3/envs/nanomix/bin/nanomix deconvolute -a atlas_immune_cell.tsv xx13_nanomix.tsv > test_output_deconv.tsv\
\
# 4. Run deconvolution\
for f in ~/*_filtered.tsv; do\
    sample=$(basename "$f" _filtered.tsv)\
    echo "Processing $sample..."\
\
    # Step 1: filter and rename columns\
    nanomix_input="$\{f/_filtered.tsv/_nanomix.tsv\}"\
    awk -F'\\t' '\
    NR==1\{for(i=1;i<=NF;i++) col[$i]=i; print "chr\\tstart\\tend\\ttotal_calls\\tmodified_calls"; next\}\
    \{print $col["chrom"]"\\t"$col["start"]"\\t"$col["end"]"\\t"$col["Nvalid_cov"]"\\t"$col["Nmod"]\}\
    ' "$f" > "$nanomix_input"\
    echo "  Filtered: $nanomix_input"\
\
    # Step 2: deconvolute\
    deconv_output="$\{f/_filtered.tsv/_deconv.tsv\}"\
    /Users/hazelmilla/miniconda3/envs/nanomix/bin/nanomix deconvolute -a atlas_immune_cell.tsv "$nanomix_input" > "$deconv_output"\
    echo "  Deconvoluted: $deconv_output"\
done\
}