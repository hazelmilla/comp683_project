{\rtf1\ansi\ansicpg1252\cocoartf2867
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 ArialMT;\f1\fnil\fcharset0 .AppleSystemUIFontMonospaced-Regular;}
{\colortbl;\red255\green255\blue255;\red97\green170\blue255;\red0\green0\blue0;\red230\green232\blue236;
\red139\green231\blue78;\red140\green232\blue80;\red139\green231\blue79;\red96\green169\blue255;\red140\green232\blue80;
\red96\green168\blue255;\red229\green232\blue236;\red140\green232\blue80;\red237\green244\blue251;}
{\*\expandedcolortbl;;\cssrgb\c44671\c73040\c100000;\cssrgb\c0\c0\c0;\cssrgb\c91987\c92775\c94114;
\cssrgb\c60296\c90813\c37808;\cssrgb\c60541\c91094\c38317;\cssrgb\c60419\c90954\c38063;\cssrgb\c44422\c72747\c100000;\cssrgb\c60663\c91234\c38571;
\cssrgb\c44173\c72453\c100000;\cssrgb\c91876\c92662\c94116;\cssrgb\c60663\c91234\c38571;\cssrgb\c94359\c96719\c98837;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 \
# Check data in tsv file\
\
\pard\pardeftab720\partightenfactor0
\cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 head\strokec3  -1 ~/ont_m.tsv \strokec4 |\strokec3  \strokec2 tr\strokec3  \strokec5 '\\t'\strokec3  \strokec5 \'91\\n\'92 \strokec6 # colnames\
\strokec7 head\strokec6  -2 ~/ont_m.tsv | \strokec7 tail\strokec6  -1 | tr '\\t' '\\n' # first row\
\strokec8 cut\strokec3  -f1 ~/ont_m.tsv \strokec4 |\strokec3  \strokec2 tail\strokec3  -n +2 # no rownames\
\
# We need a tsv file with the following columns:\
# \{chr, start, end, total_calls, modified_calls\}\
\
\strokec9 # We have chrom, start, end, \'85. \kerning1\expnd0\expndtw0 \CocoaLigature0 \outl0\strokewidth0 Nvalid_cov, \expnd0\expndtw0\kerning0
\CocoaLigature1 \outl0\strokewidth0 \strokec9 Nmod\
\
awk -F'\\t' '\
NR==1\{for(i=1;i<=NF;i++) col[$i]=i; print "chr\\tstart\\tend\\ttotal_calls\\tmodified_calls"; next\}\
\{print $col["chrom"]"\\t"$col["start"]"\\t"$col["end"]"\\t"$col["Nvalid_cov"]"\\t"$col["Nmod"]\}\
' ~/ont_m.tsv > ~/ont_m_nanomix.tsv\
\
\pard\pardeftab720\partightenfactor0
\cf0 \outl0\strokewidth0 head -1 ~/ont_m_nanomix.tsv | tr '\\t' \'91\\n\'92 # colnames\outl0\strokewidth0 \strokec9 \
\outl0\strokewidth0 head -2 ~/ont_m_nanomix.tsv | tail -1 | tr '\\t' '\\n' # first row\kerning1\expnd0\expndtw0 \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 \
# Run Nanomix to get cell type proportions\
# https://github.com/Jonbroad15/nanomix/blob/main/README.md\
\
# 1. Deactivate any active environment\
conda deactivate\
\
# 2. Verify base is clean\
echo "Current env: $CONDA_DEFAULT_ENV"  # should be empty or base\
\
# 3. Activate nanomix env\
conda activate nanomix\
\
# 4. Verify correct environment and Python version\
echo "Active env: $CONDA_DEFAULT_ENV"   # should be nanomix\
which python                             # should show miniconda3/envs/nanomix/...\
python --version                         # should be 3.10.x\
\
# 5. Set PYTHONPATH\
export PYTHONPATH="/Users/hazelmilla/nanomix/python:$PYTHONPATH"\
\
# 6. Copy .so with expected name\
cp ~/nanomix/python/nanomix/nanomix.cpython-310-darwin.so \\\
   ~/nanomix/python/nanomix/_nanomix.cpython-310-darwin.so\
\
pip install setuptools\
pip install --upgrade pyranges\
\
# 7. Verify Python can find nanomix\
python -c "import sys; print('\\n'.join(sys.path))"\
\
rm /Users/hazelmilla/miniconda3/envs/nanomix/lib/python3.10/site-packages/pkg_resources\
\
cp -r /Users/hazelmilla/miniconda3/envs/nanomix/lib/python3.10/site-packages/pip/_vendor/pkg_resources \\\
      /Users/hazelmilla/miniconda3/envs/nanomix/lib/python3.10/site-packages/pkg_resources\
/Users/hazelmilla/miniconda3/envs/nanomix/bin/python -c "import pkg_resources; print('found at:', pkg_resources.__file__)"\
\
python -c "import nanomix; print('success:', nanomix.__file__)"\
\
# 8. Run deconvolution\
nanomix deconvolute 
\f1 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec13 -
\f0 a \kerning1\expnd0\expndtw0 \outl0\strokewidth0 ~/\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec13 39Bisulfite.tsv \kerning1\expnd0\expndtw0 \outl0\strokewidth0 ~/ont_m_nanomix.tsv > ~/deconv_results.tsv\
}