
version="2026.4.2"
#conda build -c conda-forge -c bioconda -c jolespin .
#conda install --use-local fastq_preprocess=${version}
anaconda upload \
    $(conda info --base)/conda-bld/noarch/fastq_preprocessor-${version}-py_0.conda
