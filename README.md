
```

 ____________________________ _____      _____  ______ _____ ____________________________ _____  ______
 |______|_____||______   |   |   __|    |_____]|_____/|     ||      |______|______|______|     ||_____/
 |      |     |______|   |   |____\|    |      |    \_|_____||_____ |____________|______||_____||    \_
                                                                                                                                            

```
#### Description:

A fastq preprocessor for short read sequencing. Can be used as a wrapper around `fastp` for a standardized directory structure or can feed the trimmed reads into `bowtie2` if an index is provided to decontaminate sequences similar to `kneaddata` (useful in metagenomics to remove host reads). Also includes functionality to filter based on k-mer profiles and is useful for quantifying the amount of ribosomal reads.  At each stage, `seqkit stats` is run so there are read stats that can be used post hoc.  

#### About:

`__author__ = "Josh L. Espinoza"`

`__cite__ = "TBD"`

`__contact__ = "jespinoz@jcvi.org, jol.espinoz@gmail.com"`

`__developmental__ = True`

`__license__ = "BSD-3"`

`__version__ = "2022.1.13"`

#### Dependencies: 

##### Bioinformatics software:
* bbmap
* fastp
* seqkit
* bowtie2

##### Python packages: 
* pandas
* genopype
* soothsayer_utils

#### Installation: 
```
# Conda [Recommended]
conda install -c jolespin fastq_preprocessor 

# PyPI [See note below]
pip install fastq_preprocessor
``` 

Note: If `pip` is used for installation then it assumes bioinformatics packages are installed in the same `conda` environment. 

#### Usage:

```bash
fastq_preprocessor % fastq_preprocessor -h
usage: fastq_preprocessor -1 <reads_1.fq> -2 <reads_2.fq> -n <name> -o <output_directory> |Optional| -x <reference_index> -k <kmer_database>

    Running: fastq_preprocessor v2022.01.13 via Python v3.9.9 | /Users/jespinoz/anaconda3/envs/test_env/bin/python3.9

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -1 FORWARD_READS, --forward_reads FORWARD_READS
                        path/to/reads_1.fastq
  -2 REVERSE_READS, --reverse_reads REVERSE_READS
                        path/to/reads_2.fastq
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: preprocessed]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv. Must have at least 2 columns [name, executable] [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint
  -v, --version         show program's version number and exit

Fastp arguments:
  -m MINIMUM_READ_LENGTH, --minimum_read_length MINIMUM_READ_LENGTH
                        Fastp | Minimum read length [Default: 75]
  -a ADAPTERS, --adapters ADAPTERS
                        Fastp | path/to/adapters.fasta [Default: detect]
  --fastp_options FASTP_OPTIONS
                        Fastp | More options (e.g. --arg 1 ) [Default: '']

Bowtie2 arguments:
  -x CONTAMINATION_INDEX, --contamination_index CONTAMINATION_INDEX
                        Bowtie2 | path/to/contamination_index
                        (e.g., Human GRCh38 from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids//GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz)
  --retain_trimmed_reads RETAIN_TRIMMED_READS
                        Retain fastp trimmed fastq after decontamination. 0=No, 1=yes [Default: 0]
  --retain_decontaminated_reads RETAIN_DECONTAMINATED_READS
                        Retain decontaminated fastq after decontamination. 0=No, 1=yes [Default: 0]
  --bowtie2_options BOWTIE2_OPTIONS
                        Bowtie2 | More options (e.g. --arg 1 ) [Default: '']
                        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

BBDuk arguments:
  -k KMER_DATABASE, --kmer_database KMER_DATABASE
                        BBDuk | path/to/kmer_database
                        (e.g., ribokmers.fa.gz from https://drive.google.com/file/d/0B3llHR93L14wS2NqRXpXakhFaEk/view?usp=sharing)
  --kmer_size KMER_SIZE
                        BBDuk | k-mer size [Default: 31]
  --retain_kmer_hits RETAIN_KMER_HITS
                        Retain reads that map to k-mer database. 0=No, 1=yes [Default: 0]
  --retain_non_kmer_hits RETAIN_NON_KMER_HITS
                        Retain reads that do not map to k-mer database. 0=No, 1=yes [Default: 0]
  --bbduk_options BBDUK_OPTIONS
                        BBDuk | More options (e.g., --arg 1) [Default: '']

Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)

```

#### Output:

```
=====================
fastq_preprocessor.py
=====================
--------------
Configuration:
--------------
................
Name: test_kmers
................
Python version: 3.7.10 | packaged by conda-forge | (default, Feb 19 2021, 16:07:37)  [GCC 9.3.0]
Python path: /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/python
Script version: 2022.01.12
Moment: 2022-01-13 18:28:40
Directory: /local/ifs3_scratch/CORE/jespinoz/Testing
Commands:
['/home/jespinoz/Algorithms/Scripts/fastq_preprocessor.py', '-1', 'S003_R2_POST-PE-N722-S520-1_S18_R1_001.fastq.gz', '-2', 'S003_R2_POST-PE-N722-S520-1_S18_R2_001.fastq.gz', '-n', 'test_kmers', '-p', '16', '-x', '/usr/local/scratch/CORE/jespinoz/db/genomes/human/GRCh38.p13/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index', '--kmer_database', '/usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bbmap_env/ribokmers.fa.gz', '-m', '75']
------------------------------------------------------------------
Adding executables to path from the following source: CONDA_PREFIX
------------------------------------------------------------------
repair.sh --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/repair.sh
bbduk.sh --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/bbduk.sh
seqkit --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/seqkit
fastp --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/fastp
bowtie2 --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/bowtie2

===========================
. .. ... Compiling ... .. .
===========================
Step: 1, fastp | log_prefix = 1__fastp | Quality trimming and adapter removal
Step: 2, bowtie2 | log_prefix = 2__bowtie2 | Decontaminate reads based on a reference
Step: 3, bbduk | log_prefix = 3__bbduk | Decontaminate reads based on k-mer database
Step: 4, synopsis | log_prefix = 4__synopsis | Symlinking and merging relevant output files
_____________________________________________________
. .. ... fastq_preprocessor.py || test_kmers ... .. .
_____________________________________________________

=========
. fastp .
=========
Input: ['S003_R2_POST-PE-N722-S520-1_S18_R1_001.fastq.gz', 'S003_R2_POST-PE-N722-S520-1_S18_R2_001.fastq.gz']
Output: ['preprocessed/test_kmers/intermediate/1__fastp/seqkit_stats.tsv']

Command:
( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/fastp --in1 S003_R2_POST-PE-N722-S520-1_S18_R1_001.fastq.gz --in2 S003_R2_POST-PE-N722-S520-1_S18_R2_001.fastq.gz --stdout -h preprocessed/test_kmers/intermediate/1__fastp/fastp.html -j preprocessed/test_kmers/intermediate/1__fastp/fastp.json -l 75 --thread 16 --detect_adapter_for_pe  ) | ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/repair.sh in=stdin.fastq out1=preprocessed/test_kmers/intermediate/1__fastp/trimmed_1.fastq.gz out2=preprocessed/test_kmers/intermediate/1__fastp/trimmed_2.fastq.gz outs=preprocessed/test_kmers/intermediate/1__fastp/trimmed_singletons.fastq.gz overwrite=t ) && ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/seqkit stats -T -j 16 S003_R2_POST-PE-N722-S520-1_S18_R1_001.fastq.gz S003_R2_POST-PE-N722-S520-1_S18_R2_001.fastq.gz preprocessed/test_kmers/intermediate/1__fastp/trimmed_1.fastq.gz preprocessed/test_kmers/intermediate/1__fastp/trimmed_2.fastq.gz preprocessed/test_kmers/intermediate/1__fastp/trimmed_singletons.fastq.gz > preprocessed/test_kmers/intermediate/1__fastp/seqkit_stats.tsv )

Validating the following input files:
[=] File exists (2785 MB): S003_R2_POST-PE-N722-S520-1_S18_R1_001.fastq.gz
[=] File exists (2730 MB): S003_R2_POST-PE-N722-S520-1_S18_R2_001.fastq.gz

Running. .. ... .....

Log files:
preprocessed/test_kmers/log/1__fastp.*

Validating the following output files:
[=] File exists (546 bytes): preprocessed/test_kmers/intermediate/1__fastp/seqkit_stats.tsv

Duration: 00:11:18

===========
. bowtie2 .
===========
Input: ['preprocessed/test_kmers/intermediate/1__fastp/trimmed_1.fastq.gz', 'preprocessed/test_kmers/intermediate/1__fastp/trimmed_2.fastq.gz']
Output: ['preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_1.fastq.gz', 'preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_2.fastq.gz', 'preprocessed/test_kmers/intermediate/2__bowtie2/seqkit_stats.tsv']

Command:
( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/bowtie2 -x /usr/local/scratch/CORE/jespinoz/db/genomes/human/GRCh38.p13/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -p 16 -1 preprocessed/test_kmers/intermediate/1__fastp/trimmed_1.fastq.gz -2 preprocessed/test_kmers/intermediate/1__fastp/trimmed_2.fastq.gz --seed 0 --un-gz preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_singletons.fastq.gz --al-gz preprocessed/test_kmers/intermediate/2__bowtie2/contaminated_singletons.fastq.gz --un-conc preprocessed/test_kmers/intermediate/2__bowtie2/TMP__cleaned_%.fastq --al-conc preprocessed/test_kmers/intermediate/2__bowtie2/TMP__contaminated_%.fastq  > /dev/null ) && ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/repair.sh in1=preprocessed/test_kmers/intermediate/2__bowtie2/TMP__cleaned_1.fastq in2=preprocessed/test_kmers/intermediate/2__bowtie2/TMP__cleaned_2.fastq out1=preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_1.fastq.gz out2=preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_2.fastq.gz overwrite=t ) && rm preprocessed/test_kmers/intermediate/2__bowtie2/TMP__cleaned_*.fastq && ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/repair.sh in1=preprocessed/test_kmers/intermediate/2__bowtie2/TMP__contaminated_1.fastq in2=preprocessed/test_kmers/intermediate/2__bowtie2/TMP__contaminated_2.fastq out1=preprocessed/test_kmers/intermediate/2__bowtie2/contaminated_1.fastq.gz out2=preprocessed/test_kmers/intermediate/2__bowtie2/contaminated_2.fastq.gz overwrite=t ) && rm preprocessed/test_kmers/intermediate/2__bowtie2/TMP__contaminated_*.fastq && ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/seqkit stats -T -j 16 preprocessed/test_kmers/intermediate/2__bowtie2/*.fastq.gz > preprocessed/test_kmers/intermediate/2__bowtie2/seqkit_stats.tsv ) && rm -rf preprocessed/test_kmers/intermediate/1__fastp/*.fastq.gz && rm -rf preprocessed/test_kmers/intermediate/2__bowtie2/contaminated_*.fastq.gz

Running. .. ... .....

Log files:
preprocessed/test_kmers/log/2__bowtie2.*

Validating the following output files:
[=] File exists (1325 MB): preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_1.fastq.gz
[=] File exists (1303 MB): preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_2.fastq.gz
[=] File exists (693 bytes): preprocessed/test_kmers/intermediate/2__bowtie2/seqkit_stats.tsv

Duration: 00:36:48

=========
. bbduk .
=========
Input: ['preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_1.fastq.gz', 'preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_2.fastq.gz']
Output: ['preprocessed/test_kmers/intermediate/3__bbduk/seqkit_stats.tsv']

Command:
( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/bbduk.sh zl=1 overwrite=t threads=16 in1=preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_1.fastq.gz in2=preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_2.fastq.gz k=31 minlen=75 out1=preprocessed/test_kmers/intermediate/3__bbduk/non-kmer_hits_1.fastq.gz out2=preprocessed/test_kmers/intermediate/3__bbduk/non-kmer_hits_2.fastq.gz outm1=preprocessed/test_kmers/intermediate/3__bbduk/kmer_hits_1.fastq.gz outm2=preprocessed/test_kmers/intermediate/3__bbduk/kmer_hits_2.fastq.gz  ) && ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/seqkit stats -T -j 16 preprocessed/test_kmers/intermediate/3__bbduk/*.fastq.gz > preprocessed/test_kmers/intermediate/3__bbduk/seqkit_stats.tsv ) && rm -rf preprocessed/test_kmers/intermediate/3__bbduk/kmer_hits_*.fastq.gz && rm -rf preprocessed/test_kmers/intermediate/3__bbduk/non-kmer_hits_*.fastq.gz

Running. .. ... .....

Log files:
preprocessed/test_kmers/log/3__bbduk.*

Validating the following output files:
[=] File exists (448 bytes): preprocessed/test_kmers/intermediate/3__bbduk/seqkit_stats.tsv

Duration: 00:38:13

============
. synopsis .
============
Input: ['preprocessed/test_kmers/intermediate/*/*.fastq.gz']
Output: ['preprocessed/test_kmers/output/*.fastq.gz', 'preprocessed/test_kmers/output/seqkit_stats.concatenated.tsv']

Command:

python -c "import glob, pandas as pd; pd.concat(map(lambda fp: pd.read_csv(fp, sep='	', index_col=0), glob.glob('preprocessed/test_kmers/intermediate/*/seqkit_stats.tsv')), axis=0).to_csv('preprocessed/test_kmers/output/seqkit_stats.concatenated.tsv', sep='	')"
 ln -f -s /local/ifs3_scratch/CORE/jespinoz/Testing/preprocessed/test_kmers/intermediate/*/*.fastq.gz preprocessed/test_kmers/output

Validating the following input files:
[=] File exists (20 bytes): preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_singletons.fastq.gz
[=] File exists (1303 MB): preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_2.fastq.gz
[=] File exists (1325 MB): preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_1.fastq.gz

Running. .. ... .....

Log files:
preprocessed/test_kmers/log/4__synopsis.*

Validating the following output files:
[=] File exists (20 bytes): /local/ifs3_scratch/CORE/jespinoz/Testing/preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_singletons.fastq.gz
[=] File exists (1303 MB): /local/ifs3_scratch/CORE/jespinoz/Testing/preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_2.fastq.gz
[=] File exists (1325 MB): /local/ifs3_scratch/CORE/jespinoz/Testing/preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_1.fastq.gz
[=] File exists (1571 bytes): preprocessed/test_kmers/output/seqkit_stats.concatenated.tsv

Duration: 00:38:14


........................
Total duration: 00:38:14
........................

# Directory structure
test_kmers/
├── checkpoints
│   ├── 1__fastp
│   ├── 2__bowtie2
│   ├── 3__bbduk
│   └── 4__synopsis
├── commands.sh
├── intermediate
│   ├── 1__fastp
│   │   ├── fastp.html
│   │   ├── fastp.json
│   │   └── seqkit_stats.tsv
│   ├── 2__bowtie2
│   │   ├── cleaned_1.fastq.gz
│   │   ├── cleaned_2.fastq.gz
│   │   ├── cleaned_singletons.fastq.gz
│   │   └── seqkit_stats.tsv
│   └── 3__bbduk
│       └── seqkit_stats.tsv
├── log
│   ├── 1__fastp.e
│   ├── 1__fastp.o
│   ├── 1__fastp.returncode
│   ├── 2__bowtie2.e
│   ├── 2__bowtie2.o
│   ├── 2__bowtie2.returncode
│   ├── 3__bbduk.e
│   ├── 3__bbduk.o
│   ├── 3__bbduk.returncode
│   ├── 4__synopsis.e
│   ├── 4__synopsis.o
│   └── 4__synopsis.returncode
├── output
│   ├── cleaned_1.fastq.gz -> /local/ifs3_scratch/CORE/jespinoz/Testing/preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_1.fastq.gz
│   ├── cleaned_2.fastq.gz -> /local/ifs3_scratch/CORE/jespinoz/Testing/preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_2.fastq.gz
│   ├── cleaned_singletons.fastq.gz -> /local/ifs3_scratch/CORE/jespinoz/Testing/preprocessed/test_kmers/intermediate/2__bowtie2/cleaned_singletons.fastq.gz
│   └── seqkit_stats.concatenated.tsv
└── tmp

8 directories, 29 files```