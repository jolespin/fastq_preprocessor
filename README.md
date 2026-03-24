
```

 ____________________________ _____      _____  ______ _____ ____________________________ _____  ______
 |______|_____||______   |   |   __|    |_____]|_____/|     ||      |______|______|______|     ||_____/
 |      |     |______|   |   |____\|    |      |    \_|_____||_____ |____________|______||_____||    \_
                                                                                                                                            

```
#### Description:

A fastq preprocessor for short and long read sequencing. For short reads, it is a wrapper around `fastp` for a standardized directory structure and can optionally feed the trimmed reads into `bowtie2` if an index is provided to decontaminate sequences similar to `kneaddata` (useful in metagenomics to remove host reads). For long reads, it uses `fastplong` instead of `fastp` and `minimap2` instead of `bowtie2`.

Also includes functionality to filter based on k-mer profiles and is useful for quantifying the amount of ribosomal reads.  At each stage, `seqkit stats` is run so there are read stats that can be used post hoc.  

#### About:

`__developer__ = "Josh L. Espinoza"`

`__cite__ = "Espinoza JL, Dupont CL. VEBA: a modular end-to-end suite for in silico recovery, clustering, and analysis of prokaryotic, microeukaryotic, and viral genomes from metagenomes. BMC Bioinformatics. 2022 Oct 12;23(1):419. doi: 10.1186/s12859-022-04973-8. PMID: 36224545."`

`__contact__ = "jespinoz@jcvi.org, jol.espinoz@gmail.com"`

`__developmental__ = True`

`__license__ = "Apache 2.0"`

`__version__ = "2026.3.24"`

#### Dependencies: 

##### Bioinformatics software:
* bbmap
* fastp
* seqkit
* bowtie2
* minimap2
* samtools
* fastplong

##### Python packages:
* pandas
* sh
* loguru

#### Installation: 
```
# Conda [Recommended]
conda install -c jolespin fastq_preprocessor 

# PyPI [See note below]
pip install fastq_preprocessor
``` 

Note: If `pip` is used for installation then it assumes bioinformatics packages are installed in the same `conda` environment. 

#### Usage:

**Wrapper:** 

```
fastq_preprocessor -h
usage: fastq_preprocessor [-h] [-c] [-v] {short,long}

A fastq preprocessor for short and long read sequencing. Optional contamination removal.

positional arguments:
  {short,long}    `fastq_preprocessor` program for preprocessing. `short` for Illumina and `long` for ONT/PacBio.

optional arguments:
  -h, --help      show this help message and exit
  -c, --citation  Show full citation (doi: 10.1186/s12859-022-04973-8)
  -v, --version   show program's version number and exit
```

**Illumina reads (`short`)**

```
fastq_preprocessor short -h
usage: fastq_preprocessor -1 <reads_1.fq> -2 <reads_2.fq> -n <name> -o <output_directory> |Optional| -x <reference_index> -k <kmer_database>

    Running: fastq_preprocessor v2026.3.24 via Python v3.9.9 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -1 FORWARD_READS, --forward_reads FORWARD_READS
                        path/to/reads_1.fastq[.gz]
  -2 REVERSE_READS, --reverse_reads REVERSE_READS
                        path/to/reads_2.fastq[.gz]
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
  --low_complexity_filter LOW_COMPLEXITY_FILTER
                        Fastp | Enable low complexity filter. 0=No, 1=Yes [Default: 1]
  --fastp_options FASTP_OPTIONS
                        Fastp | More options (e.g. --arg 1 ) [Default: '']

Bowtie2 arguments:
  -x CONTAMINATION_INDEX, --contamination_index CONTAMINATION_INDEX
                        Bowtie2 | path/to/contamination_index
                        (e.g., Human T2T assembly from https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip)
  --retain_trimmed_reads RETAIN_TRIMMED_READS
                        Retain fastp trimmed fastq after decontamination. 0=No, 1=yes [Default: 0]
  --retain_contaminated_reads RETAIN_CONTAMINATED_READS
                        Retain contaminated fastq after decontamination. 0=No, 1=yes [Default: 0]
  --bowtie2_options BOWTIE2_OPTIONS
                        Bowtie2 | More options (e.g. --arg 1 ) [Default: '']
                        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

BBDuk arguments:
  -k KMER_DATABASE, --kmer_database KMER_DATABASE
                        BBDuk | path/to/kmer_database
                        (e.g., ribokmers.fa.gz from https://figshare.com/ndownloader/files/36220587)
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

**Oxford Nanopore and PacBio reads (`long`)**

```
fastq_preprocessor long -h
usage: fastq_preprocessor -i <reads.fq[.gz]> -n <name> -o <output_directory> |Optional| -x <reference_index> -k <kmer_database>

    Running: fastq_preprocessor v2026.3.24 via Python v3.9.9 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i READS, --reads READS
                        path/to/reads.fastq[.gz]
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

Fastplong arguments:
  -m MINIMUM_READ_LENGTH, --minimum_read_length MINIMUM_READ_LENGTH
                        Fastplong | Minimum read length [Default: 500]
  -q MINIMUM_QUALITY_SCORE, --minimum_quality_score MINIMUM_QUALITY_SCORE
                        Fastplong | Minimum quality score [Default: 10]
  --low_complexity_filter LOW_COMPLEXITY_FILTER
                        Fastplong | Enable low complexity filter. 0=No, 1=Yes [Default: 1]
  --fastplong_options FASTPLONG_OPTIONS
                        Fastplong | More options (e.g. --arg 1 ) https://github.com/OpenGene/fastplong [Default: '']

MiniMap2 arguments:
  -x CONTAMINATION_INDEX, --contamination_index CONTAMINATION_INDEX
                        MiniMap2 | path/to/contamination_index
                        (e.g., Human T2T assembly from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz)
  --minimap2_preset MINIMAP2_PRESET
                        MiniMap2 | MiniMap2 preset {map-pb, map-ont, map-hifi} [Default: map-ont]
  --retain_trimmed_reads RETAIN_TRIMMED_READS
                        Retain Fastplong trimmed fastq after decontamination. 0=No, 1=yes [Default: 0]
  --retain_contaminated_reads RETAIN_CONTAMINATED_READS
                        Retain contaminated fastq after decontamination. 0=No, 1=yes [Default: 0]
  --minimap2_options MINIMAP2_OPTIONS
                        MiniMap2 | More options (e.g. --arg 1 ) [Default: '']
                        https://github.com/lh3/minimap2

BBDuk arguments:
  -k KMER_DATABASE, --kmer_database KMER_DATABASE
                        BBDuk | path/to/kmer_database
                        (e.g., ribokmers.fa.gz from https://figshare.com/ndownloader/files/36220587)
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

