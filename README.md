
```

 ____________________________ _____      _____  ______ _____ ____________________________ _____  ______
 |______|_____||______   |   |   __|    |_____]|_____/|     ||      |______|______|______|     ||_____/
 |      |     |______|   |   |____\|    |      |    \_|_____||_____ |____________|______||_____||    \_
                                                                                                                                            

```
#### Description:

A fastq preprocessor for short and long read sequencing. For short reads, it is a wrapper around `fastp` for a standardized directory structure and can optionally feed the trimmed reads into `strobealign` if a reference FASTA or pre-built index is provided to decontaminate sequences similar to `kneaddata` (useful in metagenomics to remove host reads). For long reads, it uses `fastplong` instead of `fastp` and `minimap2` instead of `strobealign`.

Also includes functionality to filter based on k-mer profiles and is useful for quantifying the amount of ribosomal reads.  At each stage, `seqkit stats` is run so there are read stats that can be used post hoc.  

#### Citations:
Espinoza JL, Phillips A, Prentice MB, Tan GS, Kamath PL, Lloyd KG, Dupont CL. Unveiling the microbial realm with VEBA 2.0: a modular bioinformatics suite for end-to-end genome-resolved prokaryotic, (micro)eukaryotic and viral multi-omics from either short- or long-read sequencing. Nucleic Acids Res. 2024 Aug 12;52(14):e63. doi: 10.1093/nar/gkae528. PMID: 38909293; PMCID: PMC11317156.

Espinoza JL, Dupont CL. VEBA: a modular end-to-end suite for in silico recovery, clustering, and analysis of prokaryotic, microeukaryotic, and viral genomes from metagenomes. BMC Bioinformatics. 2022 Oct 12;23(1):419. doi: 10.1186/s12859-022-04973-8. PMID: 36224545."


#### Dependencies: 

##### Bioinformatics software:
* seqkit
* bbmap
* fastp
* fastplong
* strobealign
* minimap2
* samtools

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
usage: fastq_preprocessor -1 <reads_1.fq> -2 <reads_2.fq> -n <name> -o <output_directory> |Optional| -x <contamination_reference.fasta> [-i/--use_index] -k <kmer_database>

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

Strobealign arguments:
  -x CONTAMINATION_FASTA, --contamination_fasta CONTAMINATION_FASTA
                        Strobealign | path/to/contamination_reference.fasta[.gz]
                        (e.g., Human T2T assembly from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz)
  -i, --use_index       Strobealign | Use pre-built .sti index (must exist next to the FASTA). Maps to strobealign --use-index.
  --retain_trimmed_reads RETAIN_TRIMMED_READS
                        Retain fastp trimmed fastq after decontamination. 0=No, 1=yes [Default: 0]
  --retain_contaminated_reads RETAIN_CONTAMINATED_READS
                        Retain contaminated fastq after decontamination. 0=No, 1=yes [Default: 0]
  --strobealign_options STROBEALIGN_OPTIONS
                        Strobealign | More options (e.g. --arg 1 ) [Default: '']
                        https://github.com/ksahlin/strobealign

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

```

**Oxford Nanopore and PacBio reads (`long`)**

```
fastq_preprocessor long -h

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
```

**Strobealign wrapper (`strobealign_wrapper`)**

Wraps `strobealign` to split aligned output into mapped and unmapped paired-end FASTQ files and/or coordinate-sorted BAM files. All unrecognized arguments are forwarded directly to `strobealign`.

```
strobealign_wrapper -h
usage: strobealign_wrapper [strobealign_options] reference reads1 reads2 --mapped_fastq PATH --unmapped_fastq PATH --bam PATH

Required I/O arguments:
  reference             Reference FASTA for strobealign
  reads1                Forward reads FASTQ
  reads2                Reverse reads FASTQ
  --mapped_fastq        Output path for mapped reads.
                        Use % for paired: mapped_%.fastq.gz -> mapped_1.fastq.gz, mapped_2.fastq.gz
                        Without %: interleaved output. % is replaced with 1 or 2.
  --unmapped_fastq      Output path for unmapped reads.
                        Use % for paired: unmapped_%.fastq.gz -> unmapped_1.fastq.gz, unmapped_2.fastq.gz
                        Without %: interleaved output. % is replaced with 1 or 2.
  --bam                 Output path for coordinate-sorted BAM file
  --coverage            Run samtools coverage on BAM and write {bam}.coverage.tsv (requires --bam)
  --depth               Run samtools depth -a on BAM and write {bam}.depth.tsv (requires --bam)

Utility arguments:
  -t, --threads         Threads for strobealign [Default: 1]
  --samtools_threads    Threads for samtools fastq/sort [Default: 1]
  --no_repair           Disable repair.sh post-processing
  -c, --compression_level
                        Compression level [0..9] for samtools fastq bgzf output [Default: 6]
  -T, --temporary_directory
                        Temporary directory for samtools collate/sort and named pipes [Default: /tmp/strobealign_wrapper]
  -v, --version         show program's version number and exit
```

Example — split into paired files:
```
strobealign_wrapper -t 8 \
  reference.fasta.gz \
  reads_1.fastq.gz \
  reads_2.fastq.gz \
  --mapped_fastq output/mapped_%.fastq.gz \
  --unmapped_fastq output/unmapped_%.fastq.gz
```

Example — split into interleaved files:
```
strobealign_wrapper -t 8 \
  reference.fasta.gz \
  reads_1.fastq.gz \
  reads_2.fastq.gz \
  --mapped_fastq output/mapped.fastq.gz \
  --unmapped_fastq output/unmapped.fastq.gz
```

Example — BAM with coverage and depth stats:
```
strobealign_wrapper -t 8 \
  reference.fasta.gz \
  reads_1.fastq.gz \
  reads_2.fastq.gz \
  --bam output/aligned.bam \
  --coverage --depth
```

