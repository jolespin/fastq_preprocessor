**Changes:**
* [2026.4.21] - `strobealign_wrapper --coverage/--depth` now writes gzipped `coverage.tsv.gz` and `depth.tsv.gz` by piping `samtools coverage`/`samtools depth -aa -H` through `pigz` (threads = `-t/--threads`, level = `-c/--compression_level`). `coverage_breadth` accepts both `.tsv` and `.tsv.gz` transparently via pandas `compression='infer'`.
* [2026.4.21] - Made `samtools collate` conditional in `strobealign_wrapper` producer pipeline: included only when `--mapped_fastq` or `--unmapped_fastq` is requested. Skipped for BAM-only runs (and `--coverage`/`--depth`) since `samtools sort` does not benefit from name-grouped input [issue/#46](https://github.com/jolespin/fastq_preprocessor/issues/46)
* [2026.4.20] - Added `--contigs_to_genomes` to `coverage_breadth` for per-genome breadth of coverage. Accepts a 2-column TSV (`id_contig<TAB>id_genome`, no header); when provided, output becomes `id_genome<TAB>percentage_of_bases_covered` with one row per genome. Errors if a mapped contig is absent from the depth/coverage data; warns and skips data contigs not in the mapping
* [2026.4.20] - Changed default output of `coverage_breadth` (no `--contigs_to_genomes`) to a 2-column TSV with header `filepath<TAB>percentage_of_bases_covered` (was previously a single-line `percentage_of_bases_covered<TAB><value>`)
* [2026.4.17] - Added `coverage_breadth` CLI tool that computes genome-level coverage breadth (ratio of bases covered at ≥ `--minimum_depth`) from `samtools depth` (with or without `-aa`) or `samtools coverage` output; outputs a single `percentage_of_bases_covered\t<value>` row (0–100 scale) to stdout or a file
* [2026.4.17] - Added `-c/--compression_level` to `fastq_preprocessor_short` (applied to `repair.sh zl=`) and `fastq_preprocessor_long` (applied to `pigz`) for configurable gzip compression of output FASTQ files (default: 6)
* [2026.4.17] - Added random run ID (8-character hex UUID) to `strobealign_wrapper` temporary file names (FIFOs, `samtools collate`/`sort` prefixes) to prevent collisions during simultaneous runs; temporary files are removed on both success and failure via `atexit`
* [2026.4.17] - Changed `samtools depth -a` to `samtools depth -aa` in `strobealign_wrapper` to include all positions even for contigs with no mapped reads
* [2026.4.13] - Symlink `fastp` `.html`/`.json` and `fastplong` `.html`/`.json` report files to the output directory [issue/#30](https://github.com/jolespin/fastq_preprocessor/issues/30)
* [2026.4.13] - Added `--coverage` and `--depth` flags to `strobealign_wrapper` for `samtools coverage` and `samtools depth` output from BAM files [issue/#22](https://github.com/jolespin/fastq_preprocessor/issues/22)
* [2026.4.10] - Renamed `strobealign_split_reads` to `strobealign_wrapper` to better reflect broader purpose (BAM output, future coverage support) [issue/#25](https://github.com/jolespin/fastq_preprocessor/issues/25)
* [2026.4.10] - Fixed `strobealign` index support: replaced `-i/--contamination_index` with `-i/--use_index` boolean flag that maps to strobealign's native `--use-index` [issue/#24](https://github.com/jolespin/fastq_preprocessor/issues/24)
* [2026.4.10] - Added `-c/--compression_level` to `strobealign_split_reads` (default 6) applied to both `samtools fastq -c` and `repair.sh zl=` for consistent gzip compression
* [2026.4.10] - Added `--bam` flag to `strobealign_split_reads` for coordinate-sorted BAM output
* [2026.4.10] - Changed `--temporary_prefix` to `-T/--temporary_directory` in `strobealign_split_reads` with task-specific prefixes derived internally
* [2026.4.10] - Generalized `strobealign_split_reads` consumer architecture to support 1–3 concurrent outputs in any combination
* [2026.4.9] - Replaced `bowtie2` with `strobealign` in `fastq_preprocessor_short` [issue/#3](https://github.com/jolespin/fastq_preprocessor/issues/3)
* [2026.4.9] - Added `strobealign_split_reads` wrapper [issue/#13](https://github.com/jolespin/fastq_preprocessor/issues/13)
* [2026.4.2] - Simplified `__version__` imports
* [2026.4.2] - Changed `loguru` logging to write to `stdout` by default with errors also reported to `stderr`
* [2026.3.24] - Changed default output path from `veba_output/preprocess/` back to `preprocessed/`
* [2026.3.13] - Added `--low_complexity_filter` to both `fastq_preprocessor_short.py` and `fastq_preprocessor_long.py` [issue/#6](https://github.com/jolespin/fastq_preprocessor/issues/6)
* [2026.3.9] - Replaced `chopper` with `fastplong` for long-read quality trimming; fastplong handles gzipped I/O natively and generates HTML/JSON reports
* [2026.3.5] - Removed `GenoPype` dependency to now use `sh` [issue/#7](https://github.com/jolespin/fastq_preprocessor/issues/7)
* [2026.3.3] - Changed `scripts` to `entrypoints` [issue/#5](https://github.com/jolespin/fastq_preprocessor/issues/5)
* [2023.11.30] - Changed default output path from `preprocessed/` to `veba_output/preprocess/`
* [2023.11.28] - Deprecated `bin/` directory in favor of `fastq_preprocessor/`
* [2023.11.28] - Wrapped both `fastq_preprocessor_short.py` and `fastq_preprocessor_long.py` with `fastq_preprocessor`. 
* [2023.11.28] - Changed `fastq_preprocessor.py` to `fastq_preprocessor_short.py` and added `fastq_preprocessor_long.py`.
* [2023.11.28] - Added `__genopype_version__` to programs
* [2023.7.24] - Added `repair.sh` before `FastP` to account for mispaired input reads.
