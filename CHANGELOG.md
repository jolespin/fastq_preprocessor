**Changes:**
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
