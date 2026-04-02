**Changes:**
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
