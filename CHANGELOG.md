**Changes:**

* **2023.11.28**
	* Deprecated `bin/` directory in favor of `fastq_preprocessor/`
	* Wrapped both `fastq_preprocessor_short.py` and `fastq_preprocessor_long.py` with `fastq_preprocessor`. 
	* Changed `fastq_preprocessor.py` to `fastq_preprocessor_short.py` and added `fastq_preprocessor_long.py`.
	* Added `__genopype_version__` to programs


* **2023.7.24**
	* Added `repair.sh` before `FastP` to account for mispaired input reads.
