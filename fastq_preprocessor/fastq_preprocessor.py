#!/usr/bin/env python

# Built-ins
import os, sys, argparse, importlib

# Version
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.11.29"

# Accepted modules
accepted_programs = ["short", "long"]
script_directory  =  os.path.dirname(os.path.abspath( __file__ ))

# Controller
def main(argv=None):
    parser = argparse.ArgumentParser(prog="fastq_preprocessor",description="A fastq preprocessor for short and long read sequencing. Optional contamination removal.", add_help=True)
    parser.add_argument("program", choices=accepted_programs, help="`fastq_preprocessor` program for preprocessing. `short` for Illumina and `long` for ONT/PacBio.")
    parser.add_argument("-c", "--citation", action='version', help="Show full citation (doi: 10.1186/s12859-022-04973-8)", version="Espinoza JL, Dupont CL.\nVEBA: a modular end-to-end suite for in silico recovery, clustering, and analysis of prokaryotic, microeukaryotic, and viral genomes from metagenomes.\nBMC Bioinformatics. 2022 Oct 12;23(1):419. doi: 10.1186/s12859-022-04973-8. PMID: 36224545.")
    parser.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    opts = parser.parse_args(argv)
    return opts.program


# Initialize
if __name__ == "__main__":
    # Check version
    python_version = sys.version.split(" ")[0]
    condition_1 = int(python_version.split(".")[0]) == 3
    condition_2 = int(python_version.split(".")[1]) >= 6
    assert all([condition_1, condition_2]), "Python version must be >= 3.6.  You are running: {}\n{}".format(python_version, sys.executable)
    # Get the algorithm
    program = main([sys.argv[1]])
    module = importlib.import_module("fastq_preprocessor_{}".format(program))
    # module.main(sys.argv[2:])
    module.main(sys.argv[2:])

