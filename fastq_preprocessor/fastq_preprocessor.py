#!/usr/bin/env python

# Built-ins
import os, sys, argparse, importlib

# Version
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2026.3.24"

# Accepted modules
accepted_programs = ["short", "long"]
script_directory  =  os.path.dirname(os.path.abspath( __file__ ))

# Controller
def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="fastq_preprocessor",
        description="A fastq preprocessor for short and long read sequencing. Optional contamination removal.",
        add_help=False,
    )
    parser.add_argument("program", choices=accepted_programs, nargs="?", default=None, help="`fastq_preprocessor` program for preprocessing. `short` for Illumina and `long` for ONT/PacBio.")
    parser.add_argument("-h", "--help", action="store_true", default=False, help="show this help message and exit")
    parser.add_argument("-c", "--citation", action='version', help="Show full citation (doi: 10.1186/s12859-022-04973-8)", version="Espinoza JL, Dupont CL.\nVEBA: a modular end-to-end suite for in silico recovery, clustering, and analysis of prokaryotic, microeukaryotic, and viral genomes from metagenomes.\nBMC Bioinformatics. 2022 Oct 12;23(1):419. doi: 10.1186/s12859-022-04973-8. PMID: 36224545.")
    parser.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    opts, passthrough = parser.parse_known_args(argv)

    # No subcommand: show top-level help
    if opts.program is None:
        parser.print_help()
        sys.exit(0 if opts.help else 1)

    # Forward --help to subcommand
    if opts.help:
        passthrough.append("--help")

    if opts.program == "short":
        from fastq_preprocessor.fastq_preprocessor_short import main as run
    elif opts.program == "long":
        from fastq_preprocessor.fastq_preprocessor_long import main as run

    run(passthrough)