#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, atexit
from subprocess import (
    run,
    PIPE,
    DEVNULL,
)

from loguru import logger

from ._utils import (
    create_directory,
    format_path,
)

__program__ = os.path.split(sys.argv[0])[-1]
from fastq_preprocessor import __version__

# =============
# FIFO cleanup
# =============
fifo_paths = []

def cleanup_fifos():
    for p in fifo_paths:
        if os.path.exists(p):
            os.remove(p)

atexit.register(cleanup_fifos)

# =============
# Functions
# =============
def get_strobealign_help():
    try:
        result = run(
            ["strobealign", "-h"],
            capture_output=True,
            text=True,
        )
        return result.stderr or result.stdout
    except FileNotFoundError:
        return "strobealign not found in PATH"


def resolve_output_paths(path):
    if "%" in path:
        return {
            "paired": True,
            "path_1": path.replace("%", "_1"),
            "path_2": path.replace("%", "_2"),
            "path": path,
        }
    else:
        return {
            "paired": False,
            "path_1": None,
            "path_2": None,
            "path": path,
        }


def build_consumer_cmd(samtools_flag, samtools_threads, fifo_or_stdin, output_info, no_repair):
    paired = output_info["paired"]

    if paired and not no_repair:
        # Paired with repair
        cmd = (
            "samtools fastq -@ {threads} {flag} -0 /dev/null -s /dev/null -n {input}"
            " | repair.sh repair=t overwrite=t threads=1 in=stdin.fastq"
            " out1={out1} out2={out2}"
        ).format(
            threads=samtools_threads,
            flag=samtools_flag,
            input=fifo_or_stdin,
            out1=output_info["path_1"],
            out2=output_info["path_2"],
        )
    elif paired and no_repair:
        # Paired without repair
        cmd = (
            "samtools fastq -@ {threads} {flag}"
            " -1 {out1} -2 {out2}"
            " -0 /dev/null -s /dev/null -n {input}"
        ).format(
            threads=samtools_threads,
            flag=samtools_flag,
            out1=output_info["path_1"],
            out2=output_info["path_2"],
            input=fifo_or_stdin,
        )
    elif not paired and not no_repair:
        # Interleaved with repair
        cmd = (
            "samtools fastq -@ {threads} {flag} -0 /dev/null -s /dev/null -n {input}"
            " | repair.sh repair=t overwrite=t threads=1 in=stdin.fastq"
            " out={out}"
        ).format(
            threads=samtools_threads,
            flag=samtools_flag,
            input=fifo_or_stdin,
            out=output_info["path"],
        )
    else:
        # Interleaved without repair
        cmd = (
            "samtools fastq -@ {threads} {flag}"
            " -0 /dev/null -s /dev/null -n {input}"
            " > {out}"
        ).format(
            threads=samtools_threads,
            flag=samtools_flag,
            input=fifo_or_stdin,
            out=output_info["path"],
        )

    return cmd


def build_cmd(opts, strobealign_extra_args):
    parts = []

    # Resolve output paths
    mapped_info = resolve_output_paths(opts.mapped_fastq) if opts.mapped_fastq else None
    unmapped_info = resolve_output_paths(opts.unmapped_fastq) if opts.unmapped_fastq else None

    # Create output directories
    for info in [mapped_info, unmapped_info]:
        if info is not None:
            if info["paired"]:
                create_directory(os.path.dirname(info["path_1"]))
                create_directory(os.path.dirname(info["path_2"]))
            else:
                create_directory(os.path.dirname(info["path"]))

    # Build producer
    extra_args = " ".join(strobealign_extra_args) if strobealign_extra_args else ""
    producer_tokens = [
        "strobealign -t {threads}{extra} {ref} {r1} {r2}".format(
            threads=opts.threads,
            extra=" " + extra_args if extra_args else "",
            ref=opts.reference,
            r1=opts.reads1,
            r2=opts.reads2,
        ),
        "samtools view -h -",
        "samtools collate -u -O - {prefix}".format(prefix=opts.temporary_prefix),
    ]
    producer = " | ".join(producer_tokens)

    has_mapped = mapped_info is not None
    has_unmapped = unmapped_info is not None

    if has_mapped and has_unmapped:
        # Both outputs: tee + named pipes
        mapped_fifo = "{}_mapped_fifo".format(opts.temporary_prefix)
        unmapped_fifo = "{}_unmapped_fifo".format(opts.temporary_prefix)
        fifo_paths.extend([mapped_fifo, unmapped_fifo])

        parts.append("mkfifo {} {}".format(mapped_fifo, unmapped_fifo))
        parts.append("")

        mapped_consumer = build_consumer_cmd(
            "-F 12", opts.samtools_threads, mapped_fifo, mapped_info, opts.no_repair,
        )
        parts.append("{} &".format(mapped_consumer))
        parts.append("")

        unmapped_consumer = build_consumer_cmd(
            "-f 12", opts.samtools_threads, unmapped_fifo, unmapped_info, opts.no_repair,
        )
        parts.append("{} &".format(unmapped_consumer))
        parts.append("")

        parts.append("{} | tee {} > {}".format(producer, mapped_fifo, unmapped_fifo))
        parts.append("")
        parts.append("wait")
        parts.append("rm -f {} {}".format(mapped_fifo, unmapped_fifo))

    elif has_mapped:
        # Only mapped: direct pipe with -F 12
        mapped_consumer = build_consumer_cmd(
            "-F 12", opts.samtools_threads, "-", mapped_info, opts.no_repair,
        )
        parts.append("{} | {}".format(producer, mapped_consumer))

    elif has_unmapped:
        # Only unmapped: direct pipe with -f 12
        unmapped_consumer = build_consumer_cmd(
            "-f 12", opts.samtools_threads, "-", unmapped_info, opts.no_repair,
        )
        parts.append("{} | {}".format(producer, unmapped_consumer))

    return "\n".join(parts)


# =============
# Main
# =============
def main(args=None):
    # Setup
    script_directory = os.path.dirname(os.path.abspath(__file__))
    script_filename = __program__

    description = """
    Running: {} v{} via Python v{} | {}""".format(
        __program__, __version__, sys.version.split(" ")[0], sys.executable,
    )
    usage = "{} [strobealign_options] reference reads1 reads2 --mapped_fastq PATH --unmapped_fastq PATH".format(__program__)
    epilog = "Strobealign help:\n\n{}".format(get_strobealign_help())

    # Parser
    parser = argparse.ArgumentParser(
        description=description,
        usage=usage,
        epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("reference", type=str, help="Reference FASTA for strobealign")
    parser_io.add_argument("reads1", type=str, help="Forward reads FASTQ")
    parser_io.add_argument("reads2", type=str, help="Reverse reads FASTQ")
    parser_io.add_argument("--mapped_fastq", type=str, default=None,
                           help="Output path for mapped reads.\nUse %% for paired: mapped_%%.fastq.gz -> mapped_1.fastq.gz, mapped_2.fastq.gz\nWithout %%: interleaved output")
    parser_io.add_argument("--unmapped_fastq", type=str, default=None,
                           help="Output path for unmapped reads.\nUse %% for paired: unmapped_%%.fastq.gz -> unmapped_1.fastq.gz, unmapped_2.fastq.gz\nWithout %%: interleaved output")

    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("-t", "--threads", type=int, default=1, help="Threads for strobealign [Default: 1]")
    parser_utility.add_argument("--samtools_threads", type=int, default=1, help="Threads for samtools fastq gzip compression [Default: 1]")
    parser_utility.add_argument("--no_repair", action="store_true", default=False, help="Disable repair.sh post-processing")
    parser_utility.add_argument("-T", "--temporary_prefix", type=str, default="/tmp/strobealign_split_reads",
                                help="Temp prefix for samtools collate and named pipes [Default: /tmp/strobealign_split_reads]")
    parser_utility.add_argument("-v", "--version", action="version", version="{} v{}".format(__program__, __version__))

    # Parse
    opts, strobealign_extra_args = parser.parse_known_args(args)

    # Validate
    if opts.mapped_fastq is None and opts.unmapped_fastq is None:
        parser.error("At least one of --mapped_fastq or --unmapped_fastq must be provided")

    for label, path in [("--mapped_fastq", opts.mapped_fastq), ("--unmapped_fastq", opts.unmapped_fastq)]:
        if path is not None and path.count("%") > 1:
            parser.error("{} path must contain at most one '%' wildcard".format(label))

    # Resolve paths
    opts.reference = format_path(opts.reference)
    opts.reads1 = format_path(opts.reads1)
    opts.reads2 = format_path(opts.reads2)
    if opts.mapped_fastq is not None:
        opts.mapped_fastq = format_path(opts.mapped_fastq)
    if opts.unmapped_fastq is not None:
        opts.unmapped_fastq = format_path(opts.unmapped_fastq)

    assert os.path.exists(opts.reference), "Reference not found: {}".format(opts.reference)
    assert os.path.exists(opts.reads1), "Forward reads not found: {}".format(opts.reads1)
    assert os.path.exists(opts.reads2), "Reverse reads not found: {}".format(opts.reads2)

    # Configure logger
    logger.remove()
    logger.add(sys.stdout, format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | {message}")
    logger.add(sys.stderr, format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | {message}", level="ERROR")

    # Log parameters
    logger.info("=" * 60)
    logger.info("Program: {}", __program__)
    logger.info("Version: {}", __version__)
    logger.info("Reference: {}", opts.reference)
    logger.info("Reads 1: {}", opts.reads1)
    logger.info("Reads 2: {}", opts.reads2)
    logger.info("Mapped FASTQ: {}", opts.mapped_fastq)
    logger.info("Unmapped FASTQ: {}", opts.unmapped_fastq)
    logger.info("Threads: {}", opts.threads)
    logger.info("Samtools threads: {}", opts.samtools_threads)
    logger.info("Repair: {}", not opts.no_repair)
    logger.info("Temporary prefix: {}", opts.temporary_prefix)
    if strobealign_extra_args:
        logger.info("Strobealign extra args: {}", " ".join(strobealign_extra_args))
    logger.info("=" * 60)

    # Build and execute
    cmd_str = build_cmd(opts, strobealign_extra_args)
    logger.info("Command:\n{}", cmd_str)
    run(["bash", "-c", cmd_str], check=True)
    logger.success("Completed successfully")


if __name__ == "__main__":
    main()
