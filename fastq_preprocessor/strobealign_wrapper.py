#!/usr/bin/env python

from __future__ import print_function, division
import sys, os, argparse, atexit, uuid
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
run_id = uuid.uuid4().hex[:8]

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
            "path_1": path.replace("%", "1"),
            "path_2": path.replace("%", "2"),
            "path": path,
        }
    else:
        return {
            "paired": False,
            "path_1": None,
            "path_2": None,
            "path": path,
        }


def build_consumer_cmd(samtools_flag, samtools_threads, fifo_or_stdin, output_info, no_repair, compression_level=6):
    paired = output_info["paired"]

    if paired and not no_repair:
        # Paired with repair
        cmd = (
            "samtools fastq -@ {threads} -c {compression_level} {flag} -0 /dev/null -s /dev/null -n {input}"
            " | repair.sh repair=t overwrite=t zl={compression_level} threads=1 in=stdin.fastq"
            " out1={out1} out2={out2}"
        ).format(
            threads=samtools_threads,
            compression_level=compression_level,
            flag=samtools_flag,
            input=fifo_or_stdin,
            out1=output_info["path_1"],
            out2=output_info["path_2"],
        )
    elif paired and no_repair:
        # Paired without repair
        cmd = (
            "samtools fastq -@ {threads} -c {compression_level} {flag}"
            " -1 {out1} -2 {out2}"
            " -0 /dev/null -s /dev/null -n {input}"
        ).format(
            threads=samtools_threads,
            compression_level=compression_level,
            flag=samtools_flag,
            out1=output_info["path_1"],
            out2=output_info["path_2"],
            input=fifo_or_stdin,
        )
    elif not paired and not no_repair:
        # Interleaved with repair
        cmd = (
            "samtools fastq -@ {threads} -c {compression_level} {flag} -0 /dev/null -s /dev/null -n {input}"
            " | repair.sh repair=t overwrite=t zl={compression_level} threads=1 in=stdin.fastq"
            " out={out}"
        ).format(
            threads=samtools_threads,
            compression_level=compression_level,
            flag=samtools_flag,
            input=fifo_or_stdin,
            out=output_info["path"],
        )
    else:
        # Interleaved without repair
        cmd = (
            "samtools fastq -@ {threads} -c {compression_level} {flag}"
            " -0 /dev/null -s /dev/null -n {input}"
            " > {out}"
        ).format(
            threads=samtools_threads,
            compression_level=compression_level,
            flag=samtools_flag,
            input=fifo_or_stdin,
            out=output_info["path"],
        )

    return cmd


def build_cmd(opts, strobealign_extra_args, compression_level=6, coverage=False, depth=False):
    parts = []

    # Resolve output paths
    mapped_info = resolve_output_paths(opts.mapped_fastq) if opts.mapped_fastq else None
    unmapped_info = resolve_output_paths(opts.unmapped_fastq) if opts.unmapped_fastq else None
    bam_path = opts.bam if hasattr(opts, "bam") else None

    # Create output directories
    for info in [mapped_info, unmapped_info]:
        if info is not None:
            if info["paired"]:
                create_directory(os.path.dirname(info["path_1"]))
                create_directory(os.path.dirname(info["path_2"]))
            else:
                create_directory(os.path.dirname(info["path"]))
    if bam_path:
        create_directory(os.path.dirname(os.path.abspath(bam_path)))

    # Create temporary directory
    create_directory(opts.temporary_directory)

    has_mapped = mapped_info is not None
    has_unmapped = unmapped_info is not None
    has_bam = bam_path is not None

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
    ]
    # samtools collate name-groups mate pairs for samtools fastq; skip it when
    # no FASTQ output is requested since samtools sort does not benefit from it.
    if has_mapped or has_unmapped:
        producer_tokens.append(
            "samtools collate -u -O - {prefix}".format(
                prefix=os.path.join(opts.temporary_directory, "samtools_collate_{}".format(run_id)),
            )
        )
    producer = " | ".join(producer_tokens)

    # Build consumer commands
    consumers = []
    if has_mapped:
        consumers.append(("mapped", build_consumer_cmd(
            "-F 12", opts.samtools_threads, None, mapped_info, opts.no_repair, compression_level,
        )))
    if has_unmapped:
        consumers.append(("unmapped", build_consumer_cmd(
            "-f 12", opts.samtools_threads, None, unmapped_info, opts.no_repair, compression_level,
        )))
    if has_bam:
        consumers.append(("bam", "samtools sort -@ {threads} -T {tmp} -o {out} -".format(
            threads=opts.samtools_threads,
            tmp=os.path.join(opts.temporary_directory, "samtools_sort_{}".format(run_id)),
            out=bam_path,
        )))

    if len(consumers) == 0:
        # Should not happen due to validation, but just in case
        parts.append(producer + " > /dev/null")

    elif len(consumers) == 1:
        # Single consumer: direct pipe, no fifos
        label, consumer_cmd = consumers[0]
        # Replace the placeholder input with stdin
        consumer_cmd = _set_consumer_input(consumer_cmd, label, "-")
        parts.append("{} | {}".format(producer, consumer_cmd))

    else:
        # Multiple consumers: use named pipes for all but the last
        # Last consumer reads from stdout via >
        fifos = []
        for i, (label, consumer_cmd) in enumerate(consumers[:-1]):
            fifo = os.path.join(opts.temporary_directory, "{}_{}_fifo".format(label, run_id))
            fifos.append(fifo)
            fifo_paths.append(fifo)
            consumer_cmd = _set_consumer_input(consumer_cmd, label, fifo)
            parts.append("")  # blank line for readability
            parts.append("{} &".format(consumer_cmd))

        # Last consumer reads from the final fifo via >
        last_label, last_consumer_cmd = consumers[-1]
        last_fifo = os.path.join(opts.temporary_directory, "{}_{}_fifo".format(last_label, run_id))
        fifos.append(last_fifo)
        fifo_paths.append(last_fifo)
        last_consumer_cmd = _set_consumer_input(last_consumer_cmd, last_label, last_fifo)
        parts.append("")
        parts.append("{} &".format(last_consumer_cmd))

        # mkfifo at the top
        parts.insert(0, "mkfifo {}".format(" ".join(fifos)))

        # Producer tees to all fifos except last, stdout goes to last
        tee_targets = fifos[:-1]
        parts.append("")
        parts.append("{} | tee {} > {}".format(producer, " ".join(tee_targets), fifos[-1]))
        parts.append("")
        parts.append("wait")
        parts.append("rm -f {}".format(" ".join(fifos)))

    if (coverage or depth) and has_bam:
        parts.append("")
        post_cmds = ["samtools index {}".format(bam_path)]
        if coverage:
            post_cmds.append("samtools coverage {bam} > {bam}.coverage.tsv".format(bam=bam_path))
        if depth:
            post_cmds.append("samtools depth -aa -H {bam} > {bam}.depth.tsv".format(bam=bam_path))
        parts.append(" \\\n&& ".join(post_cmds))

    return "\n".join(parts)


def _set_consumer_input(consumer_cmd, label, input_path):
    """Replace the placeholder input in a consumer command.
    
    For bam consumers, append the input path.
    For fastq consumers, the input was built with None as placeholder.
    """
    if label == "bam":
        # samtools sort command ends with '- ', replace trailing '-' with actual input
        if consumer_cmd.endswith(" -"):
            return consumer_cmd[:-1] + input_path
        return consumer_cmd
    else:
        # fastq consumer commands have None as the input placeholder
        return consumer_cmd.replace("None", input_path)


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
    usage = "{} [strobealign_options] reference reads1 reads2 --mapped_fastq PATH --unmapped_fastq PATH --bam PATH".format(__program__)
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
    parser_io.add_argument("--bam", type=str, default=None,
                           help="Output path for coordinate-sorted BAM file")
    parser_io.add_argument("--coverage", action="store_true", default=False,
                           help="Run samtools coverage on BAM and write {bam}.coverage.tsv (requires --bam)")
    parser_io.add_argument("--depth", action="store_true", default=False,
                           help="Run samtools depth -a on BAM and write {bam}.depth.tsv (requires --bam)")

    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("-t", "--threads", type=int, default=1, help="Threads for strobealign [Default: 1]")
    parser_utility.add_argument("--samtools_threads", type=int, default=1, help="Threads for samtools fastq/sort [Default: 1]")
    parser_utility.add_argument("--no_repair", action="store_true", default=False, help="Disable repair.sh post-processing")
    parser_utility.add_argument("-c", "--compression_level", type=int, default=6, help="Compression level [0..9] for samtools fastq bgzf output [Default: 6]")
    parser_utility.add_argument("-T", "--temporary_directory", type=str, default="/tmp/strobealign_wrapper",
                                help="Temporary directory for samtools collate/sort and named pipes [Default: /tmp/strobealign_wrapper]")
    parser_utility.add_argument("-v", "--version", action="version", version="{} v{}".format(__program__, __version__))

    # Parse
    opts, strobealign_extra_args = parser.parse_known_args(args)

    # Validate
    if opts.mapped_fastq is None and opts.unmapped_fastq is None and opts.bam is None:
        parser.error("At least one of --mapped_fastq, --unmapped_fastq, or --bam must be provided")

    if (opts.coverage or opts.depth) and opts.bam is None:
        parser.error("--coverage and --depth require --bam to be provided")

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
    if opts.bam is not None:
        opts.bam = format_path(opts.bam)

    assert os.path.exists(opts.reference), "Reference not found: {}".format(opts.reference)
    assert os.path.exists(opts.reads1), "Forward reads not found: {}".format(opts.reads1)
    assert os.path.exists(opts.reads2), "Reverse reads not found: {}".format(opts.reads2)

    # Create temporary directory
    create_directory(opts.temporary_directory)

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
    logger.info("BAM: {}", opts.bam)
    logger.info("Coverage: {}", opts.coverage)
    logger.info("Depth: {}", opts.depth)
    logger.info("Threads: {}", opts.threads)
    logger.info("Samtools threads: {}", opts.samtools_threads)
    logger.info("Repair: {}", not opts.no_repair)
    logger.info("Compression level: {}", opts.compression_level)
    logger.info("Temporary directory: {}", opts.temporary_directory)
    if strobealign_extra_args:
        logger.info("Strobealign extra args: {}", " ".join(strobealign_extra_args))
    logger.info("=" * 60)

    # Build and execute
    cmd_str = build_cmd(opts, strobealign_extra_args, compression_level=opts.compression_level,
                        coverage=opts.coverage, depth=opts.depth)
    logger.info("Command:\n{}", cmd_str)
    run(["bash", "-c", cmd_str], check=True)
    logger.success("Completed successfully")


if __name__ == "__main__":
    main()