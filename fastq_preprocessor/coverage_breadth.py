#!/usr/bin/env python

from __future__ import print_function, division
import sys, os, argparse, gzip

import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
from fastq_preprocessor import __version__


def parse_fasta_lengths(fasta_path):
    opener = gzip.open if fasta_path.endswith('.gz') else open
    lengths = {}
    name, length = None, 0
    with opener(fasta_path, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        lengths[name] = length
    return lengths


def main(args=None):
    description = """
    Running: {} v{} via Python v{} | {}""".format(
        __program__, __version__, sys.version.split(" ")[0], sys.executable,
    )
    usage = "{} (--depth PATH | --coverage PATH) [--fasta PATH] [-o PATH] [--minimum_depth INT]".format(__program__)

    parser = argparse.ArgumentParser(
        description=description,
        usage=usage,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument('--depth', type=str, default=None,
                           help='Path to samtools depth TSV.\nIf generated with -aa, all positions are present (--fasta not needed).\nIf generated without -aa, provide --fasta for reference lengths.')
    parser_io.add_argument('--fasta', type=str, default=None,
                           help='Reference FASTA (required when samtools depth was run without -aa)')
    parser_io.add_argument('--coverage', type=str, default=None,
                           help='Path to samtools coverage TSV')
    parser_io.add_argument('-o', '--output', type=str, default='-',
                           help='Output path [Default: stdout]')

    parser_params = parser.add_argument_group('Parameter arguments')
    parser_params.add_argument('--minimum_depth', type=int, default=1,
                               help='Minimum depth to count a position as covered [Default: 1]')

    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument('-v', '--version', action='version',
                                version='{} v{}'.format(__program__, __version__))

    opts = parser.parse_args(args)

    # Validate
    if opts.depth is None and opts.coverage is None:
        parser.error('Must provide --depth or --coverage')
    if opts.depth is not None and opts.coverage is not None:
        parser.error('Cannot provide both --depth and --coverage')
    if opts.coverage is not None and opts.minimum_depth > 1:
        parser.error('--minimum_depth > 1 is not supported with --coverage; use --depth with -aa instead')

    # Compute
    if opts.coverage is not None:
        # Mode A: samtools coverage TSV
        df = pd.read_csv(opts.coverage, sep='\t')
        df.columns = [col.lstrip('#') for col in df.columns]
        total_length = int((df['endpos'] - df['startpos'] + 1).sum())
        total_covered = int(df['covbases'].sum())
    elif opts.fasta is not None:
        # Mode C: samtools depth without -aa; use FASTA lengths for total positions
        df = pd.read_csv(opts.depth, sep='\t', comment='#', header=None,
                         names=['id', 'position', 'depth'])
        fasta_lengths = parse_fasta_lengths(opts.fasta)
        total_length = sum(fasta_lengths.values())
        total_covered = int((df['depth'] >= opts.minimum_depth).sum())
    else:
        # Mode B: samtools depth with -aa; all positions present in file
        df = pd.read_csv(opts.depth, sep='\t', comment='#', header=None,
                         names=['id', 'position', 'depth'])
        total_length = len(df)
        total_covered = int((df['depth'] >= opts.minimum_depth).sum())

    percentage = total_covered / total_length * 100
    output_line = 'percentage_of_bases_covered\t{}\n'.format(percentage)

    if opts.output == '-':
        sys.stdout.write(output_line)
    else:
        with open(opts.output, 'w') as f:
            f.write(output_line)


if __name__ == '__main__':
    main()
