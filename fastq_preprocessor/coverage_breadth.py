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
                           help='Path to samtools depth TSV (.tsv or .tsv.gz).\nIf generated with -aa, all positions are present (--fasta not needed).\nIf generated without -aa, provide --fasta for reference lengths.')
    parser_io.add_argument('--fasta', type=str, default=None,
                           help='Reference FASTA (required when samtools depth was run without -aa)')
    parser_io.add_argument('--coverage', type=str, default=None,
                           help='Path to samtools coverage TSV (.tsv or .tsv.gz)')
    parser_io.add_argument('--contigs_to_genomes', type=str, default=None,
                           help='Optional 2-column TSV (id_contig<TAB>id_genome, no header).\nWhen provided, breadth is computed per genome instead of globally.')
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

    input_filepath = opts.depth if opts.depth is not None else opts.coverage

    # Build per-contig (covered, length) Series
    if opts.coverage is not None:
        # Mode A: samtools coverage TSV
        df = pd.read_csv(opts.coverage, sep='\t')
        df.columns = [col.lstrip('#') for col in df.columns]
        df = df.set_index('rname')
        contig_length = (df['endpos'] - df['startpos'] + 1).astype(int)
        contig_covered = df['covbases'].astype(int)
    elif opts.fasta is not None:
        # Mode C: samtools depth without -aa; use FASTA lengths for per-contig length
        df = pd.read_csv(opts.depth, sep='\t', comment='#', header=None,
                         names=['id', 'position', 'depth'])
        fasta_lengths = parse_fasta_lengths(opts.fasta)
        contig_length = pd.Series(fasta_lengths, dtype=int)
        contig_covered = (df[df['depth'] >= opts.minimum_depth]
                          .groupby('id').size().astype(int)
                          .reindex(contig_length.index, fill_value=0))
    else:
        # Mode B: samtools depth with -aa; all positions present in file
        df = pd.read_csv(opts.depth, sep='\t', comment='#', header=None,
                         names=['id', 'position', 'depth'])
        contig_length = df.groupby('id').size().astype(int)
        contig_covered = (df[df['depth'] >= opts.minimum_depth]
                          .groupby('id').size().astype(int)
                          .reindex(contig_length.index, fill_value=0))

    if opts.contigs_to_genomes is not None:
        mapping = pd.read_csv(opts.contigs_to_genomes, sep='\t', header=None,
                              names=['id_contig', 'id_genome'], dtype=str)

        contig_length.index = contig_length.index.astype(str)
        contig_covered.index = contig_covered.index.astype(str)
        data_contigs = set(contig_length.index)
        mapping_contigs = set(mapping['id_contig'])

        missing_from_data = mapping_contigs - data_contigs
        if missing_from_data:
            sample = ', '.join(sorted(missing_from_data)[:5])
            parser.error(
                '{} contig(s) in --contigs_to_genomes are not present in the '
                'coverage/depth data. Examples: {}'.format(len(missing_from_data), sample)
            )

        unmapped = data_contigs - mapping_contigs
        if unmapped:
            sample = ', '.join(sorted(unmapped)[:5])
            sys.stderr.write(
                'Warning: {} contig(s) present in the data are not in '
                '--contigs_to_genomes and will be skipped. Examples: {}\n'.format(
                    len(unmapped), sample)
            )

        contig_to_genome = mapping.set_index('id_contig')['id_genome']
        aligned_length = contig_length.loc[contig_to_genome.index]
        aligned_covered = contig_covered.loc[contig_to_genome.index]
        genome_length = aligned_length.groupby(contig_to_genome.values).sum()
        genome_covered = aligned_covered.groupby(contig_to_genome.values).sum()

        lines = ['id_genome\tpercentage_of_bases_covered\n']
        for genome_id in genome_length.index:
            total_length = int(genome_length[genome_id])
            if total_length == 0:
                sys.stderr.write(
                    "Warning: genome '{}' has total length 0; skipping.\n".format(genome_id)
                )
                continue
            total_covered = int(genome_covered[genome_id])
            percentage = total_covered / total_length * 100
            lines.append('{}\t{}\n'.format(genome_id, percentage))
        output_text = ''.join(lines)
    else:
        total_length = int(contig_length.sum())
        total_covered = int(contig_covered.sum())
        percentage = total_covered / total_length * 100
        output_text = 'filepath\tpercentage_of_bases_covered\n{}\t{}\n'.format(
            input_filepath, percentage
        )

    if opts.output == '-':
        sys.stdout.write(output_text)
    else:
        with open(opts.output, 'w') as f:
            f.write(output_text)


if __name__ == '__main__':
    main()
