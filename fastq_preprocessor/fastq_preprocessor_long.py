#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

pd.options.display.max_colwidth = 100

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.11.29"

# .............................................................................
# Primordial
# .............................................................................
# Fastp
def get_chopper_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
    # Command
    cmd = [
    # fastp
    "(",
    "gunzip -c" if input_filepaths[0].endswith(".gz") else "cat",
    input_filepaths[0],
    "|",
    os.environ["chopper"],
    "-l {}".format(opts.minimum_read_length),
    "-q {}".format(opts.minimum_quality_score),
    "--threads {}".format(opts.n_jobs),
    opts.chopper_options,
    "|",
    os.environ["pigz"],
    "-p {}".format(opts.n_jobs),
    ">",
    os.path.join(output_directory, "trimmed.fastq.gz"),
    ")",
    # Seqkit
    "&&",
    "(",
    os.environ["seqkit"],
    "stats",
    "-T",
    "-j {}".format(opts.n_jobs),
    input_filepaths[0],
    os.path.join(output_directory, "trimmed.fastq.gz"),
    ">",
    os.path.join(output_directory, "seqkit_stats.tsv"),
    ")",
    ]
    return cmd

# MiniMap2
def get_minimap2_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
    # Command
    cmd = [
    "(",
    os.environ["minimap2"],
    "-a",
    "-t {}".format(opts.n_jobs),
    "-x {}".format(opts.minimap2_preset),
    opts.minimap2_options,
    opts.contamination_index,
    input_filepaths[0],
    ">",
    os.path.join(directories["tmp"], "minimap2.sam"),

        "&&",

    "cat",
    os.path.join(directories["tmp"], "minimap2.sam"),
    "|",
    os.environ["samtools"],
    "fastq",
    "--threads {}".format(opts.n_jobs),
    "-f 4",
    "-",
    "|",
    os.environ["pigz"],
    "-p {}".format(opts.n_jobs),
    ">",
    os.path.join(output_directory, "contaminated.fastq.gz"),


        "&&",

    "cat",
    os.path.join(directories["tmp"], "minimap2.sam"),
    "|",
    os.environ["samtools"],
    "fastq",
    "--threads {}".format(opts.n_jobs),
    "-F 4",
    "-",
    "|",
    os.environ["pigz"],
    "-p {}".format(opts.n_jobs),
    ">",
    os.path.join(output_directory, "cleaned.fastq.gz"),
    ")",

        "&&",

    "rm -rf {}".format(os.path.join(directories["tmp"], "minimap2.sam")),

    # Seqkit
        "&&",

    os.environ["seqkit"],
    "stats",
    "-T",
    "-j {}".format(opts.n_jobs),
    os.path.join(output_directory, "*.fastq.gz"),

    ">",
    os.path.join(output_directory, "seqkit_stats.tsv"),
    ]

    # Remove trimmed reads 
    if not opts.retain_trimmed_reads:
        cmd += [
        "&&",
        "rm -rf {}".format(os.path.join( directories[("intermediate",  "1__chopper")], "trimmed.fastq.gz")),
        ]
    # Remove decontaminated reads
    if not opts.retain_contaminated_reads:
        cmd += [
        "&&",
        "rm -rf {}".format(os.path.join( output_directory, "contaminated.fastq.gz")),
        ]
            
    return cmd

# BBDuk
def get_bbduk_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
    cmd = [ 
            "(",
        os.environ["bbduk.sh"],
        "zl=1", # Most likely will delete these files
        "overwrite=t",
        "threads={}".format(opts.n_jobs),
        "in={}".format(input_filepaths[0]),
        "ref={}".format(opts.kmer_database),
        # "refstats={}".format(os.path.join(output_directory, "bbduk_refstats.txt")),
        # "stats={}".format(os.path.join(output_directory, "bbduk_stats.txt")),
        "k={}".format(opts.kmer_size),
        "minlen={}".format(opts.minimum_read_length),
        "out={}".format(os.path.join(output_directory, "non-kmer_hits.fastq.gz")),
        "outm={}".format(os.path.join(output_directory, "kmer_hits.fastq.gz")),
        opts.bbduk_options,
        ")",

        # Seqkit
        "&&",
        "(",
        os.environ["seqkit"],
        "stats",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "*.fastq.gz"),
        ">",
        os.path.join(output_directory, "seqkit_stats.tsv"),
        ")",
        ]


    # Remove bbduk kmer hits
    if not opts.retain_kmer_hits:
        cmd += [
        "&&",
        "rm -rf {}".format(os.path.join( output_directory, "kmer_hits.fastq.gz")),
        ]
    # Remove bbduk non-kmer hits
    if not opts.retain_non_kmer_hits:
        cmd += [
        "&&",
        "rm -rf {}".format(os.path.join( output_directory, "non-kmer_hits.fastq.gz")),
        ]

            
    return cmd


# Symlink
def get_synopsis_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
"""
python -c "import glob, pandas as pd; pd.concat(map(lambda fp: pd.read_csv(fp, sep='\t', index_col=0), glob.glob('{}')), axis=0).to_csv('{}', sep='\t')"
""".format(
            os.path.join(directories["intermediate"],"*/seqkit_stats.tsv"),
            os.path.join(output_directory,"seqkit_stats.concatenated.tsv"),
            ),
        ]

    cmd += [
    "DST={}; (for SRC in {}; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        output_directory,
        " ".join(input_filepaths), 
        )
    ]
        

    return cmd

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = set([])

    required_executables={
                # "repair.sh",
                "pigz",
                "samtools",
                "bbduk.sh",
                "minimap2",
                "chopper",
                "seqkit",
    } | accessory_scripts

    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in required_executables:
            executables[name] = os.path.join(os.environ["CONDA_PREFIX"], "bin", name)
    else:
        opts.path_config = format_path(opts.path_config)
        assert os.path.exists(opts.path_config), "config file does not exist.  Have you created one in the following directory?\n{}\nIf not, either create one, check this filepath:{}, or give the path to a proper config file using --path_config".format(opts.script_directory, opts.path_config)
        assert os.stat(opts.path_config).st_size > 1, "config file seems to be empty.  Please add 'name' and 'executable' columns for the following program names: {}".format(required_executables)
        df_config = pd.read_csv(opts.path_config, sep="\t")
        assert {"name", "executable"} <= set(df_config.columns), "config must have `name` and `executable` columns.  Please adjust file: {}".format(opts.path_config)
        df_config = df_config.loc[:,["name", "executable"]].dropna(how="any", axis=0).applymap(str)
        # Get executable paths
        executables = OrderedDict(zip(df_config["name"], df_config["executable"]))
        assert required_executables <= set(list(executables.keys())), "config must have the required executables for this run.  Please adjust file: {}\nIn particular, add info for the following: {}".format(opts.path_config, required_executables - set(list(executables.keys())))

    # Display
    for name in sorted(accessory_scripts):
        executables[name] = "python " + os.path.join(opts.script_directory, "scripts", name)
    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)

# Pipeline
def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name=__program__, description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    # =========
    # Fastp
    # =========
    step = 1
    program = "chopper"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Quality trimming"
    # i/o
    input_filepaths = [opts.reads]

    if all([
        bool(opts.contamination_index),
        not bool(opts.retain_trimmed_reads),
    ]):
        output_filenames = ["seqkit_stats.tsv"]
    else:
        output_filenames = ["trimmed.fastq.gz", "seqkit_stats.tsv"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    # Parameters
    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }
    # Command
    cmd = get_chopper_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )

    # =========
    # Bowtie2
    # =========
    if opts.contamination_index:
        step += 1
        program = "minimap2"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "Decontaminate reads based on a reference"
        # i/o
        input_filepaths = [
            os.path.join(os.path.join(directories["intermediate"], "1__chopper"), "trimmed.fastq.gz"),

        ]

        output_filenames = ["cleaned.fastq.gz", "seqkit_stats.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        # Parameters
        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }
        # Command
        cmd = get_minimap2_cmd(**params)
        pipeline.add_step(
                    id=program,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=False,
                    validate_outputs=True,
        )

    # =========
    # BBDuk
    # =========
    if opts.kmer_database is not None:
        step += 1
        program = "bbduk"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "Decontaminate reads based on k-mer database"

        if opts.contamination_index:
            # i/o
            input_filepaths = [
                os.path.join(os.path.join(directories["intermediate"], "2__minimap2"), "cleaned.fastq.gz"),
            ]

        else:
            # i/o
            input_filepaths = [
                os.path.join(os.path.join(directories["intermediate"], "1__chopper"), "trimmed.fastq.gz"),

            ]

        output_filenames = ["seqkit_stats.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        # Parameters
        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }
        # Command
        cmd = get_bbduk_cmd(**params)
        pipeline.add_step(
                    id=program,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=False,
                    validate_outputs=True,
        )


   

    # =============
    # Symlink
    # =============
    step += 1
    program = "synopsis"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories["output"] 
    description = "Symlinking and merging relevant output files"

    # i/o
    input_filepaths = [
        os.path.join(directories["intermediate"], "*", "*.fastq.gz"),
    ]
    
  
    output_filenames =  map(lambda fp: fp.split("/")[-1], input_filepaths)
    output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))
    output_filepaths += [ 
        os.path.join(output_directory, "seqkit_stats.concatenated.tsv"),
    ]
    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_synopsis_cmd(**params)
    pipeline.add_step(
            id=program,
            description = description,
            step=step,
            cmd=cmd,
            input_filepaths = input_filepaths,
            output_filepaths = output_filepaths,
            validate_inputs=True,
            validate_outputs=True,
    )

    return pipeline

# Configure parameters
def configure_parameters(opts, directories):

    assert_acceptable_arguments(opts.retain_trimmed_reads, {0,1})
    assert_acceptable_arguments(opts.retain_contaminated_reads, {0,1})
    assert_acceptable_arguments(opts.minimap2_preset, {"map-pb", "map-ont", "map-hifi"})
     
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <reads.fq[.gz]> -n <name> -o <output_directory> |Optional| -x <reference_index> -k <kmer_database>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--reads", type=str, help = "path/to/reads.fastq[.gz]", required=True)
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="preprocessed", help = "path/to/project_directory [Default: preprocessed]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv. Must have at least 2 columns [name, executable] [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=int, help = "Restart from a particular checkpoint")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Chopper
    parser_chopper = parser.add_argument_group('Chopper arguments')
    parser_chopper.add_argument("-m", "--minimum_read_length", type=int, default=500, help="Chopper | Minimum read length [Default: 500]")
    parser_chopper.add_argument("-q", "--minimum_quality_score", type=int, default=10, help="Chopper | Minimum quality score [Default: 10]")
    parser_chopper.add_argument("--chopper_options", type=str, default="", help="Chopper | More options (e.g. --arg 1 ) https://github.com/wdecoster/chopper [Default: '']")

    # MiniMap2
    parser_minimap2 = parser.add_argument_group('MiniMap2 arguments')
    parser_minimap2.add_argument("-x", "--contamination_index", type=str, help="MiniMap2 | path/to/contamination_index\n(e.g., Human T2T assembly from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz)")
    parser_minimap2.add_argument("--minimap2_preset", type=str, default="map-ont", help="MiniMap2 | MiniMap2 preset {map-pb, map-ont, map-hifi} [Default: map-ont]")
    parser_minimap2.add_argument("--retain_trimmed_reads", default=0, type=int, help = "Retain Chopper trimmed fastq after decontamination. 0=No, 1=yes [Default: 0]") 
    parser_minimap2.add_argument("--retain_contaminated_reads", default=0, type=int, help = "Retain contaminated fastq after decontamination. 0=No, 1=yes [Default: 0]")
    parser_minimap2.add_argument("--minimap2_options", type=str, default="", help="MiniMap2 | More options (e.g. --arg 1 ) [Default: '']\nhttp://bowtie-bio.sourceforge.net/bowtie2/manual.shtml")

    # BBDuk
    parser_bbduk = parser.add_argument_group('BBDuk arguments')
    parser_bbduk.add_argument("-k","--kmer_database", type=str,  help="BBDuk | path/to/kmer_database\n(e.g., ribokmers.fa.gz from https://figshare.com/ndownloader/files/36220587)")
    parser_bbduk.add_argument("--kmer_size", type=int, default=31, help="BBDuk | k-mer size [Default: 31]")
    parser_bbduk.add_argument("--retain_kmer_hits", default=0, type=int, help = "Retain reads that map to k-mer database. 0=No, 1=yes [Default: 0]")
    parser_bbduk.add_argument("--retain_non_kmer_hits", default=0, type=int, help = "Retain reads that do not map to k-mer database. 0=No, 1=yes [Default: 0]")
    parser_bbduk.add_argument("--bbduk_options", type=str, default="", help="BBDuk | More options (e.g., --arg 1) [Default: '']")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Threads
    if opts.n_jobs == -1:
        from multiprocessing import cpu_count 
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1.  To select all available threads, use -1."

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
    directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["sample"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))

    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print(format_header("Name: {}".format(opts.name), "."), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("GenoPype version:", genopype_version, file=sys.stdout)
    print("Script version:", __version__, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
    sys.stdout.flush()


    # Run pipeline
    with open(os.path.join(directories["sample"], "commands.sh"), "w") as f_cmds:
        pipeline = create_pipeline(
                     opts=opts,
                     directories=directories,
                     f_cmds=f_cmds,
        )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)

if __name__ == "__main__":
    main()
