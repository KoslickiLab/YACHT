import argparse
import sys

from . import (
    download_default_ref_db,
    download_demofiles,
    download_pretrained_ref_db,
    make_training_data_from_sketches,
    run_YACHT,
    sketch_ref_genomes,
    sketch_sample,
    standardize_yacht_output,
)
from .utils import __version__

# Custom help action
class CustomHelpAction(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super(CustomHelpAction, self).__init__(option_strings=option_strings, dest=dest, nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print_custom_help_exit()

def print_version():
    print(f"== This is YACHT version {utils.__version__} ==")
    print("== Please cite https://doi.org/10.1093/bioinformatics/btae047. ==")

def print_custom_help_exit():
    print_version()
    print("""
YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms 
in a metagenomic sample, based on average nucleotide identity (ANI).

Usage instructions:
    yacht download demo -h download YACHT demo files
    yacht download default_ref_db -h download default raw reference databases
    yacht download pretrained_ref_db -h download pretrained databases

    yacht sketch ref -h sketch reference genomes
    yacht sketch sample -h sketch metagenomics samples
        
    yacht train -h pre-process the reference genomes
    
    yacht run -h run the YACHT algorithm
    
    yacht convert -h convert YACHT result to other popular output formats

Options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit
""")
    sys.exit()

def main():
    parser = argparse.ArgumentParser(prog="yacht", add_help=False)
    parser.add_argument("-v", "--version", action="store_true", help="show program's version number and exit")
    parser.add_argument("-h", "--help", action=CustomHelpAction, help="show this help message and exit")
    subparsers = parser.add_subparsers(dest="command")

    # Train command
    train_parser = subparsers.add_parser(
        "train", description="Pre-process the reference genomes"
    )
    make_training_data_from_sketches.add_arguments(train_parser)
    train_parser.set_defaults(func=make_training_data_from_sketches.main)

    # Run command
    run_parser = subparsers.add_parser("run", description="Run the YACHT algorithm")
    run_YACHT.add_arguments(run_parser)
    run_parser.set_defaults(func=run_YACHT.main)

    # Convert command
    convert_parser = subparsers.add_parser(
        "convert", description="Convert YACHT result to other popular output formats"
    )
    standardize_yacht_output.add_arguments(convert_parser)
    convert_parser.set_defaults(func=standardize_yacht_output.main)

    # Download command with submodules
    download_parser = subparsers.add_parser(
        "download",
        description="Download YACHT demo files or default raw reference databases or pretrained databases",
    )
    download_subparsers = download_parser.add_subparsers(dest="download_subcommand")

    # Download demo files
    demo_parser = download_subparsers.add_parser(
        "demo", description="Download YACHT demo files"
    )
    download_demofiles.add_arguments(demo_parser)
    demo_parser.set_defaults(func=download_demofiles.main)

    # Download default raw reference databases
    default_ref_db_parser = download_subparsers.add_parser(
        "default_ref_db", description="Download default raw reference databases"
    )
    download_default_ref_db.add_arguments(default_ref_db_parser)
    default_ref_db_parser.set_defaults(func=download_default_ref_db.main)

    # Download pretrained databases
    pretrained_ref_db_parser = download_subparsers.add_parser(
        "pretrained_ref_db", description="Download pretrained databases"
    )
    download_pretrained_ref_db.add_arguments(pretrained_ref_db_parser)
    pretrained_ref_db_parser.set_defaults(func=download_pretrained_ref_db.main)

    # Sketch command with submodules
    sketch_parser = subparsers.add_parser(
        "sketch", description="Sketch reference genomes or metagenomics samples"
    )
    sketch_subparsers = sketch_parser.add_subparsers(dest="sketch_subcommand")

    # Sketch reference genomes
    sketch_ref_parser = sketch_subparsers.add_parser(
        "ref", description="Sketch reference genomes"
    )
    sketch_ref_genomes.add_arguments(sketch_ref_parser)
    sketch_ref_parser.set_defaults(func=sketch_ref_genomes.main)

    # Sketch metagenomics samples
    sketch_sample_parser = sketch_subparsers.add_parser(
        "sample", description="Sketch metagenomics samples"
    )
    sketch_sample.add_arguments(sketch_sample_parser)
    sketch_sample_parser.set_defaults(func=sketch_sample.main)

    args = parser.parse_args()
    
    if args.version:
        print_version()
        sys.exit()
    
    if "func" in args:
        args.func(args)
    else:
        print_custom_help_exit()


if __name__ == "__main__":
    main()
