import os
import argparse
import importlib
import sys

# YACHT does no multithreaded linear algebra: keep the BLAS libraries used by numpy/scipy from starting one
# thread per CPU core when they are imported (YACHT is typically run on many samples in parallel).
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

_SUBMODULES = (
    "download_default_ref_db",
    "download_demofiles",
    "download_pretrained_ref_db",
    "make_training_data_from_sketches",
    "run_YACHT",
    "sketch_ref_genomes",
    "sketch_sample",
    "standardize_yacht_output",
    "hypothesis_recovery_src",
    "utils",
    "yacht_index",
)


def __getattr__(name):
    # submodules are imported on first use, so that e.g. `yacht run` does not pay for the imports of `yacht convert`
    if name in _SUBMODULES:
        return importlib.import_module(f".{name}", __name__)
    if name == "__version__":
        return importlib.import_module(".utils", __name__).__version__
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

# Custom help action
class CustomHelpAction(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super(CustomHelpAction, self).__init__(option_strings=option_strings, dest=dest, nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print_custom_help_exit()

def print_version():
    from .utils import __version__
    print(f"== This is YACHT version {__version__} ==")
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

def _requested_command(argv):
    """The (command, subcommand) named on the command line, e.g. ("download", "demo"), or None entries."""
    words = [a for a in argv if not a.startswith("-")]
    command = words[0] if words else None
    subcommand = words[1] if command in ("download", "sketch") and len(words) > 1 else None
    return command, subcommand


def main():
    command, subcommand = _requested_command(sys.argv[1:])

    def add_command(subparsers, name, description, module_name, selected):
        sub = subparsers.add_parser(name, description=description)
        if selected:
            # only the module of the requested command is imported
            module = importlib.import_module(f".{module_name}", __name__)
            module.add_arguments(sub)
            sub.set_defaults(func=module.main)
        return sub

    parser = argparse.ArgumentParser(prog="yacht", add_help=False)
    parser.add_argument("-v", "--version", action="store_true", help="show program's version number and exit")
    parser.add_argument("-h", "--help", action=CustomHelpAction, help="show this help message and exit")
    subparsers = parser.add_subparsers(dest="command")

    # Train command
    add_command(subparsers, "train", "Pre-process the reference genomes", "make_training_data_from_sketches", command == "train")

    # Run command
    add_command(subparsers, "run", "Run the YACHT algorithm", "run_YACHT", command == "run")

    # Convert command
    add_command(
        subparsers, "convert", "Convert YACHT result to other popular output formats", "standardize_yacht_output", command == "convert"
    )

    # Download command with submodules
    download_parser = subparsers.add_parser(
        "download",
        description="Download YACHT demo files or default raw reference databases or pretrained databases",
    )
    download_subparsers = download_parser.add_subparsers(dest="download_subcommand")
    downloading = command == "download"

    # Download demo files
    add_command(download_subparsers, "demo", "Download YACHT demo files", "download_demofiles", downloading and subcommand == "demo")

    # Download default raw reference databases
    add_command(
        download_subparsers, "default_ref_db", "Download default raw reference databases", "download_default_ref_db",
        downloading and subcommand == "default_ref_db",
    )

    # Download pretrained databases
    add_command(
        download_subparsers, "pretrained_ref_db", "Download pretrained databases", "download_pretrained_ref_db",
        downloading and subcommand == "pretrained_ref_db",
    )

    # Sketch command with submodules
    sketch_parser = subparsers.add_parser(
        "sketch", description="Sketch reference genomes or metagenomics samples"
    )
    sketch_subparsers = sketch_parser.add_subparsers(dest="sketch_subcommand")
    sketching = command == "sketch"

    # Sketch reference genomes
    add_command(sketch_subparsers, "ref", "Sketch reference genomes", "sketch_ref_genomes", sketching and subcommand == "ref")

    # Sketch metagenomics samples
    add_command(sketch_subparsers, "sample", "Sketch metagenomics samples", "sketch_sample", sketching and subcommand == "sample")

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
