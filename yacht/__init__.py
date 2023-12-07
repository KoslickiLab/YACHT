import argparse

__version__ = '1.1.0'

from . import run_YACHT
from . import make_training_data_from_sketches
from . import standardize_yacht_output


def main():
    parser = argparse.ArgumentParser(prog='yacht')
    subparsers = parser.add_subparsers(dest='command')

    # Train command
    train_parser = subparsers.add_parser('train', description='Pre-process the reference genomes')
    make_training_data_from_sketches.add_arguments(train_parser)
    train_parser.set_defaults(func=make_training_data_from_sketches.main)

    # Run command
    run_parser = subparsers.add_parser('run', description='Run the YACHT algorithm')
    run_YACHT.add_arguments(run_parser)
    run_parser.set_defaults(func=run_YACHT.main)

    # Convert command
    convert_parser = subparsers.add_parser('convert', description='Convert YACHT result to other popular output formats')
    standardize_yacht_output.add_arguments(convert_parser)
    convert_parser.set_defaults(func=standardize_yacht_output.main)

    args = parser.parse_args()
    if 'func' in args:
        args.func(args)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()