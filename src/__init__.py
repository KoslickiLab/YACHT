import argparse

from . import run_YACHT
from . import make_training_data_from_sketches
from . import standardize_yacht_output


def main():
    parser = argparse.ArgumentParser(prog='yacht')
    subparsers = parser.add_subparsers(dest='command')

    # Run command
    run_parser = subparsers.add_parser('run')
    run_YACHT.add_arguments(run_parser)
    run_parser.set_defaults(func=run_YACHT.main)

    # Train command
    train_parser = subparsers.add_parser('train')
    make_training_data_from_sketches.add_arguments(train_parser)
    train_parser.set_defaults(func=make_training_data_from_sketches.main)

    # Convert command
    convert_parser = subparsers.add_parser('convert')
    standardize_yacht_output.add_arguments(convert_parser)
    convert_parser.set_defaults(func=standardize_yacht_output.main)

    args = parser.parse_args()
    if 'func' in args:
        args.func(args)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()