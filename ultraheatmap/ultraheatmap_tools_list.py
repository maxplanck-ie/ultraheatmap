import argparse
import sys
from ultraheatmap._version import __version__


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
                    [Tools]
                    computeOrderedMatrix
                    AddFeatureToMatrix
                    """)
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):
    if args is None and len(sys.argv) == 1:
        args = ["--help"]

    args = parse_arguments().parse_args(args)
    return args
