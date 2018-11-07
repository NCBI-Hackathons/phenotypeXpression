#!/usr/bin/env python

import argparse
from phenox.phenox import PhenoX


__version__ = '0.3.1'


def run(args):
    print('Input query: %s' % args.query_str)
    phenox = PhenoX(args.email, args.query_str, args.outprefix)
    phenox.subtype()


# argparse function for options
def getargs():
    parser = argparse.ArgumentParser(
        prog='run_phenox.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='Subclassification of disease states based on the '
                    'intersection of literature and expression'
    )

    parser.add_argument(
        '-o', metavar='prefix', dest='outprefix',
        default='phenox',
        help='choose an alternate prefix for outfiles'
    )

    parser.add_argument(
        "--version", action='version',
        version='\n'.join(['PhenoX v' + __version__])
    )

    parser.add_argument(
        "query_str", help="disease query term"
    )

    required_named = parser.add_argument_group('required arguments')

    required_named.add_argument(
        '-e', dest='email', required=True,
        help='NCBI requires an email for database queries'
    )

    return parser.parse_args()


if __name__ == '__main__':
    args = getargs()
    run(args)
