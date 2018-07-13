#!/usr/bin/env python

import argparse
from phenox.phenox import PhenoX
__version__ = '0.1.0'


def run(args):
    print('Input query: %s' % args.query_str)
    phenox = PhenoX(args.email, args.query_str)
    phenox.subtype()

# argparse function for options can include or delete
# def getargs():
#     parser = argparse \
#             .ArgumentParser(prog='run_phenox.py',
#                             formatter_class=argparse.RawTextHelpFormatter,
#                             description='Subclassification of disease states '
#                                         'based on the intersection of literat'
#                                         'ure and expression')
#     parser.add_argument('-o', metavar='prefix', dest='outprefix',
#                         default='phenox',
#                         help='choose an alternate prefix for outfiles')
#     parser.add_argument("--version", action='version',
#                         version='\n'.join(['PhenoX v' + __version__]))
#     parser.add_argument("query", help="disease query term")
#     requiredNamed = parser.add_argument_group('required arguments')
#     requiredNamed \
#             .add_argument('-e', dest='email', required=True,
#                           help='NCBI requires an email for database queries')
#     return parser.parse_args()

parser = argparse.ArgumentParser()

parser.add_argument('email')
parser.add_argument('query_str')
parser.set_defaults(func=run)

if __name__ == '__main__':
    # args = getargs() # if using getargs function
    args = parser.parse_args()
    args.func(args)