
import argparse

from phenox.phenox import PhenoX


def run(args):
    print('Input query: %s' % args.query_str)
    phenox = PhenoX(args.query_str)
    phenox.subtype()


parser = argparse.ArgumentParser()

parser.add_argument('query_str')
parser.set_defaults(func=run)

if __name__ == '__main__':
    args = parser.parse_args()
    args.func(args)