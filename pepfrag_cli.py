import argparse
import collections
import re
import sys
from typing import Tuple

import pepfrag


class StoreIonTypeDict(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        self._nargs = nargs
        super().__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if option_string is None:
            raise ValueError('Expected option_string')
        d = getattr(namespace, self.dest)
        if d is None:
            d = collections.defaultdict(list)
            setattr(namespace, self.dest, d)

        field = pepfrag.IonType[option_string[2:]]

        if not values:
            d[field] = []
            return

        for val in values:
            split_val = val.split(',')
            if len(split_val) == 2:
                d[field].append(tuple(split_val))
            elif len(split_val) == 1:
                d[field].append(val)
            # TODO: invalid format


def parse_args():
    """
    Parses the command line arguments.

    """
    parser = argparse.ArgumentParser(description="Peptide fragmentation.")
    parser.add_argument(
        'sequence',
        help='Peptide sequence, including modifications'
    )
    parser.add_argument(
        'charge',
        type=int,
        help='Charge state of the peptide'
    )
    parser.add_argument(
        '-r',
        '--radical',
        action='store_true'
    )
    for ion_type in pepfrag.IonType:
        parser.add_argument(
            f'--{ion_type.name}',
            nargs='*',
            dest='ion_types',
            action=StoreIonTypeDict
        )
    return parser.parse_args()


MOD_PATTERN = re.compile(r'\[(\w+)\]')


def parse_sequence(seq: str) -> Tuple[str, Tuple[pepfrag.ModSite, ...]]:
    print(MOD_PATTERN.findall(seq))
    return seq, tuple()


def construct_peptide(args: argparse.Namespace) -> pepfrag.Peptide:
    """
    Constructs a `Peptide` instance from the passed command line arguments.

    Returns:
        Peptide.

    """
    seq, mods = parse_sequence(args.sequence)
    return pepfrag.Peptide(seq, args.charge, mods, radical=args.radical)


def main():
    args = parse_args()

    print(args)

    peptide = construct_peptide(args)

    try:
        fragments = (
            peptide.fragment(ion_types=args.ion_types) if args.ion_types
            else peptide.fragment()
        )
    except KeyError as ex:
        # TODO: update pepfrag to use a more specific error - e.g. UnknownNeutralLossError
        print(ex)
        sys.exit(1)

    print(fragments)


if __name__ == '__main__':
    main()
