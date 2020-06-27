import argparse
import collections
import re
import sys
from typing import List, Optional, Sequence, Tuple

import pepfrag
import unimoddb


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
            split_val = val.split('=')
            if len(split_val) == 2:
                d[field].append((split_val[0], float(split_val[1])))
            elif len(split_val) == 1:
                d[field].append(val)
            else:
                raise ValueError(f'Invalid neutral loss format: {val}')


def parse_args(args: List[str]) -> argparse.Namespace:
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
        '--cterm',
        help='C-terminal modification'
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
    return parser.parse_args(args)


MOD_PATTERN = re.compile(r'\[(\w+)\]')


def parse_sequence(
        seq_mods: str,
        ptmdb: unimoddb.UnimodDB
) -> Tuple[str, List[pepfrag.ModSite]]:
    """
    Parses the peptide sequence to extract embedded modifications.

    """
    mods: List[pepfrag.ModSite] = []
    seq = ''
    offset = 0
    for mod in MOD_PATTERN.finditer(seq_mods):
        mod_name = mod.group(1)
        start = mod.start()
        site = 'nterm' if offset == 0 else len(seq) + start - offset
        mods.append(pepfrag.ModSite(ptmdb.get_mass(mod_name), site, mod_name))
        seq += seq_mods[offset:start]
        offset = mod.end()

    seq += seq_mods[offset:]

    return seq, mods


def construct_peptide(
        seq_mods: str,
        cterm_mod: Optional[str],
        charge: int,
        radical: bool,
        ptmdb: unimoddb.UnimodDB
) -> pepfrag.Peptide:
    """
    Constructs a `Peptide` instance from the passed command line arguments.

    Returns:
        Peptide.

    """
    seq, mods = parse_sequence(seq_mods, ptmdb)

    if cterm_mod:
        mods.append(
            pepfrag.ModSite(ptmdb.get_mass(cterm_mod), 'cterm', cterm_mod)
        )

    return pepfrag.Peptide(seq, charge, mods, radical=radical)


def write_fragments(fragments: Sequence):
    print('Ion: m/z')
    for fragment in fragments:
        print(f'{fragment[1]}: {fragment[0]:.6f}')


def main(args):
    args = parse_args(args)

    ptmdb = unimoddb.UnimodDB()

    peptide = construct_peptide(
        args.sequence,
        args.cterm,
        args.charge,
        args.radical,
        ptmdb
    )

    try:
        fragments = (
            peptide.fragment(ion_types=args.ion_types) if args.ion_types
            else peptide.fragment()
        )
    except KeyError as ex:
        # TODO: update pepfrag to use a more specific error - e.g. UnknownNeutralLossError
        print(ex)
        sys.exit(1)

    write_fragments(fragments)


if __name__ == '__main__':
    main(sys.argv[1:])
