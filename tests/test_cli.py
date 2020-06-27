import unittest.mock

from pepfrag import IonType, ModSite, Peptide
import unimoddb

from pepfrag_cli import construct_peptide, parse_args, parse_sequence


PTMDB = unimoddb.UnimodDB()


class TestIonTypeArgParser(unittest.TestCase):
    def test_no_ion_types(self):
        args = parse_args(['ACPLK', '2'])
        self.assertEqual(None, args.ion_types)

    def test_ion_types(self):
        args = parse_args([
            'ACPLK', '2',
            '--b', 'NH3',
            '--y', 'TEST=12.01', 'CO',
            '--precursor'
        ])
        self.assertEqual(
            {
                IonType.b: ['NH3'],
                IonType.y: [('TEST', 12.01), 'CO'],
                IonType.precursor: []
            },
            args.ion_types
        )


class TestParseSequence(unittest.TestCase):
    def test_no_modifications(self):
        self.assertEqual(('ACPLK', []), parse_sequence('ACPLK', PTMDB))

    def test_modifications(self):
        self.assertEqual(
            ('ACPLKYMR', [
                ModSite(304.20536, 'nterm', 'iTRAQ8plex'),
                ModSite(15.994915, 3, 'Oxidation'),
                ModSite(44.985078, 6, 'Nitro'),
            ]),
            parse_sequence('[iTRAQ8plex]ACP[Oxidation]LKY[Nitro]MR', PTMDB)
        )


@unittest.mock.patch('pepfrag_cli.parse_sequence', autospec=True)
class TestPeptideConstruction(unittest.TestCase):
    def test_no_cterm(self, mock_parse):
        seq = 'TLPMWYK'
        mods = [
            ModSite(304.20536, 'nterm', 'iTRAQ8plex'),
            ModSite(15.994915, 5, 'Oxidation')
        ]
        mock_parse.return_value = (seq, mods)
        self.assertEqual(
            Peptide(seq, 3, mods),
            construct_peptide(
                '[iTRAQ8plex]TLPMW[Oxidation]YK', None, 3, False, PTMDB
            )
        )

    def test_cterm(self, mock_parse):
        seq = 'TLPMWYK'
        mods = [
            ModSite(304.20536, 'nterm', 'iTRAQ8plex'),
            ModSite(15.994915, 5, 'Oxidation'),
            ModSite(37.955882, 'cterm', 'Cation:K'),
        ]
        mock_parse.return_value = (seq, mods)
        self.assertEqual(
            Peptide(seq, 3, mods),
            construct_peptide(
                '[iTRAQ8plex]]TLPMW[Oxidation]YK', 'Cation:K', 3, False, PTMDB
            )
        )

    def test_radical(self, mock_parse):
        seq = 'TLPMWYK'
        mods = [
            ModSite(304.20536, 'nterm', 'iTRAQ8plex'),
            ModSite(15.994915, 5, 'Oxidation')
        ]
        mock_parse.return_value = (seq, mods)
        self.assertEqual(
            Peptide(seq, 2, mods, radical=True),
            construct_peptide(
                '[iTRAQ8plex]TLPMW[Oxidation]YK', None, 2, True, PTMDB
            )
        )


if __name__ == '__main__':
    unittest.main()
