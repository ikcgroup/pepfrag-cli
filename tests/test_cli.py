import unittest

from pepfrag import IonType, ModSite, Peptide

from pepfrag_cli import parse_sequence


class TestParseSequence(unittest.TestCase):
    def test_no_modifications(self):
        self.assertEqual(('ACPLK', tuple()), parse_sequence('ACPLK'))

    def test_nterm_modification(self):
        self.assertEqual(
            ('ACPLK', (ModSite()))
        )


if __name__ == '__main__':
    unittest.main()
