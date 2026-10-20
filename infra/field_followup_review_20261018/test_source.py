import sys
import unittest
from pathlib import Path

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import field_followup_20261018 as f


class FollowupTests(unittest.TestCase):
    def test_hidden_counts_cannot_change_gene_inputs(self):
        x = np.arange(60).reshape(10, 6)
        genes = np.array(list('ABCDEF')); vocab = list('CDEF')
        z, qc = f.gene_input(sparse.csr_matrix(x), genes, vocab)
        x[:, :2] += 100000
        zz, qq = f.gene_input(sparse.csr_matrix(x), genes, vocab)
        np.testing.assert_array_equal(z, zz)
        np.testing.assert_array_equal(qc, qq)

    def test_common_vocabulary_retains_absence_and_excludes_target(self):
        self.assertEqual(f.common_vocab({'A', 'B', 'C'}, [{'A', 'B'}, {'B', 'C'}], {'A'}), ['B'])

    def test_duplicate_unused_gene_does_not_invalidate_visible_input(self):
        x = sparse.csr_matrix(np.ones((4,4)))
        f.gene_input(x, np.array(['A','A','B','C']), ['B','C'])
        with self.assertRaisesRegex(ValueError, 'exactly once'):
            f.gene_input(x, np.array(['A','A','B','C']), ['A','B'])

    def test_dictionary_sampling_balances_units_with_two_sections(self):
        rows = [{'section_id': 'a', 'unit_id': 'p', 'n_spots': 1000},
                {'section_id': 'b', 'unit_id': 'p', 'n_spots': 1000},
                {'section_id': 'c', 'unit_id': 'q', 'n_spots': 1000}]
        idx = f.dictionary_indices(rows, 128, 17)
        self.assertEqual(len(idx['a']) + len(idx['b']), len(idx['c']))

    def test_inner_partition_is_unit_based(self):
        rows = [{'section_id': 'a', 'unit_id': 'p'}, {'section_id': 'b', 'unit_id': 'p'}]
        self.assertEqual(f.inner_fold(rows[0]['unit_id']), f.inner_fold(rows[1]['unit_id']))


if __name__ == '__main__':
    unittest.main()
