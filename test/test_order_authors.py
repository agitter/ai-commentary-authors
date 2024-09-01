"""
Test functions for order_authors.py
"""

import math
import pytest
from pathlib import Path

import order_authors

GB1_SEQ = 'MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE'
PDB_DIR = Path('test/pdbs/')


class TestOrderAuthors:
    @pytest.mark.parametrize('name, copies, seq',
                             [("Im A.N. E'rror", 1, 'IMANERRR'),
                              ("Im A.N. E'rror", 2, 'IMANERRRIMANERRR'),
                              (GB1_SEQ, 1,
                               GB1_SEQ),
                              ('~HAL 9000~', 10, 'HALHALHALHALHALHALHALHALHALHAL')])
    def test_name_to_aa_valid(self, name: str, copies: int, seq: str):
        assert order_authors.name_to_aa(name, copies) == seq

    def test_name_to_aa_invalid(self):
        with pytest.raises(ValueError):
            assert order_authors.name_to_aa('~HAL 9000~', 0)
        with pytest.raises(ValueError):
            assert order_authors.name_to_aa('~HAL 9000~', -1)

    def test_write_esmfold_pdb_existing(self):
        assert order_authors.write_esmfold_pdb(GB1_SEQ, PDB_DIR).exists()

    def test_write_esmfold_pdb_new(self):
        seq = 'HITHERE'
        pdb_path = Path(PDB_DIR, f'{seq}.pdb')
        pdb_path.unlink(missing_ok=True)
        assert order_authors.write_esmfold_pdb(seq, PDB_DIR).exists()

    def test_extract_plddt(self):
        pdb_path = Path(PDB_DIR, f'{GB1_SEQ}.pdb')
        plddt = order_authors.extract_plddt(pdb_path)
        assert math.isclose(plddt, 0.814977)


