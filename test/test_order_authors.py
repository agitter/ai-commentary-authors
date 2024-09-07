"""
Test functions for order_authors.py
"""

import math
from pathlib import Path

import pandas as pd
import pytest

import order_authors

GB1_SEQ = 'MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE'
PDB_DIR = Path('test/pdbs/')


class TestOrderAuthors:
    @pytest.mark.parametrize('name, seq',
                             [("Im A.N. E'rror", 'IMANERRR'),
                              (GB1_SEQ, GB1_SEQ),
                              ('~HAL 9000~', 'HAL')])
    def test_name_to_aa_valid(self, name: str, seq: str):
        assert order_authors.name_to_aa(name) == seq

    def test_write_esmfold_pdb_existing(self):
        assert order_authors.write_esmfold_pdb(GB1_SEQ, PDB_DIR, 1).exists()

    def test_write_esmfold_pdb_new(self):
        seq = 'HITHERE'
        copies = 1
        pdb_path = Path(PDB_DIR, f'{seq}-{copies}.pdb')
        pdb_path.unlink(missing_ok=True)
        assert order_authors.write_esmfold_pdb(seq, PDB_DIR, copies).exists()
        pdb_path.unlink(missing_ok=True)

    def test_write_esmfold_pdb_invalid(self):
        with pytest.raises(ValueError):
            order_authors.write_esmfold_pdb(GB1_SEQ, PDB_DIR, 0)
        with pytest.raises(ValueError):
            order_authors.write_esmfold_pdb(GB1_SEQ, PDB_DIR, -1)

    def test_extract_plddt(self):
        pdb_path = Path(PDB_DIR, f'{GB1_SEQ}-1.pdb')
        plddt = order_authors.extract_plddt(pdb_path)
        assert math.isclose(plddt, 0.849565, abs_tol=1e-6)

        pdb_path = Path(PDB_DIR, f'{GB1_SEQ}-2.pdb')
        plddt = order_authors.extract_plddt(pdb_path)
        assert math.isclose(plddt, 0.814977, abs_tol=1e-6)

    def test_fetch_pdbs(self):
        author_df = order_authors.fetch_pdbs(Path('test/example_authors.txt'), PDB_DIR, 2)
        expected_author_df = pd.read_csv('test/unsorted_author_df.tsv', sep='\t')

        # Consider operating system path differences and floating point imprecision before comparing
        author_df['PDBs'] = author_df['PDBs'].map(convert_path)
        expected_author_df['PDBs'] = expected_author_df['PDBs'].map(convert_path)
        author_df['pLDDTs'] = author_df['pLDDTs'].round(5)
        expected_author_df['pLDDTs'] = expected_author_df['pLDDTs'].round(5)

        assert author_df.equals(expected_author_df)

    def test_write_ordered_authors(self):
        author_df_path = Path('test/sorted_authors.tsv')
        author_df = pd.read_csv('test/unsorted_author_df.tsv', sep='\t')
        order_authors.write_ordered_authors(author_df, author_df_path)

        # Read back in the sorted author list from the file
        # Only test the dataframe, not the standalone list
        author_df = pd.read_csv(author_df_path, sep='\t')

        author_df_path.unlink(missing_ok=True)
        author_df_path.with_suffix('.txt').unlink(missing_ok=True)

        expected_author_df = pd.read_csv('test/sorted_author_df.tsv', sep='\t')

        # Consider operating system path differences and floating point imprecision before comparing
        author_df['PDBs'] = author_df['PDBs'].map(convert_path)
        expected_author_df['PDBs'] = expected_author_df['PDBs'].map(convert_path)
        author_df['pLDDTs'] = author_df['pLDDTs'].round(5)
        expected_author_df['pLDDTs'] = expected_author_df['pLDDTs'].round(5)

        assert author_df.equals(expected_author_df)


def convert_path(file_path: Path) -> str:
    """
    File paths have to be converted for the stored expected output files because otherwise the dataframes may not
    match if the test is run on a different operating system than the one used when the expected output was generated
    due to Linux versus Windows file path conventions
    :param file_path: input file Path
    :return: string representation of the converted Path
    """
    return str(Path(file_path))
