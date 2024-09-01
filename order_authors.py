"""
Order authors in a file
- Remove all characters besides the 20 canonical amino acids
- Concatenate 10 copies of the amino acid sequence
- Use the ESMFold API to predict a protein structure from the sequence
- Extract the pLDDT from the structure
- Sort the names by the corresponding pLDDT in descending order
"""

import argparse
import re
import requests
import urllib3
import biotite.structure.io as bsio
from pathlib import Path
import pandas as pd

# Disable the expected warning from write_esmfold_pdb
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


def parse_arguments():
    """
    Process command line arguments.
    @return arguments
    """
    parser = argparse.ArgumentParser(
        description='Order author names by pLDDT of corresponding predicted protein structures from ESMFold'
    )
    parser.add_argument('--input', type=Path, default=Path('authors.txt'), help='Input unordered author list')
    parser.add_argument('--output', type=Path, default=Path('ordered_authors.txt'), help='Output sorted author list '
                                                                                         'without a file extension. '
                                                                                         '.txt and .tsv outputs will '
                                                                                         'be created from this path.')
    parser.add_argument('--pdb_dir', type=Path, default=Path('pdbs'), help='Path to cache predicted structures')

    return parser.parse_args()


def name_to_aa(name: str, copies: int = 1) -> str:
    """
    Covert an author name to an amino acid sequence. String all characters besides the 20 canonical amino acids.
    Concatenate the specified number of copies of the amino acid sequence under the assumption that most author names
    are very short.
    :param name: The author name
    :param copies: The number of copies to concatenate
    :return: The converted, concatenated author name as an amino acid sequence
    """
    if copies < 1:
        raise ValueError(f"copies must be a positive integer but was {copies}")

    name = name.upper()
    name = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', name)
    return name * copies


def write_esmfold_pdb(seq: str, pdb_dir: Path) -> Path:
    """
    Use the ESM Metagenomic Atlas API to run ESMFold and predict a protein structure from the amino acid sequence.
    Store the sequence as a .pdb file. If the PDB file already exists in the pdb_dir, the API is not called again.
    :param seq: The amino acid sequence
    :param pdb_dir: The directory in which to store the output .pdb file. The filename will be the sequence with .pdb
    added.
    :return: Path to the PDB file that was written or already exists
    """
    pdb_path = Path(pdb_dir, f'{seq}.pdb')
    if pdb_path.exists():
        print(f'{seq}.pdb already exists')
    else:
        # https://curlconverter.com/ to help create the request from the original curl command
        headers = {
            'Content-Type': 'application/x-www-form-urlencoded',
        }

        # ESM Metagenomic Atlas API https://esmatlas.com/about#api
        # Encountered certificate issues described at https://github.com/facebookresearch/esm/discussions/665
        # so set verify=False
        response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=seq, verify=False)

        with open(pdb_path, 'w') as f:
            f.write(response.text)
        print(f'Wrote {seq}.pdb')

    return pdb_path


def extract_plddt(pdb_path: Path) -> float:
    """
    Return the pLDDT from the predicted PDB structure
    :param pdb_path: Path to the .pdb file
    :return: pLDDT value
    """
    # pLDDT calculation from https://github.com/facebookresearch/esm?tab=readme-ov-file#esmfold-structure-prediction-
    # Confirmed by discussion at https://github.com/facebookresearch/esm/discussions/608
    struct = bsio.load_structure(str(pdb_path), extra_fields=['b_factor'])
    return struct.b_factor.mean()


def fetch_pdbs(authors: Path, pdb_dir: Path) -> pd.DataFrame:
    """
    For each author name, convert to an amino acid sequence, predict the protein structure with ESMFold, and extract
    the pLDDT of the predicted structure.
    :param authors: a list of authors, one per line, with no header
    :param pdb_dir: The directory in which to store the output .pdb files
    :return:
    """
    if not pdb_dir.exists():
        pdb_dir.mkdir(exist_ok=True)

    author_df = pd.read_csv(authors, header=None, names=['Authors'])
    author_df['Seqs'] = author_df['Authors'].apply(name_to_aa, copies=1)
    author_df['PDBs'] = author_df['Seqs'].apply(write_esmfold_pdb, pdb_dir=pdb_dir)
    author_df['pLDDTs'] = author_df['PDBs'].apply(extract_plddt)

    return author_df


def write_ordered_authors(pdb_df: pd.DataFrame, output: Path) -> None:
    """
    Sort the author dataframe by pLDDT and write both output files
    :param pdb_df: the author dataframe with 'Authors' and 'pLDDT' columns
    :param output: the Path to the output without extension that will be used to create a .txt file with only the
    sorted author names and a .tsv file with the entire sorted dataframe
    :return: None
    """
    out_txt = output.with_suffix('.txt')
    out_tsv = output.with_suffix('.tsv')

    pdb_df = pdb_df.sort_values(by='pLDDTs', ascending=False)
    pdb_df['Authors'].to_csv(out_txt, header=False, index=False)
    pdb_df.to_csv(out_tsv, header=True, index=False, sep='\t')


def main():
    """
    Parse arguments and author ordering
    """
    args = parse_arguments()

    pdb_df = fetch_pdbs(args.input, args.pdb_dir)
    write_ordered_authors(pdb_df, args.output)

    # write_esmfold_pdb('MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE', args.pdb_dir)


if __name__ == '__main__''':
    main()
