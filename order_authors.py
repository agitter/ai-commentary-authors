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
from pathlib import Path
import pandas as pd


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


def name_to_aa(name: str) -> str:
    name = name.upper()
    return re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', name)


def write_esmfold_pdb(seq: str, pdb_dir: Path) -> None:
    # https://curlconverter.com/ to help create the request from the original curl command
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }

    # ESM Metagenomic Atlas API https://esmatlas.com/about#api
    # Encountered certificate issues described at https://github.com/facebookresearch/esm/discussions/665
    # so set verify=False
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=seq, verify=False)

    pdb_path = Path(pdb_dir, f'{seq}.pdb')
    with open(pdb_path, 'w') as f:
        f.write(response.text)
    print(f'Wrote {seq}.pdb')


def fetch_pdbs(authors: Path, pdb_dir: Path) -> pd.DataFrame:
    if not pdb_dir.exists():
        pdb_dir.mkdir(exist_ok=True)

    author_df = pd.read_csv(authors, header=None, names=['Authors'])
    author_df['Seqs'] = author_df['Authors'].apply(name_to_aa)

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

    pdb_df = pdb_df.sort_values(by='Authors')
    pdb_df['Authors'].to_csv(out_txt, header=False, index=False)
    pdb_df.to_csv(out_tsv, header=True, index=False, sep='\t')


def main():
    """
    Parse arguments and author ordering
    """
    args = parse_arguments()

    pdb_df = fetch_pdbs(args.input, args.pdb_dir)
    # pdb_df = extract_plddt(pdb_df)
    write_ordered_authors(pdb_df, args.output)

    write_esmfold_pdb('MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE', args.pdb_dir)


if __name__ == '__main__''':
    main()
