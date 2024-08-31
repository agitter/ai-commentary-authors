"""
Order authors in a file
- Remove all characters besides the 20 canonical amino acids
- Concatenate 10 copies of the amino acid sequence
- Use the ESMFold API to predict a protein structure from the sequence
- Extract the pLDDT from the structure
- Sort the names by the corresponding pLDDT in descending order
"""

import argparse
from pathlib import Path


def parse_arguments():
    """
    Process command line arguments.
    @return arguments
    """
    parser = argparse.ArgumentParser(
        description="Order author names by pLDDT of corresponding predicted protein structures from ESMFold"
    )
    parser.add_argument("--input", type=Path, default=Path('authors.txt'), help="Input unordered author list")
    parser.add_argument("--output", type=Path, default=Path('ordered_authors.txt'), help="Output sorted author list")
    parser.add_argument("--pdb_dir", type=Path, default=Path('pdbs'), help="Path to cache predicted structures")

    return parser.parse_args()


def main():
    """
    Parse arguments and author ordering
    """
    args = parse_arguments()

    pdb_df = fetch_pdbs(args.input, args.pdb_dir)
    pdb_df = extract_plddt(pdb_df)
    write_ordered_authors(pdb_df, args.output)


if __name__ == "__main__":
    main()
