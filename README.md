# "A renewed call for open artificial intelligence in biomedicine" author ordering
[![GitHub Actions](https://github.com/agitter/ai-commentary-authors/actions/workflows/test.yml/badge.svg)](https://github.com/agitter/ai-commentary-authors/actions/workflows/test.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13723236.svg)](https://doi.org/10.5281/zenodo.13723236)

This repository is used to order authors (except the first author) in the commentary manuscript "A renewed call for open artificial intelligence in biomedicine".
The algorithm was designed in advance so that the final author order would be unknown to any of the authors, who were only told the "author list will be sorted by some unexpected, non-alphabetical key".
It can accommodate adding new authors during revisions.

Ordering algorithm:
- Remove all characters besides the 20 canonical amino acids
- Concatenate 10 copies of the amino acid sequence
- Use the [ESMFold](https://doi.org/10.1126/science.ade2574) [API](https://esmatlas.com/about#api) to predict a protein structure from the sequence
- Extract the [pLDDT](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/) from the structure
- Sort the names by the corresponding pLDDT in descending order

## Running the script
The script `order_authors.py` orders authors from an input text file with no header and one author per line.

Example using data available in the repository:
```bash
python order_authors.py --input test/example_authors.txt --output ordered_authors --pdb_dir test/pdbs --copies 2
```

Usage:
```
usage: order_authors.py [-h] [--input INPUT] [--output OUTPUT]
                        [--pdb_dir PDB_DIR] [--copies COPIES]

Order author names by pLDDT of corresponding predicted protein structures from
ESMFold

options:
  -h, --help         show this help message and exit
  --input INPUT      Input unordered author list
  --output OUTPUT    Output sorted author list without a file extension. .txt
                     and .tsv outputs will be created from this path.
  --pdb_dir PDB_DIR  Path to cache predicted structures
  --copies COPIES    Number of copies of each author name to concatenate
                     (default 10)
```
