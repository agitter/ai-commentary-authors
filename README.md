# "A renewed call for open artificial intelligence in biomedicine" author ordering
This repository is used to order authors in the commentary manuscript "A renewed call for open artificial intelligence in biomedicine".
The algorithm was designed in advance so that the final author order would be unknown to any of the authors, who were only told the "author list will be sorted by some unexpected, non-alphabetical key".
It can accommodate adding new authors during revisions.

Ordering algorithm:
- Remove all characters besides the 20 canonical amino acids
- Concatenate 10 copies of the amino acid sequence
- Use the ESMFold [API](https://esmatlas.com/about#api) to predict a protein structure from the sequence
- Extract the [pLDDT](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/) from the structure
- Sort the names by the corresponding pLDDT in descending order
