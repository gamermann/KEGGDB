# KEGGDB

**A Python package for accessing, processing, and analyzing KEGG database information.**

## Overview

KEGGDB is a Python package designed to work with data from the Kyoto Encyclopedia of Genes and Genomes (KEGG) database. It was developed specifically to link metabolites through chemical reactions, with the goal of reconstructing abstract graphs representing organism metabolism.

## Features

- Access to KEGG data through its REST API.
- Reconstruction of chemical formulas for compounds and glycans.
- Identification and comparison of chemical compounds.
- Linking metabolites through the transformations they undergo in biochemical reactions.
- Construction of abstract graphs representing metabolic networks.

## Installation

Install KEGGDB using pip:

```bash
pip install KEGGDB
```

## Requirements

- Python >= 3.10

## Usage

The following example illustrates a basic use of KEGGDB:

```python
from KEGGDB import *

kegg = KEGGDB()

water = kegg.giveMe('C00001')
formula = water.getChemDic()
print(formula.toFormula())

reaction = ChemReac(kegg.giveMe('R00001').EQUATION, name='R00001', db=kegg)
print(reaction.latex())
scheme = printReac(reaction)
print(scheme)
pairs, tab = getPairsTab(reaction)
kegg.showPairs(pairs)
```

For more detailed examples and explanations, see the documentation.

## Documentation

The complete manual is available in the [KEGGDB Manual](https://github.com/gamermann/KEGGDB/blob/main/docs/Manual.pdf).

## KEGG Database Usage

KEGGDB retrieves data from the KEGG database through the
KEGG REST API and stores the retrieved data locally.

KEGG data is copyrighted by Kanehisa Laboratories.
Non-academic use of KEGG data requires a license.

Users are responsible for ensuring that their use of KEGG
data complies with the applicable licensing conditions.

KEGGDB is an independent software package and is not
affiliated with or endorsed by Kanehisa Laboratories.

For further information, please refer to the official KEGG
documentation:

- [KEGG REST API](https://www.genome.jp/kegg/rest/)
- [KEGG Copyright and Disclaimer](https://www.kegg.jp/kegg/legal.html)

## Acknowledgements

Daniel Gamermann acknowledges the financial support provided by the Brazilian National Council for Scientific and Technological Development (CNPq) through the Universal Call grant No. 402487/2023-0 and the Research Productivity Fellowship grant No. 305872/2022-2.

## License

KEGGDB is distributed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).

This project is distributed under the terms of the [LICENSE](LICENSE) file.

This license applies to the KEGGDB source code only. It does not
grant rights to KEGG database content or services.


## Third-Party Data

KEGGDB includes two auxiliary text files containing chemical information used by the package:

- **`AminoAcids.txt`**: contains amino acid names, abbreviations, chemical formulas, molecular weights, and physicochemical properties. The data are consistent with the [Amino Acids Reference Chart](https://www.sigmaaldrich.com/CH/en/technical-documents/technical-article/protein-biology/protein-structural-analysis/amino-acid-reference-chart) published by Sigma-Aldrich, which cites *Handbook of Chemistry and Physics*, 72nd Edition, CRC Press (1991).

- **`monosachs.txt`**: contains monosaccharide nomenclature, abbreviations, and classifications based on the [Symbol Nomenclature for Glycans (SNFG)](https://www.ncbi.nlm.nih.gov/glycans/snfg2_0_3.html), published by the National Center for Biotechnology Information (NCBI).

These files are included as auxiliary reference data and are distinct from the KEGG database content retrieved through the KEGG REST API.

The GPL-3.0-or-later license applies to the KEGGDB source code. It does not imply that third-party data are covered by the same license. Users should consult the original sources for applicable terms of use.

## Authors

See [AUTHORS](AUTHORS) for information about the contributors.

## Citation

If you use KEGGDB in your research, please cite:

Pellizzaro, J. A., Gamermann, D., and Triana Dopico, J. (2026).
*Metabolic Network Properties: Comprehensive Analysis Across Domains.*
arXiv:2608.21168.

[https://arxiv.org/abs/2608.21168](https://arxiv.org/abs/2608.21168)

When using KEGG data, please also cite the
appropriate KEGG publications:

[https://www.genome.jp/kegg/kegg1.html](https://www.genome.jp/kegg/kegg1.html)
