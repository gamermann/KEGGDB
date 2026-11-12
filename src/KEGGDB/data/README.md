## Third-Party Data

KEGGDB includes two auxiliary text files containing chemical information used by the package:

- **`AminoAcids.txt`**: contains amino acid names, abbreviations, chemical formulas, molecular weights, and physicochemical properties. The data are consistent with the [Amino Acids Reference Chart](https://www.sigmaaldrich.com/CH/en/technical-documents/technical-article/protein-biology/protein-structural-analysis/amino-acid-reference-chart) published by Sigma-Aldrich, which cites *Handbook of Chemistry and Physics*, 72nd Edition, CRC Press (1991).

- **`monosachs.txt`**: contains monosaccharide nomenclature, abbreviations, and classifications based on the [Symbol Nomenclature for Glycans (SNFG)](https://www.ncbi.nlm.nih.gov/glycans/snfg2_0_3.html), published by the National Center for Biotechnology Information (NCBI).

These files are included as auxiliary reference data and are distinct from the KEGG database content retrieved through the KEGG REST API.

The GPL-3.0-or-later license applies to the KEGGDB source code. It does not imply that third-party data are covered by the same license. Users should consult the original sources for applicable terms of use.
