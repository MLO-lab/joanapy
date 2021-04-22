# joanapy
Joint continuous multi-level Ontology enrichment ANAlysis (JOANA).

## Build and install joanapy
Make sure that your current directory is
```
../joanapy/src
```
and then run the following command to build and install joanapy:
```
python3 -m build
python3 -m pip install dist/joanapy-0.0.1.tar.gz --no-deps joanapy
```
The package can be uninstalled with the following command:
```
pip uninstall -y joanapy
```

## Example
An example for how to run multi-level ontology enrichment analysis with joanapy can be found in the joanapy/tests folder.

## Fitting a mixture of Beta distributions
Code was adapted from Schröder C, Rahmann S. A hybrid parameter estimation algorithm for beta mixtures and applications to methylation state classification. Algorithms Mol Biol. 2017 Aug 18;12:21. doi: 10.1186/s13015-017-0112-1. PMID: 28828033; PMCID: PMC5563068 (https://bitbucket.org/genomeinformatics/betamix/src/master/).