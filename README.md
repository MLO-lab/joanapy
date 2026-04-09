# joanapy

Joint continuous multi-Omics enrichment ANAlysis (JOANA) from MLO Lab.

## Install joanapy
We recommend utilizing conda for environment management, and pip for installing joanapy as a Python package. Follow these instructions to set up and activate joanapy

```
conda create -n joana python=3.11
conda activate joana
```
Before installing joanapy try to install mono on your conda environment by the following command

```
conda install conda-forge::mono
```

Make sure that your working directory is JOANA-main which you have downloaded. 
Use pip to install joanapy on your conda environment. 

```
pip install .
```

after installing joanapy on joana environment run JOANA through run-joana function.

```
run-joana -o <omics1.txt> [-o2 <omics2.txt>] -p <pathwayfile.gmt> -d <output_directory> [-m <min_num_genes>]

-o <omics1.txt> 
Path to the primary omics input file.
The file must be a two-column tab-delimited file:
1. Gene name
2. Significance score (e.g., q-value) corresponding to the gene

-o2 <omics2.txt> (optional)
Path to the second omics input file for multi-omics analysis.
Format must be the same as -o.

-p <pathwayfile.gmt>
Path to the pathway file in GMT format, containing biological pathways to be tested for enrichment.

-m <min_num_genes> (optional)
A value in the range [0, 1) (default: 0).
Defines the minimum proportion of genes in a pathway that must have measurements.
Example:
-m 0.5 → Only pathways where at least 50% of genes have measurements will be considered.

-d <output_directory>
Path to the directory where JOANA results will be saved.

```
Note:
Use full paths (e.g., /home/user/data/file.txt) or relative paths (e.g., data/file.txt) depending on your working directory.
The -o2 parameter is only required for multi-omics analysis.


Input file format (-o and -o2)
The input files specified by -o and -o2 must contain two columns with the following structure:
1. Gene identifier (e.g., gene symbol)
2. Numeric score (e.g., q-value or p-value)

⚠️ Important:
- Files must not contain a header row
- Only the first two columns are used
- The second column must contain numeric values

Supported file types
JOANA supports the following formats:
- .txt → whitespace-separated (spaces or tabs)
- .tsv → tab-separated


Example (TXT / TSV format)
```
A2ML1  0.025202476125022
A3GALT2  0.878666355638669
A4GALT  0.983155339235838
A4GNT  0.971337673847852
AAAS  0.0863723498889275
AACS  0.230709278931887
AADAC  0.881216487254285
```

The 'gmt' file could be downloaded from msigDB or any other desired biological pathway file with gmt format.



```
run-joana -o /path/to/omics1.txt -p /path/to/pathway.gmt -d /path/to/dirOutputs/

```
And to execute JAOAN on multi-omics data the command line would be:

```
run-joana -o /path/to/omics1.txt -o2 /path/to/omics2.txt -p /path/to/pathway.gmt -m 0.7 -d /path/to/dirOutputs/

```
Note:
When dealing with multi-omics data, the '-o' input file serves as the reference file, and missing values in the second modality '-o2' are handled based on the reference data-modality. It's crucial to select the file with more gene measurements as the reference, as this provides better data integrity and completeness.

## Example with Sample Data
You can quickly test the tool using the included sample data:

```
run-joana -o ./sample_data/rna.txt -o2 ./sample_data/prot.txt -p ./sample_data/h.all.v6.2.symbols.gmt -m 0.7 -d ./dirOutputs/

```
## Original Data
For evaluation of JOANA on real data we used supplementary Excel files from:

[Lung adenocarcinoma dataset](https://www.sciencedirect.com/science/article/pii/S0092867420307443?via%3Dihub#app2)(Table S2)<br>

[Hot tumor dataset](https://www.sciencedirect.com/science/article/pii/S0092867420314513?via%3Dihub#app2)(Table S2)

[Myeloma (Single-cell transcriptomics dataset)](https://www.nature.com/articles/s41591-018-0269-2#data-availability)(GSE117156)

[coding and non-coding mutations](http://docs.icgc.org/pcawg) unfortunatly it is not accessible anymore

[LDC mouse hepatocyte dataset]()

## Uninstall joanapy
The package can be uninstalled with the following command:

```
pip uninstall joanapy
```


## Fitting a mixture of Beta distributions
Code was adapted from Schröder C, Rahmann S. A hybrid parameter estimation algorithm for beta mixtures and applications to methylation state classification. Algorithms Mol Biol. 2017 Aug 18;12:21. doi: 10.1186/s13015-017-0112-1. PMID: 28828033; PMCID: PMC5563068 (https://bitbucket.org/genomeinformatics/betamix/src/master/).
