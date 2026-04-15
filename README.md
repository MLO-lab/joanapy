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
## Data Sources

For the evaluation of JOANA on real data, we used supplementary Excel files from the following studies:

[Lung adenocarcinoma dataset](https://www.sciencedirect.com/science/article/pii/S0092867420307443?via%3Dihub#app2) (Table S2, Table S3)<br>
[Hot tumor dataset](https://www.sciencedirect.com/science/article/pii/S0092867420314513?via%3Dihub#app2) (Table S2)<br>
[Myeloma (single-cell transcriptomics dataset)](https://www.nature.com/articles/s41591-018-0269-2#data-availability) (GSE117156)<br>
[Coding and non-coding mutations dataset](http://docs.icgc.org/pcawg) *(currently not accessible)*<br>
[LDC mouse hepatocyte dataset](https://pmc.ncbi.nlm.nih.gov/articles/PMC10324225/#sec27) (13072_2023_504_MOESM3_ESM.bed,13072_2023_504_MOESM2_ESM.xlsx)

---

## Data Preprocessing

We provide the scripts used for preprocessing each dataset in the `scripts/` folder.

Due to licensing restrictions, we do not redistribute the original datasets.  
Please download the data from the sources above and run the preprocessing scripts to reproduce the inputs used in JOANA.

## Troubleshooting

If you encounter issues while installing or running JOANA, check the following common problems and solutions:

### 1. run-joana: command not found
Cause: The package was not installed correctly or the environment is not activated.

#### Solution:

Ensure your conda environment is activated:
  ```
  conda activate joana
  ```
Reinstall the package:
  ```
  pip install .
  ```
Verify installation:
  ```
  pip list | grep joana
  ```

### 2. Error related to mono not found

Cause: JOANA depends on Mono, which may not be installed properly.

#### Solution:
```
conda install -c conda-forge mono
```
Then verify:
  ```
  mono --version
  ```
### 3. Input file format errors

Cause: Input file does not follow the required two-column structure.
#### Common error:
```
ValueError: Error reading file 'omics1.txt'. Ensure it is a properly formatted two-column file (geneSymbol, q-value) without extra columns.
Original error: Error tokenizing data. C error: Expected 2 fields in line 2, saw 3
```
#### Explanation:
- JOANA expects exactly 2 columns only:
    1. Gene name
    2. Numeric score
- This error usually happens when:
    A header AND row index are both included → creates a third column
- Files with only a header or only row names may still work, but both together will cause failure

#### Solution:
- Ensure the file has exactly two columns only
- Remove any extra index column when exporting (e.g., from Excel or pandas)
- Save as .txt or .tsv
- Correct format example:
  ```
  A2ML1  0.025202476125022
  A3GALT2  0.878666355638669
  A4GALT  0.983155339235838
  ```
### 4. Unsupported file type
Cause: Input file is not in a supported format (e.g., .csv).
#### Common error:
```
Unsupported file type '.csv' for file 'prot1.csv'. Only '.txt' and '.tsv' files are supported
```
#### Solution:
- Convert .csv to .tsv or .txt
  
### 5. File/path errors (No such file or directory)

Cause: Cause: Incorrect file path or missing file (applies to input files and .gmt pathway file)..

#### Solution:
- Ensure the file exists and path is correct:
  ```
  ls /path/to/pathway.gmt
  ```
- Use absolute paths:
  ```
  /home/user/data/file.txt
  ```
- Or verify your current working directory for relative paths:
  ```
  pwd
  ```
### 6. Empty output or pandas.errors.EmptyDataError

Cause:
- -m parameter too strict (filters out all pathways)
- Input data does not meet minimum gene coverage

#### Common error:
```
  pandas.errors.EmptyDataError: No columns to parse from file
```
#### Solution:
- Try lowering the -m threshold:
  ```
  -m 0.3
  ```
- Check that enough genes in your dataset overlap with pathway genes
- Ensure input files are not empty and properly formatted

### 7. Unrecognized arguments error
Cause: Incorrect command syntax (using -o1 instead of -o).

#### Common error:
```
run-joana: error: unrecognized arguments: omics1.txt
```
#### Explanation:
- The correct flag for the primary omics file is -o, not -o1
- Using -o1 causes the argument parser to misinterpret the command
#### Solution:
- Use the correct command format:
  ```
  run-joana -o omics1.txt -o2 omics2.txt -p h.all.v6.2.symbols.gmt -d ./output
  ```
- General syntax:
  ```
  run-joana -o <omics1.txt> [-o2 <omics2.txt>] -p <pathway.gmt> -d <output_directory>
  ```
### 8. TypeError: expected str, bytes or os.PathLike object, not NoneType
Cause: Missing required -d (output directory) argument.
#### Common error:
```
TypeError: expected str, bytes or os.PathLike object, not NoneType
```

Example of incorrect command:
```
run-joana -o omics1.txt -o2 omics2.txt -p h.all.v6.2.symbols.gmt
```

#### Solution:
Always include the -d argument:
```
run-joana -o omics1.txt -o2 omics2.txt -p h.all.v6.2.symbols.gmt -d ./output
```

#### Explanation:
- The -d parameter (output directory) is required
- If -d is not provided, JOANA receives None instead of a path, causing this error
### 9. Permission errors when writing output

Cause: No write access to output directory.

#### Solution:
- Use a directory you own:
  ```
  mkdir -p ./dirOutputs
  ```
- Or change permissions:
  ```
  chmod u+w /path/to/output_directory
  ```

### 10. Python or dependency issues

Cause: Incompatible Python version or missing dependencies.

#### Solution:
- Ensure Python 3.11 is used:
  ```
  python --version
  ```
- Recreate environment:
  ```
  conda remove -n joana --all
  conda create -n joana python=3.11
  conda activate joana
  ```
### 11. Unexpected crashes or errors
If you encounter an issue not listed here:
- Double-check all input formats and parameters
- Run with minimal example:
  ```
  run-joana -o ./sample_data/rna.txt -p ./sample_data/h.all.v6.2.symbols.gmt -d ./dirOutputs/
  ```
- Ensure all dependencies are installed

## Uninstall joanapy
The package can be uninstalled with the following command:

```
pip uninstall joanapy
```


## Fitting a mixture of Beta distributions
Code was adapted from Schröder C, Rahmann S. A hybrid parameter estimation algorithm for beta mixtures and applications to methylation state classification. Algorithms Mol Biol. 2017 Aug 18;12:21. doi: 10.1186/s13015-017-0112-1. PMID: 28828033; PMCID: PMC5563068 (https://bitbucket.org/genomeinformatics/betamix/src/master/).
