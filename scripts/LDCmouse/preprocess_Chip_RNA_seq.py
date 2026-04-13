import pandas as pd
import statsmodels.stats
from statsmodels.stats import multitest
import biomart
import mygene
from pybedtools import BedTool
import GTF_Processing



"""Read the Additional File 3 from doi.org/10.1186/s13072-023-00504-8 and link the regions to the nearest gene."""
######## Reading downloaded .bed file from supp paper 
chip_table = pd.read_table('pathto/rawdata/LDCmouse/Supp.bed', sep='\t', skiprows=1,
                           header=None)
chip_regions = BedTool('\n'.join([chip_table[0] + '\t' + chip_table[1].astype(str) + '\t' + chip_table[2].astype(str)][0]), from_string=True)

########### 
# Please download the mouse GENCODE annotation file (vM21) and place it in:
# pathto/rawdata/LDCmouse/
# 
#
# Update the path below to point to the downloaded file.
#  
annotation = 'pathto/rawdata/LDCmouse/gencode.vM21.annotation.gtf.gz'
gene_tss = GTF_Processing.gene_window_bed(annotation, extend=0, tss_type='5')
bed_closest = chip_regions.sort().closest(gene_tss.sort(), t='first')
chip_closest = pd.DataFrame([x.fields for x in bed_closest])
chip_closest.index = chip_closest[0] + '\t' + chip_closest[1].astype(str) + '\t' + chip_closest[2].astype(str)

# Merge the DataFrames to have the nearest gene in the supplementary table.
chip_table.index = chip_table[0] + '\t' + chip_table[1].astype(str) + '\t' + chip_table[2].astype(str)
chip_table = chip_table.join(chip_closest[[3, 4, 5, 6]], how='left', rsuffix='intersect_col')
chip_table.to_csv('pathto/rawData/LDCmouse/Supp. File 3_closest5TSS.bed',
                  header=False, index=False, sep='\t')

###################  
Chip_seq = pd.read_csv(
    'pathto/rawdata/LDCmouse/Supp. File 3_closest5TSS.bed', sep='\t', header=None)
p_corrected = statsmodels.stats.multitest.fdrcorrection(Chip_seq.iloc[:, 5])[1]
Chip_seq['p_corrected'] = p_corrected
Chip_seq = Chip_seq[Chip_seq.iloc[:, 21] != '.']
Chip_seq = Chip_seq[[21, 'p_corrected']]
Chip_seq[21] = Chip_seq[21].apply(lambda x: x.split('.')[0])
Chip_seq.drop_duplicates(inplace=True)

Chip_seq = Chip_seq.groupby(21)['p_corrected'].apply(lambda x: x.min())

Chip_seq.to_csv(
    'pathto/rawdata/LDCmouse/Chip_seq_data_joana_logFC.tsv', sep='\t', header=False)

RNA_seq = pd.read_csv(
    'pathto/rawdata/LDCmouse/ShvsCt_deseq2_diff_expressed_genes.txt', sep='\t', index_col=0)

RNA_seq['row'] = RNA_seq['row'].apply(lambda x: x.split('.')[0])
RNA_seq['abs_log2FC'] = RNA_seq['log2FoldChange'].abs()
RNA_seq = RNA_seq.sort_values(by='abs_log2FC', ascending=False)
RNA_seq_keep = RNA_seq[RNA_seq['row'].isin(CHIP_seq[0].to_list())]
#RNA_seq = RNA_seq.iloc[:10000, :]
# RNA_seq = RNA_seq[RNA_seq['abs_log2FC'] > 1]
RNA_seq = pd.concat([RNA_seq, RNA_seq_keep])
RNA_seq.drop_duplicates(inplace=True)

p_corrected = statsmodels.stats.multitest.fdrcorrection(RNA_seq['pvalue'])[1]
RNA_seq['p_corrected'] = p_corrected
#plot = sns.histplot(RNA_seq['p_corrected'])
#fig = plot.get_figure()
#fig.savefig("RNA_qval_hist_10000.png")

RNA_seq = RNA_seq[['row', 'p_corrected','log2FoldChange']]
RNA_seq.shape
RNA_seq.drop_duplicates(inplace=True)
RNA_seq.shape

RNA_seq.to_csv(
    'pathto/rawdata/LDCmouse/RNA_seq_data_joana_logFC.tsv', sep='\t', index=False, header=False)
####### mapping Entrez ID for rna data
mg = mygene.MyGeneInfo()

# Query mygene for mouse
result = mg.querymany(
    RNA_seq['ensembl'].tolist(),
    scopes='ensembl.gene',
    fields='entrezgene',
    species='mouse'
)

# Convert result to DataFrame
mapping = pd.DataFrame(result)
mapping = mapping[['query', 'entrezgene']].rename(columns={'query':'ensembl', 'entrezgene':'entrez'})
mapping.dropna(inplace=True)

# Merge with RNA-seq data
RNA_seq_mapped = RNA_seq.merge(mapping, on='ensembl', how='left')

# Save final table
RNA_seq_mapped.to_csv('./RNA_seq_mouse_with_entrez.tsv', sep='\t', index=False)

####### mapping Entrez ID for chip data

mg = mygene.MyGeneInfo()

# Query mygene for mouse
result = mg.querymany(
    Chip_seq['ensembl'].tolist(),
    scopes='ensembl.gene',
    fields='entrezgene',
    species='mouse'
)

# Convert result to DataFrame
mapping = pd.DataFrame(result)
mapping = mapping[['query', 'entrezgene']].rename(columns={'query':'ensembl', 'entrezgene':'entrez'})
mapping.dropna(inplace=True)

# Merge with RNA-seq data
chip_seq_mapped = Chip_seq.merge(mapping, on='ensembl', how='left')

# Save final table
chip_seq_mapped.to_csv('pathto/rawdata/LDCmouse/chip_seq_mouse_with_entrez.tsv', sep='\t', index=False)
