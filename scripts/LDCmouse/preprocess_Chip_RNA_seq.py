import pandas as pd
import statsmodels.stats
from statsmodels.stats import multitest
import biomart
import mygene

Chip_seq = pd.read_csv(
    'Supp. File 3_closest5TSS.bed', sep='\t', header=None)
p_corrected = statsmodels.stats.multitest.fdrcorrection(Chip_seq.iloc[:, 5])[1]
Chip_seq['p_corrected'] = p_corrected
Chip_seq = Chip_seq[Chip_seq.iloc[:, 21] != '.']
Chip_seq = Chip_seq[[21, 'p_corrected']]
Chip_seq[21] = Chip_seq[21].apply(lambda x: x.split('.')[0])
Chip_seq.drop_duplicates(inplace=True)

Chip_seq = Chip_seq.groupby(21)['p_corrected'].apply(lambda x: x.min())

Chip_seq.to_csv(
    './Chip_seq_data_joana_logFC.tsv', sep='\t', header=False)

RNA_seq = pd.read_csv(
    'ShvsCt_deseq2_diff_expressed_genes.txt', sep='\t', index_col=0)

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
    './RNA_seq_data_joana_logFC.tsv', sep='\t', index=False, header=False)
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
chip_seq_mapped.to_csv('./chip_seq_mouse_with_entrez.tsv', sep='\t', index=False)
