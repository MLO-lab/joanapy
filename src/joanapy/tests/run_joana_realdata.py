from pandas.core.dtypes import missing
from joanapy.enrichment_obj import ENRICHMENT_OBJ
from joanapy.joana import JOANA
from joanapy.joana_helper import transform_gmt_assignment_matrix
import os
import pandas as pd

filename_output_cooperative = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/joanapy_test/output_joana_cooperative_single_spescies_test_min_term_utilization1.csv'
filename_assignment_matrix_gmt = '/Users/akopf/Documents/Projects/JOANA/data/term_annotations/Human_GOALL_with_GO_iea_March_01_2021_symbol.gmt'
# filename_assignment_matrix_gmt = '/Users/akopf/Documents/Projects/JOANA/data/term_annotations/Human_GOALL_with_GO_iea_March_01_2021_symbol_subset_for_tests.gmt'
filename_gene_ids = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/gene_symbols.txt'
filename_qvalues_first = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/qvalues_first_species_symptomatic.txt'
filename_qvalues_first_filtered = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/joanapy_test/qvalues_first_species_symptomatic_filtered.txt'
filename_qvalues_second = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/qvalues_second_species_asymptomatic.txt'
filename_qvalues_second_filtered = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/joanapy_test/qvalues_second_species_asymptomatic_filtered.txt'
filename_missing_first = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/missing_first_species_symptomatic.txt'
filename_missing_first_filtered = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/joanapy_test/missing_first_species_symptomatic_filtered.txt'
filename_missing_second = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/missing_second_species_asymptomatic.txt'
filename_missing_second_filtered = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/joanapy_test/missing_second_species_asymptomatic_filtered.txt'
min_term_size = 5
max_term_size = 1000
prior_pA = 1
goodness_of_fit = True
plot_components = True
signif_threshold=0.1
quantile=0.95
min_term_utilization = 1.


# read in first species
gene_ids = pd.read_table(filename_gene_ids, header=None)
filename_terms = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/joanapy_test/Human_GOALL_with_GO_iea_March_01_2021_symbol_terms_utilization%f.txt' %(min_term_size, max_term_size, min_term_utilization)
filename_assignment_matrix = '/Users/akopf/Documents/Projects/JOANA/data/RealData/myeloma/compare_asymptomatic_vs_symptomatic/BMPC#/joanapy_test/Human_GOALL_with_GO_iea_March_01_2021_symbol_utilization%f.txt' %(min_term_size, max_term_size, min_term_utilization)
if not os.path.exists(os.path.dirname(filename_assignment_matrix)):
    os.makedirs(os.path.dirname(filename_assignment_matrix))
ind_filtered_qvalues = transform_gmt_assignment_matrix(filename_assignment_matrix_gmt, gene_ids.values.squeeze().tolist(), filename_assignment_matrix, filename_terms, min_term_utilization=min_term_utilization)

# filtere qvalues
qvalues_first = pd.read_table(filename_qvalues_first, header=None)
qvalues_first_filtered = qvalues_first.iloc[ind_filtered_qvalues]
qvalues_first_filtered.to_csv(filename_qvalues_first_filtered, header=False, index=False)
qvalues_second = pd.read_table(filename_qvalues_second, header=None)
qvalues_second_filtered = qvalues_second.iloc[ind_filtered_qvalues]
qvalues_second_filtered.to_csv(filename_qvalues_second_filtered, header=False, index=False)

missing_first = pd.read_table(filename_missing_first, header=None)
missing_first_filtered = missing_first.iloc[ind_filtered_qvalues]
missing_first_filtered.to_csv(filename_missing_first_filtered, header=False, index=False)
missing_second = pd.read_table(filename_missing_second, header=None)
missing_second_filtered = missing_second.iloc[ind_filtered_qvalues]
missing_second_filtered.to_csv(filename_missing_second_filtered, header=False, index=False)


# run JOANA cooperative
enrichment_obj_cooperative = ENRICHMENT_OBJ(filename_assignment_matrix, filename_terms,
                                            filename_qvalues_first_filtered,
                                            filename_qvalues_second=filename_qvalues_second_filtered,
                                            filename_missing_first=filename_missing_first_filtered,
                                            filename_missing_second=filename_missing_second_filtered)
joana_cooperative = JOANA(enrichment_obj_cooperative)
joana_cooperative.run(filename_output_cooperative, goodness_of_fit=goodness_of_fit, plot_components=plot_components, prior_pA=prior_pA, min_term_size=min_term_size, max_term_size=max_term_size)
joana_cooperative.plot_hidden_significant_genes(os.path.dirname(filename_output_cooperative), signif_threshold=signif_threshold, quantile=quantile) 
# joana_single_species.plot(filename_output_cooperative.replace('.csv', '.pdf'))
print(enrichment_obj_cooperative.joana_output.sort_values(by='Cooperative', ascending=False).head())