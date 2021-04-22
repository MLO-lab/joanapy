from joanapy.enrichment_obj import ENRICHMENT_OBJ
from joanapy.joana import JOANA
import os

filename_output_single_species = '/Users/akopf/Documents/Projects/JOANA/data/CooperativeToyData/ExpUniformBeta/PValues/TestData_Alpha1_Beta4_Alpha1_Beta4/50/output_cont_single-species_JOANAConsoleApp_test.csv'
filename_output_cooperative = '/Users/akopf/Documents/Projects/JOANA/data/CooperativeToyData/ExpUniformBeta/PValues/TestData_Alpha1_Beta4_Alpha1_Beta4/50/output_cont_cooperative_JOANAConsoleApp_test.csv'
filename_assignment_matrix = '/Users/akopf/Documents/Projects/JOANA/data/CooperativeToyData/ExpUniformBeta/PValues/assignmentMatrix.txt'
filename_terms = '/Users/akopf/Documents/Projects/JOANA/data/CooperativeToyData/ExpUniformBeta/PValues/terms.txt'
filename_qvalues_first = '/Users/akopf/Documents/Projects/JOANA/data/CooperativeToyData/ExpUniformBeta/PValues/TestData_Alpha1_Beta4_Alpha1_Beta4/50/qvalues_first_species.txt'
filename_qvalues_second = '/Users/akopf/Documents/Projects/JOANA/data/CooperativeToyData/ExpUniformBeta/PValues/TestData_Alpha1_Beta4_Alpha1_Beta4/50/qvalues_second_species.txt'
filename_missing_first = '/Users/akopf/Documents/Projects/JOANA/data/CooperativeToyData/ExpUniformBeta/PValues/TestData_Alpha1_Beta4_Alpha1_Beta4/50/missing_first_species.txt'
filename_missing_second = '/Users/akopf/Documents/Projects/JOANA/data/CooperativeToyData/ExpUniformBeta/PValues/TestData_Alpha1_Beta4_Alpha1_Beta4/50/missing_second_species.txt'
min_term_size = 75
max_term_size = 3000
prior_pA = 1
goodness_of_fit = True
plot_components = True
signif_threshold=0.1
quantile=0.95

# run single species model
enrichment_obj_single_species = ENRICHMENT_OBJ(filename_assignment_matrix, filename_terms, filename_qvalues_first)
joana_single_species = JOANA(enrichment_obj_single_species)
joana_single_species.run(filename_output_single_species, goodness_of_fit=goodness_of_fit, plot_components=plot_components, prior_pA=prior_pA, min_term_size=min_term_size, max_term_size=max_term_size)
joana_single_species.plot_hidden_significant_genes(os.path.dirname(filename_output_single_species), signif_threshold=signif_threshold, quantile=quantile)
print(enrichment_obj_single_species.joana_output.sort_values(by='Single-Species', ascending=False).head())
# joana_single_species.plot(filename_output_single_species.replace('.csv', '.pdf'))



# run cooperative model
enrichment_obj_cooperative = ENRICHMENT_OBJ(filename_assignment_matrix, filename_terms,
                                            filename_qvalues_first,
                                            filename_qvalues_second=filename_qvalues_second,
                                            filename_missing_first=filename_missing_first,
                                            filename_missing_second=filename_missing_second)
joana_cooperative = JOANA(enrichment_obj_cooperative)
joana_cooperative.run(filename_output_cooperative, goodness_of_fit=goodness_of_fit, plot_components=plot_components, prior_pA=prior_pA, min_term_size=min_term_size, max_term_size=max_term_size)
joana_cooperative.plot_hidden_significant_genes(os.path.dirname(filename_output_cooperative), signif_threshold=signif_threshold, quantile=quantile)
# joana_single_species.plot(filename_output_cooperative.replace('.csv', '.pdf'))
print(enrichment_obj_cooperative.joana_output.sort_values(by='Cooperative', ascending=False).head())


# load objects
enrichment_obj_single_species_reloaded = ENRICHMENT_OBJ.load(filename_output_single_species.replace('.csv', '.pickle'))
joana_single_species_plotting = JOANA(enrichment_obj_single_species_reloaded)
joana_single_species_plotting.plot(filename_output_single_species.replace('.csv', '.pdf'))
joana_single_species_plotting.plot_barplot(filename_output_single_species.replace('.csv', '_barplot.pdf'))


enrichment_obj_cooperative_reloaded = ENRICHMENT_OBJ.load(filename_output_cooperative.replace('.csv', '.pickle'))
joana_cooperative_plotting = JOANA(enrichment_obj_cooperative_reloaded)
joana_cooperative_plotting.plot(filename_output_cooperative.replace('.csv', '.pdf'), top_percentile_edges=0.9, n_top_terms=10)
joana_cooperative_plotting.plot_barplot(filename_output_cooperative.replace('.csv', '_barplot.pdf'))

