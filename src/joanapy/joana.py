from numpy.lib import index_tricks
from numpy.lib.function_base import quantile
from seaborn.palettes import color_palette
from joanapy.moment_fitting_core import MOMENT_FITTING
import os
from subprocess import call
from joanapy.enrichment_obj import ENRICHMENT_OBJ
import pandas as pd
import numpy as np
import psutil

import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.artist import Artist
from igraph import BoundingBox, Graph, palettes
import seaborn as sns

FILEPATH_JOANA_CLASS = os.path.realpath(__file__)
FILENAME_JOANA = os.path.join(os.path.dirname(FILEPATH_JOANA_CLASS), 'joana_app', 'MonaConsoleApp.exe')


class JOANA:

    def __init__(self, enrichment_obj: ENRICHMENT_OBJ):
        self.enrichment_obj = enrichment_obj


    def run(self, filename_output, filename_moment_fit_first=None, filename_moment_fit_second=None, tolerance_fitting=1E-5, steps_fitting=30000, verbose=True, goodness_of_fit=False, plot_components=False, prior_pA=1, min_term_size=0, max_term_size=100000, save_enrichment_obj=True):
        
        if not filename_output.endswith('.csv'):
            raise ValueError('The output filename needs to be a .csv file.')

        self.enrichment_obj.filename_joana_output = filename_output
        filename_output_ending = filename_output.split('.')[-1]
        dir_output = os.path.dirname(filename_output)
        if not os.path.exists(dir_output):
            os.makedirs(dir_output)

        # get term sizes
        term_sizes = np.zeros((len(self.enrichment_obj.terms), ))
        for i in range(len(self.enrichment_obj.assignment_matrix)):
            term_sizes[np.asarray(self.enrichment_obj.assignment_matrix[i])] += 1
        pd.DataFrame({'terms' : self.enrichment_obj.terms.values.flatten(), 'size' : term_sizes}).to_csv(self.enrichment_obj.filename_terms.replace('.txt', '_sizes.txt'), index=False, header=True, sep='&')

        # filter terms
        ind_keep = np.intersect1d(np.where(min_term_size <= term_sizes)[0], np.where(term_sizes <= max_term_size)[0])
        index_mapping = dict()
        for i in range(len(ind_keep)):
            index_mapping[ind_keep[i]] = i
        
        # filter assignment matrix
        assignment_matrix_filtered = list()
        for i in range(len(self.enrichment_obj.assignment_matrix)):
            temp_gene_assignemnt = list(self.enrichment_obj.assignment_matrix[i])

            # remove filtered terms
            list_remove_ind = list()
            for ind in temp_gene_assignemnt:
                if not ind in ind_keep:
                    list_remove_ind.append(ind)
            [temp_gene_assignemnt.remove(x) for x in list_remove_ind]

            # remap indices
            assignment_matrix_filtered.append(','.join(np.asarray([str(index_mapping[x]) for x in temp_gene_assignemnt])))

        # TODO potentially filter for genes without assignment
        
        self.enrichment_obj.filename_assignment_matrix = os.path.join(self.enrichment_obj.filename_assignment_matrix.replace('.txt', '_filtered_min%d_max%d.txt' %(min_term_size, max_term_size)))
        with open(self.enrichment_obj.filename_assignment_matrix, "w") as myfile:
            for line in assignment_matrix_filtered:
                myfile.write(str(line) + '\n')

        # filter terms
        terms_filtered = self.enrichment_obj.terms.iloc[ind_keep]
        self.enrichment_obj.filename_terms = self.enrichment_obj.filename_terms.replace('.txt', '_filtered_min%d_max%d.txt' %(min_term_size, max_term_size))
        terms_filtered.to_csv(self.enrichment_obj.filename_terms, index=False, header=False)


        if self.enrichment_obj.type == 'single-species':
            if filename_moment_fit_first is None:
                # first species
                filename_moment_fit_first = os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting.txt')
                moment_fitting_first = MOMENT_FITTING(self.enrichment_obj.qvalues_first,
                                                    filename_moment_fit_first,
                                                    tolerance=tolerance_fitting, steps=steps_fitting)
                moment_fitting_first.run(verbose=verbose)
                if goodness_of_fit:
                    moment_fitting_first.goodness_of_fit(os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting_gof.txt'), plot_histograms=True)
                if plot_components:
                    moment_fitting_first.plot_components_density(os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting_components.pdf'))
                self.enrichment_obj.moment_fit_first = moment_fitting_first

            while self.__has_handle(FILENAME_JOANA):
                continue

            call(['mono',
                  FILENAME_JOANA,
                  '1',
                  self.enrichment_obj.filename_assignment_matrix,
                  self.enrichment_obj.filename_terms,
                  self.enrichment_obj.filename_qvalues_first,
                  filename_moment_fit_first,
                  filename_output,
                  str(int(prior_pA))])


        if self.enrichment_obj.type == 'cooperative':
            if filename_moment_fit_first is None:
                # first species
                filename_moment_fit_first = os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting.txt')
                ind_not_missing_first = self.enrichment_obj.missing_first.values == 0
                moment_fitting_first = MOMENT_FITTING(self.enrichment_obj.qvalues_first[ind_not_missing_first.flatten()],
                                                    filename_moment_fit_first,
                                                    tolerance=tolerance_fitting, steps=steps_fitting)
                moment_fitting_first.run(verbose=verbose)
                if goodness_of_fit:
                    moment_fitting_first.goodness_of_fit(os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting_gof.txt'), plot_histograms=True)
                if plot_components:
                    moment_fitting_first.plot_components_density(os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting_components.pdf'))
                self.enrichment_obj.moment_fit_first = moment_fitting_first

            if filename_moment_fit_second is None:
                # second species
                filename_moment_fit_second = os.path.join(os.path.dirname(filename_output), 'qvalues_second_moment_fitting.txt')
                ind_not_missing_second = self.enrichment_obj.missing_second.values == 0
                moment_fitting_second = MOMENT_FITTING(self.enrichment_obj.qvalues_second[ind_not_missing_second.flatten()], filename_moment_fit_second,
                                                    tolerance=tolerance_fitting, steps=steps_fitting)
                moment_fitting_second.run(verbose=verbose)
                if goodness_of_fit:
                    moment_fitting_second.goodness_of_fit(os.path.join(os.path.dirname(filename_output),'qvalues_second_moment_fitting_gof.txt'), plot_histograms=True)
                if plot_components:
                    moment_fitting_second.plot_components_density(os.path.join(os.path.dirname(filename_output),'qvalues_second_moment_fitting_components.pdf'))
                self.enrichment_obj.moment_fit_second = moment_fitting_second

            while self.__has_handle(FILENAME_JOANA):
                continue

            call(['mono',
                FILENAME_JOANA,
                '0',
                self.enrichment_obj.filename_assignment_matrix,
                self.enrichment_obj.filename_terms,
                self.enrichment_obj.filename_qvalues_first,
                self.enrichment_obj.filename_qvalues_second,
                self.enrichment_obj.filename_missing_first,
                self.enrichment_obj.filename_missing_second,
                filename_moment_fit_first,
                filename_moment_fit_second,
                filename_output,
                str(int(prior_pA))])

        # load joana output in enrichment object
        if os.path.exists(filename_output):
            self.enrichment_obj.joana_output = pd.read_csv(filename_output, index_col='Terms')

        # save enrichment_obj after JOANA run
        if save_enrichment_obj:
            self.enrichment_obj.save(filename_output.replace(filename_output_ending, 'pickle'))


    def __has_handle(self, fpath):
        for proc in psutil.process_iter():
            try:
                for item in proc.open_files():
                    if fpath == item.path:
                        return True
            except Exception:
                pass

        return False

    # def __is_locked(self, filepath):
    #     """Checks if a file is locked by opening it in append mode.
    #     If no exception thrown, then the file is not locked.
    #     """
    #     locked = None
    #     file_object = None
    #     if os.path.exists(filepath):
    #         try:
    #             print("Trying to open %s." % filepath)
    #             buffer_size = 8
    #             # Opening file in append mode and read the first 8 characters.
    #             file_object = open(filepath, 'a', buffer_size)
    #             if file_object:
    #                 print("%s is not locked." % filepath)
    #                 locked = False
    #         except IOError:
    #             print("File is locked (unable to open in append mode). %s." %message)
    #             locked = True
    #         finally:
    #             if file_object:
    #                 file_object.close()
    #                 print("%s closed." % filepath)
    #     else:
    #         print("%s not found." % filepath)
    #     return locked

    def plot_hidden_significant_genes(self, dir_output, signif_threshold=0.1, quantile=0.95, output_filetag=''):
        if self.enrichment_obj.joana_output is None:
            print('Please run JOANA model at first or load a JOANA output.')
        else:
            if self.enrichment_obj.type == 'single-species':
                self._plot_hidden_significant_genes_single_species(dir_output, signif_threshold, quantile, output_filetag)
            elif self.enrichment_obj.type == 'cooperative':
                self.__plot_hidden_significant_genes_cooperative(dir_output, signif_threshold, quantile, output_filetag)


    def _plot_hidden_significant_genes_single_species(self, dir_output, signif_threshold=0.1, quantile=0.95, output_filetag=''):
        # define 0.95 quantile of signif beta
        a_first = self.enrichment_obj.moment_fit_first.ab[0][0]
        b_first = self.enrichment_obj.moment_fit_first.ab[0][1]
        empirical_quantile = np.round(np.quantile(np.random.beta(a_first, b_first, (10000000,)), quantile), 3)

        col = 'Single-Species'
        joana_output = self.enrichment_obj.joana_output.sort_values(by=col, ascending=False)
        joana_output_enriched = joana_output.iloc[joana_output[col].values > 0.5]

        qvalues_enriched_terms = list()
        genes_beta_signif_percentage = list()
        genes_beta_signif_absolut_numbers = list()
        for i in range(joana_output_enriched.shape[0]):
            temp_ind_term = np.where(self.enrichment_obj.terms[0].values == joana_output_enriched.index.values[i])[0]

            # find assigned genes
            ind_genes_assigned = list()
            for j in range(len(self.enrichment_obj.assignment_matrix)):
                if temp_ind_term in self.enrichment_obj.assignment_matrix[j]:
                    ind_genes_assigned.append(j)
            ind_genes_assigned = np.asarray(ind_genes_assigned)

            qvalues_enriched_terms.append(self.enrichment_obj.qvalues_first[ind_genes_assigned].flatten())
            genes_beta_signif_absolut_numbers.append(np.sum(np.logical_and(qvalues_enriched_terms[-1] > signif_threshold, qvalues_enriched_terms[-1] < empirical_quantile)))
            genes_beta_signif_percentage.append(genes_beta_signif_absolut_numbers[-1]/len(ind_genes_assigned))

        # histogram of percentage genes beta significant
        plt.figure()
        plt.hist(genes_beta_signif_percentage, bins=20)
        plt.xlabel('Percentage %s < genes < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
        plt.ylabel('Frequency', fontsize=18)
        plt.xlim([-0.05, 1.05])
        if output_filetag == '':
            plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
        plt.close()

        # histogram of number of genes beta significant
        plt.figure()
        plt.hist(genes_beta_signif_absolut_numbers, bins=20)
        plt.xlabel('#genes %s < gene < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
        plt.ylabel('Frequency', fontsize=18)
        if output_filetag == '':
            plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
        plt.close()

        # scatter percentage of genes beta significant vs. probability
        plt.figure()
        plt.scatter(joana_output_enriched[col].values.flatten(), genes_beta_signif_percentage)
        plt.xlabel('Probability', fontsize=18)
        plt.xlim([0.45, 1.05])
        plt.ylabel('Percentage %s < genes < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
        plt.ylim([-0.05, 1.05])
        if output_filetag == '':
            plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_vs_probability_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_vs_probability_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
        plt.close()

        # scatter number of genes beta significant vs. probability
        plt.figure()
        plt.scatter(joana_output_enriched[col].values.flatten(), genes_beta_signif_absolut_numbers)
        plt.xlabel('Probability', fontsize=18)
        plt.xlim([0.45, 1.05])
        plt.ylabel('#genes %s < gene < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
        if output_filetag == '':
            plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_vs_probability_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_vs_probability_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
        plt.close()








    def __plot_hidden_significant_genes_cooperative(self, dir_output, signif_threshold=0.1, quantile=0.95, output_filetag=''):

        # define 0.95 quantile of signif beta
        a_first = self.enrichment_obj.moment_fit_first.ab[0][0]
        b_first = self.enrichment_obj.moment_fit_first.ab[0][1]
        empirical_quantile_first = np.round(np.quantile(np.random.beta(a_first, b_first, (10000000,)), quantile), 3)

        # define 0.95 quantile of signif beta
        a_second = self.enrichment_obj.moment_fit_second.ab[0][0]
        b_second = self.enrichment_obj.moment_fit_second.ab[0][1]
        empirical_quantile_second = np.round(np.quantile(np.random.beta(a_second, b_second, (10000000,)), quantile), 3)

        for col in ['Single-Species-1', 'Single-Species-2']:
            joana_output = self.enrichment_obj.joana_output.sort_values(by=col, ascending=False)
            joana_output_enriched = joana_output.iloc[joana_output[col].values > 0.5]

            qvalues_enriched_terms = list()
            genes_beta_signif_percentage = list()
            genes_beta_signif_absolut_numbers = list()
            for i in range(joana_output_enriched.shape[0]):
                temp_ind_term = np.where(self.enrichment_obj.terms[0].values == joana_output_enriched.index.values[i])[0]

                # find assigned genes
                ind_genes_assigned = list()
                for j in range(len(self.enrichment_obj.assignment_matrix)):
                    if temp_ind_term in self.enrichment_obj.assignment_matrix[j]:
                        if col == 'Single-Species-1':
                            if self.enrichment_obj.missing_first.values[j][0] == 0:
                                ind_genes_assigned.append(j)
                        elif col == 'Single-Species-2':
                            if self.enrichment_obj.missing_second.values[j][0] == 0:
                                ind_genes_assigned.append(j)
                ind_genes_assigned = np.asarray(ind_genes_assigned)

                if col == 'Single-Species-1':
                    qvalues_enriched_terms.append(self.enrichment_obj.qvalues_first[ind_genes_assigned].flatten())
                    genes_beta_signif_absolut_numbers.append(np.sum(np.logical_and(qvalues_enriched_terms[-1] > signif_threshold, qvalues_enriched_terms[-1] < empirical_quantile_first)))
                elif col == 'Single-Species-2':
                    qvalues_enriched_terms.append(self.enrichment_obj.qvalues_second[ind_genes_assigned].flatten())
                    genes_beta_signif_absolut_numbers.append(np.sum(np.logical_and(qvalues_enriched_terms[-1] > signif_threshold, qvalues_enriched_terms[-1] < empirical_quantile_second)))
                genes_beta_signif_percentage.append(genes_beta_signif_absolut_numbers[-1] / len(ind_genes_assigned))

            # histogram of percentage genes beta significant
            plt.figure()
            plt.hist(genes_beta_signif_percentage, bins=20)
            plt.xlabel('Percentage %s < genes < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
            plt.ylabel('Frequency', fontsize=18)
            plt.xlim([-0.05, 1.05])
            if output_filetag == '':
                plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
            else:
                plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
            plt.close()

            # histogram of number of genes beta significant
            plt.figure()
            plt.hist(genes_beta_signif_absolut_numbers, bins=20)
            plt.xlabel('#genes %s < gene < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
            plt.ylabel('Frequency', fontsize=18)
            if output_filetag == '':
                plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
            else:
                plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
            plt.close()

            # scatter percentage of genes beta significant vs. probability
            plt.figure()
            plt.scatter(joana_output_enriched[col].values.flatten(), genes_beta_signif_percentage)
            plt.xlabel('Probability', fontsize=18)
            plt.xlim([0.45, 1.05])
            plt.ylabel('Percentage %s < genes < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
            plt.ylim([-0.05, 1.05])
            if output_filetag == '':
                plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_vs_probability_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
            else:
                plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_vs_probability_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
            plt.close()

            # scatter number of genes beta significant vs. probability
            plt.figure()
            plt.scatter(joana_output_enriched[col].values.flatten(), genes_beta_signif_absolut_numbers)
            plt.xlabel('Probability', fontsize=18)
            plt.xlim([0.45, 1.05])
            plt.ylabel('#genes %s < gene < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
            if output_filetag == '':
                plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_vs_probability_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
            else:
                plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_vs_probability_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
            plt.close()


        ## cooperative
        col = 'Cooperative'
        joana_output = self.enrichment_obj.joana_output.sort_values(by=col, ascending=False)
        joana_output_enriched = joana_output.iloc[joana_output[col].values > 0.5]

        qvalues_enriched_terms_first = list()
        genes_beta_signif_percentage_first = list()
        genes_beta_signif_absolut_numbers_first = list()
        qvalues_enriched_terms_second = list()
        genes_beta_signif_percentage_second = list()
        genes_beta_signif_absolut_numbers_second = list()
        for i in range(joana_output_enriched.shape[0]):
            temp_ind_term = np.where(self.enrichment_obj.terms[0].values == joana_output_enriched.index.values[i])[0]

            # find assigned genes
            ind_genes_assigned_first = list()
            ind_genes_assigned_second = list()
            for j in range(len(self.enrichment_obj.assignment_matrix)):
                if temp_ind_term in self.enrichment_obj.assignment_matrix[j]:
                    if self.enrichment_obj.missing_first.values[j][0] == 0:
                        ind_genes_assigned_first.append(j)
                    if self.enrichment_obj.missing_second.values[j][0] == 0:
                        ind_genes_assigned_second.append(j)
            ind_genes_assigned_first = np.asarray(ind_genes_assigned_first)
            ind_genes_assigned_second = np.asarray(ind_genes_assigned_second)

            if len(ind_genes_assigned_first) > 0:
                qvalues_enriched_terms_first.append(self.enrichment_obj.qvalues_first[ind_genes_assigned_first].flatten())
                genes_beta_signif_absolut_numbers_first.append(np.sum(np.logical_and(qvalues_enriched_terms_first[-1] > signif_threshold,
                                                                            qvalues_enriched_terms_first[
                                                                                -1] < empirical_quantile_first)))
                genes_beta_signif_percentage_first.append(genes_beta_signif_absolut_numbers_first[-1] / len(ind_genes_assigned_first))
            else:
                qvalues_enriched_terms_first.append(np.nan)
                genes_beta_signif_absolut_numbers_first.append(np.nan)
                genes_beta_signif_percentage_first.append(np.nan)

            if len(ind_genes_assigned_second) > 0:
                qvalues_enriched_terms_second.append(self.enrichment_obj.qvalues_second[ind_genes_assigned_second].flatten())
                genes_beta_signif_absolut_numbers_second.append(np.sum(np.logical_and(qvalues_enriched_terms_second[-1] > signif_threshold,
                                                                            qvalues_enriched_terms_second[
                                                                                -1] < empirical_quantile_second)))
                genes_beta_signif_percentage_second.append(genes_beta_signif_absolut_numbers_second[-1] / len(ind_genes_assigned_second))
            else:
                qvalues_enriched_terms_second.append(np.nan)
                genes_beta_signif_absolut_numbers_second.append(np.nan)
                genes_beta_signif_percentage_second.append(np.nan)


        # histogram of percentage genes beta significant
        plt.figure()
        plt.hist([genes_beta_signif_percentage_first, genes_beta_signif_percentage_second], color=['blue', 'red'], label=['genes', 'proteins'], bins=20)
        plt.xlabel('Percentage %s < genes < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
        plt.ylabel('Frequency', fontsize=18)
        plt.xlim([-0.05, 1.05])
        leg = plt.legend(loc='best', fancybox=True, fontsize=18)
        leg.get_frame().set_alpha(0.5)
        if output_filetag == '':
            plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
        plt.close()

        # histogram of number of genes beta significant
        plt.figure()
        plt.hist([genes_beta_signif_absolut_numbers_first, genes_beta_signif_absolut_numbers_second], color=['blue', 'red'], label=['genes', 'proteins'], bins=20)
        plt.xlabel('#genes %s < gene < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
        plt.ylabel('Frequency', fontsize=18)
        leg = plt.legend(loc='best', fancybox=True, fontsize=18)
        leg.get_frame().set_alpha(0.5)
        if output_filetag == '':
            plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
        plt.close()

        # scatter percentage of genes beta significant vs. probability
        plt.figure()
        plt.scatter(joana_output_enriched[col].values.flatten(), genes_beta_signif_percentage_first, c='blue', label='genes')
        plt.scatter(joana_output_enriched[col].values.flatten(), genes_beta_signif_percentage_second, c='red', label='proteins', alpha=0.7)
        plt.xlabel('Probability', fontsize=18)
        plt.xlim([0.45, 1.05])
        plt.ylabel('Percentage %s < genes < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
        plt.ylim([-0.05, 1.05])
        leg = plt.legend(loc='best', fancybox=True, fontsize=18)
        leg.get_frame().set_alpha(0.5)
        if output_filetag == '':
            plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_vs_probability_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_vs_probability_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
        plt.close()

        # scatter number of genes beta significant vs. probability
        plt.figure()
        plt.scatter(joana_output_enriched[col].values.flatten(), genes_beta_signif_absolut_numbers_first, c='blue', label='genes')
        plt.scatter(joana_output_enriched[col].values.flatten(), genes_beta_signif_absolut_numbers_second, c='red', label='proteins', alpha=0.7)
        plt.xlabel('Probability', fontsize=18)
        plt.xlim([0.45, 1.05])
        plt.ylabel('#genes %s < gene < Q95' % str(np.round(signif_threshold, 3)), fontsize=18)
        leg = plt.legend(loc='best', fancybox=True, fontsize=18)
        leg.get_frame().set_alpha(0.5)
        if output_filetag == '':
            plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_vs_probability_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(dir_output, 'absolut_beta_significant_genes_vs_probability_%s_%s_%s.pdf' % (str(np.round(signif_threshold, 3)), col, output_filetag)), bbox_inches='tight')
        plt.close()

    
    def plot_barplot(self, filename, threshold=0.5, fig_size=(8,6)):        
        if self.enrichment_obj.joana_output is None:
            print('Please run JOANA model at first or load a JOANA output.')
        if self.enrichment_obj.type == 'single-species':
            self.__plot_barplot_single_species(filename, threshold=threshold, fig_size=fig_size)
        elif self.enrichment_obj.type == 'cooperative':
            self.__plot_barplot_cooperative(filename, threshold=threshold, fig_size=fig_size)

    def __plot_barplot_single_species(self, filename, threshold=0.5, fig_size=(8,6)):
        joana_result = self.enrichment_obj.joana_output.sort_values(by='Single-Species', ascending=False)
        joana_result_filtered = joana_result[joana_result['Single-Species'] > threshold]
        joana_result_filtered.reset_index(inplace=True)

        plt.figure(figsize=fig_size)
        ax = sns.barplot(x='Single-Species', y='Terms', data=joana_result_filtered, color='forestgreen')
        plt.xlim([threshold, 1.])
        plt.savefig(filename, bbox_inches='tight')
        plt.close()

    def __plot_barplot_cooperative(self, filename, threshold=0.5, fig_size=(8,6)):
        joana_result = self.enrichment_obj.joana_output

        joana_cooperative = joana_result[joana_result['Cooperative'] > threshold]
        joana_cooperative_merge = pd.DataFrame({'Terms' : joana_cooperative.index.values,'Probability' : joana_cooperative['Cooperative'], 'Model' : ['Cooperative'] * joana_cooperative.shape[0]})
        joana_cooperative_merge = joana_cooperative_merge.sort_values(by='Probability', ascending=False).reset_index(drop=True)
        joana_single_species_1 = joana_result[joana_result['Single-Species-1'] > threshold]
        joana_single_species_1 = joana_single_species_1[joana_single_species_1['Cooperative'] < threshold]
        joana_single_species_1 = joana_single_species_1[joana_single_species_1['Single-Species-2'] < threshold]
        joana_single_species_1_merge = pd.DataFrame({'Terms' : joana_single_species_1.index.values,'Probability' : joana_single_species_1['Single-Species-1'], 'Model' : ['Single-Species-1'] * joana_single_species_1.shape[0]})
        joana_single_species_1_merge = joana_single_species_1_merge.sort_values(by='Probability', ascending=False).reset_index(drop=True)
        joana_single_species_2 = joana_result[joana_result['Single-Species-2'] > threshold]
        joana_single_species_2 = joana_single_species_2[joana_single_species_2['Cooperative'] < threshold]
        joana_single_species_2 = joana_single_species_2[joana_single_species_2['Single-Species-1'] < threshold]
        joana_single_species_2_merge = pd.DataFrame({'Terms' : joana_single_species_2.index.values,'Probability' : joana_single_species_2['Single-Species-2'], 'Model' : ['Single-Species-2'] * joana_single_species_2.shape[0]})
        joana_single_species_2_merge = joana_single_species_2_merge.sort_values(by='Probability', ascending=False).reset_index(drop=True)

        # iterative assignemnt of how many subplots we need based on number of active modeling results from the different species        
        active_models = int(joana_cooperative_merge.shape[0] > 0) + int(joana_single_species_1_merge.shape[0] > 0) + int(joana_single_species_2_merge.shape[0] > 0)
        if active_models == 3:
            subplots_number = 311
        elif active_models == 2:
            subplots_number = 211
        else:
            subplots_number = None

        color = {'Cooperative' : 'lime', 'Single-Species-1' : 'orange', 'Single-Species-2' : 'khaki'}
        fig = plt.figure(figsize=fig_size)
        if joana_cooperative_merge.shape[0] > 0:
            if not subplots_number is None:
                plt.subplot(subplots_number)
                subplots_number += 1
            ax = sns.barplot(y='Terms', x='Probability', hue='Model', data=joana_cooperative_merge, palette=color)
            plt.legend([],[], frameon=False)
            plt.xlabel('Probability', fontsize=18)
            plt.ylabel('Terms', fontsize=18)
            plt.title('Cooperative', fontsize=18)
            plt.xlim([threshold, 1.])
        if joana_single_species_1_merge.shape[0] > 0:
            if not subplots_number is None:
                plt.subplot(subplots_number)
                subplots_number += 1
            ax = sns.barplot(y='Terms', x='Probability', hue='Model', data=joana_single_species_1_merge, palette=color)
            plt.legend([],[], frameon=False)
            plt.xlabel('Probability', fontsize=18)
            plt.ylabel('Terms', fontsize=18)
            plt.title('Single-Species-1', fontsize=20)
            plt.xlim([threshold, 1.])
        if joana_single_species_2_merge.shape[0] > 0:
            if not subplots_number is None:
                plt.subplot(subplots_number)
                subplots_number += 1
            ax = sns.barplot(y='Terms', x='Probability', hue='Model', data=joana_single_species_2_merge, palette=color)
            plt.legend([],[], frameon=False)
            plt.xlabel('Probability', fontsize=18)
            plt.ylabel('Terms', fontsize=18)
            plt.title('Single-Species-2', fontsize=20)
            plt.xlim([threshold, 1.])
        plt.savefig(filename, bbox_inches='tight')
        plt.close()



    def plot(self, filename, top_percentile_edges=0, scaling_factor=5, fig_heigth=12, fig_width=15,
             bbox=(100, 100, 600, 600), layout='kk', verbose_n_top_terms=-1, second_order='detailed', n_top_terms=10):
        if self.enrichment_obj.joana_output is None:
            print('Please run JOANA model at first or load a JOANA output.')
        else:
            if self.enrichment_obj.type == 'single-species':
                self.__plot_single_species(filename, top_percentile_edges, scaling_factor, fig_heigth, fig_width,
                                          bbox, layout, verbose_n_top_terms, second_order, n_top_terms)
            elif self.enrichment_obj.type == 'cooperative':
                self.__plot_cooperative(filename, top_percentile_edges, scaling_factor, fig_heigth, fig_width,
                                          bbox, layout, verbose_n_top_terms, second_order, n_top_terms)



    def __plot_single_species(self, filename, top_percentile_edges=0, scaling_factor=5, fig_heigth=12, fig_width=15,
             bbox=(100, 100, 600, 600), layout='kk', verbose_n_top_terms=-1, second_order='detailed', n_top_terms=10):

        # get top n terms
        output_sorted = self.enrichment_obj.joana_output.sort_values('Single-Species', ascending=False)
        prob_nth_term = output_sorted['Single-Species'].values[n_top_terms]
        argmax_prob_nth_term = np.max(np.where(output_sorted['Single-Species'].values == prob_nth_term)[0])
        if n_top_terms < argmax_prob_nth_term:
            top_single_species = self.enrichment_obj.joana_output.sort_values('Single-Species', ascending=False).iloc[:argmax_prob_nth_term]
        else:
            top_single_species = self.enrichment_obj.joana_output.sort_values('Single-Species', ascending=False).iloc[:n_top_terms]
        ind_prob = np.where(top_single_species['Single-Species'].values > 0.5)[0]
        top_cooperative = top_single_species.iloc[ind_prob]

        prob_terms = top_cooperative
        prob_terms_unique = prob_terms.index.unique()
        ind_unique = np.asarray([np.where(prob_terms.index.values == x)[0][0] for x in prob_terms_unique])
        prob_terms = prob_terms.iloc[ind_unique]
        n_prob_terms = len(prob_terms_unique)

        # # get term size and term annotations
        # ind_terms = np.asarray([np.where(np.asarray(self.enrichment_obj.terms) == x)[0] for x in prob_terms_unique])
        # term_sizes = np.zeros((n_prob_terms,))
        # term_annotations = [list() for i in range(n_prob_terms)]
        # for i in range(len(self.enrichment_obj.assignment_matrix)):
        #     for j in range(n_prob_terms):
        #         if ind_terms[j] in self.enrichment_obj.assignment_matrix[i]:
        #             term_sizes[j] += 1
        #             term_annotations[j].append(i)

        # select terms based on ordering detailed or bigger picture with big terms
        prob_terms_size = pd.DataFrame.copy(prob_terms)
        # prob_terms_size['size'] = term_sizes
        if second_order == 'detailed':
            prob_terms_size = prob_terms_size.sort_values(by=['Single-Species', 'size'], ascending=[False, True])
        else:
            prob_terms_size = prob_terms_size.sort_values(by=['Single-Species', 'size'], ascending=[False, False])
        prob_terms = prob_terms_size[['Single-Species']].iloc[:n_top_terms]
        prob_terms_unique = prob_terms.index.unique()
        ind_unique = np.asarray([np.where(prob_terms.index.values == x)[0][0] for x in prob_terms_unique])
        prob_terms = prob_terms.iloc[ind_unique]
        n_prob_terms = len(prob_terms_unique)

        # how many labels to show on graph
        if verbose_n_top_terms == -1:
            verbose_n_top_terms = n_prob_terms

        # get term size and term annotations
        ind_terms = np.asarray([np.where(np.asarray(self.enrichment_obj.terms) == x)[0] for x in prob_terms_unique])
        term_sizes = np.zeros((n_prob_terms,))
        term_annotations = [list() for i in range(n_prob_terms)]
        for i in range(len(self.enrichment_obj.assignment_matrix)):
            for j in range(n_prob_terms):
                if ind_terms[j] in self.enrichment_obj.assignment_matrix[i]:
                    term_sizes[j] += 1
                    term_annotations[j].append(i)

        # compute similarity based on term overlap
        similarity_matrix_terms = np.zeros((n_prob_terms, n_prob_terms))
        for i in range(n_prob_terms):
            for j in range(n_prob_terms):
                if i != j:
                    similarity_matrix_terms[i, j] = len(np.intersect1d(term_annotations[i], term_annotations[j]))
        similarity_matrix_terms_log = np.log(1 + similarity_matrix_terms)

        # show only top percentile of connections
        median_similarity = np.percentile(similarity_matrix_terms_log.flatten(), top_percentile_edges)
        similarity_matrix_terms_log[similarity_matrix_terms_log < median_similarity] = 0

        # normalize weights between 0-1
        similarity_matrix_terms_log = (similarity_matrix_terms_log - np.min(similarity_matrix_terms_log)) / (
                    np.max(similarity_matrix_terms_log) - np.min(similarity_matrix_terms_log))
        similarity_matrix_terms_log *= scaling_factor

        # get the row, col indices of the non-zero elements in your adjacency matrix
        conn_indices = np.where(np.triu(similarity_matrix_terms_log))

        # get the weights corresponding to these indices
        weights = similarity_matrix_terms_log[conn_indices]

        # a sequence of (i, j) tuples, each corresponding to an edge from i -> j
        edges = zip(*conn_indices)

        # initialize the graph from the edge sequence
        G = Graph(edges=edges, directed=False)

        # assign node names and weights to be attributes of the vertices and edges
        # respectively
        G.es['weight'] = weights

        # I will also assign the weights to the 'width' attribute of the edges. this
        # means that igraph.plot will set the line thicknesses according to the edge
        # weights
        G.es['width'] = weights

        # vertex size based on term size
        term_sizes_logged = 1 + np.log(term_sizes) - np.min(np.log(term_sizes))
        # term_sizes_logged_scaled = 1 + (term_sizes_logged - np.min(term_sizes_logged))/(np.max(term_sizes_logged) - np.min(term_sizes_logged))
        min_term_sizes_scaled = np.percentile(term_sizes_logged, 25)
        factor_term_size = 15 / min_term_sizes_scaled
        final_vertex_sizes = 10 + term_sizes_logged * factor_term_size
        G.vs['size'] = final_vertex_sizes

        col = 'Single-Species'
        # only write significant terms
        temp_label = [''] * n_prob_terms
        for i in range(n_prob_terms):
            if prob_terms[col].values[i] > 0.5 and i < verbose_n_top_terms:
                temp_label[i] = self.enrichment_obj.terms.values[ind_terms[i][0]]

        G.vs['label'] = temp_label

        cmap = cm.get_cmap('Reds')
        color_vertex = cmap(prob_terms[col].values)
        G.vs['color'] = [tuple(color_vertex[i]) for i in range(n_prob_terms)]

        # Create the figure
        fig = plt.figure()
        fig.set_figheight(fig_heigth)
        fig.set_figwidth(fig_width)

        # Create a basic plot
        axes = fig.add_subplot(111)

        # Draw the graph over the plot
        # Two points to note here:
        # 1) we add the graph to the axes, not to the figure. This is because
        #    the axes are always drawn on top of everything in a matplotlib
        #    figure, and we want the graph to be on top of the axes.
        # 2) we set the z-order of the graph to infinity to ensure that it is
        #    drawn above all the curves drawn by the axes object itself.
        graph_artist = GraphArtist(G, bbox=bbox, layout=layout)
        graph_artist.set_zorder(float('inf'))
        axes.artists.append(graph_artist)
        axes.axis('off')

        # plot colorbar for probabilities
        img = axes.imshow(np.array([[0, 1]]), cmap=cmap)
        img.set_visible(False)
        cbar = plt.colorbar(img, ax=axes)
        cbar.ax.tick_params(labelsize=20)

        # Save the figure
        fig.savefig(filename.replace('.pdf', '_%d.pdf' % verbose_n_top_terms), bbox_inches='tight')
        plt.close()


    def __plot_cooperative(self, filename, top_percentile_edges=0, scaling_factor=5, fig_heigth=12, fig_width=15,
             bbox=(100, 100, 600, 600), layout='kk', verbose_n_top_terms=-1, second_order='detailed', n_top_terms=10):

        all_prob_terms = list()
        for col in ['Cooperative', 'Single-Species-1', 'Single-Species-2']:
            # get top n terms
            output_sorted = self.enrichment_obj.joana_output.sort_values(col, ascending=False)
            prob_nth_term = output_sorted[col].values[n_top_terms]
            argmax_prob_nth_term = np.max(np.where(output_sorted[col].values == prob_nth_term)[0])
            if n_top_terms < argmax_prob_nth_term:
                top_single_species = self.enrichment_obj.joana_output.sort_values(col, ascending=False).iloc[:argmax_prob_nth_term]
            else:
                top_single_species = self.enrichment_obj.joana_output.sort_values(col, ascending=False).iloc[:n_top_terms]
            ind_prob = np.where(top_single_species[col].values > 0.5)[0]
            top_cooperative = top_single_species.iloc[ind_prob]

            prob_terms = top_cooperative
            prob_terms_unique = prob_terms.index.unique()
            ind_unique = np.asarray([np.where(prob_terms.index.values == x)[0][0] for x in prob_terms_unique])
            prob_terms = prob_terms.iloc[ind_unique]
            n_prob_terms = len(prob_terms_unique)

            # get term size and term annotations
            ind_terms = np.asarray([np.where(np.asarray(self.enrichment_obj.terms) == x)[0] for x in prob_terms_unique])
            term_sizes = np.zeros((n_prob_terms,))
            term_annotations = [list() for i in range(n_prob_terms)]
            for i in range(len(self.enrichment_obj.assignment_matrix)):
                for j in range(n_prob_terms):
                    if ind_terms[j] in self.enrichment_obj.assignment_matrix[i]:
                        term_sizes[j] += 1
                        term_annotations[j].append(i)

            # select terms based on ordering detailed or bigger picture with big terms
            prob_terms_size = pd.DataFrame.copy(prob_terms)
            prob_terms_size['size'] = term_sizes
            if second_order == 'detailed':
                prob_terms_size = prob_terms_size.sort_values(by=[col, 'size'], ascending=[False, True])
            else:
                prob_terms_size = prob_terms_size.sort_values(by=[col, 'size'], ascending=[False, False])
            prob_terms = prob_terms_size.iloc[:n_top_terms]
            prob_terms_unique = prob_terms.index.unique()
            ind_unique = np.asarray([np.where(prob_terms.index.values == x)[0][0] for x in prob_terms_unique])
            prob_terms = prob_terms.iloc[ind_unique]
            n_prob_terms = len(prob_terms_unique)

            all_prob_terms.append(prob_terms)

        prob_terms = pd.concat(all_prob_terms)
        prob_terms_unique = prob_terms.index.unique()
        ind_prob_terms_unique = list()
        for term in prob_terms_unique:
            temp_ind_term = np.where(prob_terms.index.values == term)[0]
            if len(temp_ind_term) > 1:
                ind_prob_terms_unique.append(np.asarray([temp_ind_term[0]]))
            else:
                ind_prob_terms_unique.append(temp_ind_term)
        prob_terms = prob_terms.iloc[np.concatenate(ind_prob_terms_unique)]
        n_prob_terms = len(prob_terms_unique)

        # # get term size and term annotations
        # ind_terms = np.asarray([np.where(np.asarray(self.enrichment_obj.terms) == x)[0] for x in prob_terms_unique])
        # term_sizes = np.zeros((n_prob_terms,))
        # term_annotations = [list() for i in range(n_prob_terms)]
        # for i in range(len(self.enrichment_obj.assignment_matrix)):
        #     for j in range(n_prob_terms):
        #         if ind_terms[j] in self.enrichment_obj.assignment_matrix[i]:
        #             term_sizes[j] += 1
        #             term_annotations[j].append(i)

        # select terms based on ordering detailed or bigger picture with big terms
        # prob_terms_size = pd.DataFrame.copy(prob_terms)
        # # prob_terms_size['size'] = term_sizes
        # if second_order == 'detailed':
        #     prob_terms_size = prob_terms_size.sort_values(by=['Cooperative', 'size'], ascending=[False, True])
        # else:
        #     prob_terms_size = prob_terms_size.sort_values(by=['Cooperative', 'size'], ascending=[False, False])
        # prob_terms = prob_terms_size.iloc[:n_top_terms]
        # prob_terms_unique = prob_terms.index.unique()
        # ind_unique = np.asarray([np.where(prob_terms.index.values == x)[0][0] for x in prob_terms_unique])
        # prob_terms = prob_terms.iloc[ind_unique]
        # n_prob_terms = len(prob_terms_unique)

        # how many labels to show on graph
        if verbose_n_top_terms == -1:
            verbose_n_top_terms = n_prob_terms

        # get term size and term annotations
        ind_terms = np.asarray([np.where(np.asarray(self.enrichment_obj.terms) == x)[0] for x in prob_terms_unique])
        term_annotations = [list() for i in range(n_prob_terms)]
        for i in range(len(self.enrichment_obj.assignment_matrix)):
            for j in range(n_prob_terms):
                if ind_terms[j] in self.enrichment_obj.assignment_matrix[i]:
                    term_annotations[j].append(i)

        # compute similarity based on term overlap
        similarity_matrix_terms = np.zeros((n_prob_terms, n_prob_terms))
        for i in range(n_prob_terms):
            for j in range(n_prob_terms):
                if i != j:
                    similarity_matrix_terms[i, j] = len(np.intersect1d(term_annotations[i], term_annotations[j]))
        similarity_matrix_terms_log = np.log(1 + similarity_matrix_terms)

        # show only top percentile of connections
        median_similarity = np.percentile(similarity_matrix_terms_log.flatten(), top_percentile_edges)
        similarity_matrix_terms_log[similarity_matrix_terms_log < median_similarity] = 0

        # normalize weights between 0-1
        similarity_matrix_terms_log = (similarity_matrix_terms_log - np.min(similarity_matrix_terms_log)) / (
                    np.max(similarity_matrix_terms_log) - np.min(similarity_matrix_terms_log))
        similarity_matrix_terms_log *= scaling_factor

        # get the row, col indices of the non-zero elements in your adjacency matrix
        conn_indices = np.where(np.triu(similarity_matrix_terms_log))

        # get the weights corresponding to these indices
        weights = similarity_matrix_terms_log[conn_indices]

        # a sequence of (i, j) tuples, each corresponding to an edge from i -> j
        edges = zip(*conn_indices)

        # initialize the graph from the edge sequence
        G = Graph(edges=edges, directed=False)

        # assign node names and weights to be attributes of the vertices and edges
        # respectively
        G.es['weight'] = weights

        # I will also assign the weights to the 'width' attribute of the edges. this
        # means that igraph.plot will set the line thicknesses according to the edge
        # weights
        G.es['width'] = weights

        # vertex size based on term size
        term_sizes_logged = 1 + np.log(term_sizes) - np.min(np.log(term_sizes))
        # term_sizes_logged_scaled = 1 + (term_sizes_logged - np.min(term_sizes_logged))/(np.max(term_sizes_logged) - np.min(term_sizes_logged))
        min_term_sizes_scaled = np.percentile(term_sizes_logged, 25)
        factor_term_size = 15 / min_term_sizes_scaled
        final_vertex_sizes = 10 + term_sizes_logged * factor_term_size
        G.vs['size'] = final_vertex_sizes

        for col in ['Cooperative', 'Single-Species-1', 'Single-Species-2']:
            # only write significant terms
            temp_label = [''] * n_prob_terms
            counter = 0
            for i in range(n_prob_terms):
                if prob_terms[col].values[i] > 0.5 and counter < verbose_n_top_terms:
                    temp_label[i] = self.enrichment_obj.terms.values[ind_terms[i][0]]
                    counter += 1

            G.vs['label'] = temp_label

            cmap = cm.get_cmap('Reds')
            color_vertex = cmap(prob_terms[col].values)
            G.vs['color'] = [tuple(color_vertex[i]) for i in range(n_prob_terms)]

            # Create the figure
            fig = plt.figure()
            fig.set_figheight(fig_heigth)
            fig.set_figwidth(fig_width)

            # Create a basic plot
            axes = fig.add_subplot(111)

            # Draw the graph over the plot
            # Two points to note here:
            # 1) we add the graph to the axes, not to the figure. This is because
            #    the axes are always drawn on top of everything in a matplotlib
            #    figure, and we want the graph to be on top of the axes.
            # 2) we set the z-order of the graph to infinity to ensure that it is
            #    drawn above all the curves drawn by the axes object itself.
            graph_artist = GraphArtist(G, bbox=bbox, layout=layout)
            graph_artist.set_zorder(float('inf'))
            axes.artists.append(graph_artist)
            axes.axis('off')

            # plot colorbar for probabilities
            img = axes.imshow(np.array([[0, 1]]), cmap=cmap)
            img.set_visible(False)
            cbar = plt.colorbar(img, ax=axes)
            cbar.ax.tick_params(labelsize=20)

            # Save the figure
            fig.savefig(filename.replace('.pdf', '_%d_%s.pdf' % (n_top_terms, col)), bbox_inches='tight')
            plt.close()


class GraphArtist(Artist):
    """Matplotlib artist class that draws igraph graphs.

    Only Cairo-based backends are supported.
    """

    def __init__(self, graph, bbox, palette=None, *args, **kwds):
        """Constructs a graph artist that draws the given graph within
        the given bounding box.

        `graph` must be an instance of `igraph.Graph`.
        `bbox` must either be an instance of `igraph.drawing.BoundingBox`
        or a 4-tuple (`left`, `top`, `width`, `height`). The tuple
        will be passed on to the constructor of `BoundingBox`.
        `palette` is an igraph palette that is used to transform
        numeric color IDs to RGB values. If `None`, a default grayscale
        palette is used from igraph.

        All the remaining positional and keyword arguments are passed
        on intact to `igraph.Graph.__plot__`.
        """
        Artist.__init__(self)

        if not isinstance(graph, Graph):
            raise TypeError("expected igraph.Graph, got %r" % type(graph))

        self.graph = graph
        self.palette = palette or palettes["gray"]
        self.bbox = BoundingBox(bbox)
        self.args = args
        self.kwds = kwds

    def draw(self, renderer):
        from matplotlib.backends.backend_cairo import RendererCairo
        if not isinstance(renderer, RendererCairo):
            raise TypeError("graph plotting is supported only on Cairo backends")
        self.graph.__plot__(renderer.gc.ctx, self.bbox, self.palette, *self.args, **self.kwds)
