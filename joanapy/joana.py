from numpy.lib import index_tricks
from numpy.lib.function_base import quantile
from pandas.io.stata import ValueLabelTypeMismatch
from seaborn.palettes import color_palette
from joanapy.moment_fitting_core import MOMENT_FITTING
from joanapy.enrichment_obj import ENRICHMENT_OBJ
from joanapy.joana_helper import load_assignmentMatrix
import networkx as nx
import os
from subprocess import call
import pandas as pd
import numpy as np
import psutil
import sys
import matplotlib
#matplotlib.use("Cairo")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.artist import Artist
from matplotlib import gridspec
import seaborn as sns
import time
import igraph as ig
from matplotlib.cm import get_cmap
from igraph import Graph, plot
from igraph.drawing import BoundingBox
from igraph.drawing.colors import palettes
from matplotlib.colors import LinearSegmentedColormap


FILEPATH_JOANA_CLASS = os.path.realpath(__file__)
FILENAME_JOANA = os.path.join(os.path.dirname(FILEPATH_JOANA_CLASS), 'joana_app', 'MonaConsoleApp.exe')


class JOANA:
    
    def __init__(self, enrichment_obj: ENRICHMENT_OBJ):
        self.enrichment_obj = enrichment_obj


    ########################################################################################################################################
    ################################################################# RUN JOANA ############################################################
    ########################################################################################################################################

    def run(self, filename_output, filename_moment_fit_first=None, filename_moment_fit_second=None, tolerance_fitting=1E-5, steps_fitting=30000, verbose=True, goodness_of_fit=False, plot_components=False, prior_pA=1, min_term_size=0, max_term_size=100000, save_enrichment_obj=True, init='moment', second_comp_uniform=False):
        
        if not filename_output.endswith('.csv'):
            raise ValueError('The output filename needs to be a .csv file.')

        self.enrichment_obj.filename_joana_output = filename_output
        filename_output_ending = filename_output.split('.')[-1]
        dir_output = os.path.dirname(filename_output)
        if not os.path.exists(dir_output):
            os.makedirs(dir_output)

        # get term sizes
        term_sizes = np.zeros((len(self.enrichment_obj.terms), ))
        #for i in range(len(self.enrichment_obj.assignment_matrix)):
        #    term_sizes[np.asarray(self.enrichment_obj.assignment_matrix[i])] += 1
        #pd.DataFrame({'terms' : self.enrichment_obj.terms.values.flatten(), 'size' : term_sizes}).to_csv(self.enrichment_obj.filename_terms.replace('.txt', '_sizes.txt'), index=False, header=True, sep='&')
        term_sizes=self.enrichment_obj.term_size
        pd.DataFrame({'terms' : self.enrichment_obj.terms.values.flatten(), 'size' : term_sizes}).to_csv(self.enrichment_obj.filename_terms.replace('.txt', '_sizes.txt'), index=False, header=True, sep='&')

        #eq=np.array_equal(term_size2, term_sizes)
        # filter terms
        ind_keep = np.intersect1d(np.where(min_term_size <= term_sizes)[0], np.where(term_sizes <= max_term_size)[0])
        index_mapping = dict()
        print(range(len(ind_keep)))
        for i in range(len(ind_keep)):
            index_mapping[ind_keep[i]] = i
        
        # filter assignment matrix
        assignment_matrix_filtered = list()
        print(range(len(self.enrichment_obj.assignment_matrix)))
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
        
        # reload assignment matrix in enrichment objed
        self.enrichment_obj.assignment_matrix = load_assignmentMatrix(self.enrichment_obj.filename_assignment_matrix)

        # filter terms
        terms_filtered = self.enrichment_obj.terms.iloc[ind_keep]
        self.enrichment_obj.filename_terms = self.enrichment_obj.filename_terms.replace('.txt', '_filtered_min%d_max%d.txt' %(min_term_size, max_term_size))
        terms_filtered.to_csv(self.enrichment_obj.filename_terms, index=False, header=False)
        self.enrichment_obj.terms = terms_filtered

        start_time = time.time()

        if self.enrichment_obj.type == 'single-species':
            if filename_moment_fit_first is None:
                # first species
                filename_moment_fit_first = os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting.txt')
                moment_fitting_first = MOMENT_FITTING(self.enrichment_obj.qvalues_first,
                                                    filename_moment_fit_first,
                                                    tolerance=tolerance_fitting, steps=steps_fitting)
                moment_fitting_first.run(verbose=verbose, init=init, second_comp_uniform=second_comp_uniform)
                if goodness_of_fit:
                    moment_fitting_first.goodness_of_fit(os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting_gof.txt'), plot_histograms=True)
                if plot_components:
                    moment_fitting_first.plot_components_density(os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting_components.pdf'))
                    moment_fitting_first.plot_mixture_pdfs(os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting_mixture_pdf.pdf'))
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
                moment_fitting_first.run(verbose=verbose, init=init, second_comp_uniform=second_comp_uniform)
                if goodness_of_fit:
                    moment_fitting_first.goodness_of_fit(os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting_gof.txt'), plot_histograms=True)
                if plot_components:
                    moment_fitting_first.plot_components_density(os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting_components.pdf'))
                    moment_fitting_first.plot_mixture_pdfs(os.path.join(os.path.dirname(filename_output),'qvalues_first_moment_fitting_mixture_pdf.pdf'))
                self.enrichment_obj.moment_fit_first = moment_fitting_first

            if filename_moment_fit_second is None:
                # second species
                filename_moment_fit_second = os.path.join(os.path.dirname(filename_output), 'qvalues_second_moment_fitting.txt')
                ind_not_missing_second = self.enrichment_obj.missing_second.values == 0
                moment_fitting_second = MOMENT_FITTING(self.enrichment_obj.qvalues_second[ind_not_missing_second.flatten()], filename_moment_fit_second,
                                                    tolerance=tolerance_fitting, steps=steps_fitting)
                moment_fitting_second.run(verbose=verbose, init=init, second_comp_uniform=second_comp_uniform)
                if goodness_of_fit:
                    moment_fitting_second.goodness_of_fit(os.path.join(os.path.dirname(filename_output),'qvalues_second_moment_fitting_gof.txt'), plot_histograms=True)
                if plot_components:
                    moment_fitting_second.plot_components_density(os.path.join(os.path.dirname(filename_output),'qvalues_second_moment_fitting_components.pdf'))
                    moment_fitting_second.plot_mixture_pdfs(os.path.join(os.path.dirname(filename_output),'qvalues_second_moment_fitting_mixture_pdf.pdf'))
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

        time_needed = time.time() - start_time

        # load joana output in enrichment object
        if os.path.exists(filename_output):
            self.enrichment_obj.joana_output = pd.read_csv(filename_output, index_col='Terms')

        # save enrichment_obj after JOANA run
        if save_enrichment_obj:
            self.enrichment_obj.save(filename_output.replace(filename_output_ending, 'pickle'))

        return time_needed
    



    ########################################################################################################################################
    ################################################################# RUN JOANA DIRECTIVE ##################################################
    ########################################################################################################################################

    def run_directive(self, dir_output, statistical_direction_first, statistical_direction_second=None, filename_moment_fit_first=None, filename_moment_fit_second=None, tolerance_fitting=1E-5, steps_fitting=30000, verbose=True, goodness_of_fit=False, plot_components=False, prior_pA=1, min_term_size=0, max_term_size=100000, save_enrichment_obj=True, init='moment', second_comp_uniform=False):
        
        if self.enrichment_obj.type == 'cooperative' and statistical_direction_second is None:
            raise ValueError('Please pass a statistical direction for the second species.')

        # store in enrichment object
        # setattr(self.enrichment_obj, 'statistical_direction', statistical_direction)
        self.enrichment_obj.statistical_direction_first = statistical_direction_first
        ind_upregulated_first = np.where(statistical_direction_first > 0)[0]
        ind_downregulated_first = np.where(statistical_direction_first < 0)[0]

        if not statistical_direction_second is None:
            self.enrichment_obj.statistical_direction_second = statistical_direction_second
            ind_upregulated_second = np.where(statistical_direction_second > 0)[0]
            ind_downregulated_second = np.where(statistical_direction_second < 0)[0]

        # create output directories for both directions
        if not os.path.exists(dir_output):
            os.makedirs(dir_output)

        dir_output_up = os.path.join(dir_output, 'up_regulated')
        if not os.path.exists(dir_output_up):
            os.makedirs(dir_output_up)

        dir_output_down = os.path.join(dir_output, 'down_regulated')
        if not os.path.exists(dir_output_down):
            os.makedirs(dir_output_down)


        # get term sizes
        #term_sizes = np.zeros((len(self.enrichment_obj.terms), ))
        #for i in range(len(self.enrichment_obj.assignment_matrix)):
        #    term_sizes[np.asarray(self.enrichment_obj.assignment_matrix[i])] += 1
        term_sizes=self.enrichment_obj.term_size
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
            
        self.enrichment_obj.assignment_matrix = load_assignmentMatrix(self.enrichment_obj.filename_assignment_matrix)

        # filter terms
        terms_filtered = self.enrichment_obj.terms.iloc[ind_keep]
        self.enrichment_obj.filename_terms = self.enrichment_obj.filename_terms.replace('.txt', '_filtered_min%d_max%d.txt' %(min_term_size, max_term_size))
        terms_filtered.to_csv(self.enrichment_obj.filename_terms, index=False, header=False)
        self.enrichment_obj.terms = terms_filtered

        start_time = time.time()
        #######################################################
        ################## Moment fitting #####################
        #######################################################

        if filename_moment_fit_first is None:
            # first species
            filename_moment_fit_first = os.path.join(dir_output,'qvalues_first_moment_fitting.txt')
            moment_fitting_first = MOMENT_FITTING(self.enrichment_obj.qvalues_first, filename_moment_fit_first, tolerance=tolerance_fitting, steps=steps_fitting)
            moment_fitting_first.run(verbose=verbose, init=init, second_comp_uniform=second_comp_uniform)
            if goodness_of_fit:
                moment_fitting_first.goodness_of_fit(os.path.join(dir_output,'qvalues_first_moment_fitting_gof.txt'), plot_histograms=True)
            if plot_components:
                moment_fitting_first.plot_components_density(os.path.join(dir_output,'qvalues_first_moment_fitting_components.pdf'))
                moment_fitting_first.plot_mixture_pdfs(os.path.join(dir_output,'qvalues_first_moment_fitting_mixture_pdf.pdf'))
            self.enrichment_obj.moment_fit_first = moment_fitting_first

        if self.enrichment_obj.type == 'cooperative':
            if filename_moment_fit_second is None:
                # second species
                filename_moment_fit_second = os.path.join(dir_output, 'qvalues_second_moment_fitting.txt')
                ind_not_missing_second = self.enrichment_obj.missing_second.values == 0
                moment_fitting_second = MOMENT_FITTING(self.enrichment_obj.qvalues_second[ind_not_missing_second.flatten()], filename_moment_fit_second, tolerance=tolerance_fitting, steps=steps_fitting)
                moment_fitting_second.run(verbose=verbose, init=init, second_comp_uniform=second_comp_uniform)
                if goodness_of_fit:
                    moment_fitting_second.goodness_of_fit(os.path.join(dir_output,'qvalues_second_moment_fitting_gof.txt'), plot_histograms=True)
                if plot_components:
                    moment_fitting_second.plot_components_density(os.path.join(dir_output,'qvalues_second_moment_fitting_components.pdf'))
                    moment_fitting_second.plot_mixture_pdfs(os.path.join(dir_output,'qvalues_second_moment_fitting_mixture_pdf.pdf'))
                self.enrichment_obj.moment_fit_second = moment_fitting_second


        #######################################################
        ################### UPREGULATED #######################
        #######################################################

        # create assignment matrix for upregulated genes
        basename_assignment = os.path.basename(self.enrichment_obj.filename_assignment_matrix)
        with open(os.path.join(dir_output_up, basename_assignment), "w") as myfile:
            for i, line in enumerate(assignment_matrix_filtered):
                if i in ind_upregulated_first:
                    myfile.write(str(line) + '\n')
        

        # create qvalues file with upregulated qvalues for first species
        filename_qvalues_first_upregulated = os.path.join(dir_output_up, 'qvalues_first_upregulated.txt')
        pd.DataFrame(self.enrichment_obj.qvalues_first[ind_upregulated_first]).to_csv(filename_qvalues_first_upregulated, header=False, index=False)


        # prepare files for UPREGULATED
        if self.enrichment_obj.type == 'single-species':
            
            while self.__has_handle(FILENAME_JOANA):
                continue

            filename_output_upregulated = os.path.join(dir_output_up, 'JOANA_single-species_upregulated.csv')
            call(['mono',
                  FILENAME_JOANA,
                  '1',
                  os.path.join(dir_output_up, basename_assignment),
                  self.enrichment_obj.filename_terms,
                  filename_qvalues_first_upregulated,
                  filename_moment_fit_first,
                  filename_output_upregulated,
                  str(int(prior_pA))])


        if self.enrichment_obj.type == 'cooperative':
            # find qvalues from second species which are overlapping upregulated with the one from first species
            qvalues_second_upregulated = np.ones((len(ind_upregulated_first), ))
            missing_second_upregulated = np.zeros((len(ind_upregulated_first), ), dtype=int)
            missing_first_upregulated = np.zeros((len(ind_upregulated_first), ), dtype=int)
            counter = 0
            ind_second_not_missing = list()
            for i in range(len(self.enrichment_obj.qvalues_first)):
                if i in ind_upregulated_first and i in ind_upregulated_second:
                    qvalues_second_upregulated[counter] = self.enrichment_obj.qvalues_second[i]
                    counter += 1
                    ind_second_not_missing.append(i)
                elif i in ind_upregulated_first and not i in ind_upregulated_second:
                    missing_second_upregulated[counter] = 1
                    counter += 1
            filename_qvalues_second_upregulated = os.path.join(dir_output_up, 'qvalues_second_upregulated.txt')
            pd.DataFrame(qvalues_second_upregulated).to_csv(filename_qvalues_second_upregulated, header=False, index=False)
            ind_second_not_missing = np.asarray(ind_second_not_missing)
            filename_missing_second_upregulated = os.path.join(dir_output_up, 'missing_second_upregulated.txt')
            pd.DataFrame(missing_second_upregulated).to_csv(filename_missing_second_upregulated, header=False, index=False)
            filename_missing_first_upregulated = os.path.join(dir_output_up, 'missing_first_upregulated.txt')
            pd.DataFrame(missing_first_upregulated).to_csv(filename_missing_first_upregulated, header=False, index=False)

            while self.__has_handle(FILENAME_JOANA):
                continue
            
            filename_output_upregulated = os.path.join(dir_output_up, 'JOANA_cooperative-single-species_upregulated.csv')
            call(['mono',
                FILENAME_JOANA,
                '0',
                os.path.join(dir_output_up, basename_assignment),
                self.enrichment_obj.filename_terms,
                filename_qvalues_first_upregulated,
                filename_qvalues_second_upregulated,
                filename_missing_first_upregulated,
                filename_missing_second_upregulated,
                filename_moment_fit_first,
                filename_moment_fit_second,
                filename_output_upregulated,
                str(int(prior_pA))])

        # load joana output in enrichment object
        if os.path.exists(filename_output_upregulated):
            self.enrichment_obj.joana_output_upregulated = pd.read_csv(filename_output_upregulated, index_col='Terms')

        #######################################################
        ################## DOWNREGULATED ######################
        #######################################################

        # define the standard filename for terms, if we filter for upregulated terms this will we reassigned
        filename_terms_downregulation = self.enrichment_obj.filename_terms
        # filter terms which are predicted to be upreagulated from assignment matrix
        if not self.enrichment_obj.joana_output_upregulated is None:
            joana_output_up = self.enrichment_obj.joana_output_upregulated.copy()

            # filter results data frame for active terms
            if self.enrichment_obj.type == 'single-species':
                joana_output_up = joana_output_up[joana_output_up['Single-Species'] > 0.5]
            if self.enrichment_obj.type == 'cooperative':
                joana_output_up_cooperative = joana_output_up[joana_output_up['Cooperative'] > 0.5]
                joana_output_up_single_species_1 = joana_output_up[joana_output_up['Single-Species-1'] > 0.5]
                joana_output_up_single_species_2 = joana_output_up[joana_output_up['Single-Species-2'] > 0.5]
                joana_output_up = pd.concat([joana_output_up_cooperative, joana_output_up_single_species_1, joana_output_up_single_species_2])
            terms_to_filter_out = np.unique(joana_output_up.index.values)
            
            ind_keep = np.asarray([i for i in range(len(self.enrichment_obj.terms)) if not self.enrichment_obj.terms.values[i] in terms_to_filter_out])
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

            # filter terms
            terms_filtered = self.enrichment_obj.terms.iloc[ind_keep]
            filename_terms_downregulation = os.path.join(dir_output_down, os.path.basename(self.enrichment_obj.filename_terms).replace('.txt', '_filtered_upregulated.txt'))
            terms_filtered.to_csv(filename_terms_downregulation, index=False, header=False)



        # create assignment matrix for downregulated genes
        basename_assignment = os.path.basename(self.enrichment_obj.filename_assignment_matrix)
        with open(os.path.join(dir_output_down, basename_assignment), "w") as myfile:
            for i, line in enumerate(assignment_matrix_filtered):
                if i in ind_downregulated_first:
                    myfile.write(str(line) + '\n')
        

        # create qvalues file with downregulated qvalues for first species
        filename_qvalues_first_downregulated = os.path.join(dir_output_down, 'qvalues_first_downregulated.txt')
        pd.DataFrame(self.enrichment_obj.qvalues_first[ind_downregulated_first]).to_csv(filename_qvalues_first_downregulated, header=False, index=False)

        # prepare files for DOWNREGULATED
        if self.enrichment_obj.type == 'single-species':

            while self.__has_handle(FILENAME_JOANA):
                continue
            
            filename_output_downregulated = os.path.join(dir_output_down, 'JOANA_single-species_downregulated.csv')
            call(['mono',
                  FILENAME_JOANA,
                  '1',
                  os.path.join(dir_output_down, basename_assignment),
                  filename_terms_downregulation,
                  filename_qvalues_first_downregulated,
                  filename_moment_fit_first,
                  filename_output_downregulated,
                  str(int(prior_pA))])


        if self.enrichment_obj.type == 'cooperative':
            # find qvalues from second species which are overlapping downregulated with the one from first species
            qvalues_second_downregulated = np.ones((len(ind_downregulated_first), ))
            missing_second_downregulated = np.zeros((len(ind_downregulated_first), ), dtype=int)
            missing_first_downregulated = np.zeros((len(ind_downregulated_first), ), dtype=int)
            counter = 0
            ind_second_not_missing = list()
            for i in range(len(self.enrichment_obj.qvalues_first)):
                if i in ind_downregulated_first and i in ind_downregulated_second:
                    qvalues_second_downregulated[counter] = self.enrichment_obj.qvalues_second[i]
                    counter += 1
                    ind_second_not_missing.append(i)
                elif i in ind_downregulated_first and not i in ind_downregulated_second:
                    missing_second_downregulated[counter] = 1
                    counter += 1
            filename_qvalues_second_downregulated = os.path.join(dir_output_down, 'qvalues_second_downregulated.txt')
            pd.DataFrame(qvalues_second_downregulated).to_csv(filename_qvalues_second_downregulated, header=False, index=False)
            ind_second_not_missing = np.asarray(ind_second_not_missing)
            filename_missing_second_downregulated = os.path.join(dir_output_down, 'missing_second_downregulated.txt')
            pd.DataFrame(missing_second_downregulated).to_csv(filename_missing_second_downregulated, header=False, index=False)
            filename_missing_first_downregulated = os.path.join(dir_output_down, 'missing_first_downregulated.txt')
            pd.DataFrame(missing_first_downregulated).to_csv(filename_missing_first_downregulated, header=False, index=False)
            
            while self.__has_handle(FILENAME_JOANA):
                continue

            filename_output_downregulated = os.path.join(dir_output_down, 'JOANA_cooperative-single-species_downregulated.csv')
            call(['mono',
                FILENAME_JOANA,
                '0',
                os.path.join(dir_output_down, basename_assignment),
                filename_terms_downregulation,
                filename_qvalues_first_downregulated,
                filename_qvalues_second_downregulated,
                filename_missing_first_downregulated,
                filename_missing_second_downregulated,
                filename_moment_fit_first,
                filename_moment_fit_second,
                filename_output_downregulated,
                str(int(prior_pA))])
        
        time_needed = time.time() - start_time

        # load joana output in enrichment object
        if os.path.exists(filename_output_downregulated):
            self.enrichment_obj.joana_output_downregulated = pd.read_csv(filename_output_downregulated, index_col='Terms')

        # save enrichment_obj after JOANA run
        if save_enrichment_obj:
            self.enrichment_obj.save(os.path.join(dir_output, 'enrichment_obj_up_down_min%d_max%d.pickle' %(min_term_size, max_term_size)))

        return time_needed

        

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





    ########################################################################################################################################
    ############################################# HIDDEN SIGNIFICANT OF JOANA RESULTS ######################################################
    ########################################################################################################################################


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
        ind_terms = np.asarray([np.where(np.asarray(self.enrichment_obj.terms) == x)[0] for x in joana_output_enriched.index])
        #term_annotations = [list() for i in range(joana_output_enriched.shape[0])]
        term_sizes=self.enrichment_obj.term_size[ind_terms.flatten()]
        term_annotations = [self.enrichment_obj.term_annotation[i] for i in ind_terms.flatten()]
        qvalues_enriched_terms2 = [self.enrichment_obj.qvalues_first[i] for i in term_annotations]
        genes_beta_signif_absolut_numbers2 = [np.sum(np.logical_and(x > signif_threshold, x < empirical_quantile)) for x in  qvalues_enriched_terms2]
        genes_beta_signif_percentage2=genes_beta_signif_absolut_numbers2/term_sizes
        genes_beta_signif_percentage2= np.asarray(genes_beta_signif_percentage2)*100
        

        if genes_beta_signif_percentage2.size==0:
            print( "Notice: Skipping plot for non-significant active genes because the number of 'non-significant active (hidden-active) genes' is zero")
            return  # Exit the function early
        
        
        plt.figure()
        sns.swarmplot(data=genes_beta_signif_percentage2)

        plt.xlabel('Single-Species')
        plt.ylabel('percentages')
        plt.title('Percentages of non-significant active genes in Enriched pathways')
        plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%sBeesWarm.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        plt.close()

        

        # Set labels and title
        plt.figure()
        plt.boxplot(genes_beta_signif_percentage2)
        plt.xlabel('Single-Species')
        plt.ylabel('Percentages')
        plt.title('Percentages of non-significant active genes in Enriched pathways')
        plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%sBoxplot.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        plt.close()

        # Show the plot
        


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
            genes_beta_signif_percentage = np.asarray(genes_beta_signif_percentage) * 100
            
            if genes_beta_signif_percentage.size==0:
                print( "Notice: Skipping plot for non-significant active genes in %s because the number of 'non-significant active (hidden-active) genes' is zero" % str(col))
                continue  # Exit the function early
            
            plt.figure()
            sns.swarmplot(data=genes_beta_signif_percentage)

            plt.xlabel(col)
            plt.ylabel('percentages')
            plt.title('Percentages of non-significant active genes in Enriched pathways')
            plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%sBeesWarm.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
            plt.close()

            

            # Set labels and title
            plt.figure()
            plt.boxplot(genes_beta_signif_percentage)
            plt.xlabel(col)
            plt.ylabel('Percentages')
            plt.title('Percentages of non-significant active genes in Enriched pathways')
            plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%sBoxplot.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
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

        genes_beta_signif_percentage_first=np.asarray(genes_beta_signif_percentage_first) * 100
        genes_beta_signif_percentage_second=np.asarray(genes_beta_signif_percentage_second) * 100

        if genes_beta_signif_percentage_first.size == 0 and genes_beta_signif_percentage_second.size == 0:
                print("Notice: Skipping plot for non-significant active genes because the number of 'non-significant active (hidden-active) genes' for both modalities is zero.")
                return  # Exit the function early
        
        plt.figure()
        sns.swarmplot(data=[genes_beta_signif_percentage_first,genes_beta_signif_percentage_second],palette=['blue', 'red'], label=['Single-Species-1', 'Single-Species-2'])
        plt.legend()

        plt.xlabel(col)
        plt.ylabel('percentages')
        plt.title('Percentages of non-significant active genes in Enriched pathways')
        plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%sBeesWarm.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        plt.close()



        plt.boxplot(genes_beta_signif_percentage_first, 
        positions=[1],  # Position of the boxplot on the x-axis
        patch_artist=True,  # Fill with color
        boxprops=dict(facecolor='blue', color='black'),  # Customize box color
        whiskerprops=dict(color='black'),  # Customize whisker color
        capprops=dict(color='black'),  # Customize cap color
        medianprops=dict(color='black'))  # Customize median color

        # Create a boxplot for the second dataset (red)
        plt.boxplot(genes_beta_signif_percentage_second, 
        positions=[2],  # Position of the boxplot on the x-axis
        patch_artist=True,  # Fill with color
        boxprops=dict(facecolor='red', color='black'),  # Customize box color
        whiskerprops=dict(color='black'),  # Customize whisker color
        capprops=dict(color='black'),  # Customize cap color
        medianprops=dict(color='black'))  # Customize median color

        plt.xlabel('Groups')
        plt.ylabel('Values')
        plt.title('Percentages of non-significant active genes in Enriched pathways')

        # Set x-axis tick labels
        plt.xticks([1, 2], ['Single-Species-1', 'Single-Species-2'])
        plt.savefig(os.path.join(dir_output, 'percentage_beta_significant_genes_%s_%sBoxplot.pdf' % (str(np.round(signif_threshold, 3)), col)), bbox_inches='tight')
        plt.close()


        

            

            # Set labels and title
        

        

        

    ########################################################################################################################################
    ######################################################## BARPLOT OF JOANA RESULTS ######################################################
    ########################################################################################################################################

    
    def plot_barplot(self, filename, joana_output='joana_output', threshold=0.5):        
        if getattr(self.enrichment_obj, joana_output) is None:
            print('Please run JOANA model at first or load a JOANA output.')
        if self.enrichment_obj.type == 'single-species':
            self.__plot_barplot_single_species(filename, joana_output=joana_output, threshold=threshold,)
        elif self.enrichment_obj.type == 'cooperative':
            self.__plot_barplot_cooperative(filename, joana_output=joana_output, threshold=threshold)

    def __plot_barplot_single_species(self, filename, joana_output='joana_output', threshold=0.5):
        joana_result = getattr(self.enrichment_obj, joana_output).sort_values(by='Single-Species', ascending=False)
        joana_result_filtered = joana_result[joana_result['Single-Species'] > threshold]
        joana_result_filtered.reset_index(inplace=True)

        n_active_terms = joana_result_filtered.shape[0]
        fig_heigth_factor = (n_active_terms / 5) - 1
        if fig_heigth_factor < 0:
            fig_heigth_factor = 0
        fig_height = 1 + 2.3*fig_heigth_factor

        plt.figure(figsize=(8, fig_height))
        ax = sns.barplot(x='Single-Species', y='Terms', data=joana_result_filtered, color='forestgreen')
        plt.xlim([threshold, 1.])
        fname_ending = filename.split('.')[-1]
        plt.savefig(filename.replace('.'+fname_ending, '_' + joana_output + '.' + fname_ending), bbox_inches='tight')
        plt.close()

    def __plot_barplot_cooperative(self, filename, joana_output='joana_output', threshold=0.5):
        # automatic fig_size.. 
        # 1 inch = 177 pixel heigth
        # adding an inch is additional 77 pixel in height
        
        joana_result = getattr(self.enrichment_obj, joana_output)
        
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

        # define fig size
        n_active_terms = joana_cooperative_merge.shape[0] + joana_single_species_1_merge.shape[0] + joana_single_species_2_merge.shape[0]
        fig_heigth_factor = (n_active_terms / 5) - 1
        if fig_heigth_factor < 0:
            fig_heigth_factor = 0
        fig_height = 1 + 2.3*fig_heigth_factor

        # iterative assignemnt of how many subplots we need based on number of active modeling results from the different species     
        active_models = int(joana_cooperative_merge.shape[0] > 0) + int(joana_single_species_1_merge.shape[0] > 0) + int(joana_single_species_2_merge.shape[0] > 0)
        if active_models == 3:
            subplots_number = 311
            gs = gridspec.GridSpec(3, 1, height_ratios=[joana_cooperative_merge.shape[0], joana_single_species_1_merge.shape[0], joana_single_species_2_merge.shape[0]])    
        elif active_models == 2:
            subplots_number = 211
            if joana_cooperative_merge.shape[0] > 0 and joana_single_species_1_merge.shape[0] > 0:
                gs = gridspec.GridSpec(2, 1, height_ratios=[joana_cooperative_merge.shape[0], joana_single_species_1_merge.shape[0]])
            elif joana_cooperative_merge.shape[0] > 0 and joana_single_species_2_merge.shape[0] > 0:
                gs = gridspec.GridSpec(2, 1, height_ratios=[joana_cooperative_merge.shape[0], joana_single_species_2_merge.shape[0]])    
            else:
                gs = gridspec.GridSpec(2, 1, height_ratios=[joana_single_species_1_merge.shape[0], joana_single_species_2_merge.shape[0]])    
        else:
            subplots_number = None

        gs_counter = 0
        color = {'Cooperative' : 'lime', 'Single-Species-1' : 'orange', 'Single-Species-2' : 'khaki'}
        fig = plt.figure(figsize=(8, fig_height))
        plt.subplots_adjust(hspace=.5)
        if joana_cooperative_merge.shape[0] > 0:
            if not subplots_number is None:
                ax_sub = plt.subplot(gs[gs_counter])
                gs_counter += 1
                subplots_number += 1
            ax = sns.barplot(y='Terms', x='Probability', hue='Model', data=joana_cooperative_merge, palette=color)
            plt.legend([],[], frameon=False)
            ax.set_xlabel('Probability', fontsize=18)
            ax.set_ylabel('Terms', fontsize=18)
            ax.set_title('Cooperative', fontsize=18)
            plt.xlim([threshold, 1.])
        if joana_single_species_1_merge.shape[0] > 0:
            if not subplots_number is None:
                ax_sub = plt.subplot(gs[gs_counter])
                gs_counter += 1
                subplots_number += 1
            ax = sns.barplot(y='Terms', x='Probability', hue='Model', data=joana_single_species_1_merge, palette=color)
            plt.legend([],[], frameon=False)
            ax.set_xlabel('Probability', fontsize=18)
            ax.set_ylabel('Terms', fontsize=18)
            ax.set_title('Single-Species-1', fontsize=20)
            plt.xlim([threshold, 1.])
        if joana_single_species_2_merge.shape[0] > 0:
            if not subplots_number is None:
                ax_sub = plt.subplot(gs[gs_counter])
                gs_counter += 1
                subplots_number += 1
            ax = sns.barplot(y='Terms', x='Probability', hue='Model', data=joana_single_species_2_merge, palette=color)
            plt.legend([],[], frameon=False)
            ax.set_xlabel('Probability', fontsize=18)
            ax.set_ylabel('Terms', fontsize=18)
            ax.set_title('Single-Species-2', fontsize=20)
            plt.xlim([threshold, 1.])
        fname_ending = filename.split('.')[-1]
        plt.savefig(filename.replace('.'+fname_ending, '_' + joana_output + '.' + fname_ending), bbox_inches='tight')
        plt.close()




    ########################################################################################################################################
    ######################################################## GRAPH PLOT OF JOANA RESULTS ###################################################
    ########################################################################################################################################


    def plot(self, filename, top_percentile_edges=0, scaling_factor=10, fig_heigth=12, fig_width=15,
             bbox=(100, 100, 600, 600), layout='kk', verbose_n_top_terms=-1, second_order='detailed', n_top_terms=10,max_term_label=20):
        if self.enrichment_obj.joana_output is None:
            print('Please run JOANA model at first or load a JOANA output.')
        else:
            if self.enrichment_obj.type == 'single-species':
                self.__plot_single_species(filename, top_percentile_edges, scaling_factor, fig_heigth, fig_width,
                                          bbox, layout, verbose_n_top_terms, second_order, n_top_terms,max_term_label)
            elif self.enrichment_obj.type == 'cooperative':
                self.__plot_cooperative(filename, top_percentile_edges, scaling_factor, fig_heigth, fig_width,
                                          bbox, layout, verbose_n_top_terms, second_order, n_top_terms,max_term_label)



    def __plot_single_species(self, filename, top_percentile_edges=0, scaling_factor=10, fig_heigth=12, fig_width=15,
             bbox=(100, 100, 600, 600), layout='kk', verbose_n_top_terms=-1, second_order='detailed', n_top_terms=10,max_term_label=20):

        # get top n terms
        output_sorted = self.enrichment_obj.joana_output.sort_values('Single-Species', ascending=False)
        
        ind_prob = len(np.where(output_sorted['Single-Species'].values > 0.5)[0])
        if ind_prob <= 1:
                print("Notice: Network plotting is skipped because the number of enriched pathways is 1 or less.")
                return  # Exit the function early

        if(ind_prob<n_top_terms):
            n_top_terms=ind_prob
        top_single_species=output_sorted.iloc[:(n_top_terms)]
        
        prob_terms= top_single_species
        n_prob_terms= len(prob_terms)
        
        ind_terms = np.asarray([np.where(np.asarray(self.enrichment_obj.terms) == x)[0] for x in prob_terms.index])
        term_sizes = np.zeros((n_prob_terms,))
        term_annotations = [list() for i in range(n_prob_terms)]
        term_sizes=self.enrichment_obj.term_size[ind_terms.flatten()]
        term_annotations = [self.enrichment_obj.term_annotation[i] for i in ind_terms.flatten()]
        prob_terms_size_annotation = pd.DataFrame.copy(prob_terms)
        prob_terms_size_annotation['size'] = term_sizes
        prob_terms_size_annotation['annotation']= term_annotations
        
        
        if second_order == 'detailed':
            prob_terms_size_annotation = prob_terms_size_annotation.sort_values(by=['Single-Species', 'size'], ascending=[False, True])
            
        else:
            prob_terms_size_annotation = prob_terms_size_annotation.sort_values(by=['Single-Species', 'size'], ascending=[False, False])
            
        

        # how many labels to show on graph
        if verbose_n_top_terms == -1:
            verbose_n_top_terms = n_prob_terms
        
        

        # compute similarity based on term overlap
        similarity_matrix_terms = np.zeros((n_prob_terms, n_prob_terms))
        N_similarity_matrix_terms = np.zeros((n_prob_terms, n_prob_terms))
        for i in range(n_prob_terms):
            for j in range(n_prob_terms):
                if i != j:
                    similarity_matrix_terms[i, j] = len(np.intersect1d(prob_terms_size_annotation.iloc[i,2], prob_terms_size_annotation.iloc[j,2]))
                    N_similarity_matrix_terms[i, j]= similarity_matrix_terms[i, j]/min(prob_terms_size_annotation.iloc[i,1],prob_terms_size_annotation.iloc[j,1])
        
        N_similarity_matrix_terms_log=-np.log(1/(1+N_similarity_matrix_terms))
        N_similarity_matrix_terms_log=abs(N_similarity_matrix_terms_log)
        similarity_matrix_terms_log = np.log(1 + similarity_matrix_terms)
        
        

        # show only top percentile of connections
        median_similarity = np.percentile(N_similarity_matrix_terms_log.flatten(), top_percentile_edges)
        N_similarity_matrix_terms_log[N_similarity_matrix_terms_log < median_similarity] = 0
        

        # normalize weights between 0-1
        similarity_matrix_terms_log = (similarity_matrix_terms_log - np.min(similarity_matrix_terms_log)) / (
                    np.max(similarity_matrix_terms_log) - np.min(similarity_matrix_terms_log))
        similarity_matrix_terms_log *= scaling_factor
        N_similarity_matrix_terms_log *= scaling_factor
        

        

        # get the row, col indices of the non-zero elements in your adjacency matrix
        conn_indices = np.where(np.triu(N_similarity_matrix_terms_log))

        # get the weights corresponding to these indices
        weights = N_similarity_matrix_terms_log[conn_indices]
    
        # vertex size based on term size
        term_sizes_logged = 1 + np.log(prob_terms_size_annotation.iloc[:,1]) - np.min(np.log(prob_terms_size_annotation.iloc[:,1]))
    
        min_term_sizes_scaled = np.percentile(term_sizes_logged, 25)
        factor_term_size = 15 / min_term_sizes_scaled
        final_vertex_sizes = 10 + term_sizes_logged * factor_term_size * 500
        value_sizes=prob_terms_size_annotation.iloc[:,1]/max(prob_terms_size_annotation.iloc[:,1])
        final_vertex_sizes = value_sizes*1000
  

        col = 'Single-Species'

        temp_label=prob_terms_size_annotation.index
        strings_without_prefix = [s.split('_', 1)[-1] for s in temp_label]
        modified_strings = [s[0] + s[1:max_term_label].lower() + ("..." if len(s) > max_term_label else "") for s in strings_without_prefix]
        temp_label= modified_strings
       


        edges_list = list(zip(*conn_indices))

        # Create an igraph Graph object
        G_igraph = Graph(edges=edges_list, directed=False)

        # Assign weights to the edges
        G_igraph.es['weight'] = weights 

        # Vertex sizes based on term size
        # Assuming you have already computed 'final_vertex_sizes'

        # Create a NetworkX graph
        G_nx = nx.Graph()

        # Add nodes with attributes to the NetworkX graph
        for v_idx, label in enumerate(temp_label):
                G_nx.add_node(v_idx, size=final_vertex_sizes[v_idx], label=label)

        # Add edges with weights to the NetworkX graph
        for e in G_igraph.es:
            source = e.source
            target = e.target
            weight = e['weight']
            G_nx.add_edge(source, target, weight=weight)

        # Now G_nx is a NetworkX graph equivalent to G_igraph

        # Define a colormap
        np.random.seed(42)
        colors = [(1, 1, 1), (1, 0, 0)]  # Red to white
        cm = LinearSegmentedColormap.from_list("my_colormap", colors, N=256)


        node_sizes = [G_nx.nodes[node]['size'] for node in G_nx.nodes()]
        node_labels = {node: G_nx.nodes[node]['label'] for node in G_nx.nodes()}

        # Get edge weights from the NetworkX graph
        edge_weights = [G_nx.edges[edge]['weight'] for edge in G_nx.edges()]

        # Create the figure and axis objects
        plt.figure(figsize=(8, 6))
        ax = plt.gca()
        prob=prob_terms_size_annotation.iloc[:,0].values
        
        
        # Draw the graph
        pos = nx.spring_layout(G_nx, k=0.9)  # You may need to adjust the layout algorithm
        nx.draw(G_nx, pos, ax=ax, with_labels=False, node_size=node_sizes, node_color=[cm(value) for value in prob], width=edge_weights)

        # Draw node labels above the nodes
        for node, (x, y) in pos.items():
            plt.text(x, y + 0.1, node_labels[node], horizontalalignment='center', color='blue', fontsize=10)

        # Draw color bar
        sm = plt.cm.ScalarMappable(cmap=cm)
        sm.set_array([])
        #cbar = plt.colorbar(sm, ticks=[0, 0.5, 1])
        plt.colorbar(sm, label='Probability of activity',
                     cax=ax.inset_axes([0.95, 0.1, 0.05, 0.8]),)
        #cbar.set_label('Node Values')

        plt.title('Top enriched gene-sets')

        # Save the plot to a PDF file
        plt.savefig(filename.replace('.pdf', '_%d_topTerms.pdf' % verbose_n_top_terms))

        # Close the plot to avoid displaying it
        plt.close()



    def __plot_cooperative(self, filename, top_percentile_edges=0, scaling_factor=5, fig_heigth=12, fig_width=15,
             bbox=(100, 100, 600, 600), layout='kk', verbose_n_top_terms=-1, second_order='detailed', n_top_terms=10,max_term_label=20):
    
        # get top n terms

        for col in ['Cooperative', 'Single-Species-1', 'Single-Species-2']:
            output_sorted = self.enrichment_obj.joana_output.sort_values(col, ascending=False)
            
            ind_prob = len(np.where(output_sorted[col].values > 0.5)[0])
            
            if ind_prob <= 1:
                print("Notice: Network plotting is skipped for %s because the number of enriched pathways is 1 or less." % col)
                return  # Exit the function early
    

            n_terms=n_top_terms
            if(ind_prob<n_top_terms):
                n_terms=ind_prob
            top_single_species=output_sorted.iloc[:(n_terms)]
            
            prob_terms= top_single_species
            n_prob_terms= len(prob_terms)
            
            ind_terms = np.asarray([np.where(np.asarray(self.enrichment_obj.terms) == x)[0] for x in prob_terms.index])
            term_sizes = np.zeros((n_prob_terms,))
            term_annotations = [list() for i in range(n_prob_terms)]
            term_sizes=self.enrichment_obj.term_size[ind_terms.flatten()]
            term_annotations = [self.enrichment_obj.term_annotation[i] for i in ind_terms.flatten()]
            prob_terms_size_annotation = pd.DataFrame.copy(prob_terms)
            prob_terms_size_annotation = prob_terms_size_annotation[[col]]
            prob_terms_size_annotation['size'] = term_sizes
            prob_terms_size_annotation['annotation']= term_annotations
            
            
            if second_order == 'detailed':
                prob_terms_size_annotation = prob_terms_size_annotation.sort_values(by=[col, 'size'], ascending=[False, True])
                
            else:
                prob_terms_size_annotation = prob_terms_size_annotation.sort_values(by=[col, 'size'], ascending=[False, False])
                
            

            # how many labels to show on graph
            if verbose_n_top_terms == -1:
                verbose_n_top_terms = n_terms
                

            # compute similarity based on term overlap
            similarity_matrix_terms = np.zeros((n_prob_terms, n_prob_terms))
            N_similarity_matrix_terms = np.zeros((n_prob_terms, n_prob_terms))
            for i in range(n_prob_terms):
                for j in range(n_prob_terms):
                    if i != j:
                        similarity_matrix_terms[i, j] = len(np.intersect1d(prob_terms_size_annotation.iloc[i,2], prob_terms_size_annotation.iloc[j,2]))
                        N_similarity_matrix_terms[i, j]= similarity_matrix_terms[i, j]/min(prob_terms_size_annotation.iloc[i,1],prob_terms_size_annotation.iloc[j,1])
            
            N_similarity_matrix_terms_log=-np.log(1/(1+N_similarity_matrix_terms))
            N_similarity_matrix_terms_log=abs(N_similarity_matrix_terms_log)
            similarity_matrix_terms_log = np.log(1 + similarity_matrix_terms)
            
            

            # show only top percentile of connections
            median_similarity = np.percentile(N_similarity_matrix_terms_log.flatten(), top_percentile_edges)
            N_similarity_matrix_terms_log[N_similarity_matrix_terms_log < median_similarity] = 0
            

            # normalize weights between 0-1
            similarity_matrix_terms_log = (similarity_matrix_terms_log - np.min(similarity_matrix_terms_log)) / (
                        np.max(similarity_matrix_terms_log) - np.min(similarity_matrix_terms_log))
            similarity_matrix_terms_log *= scaling_factor
            N_similarity_matrix_terms_log *= scaling_factor
            

            

            # get the row, col indices of the non-zero elements in your adjacency matrix
            conn_indices = np.where(np.triu(N_similarity_matrix_terms_log))

            # get the weights corresponding to these indices
            weights = N_similarity_matrix_terms_log[conn_indices]
            
            # vertex size based on term size
            term_sizes_logged = 1 + np.log(prob_terms_size_annotation.iloc[:,1]) - np.min(np.log(prob_terms_size_annotation.iloc[:,1]))
            
            min_term_sizes_scaled = np.percentile(term_sizes_logged, 25)
            factor_term_size = 15 / min_term_sizes_scaled
            final_vertex_sizes = 10 + term_sizes_logged * factor_term_size * 500
            value_sizes=prob_terms_size_annotation.iloc[:,1]/max(prob_terms_size_annotation.iloc[:,1])
            final_vertex_sizes = value_sizes*1000
            

            

            temp_label=prob_terms_size_annotation.index
            strings_without_prefix = [s.split('_', 1)[-1] for s in temp_label]
            modified_strings = [s[0] + s[1:max_term_label].lower() + ("..." if len(s) > max_term_label else "") for s in strings_without_prefix]
            temp_label= modified_strings
        
            

            edges_list = list(zip(*conn_indices))

            # Create an igraph Graph object
            G_igraph = Graph(edges=edges_list, directed=False)

            # Assign weights to the edges
            G_igraph.es['weight'] = weights

            # Vertex sizes based on term size
            # Assuming you have already computed 'final_vertex_sizes'

            # Create a NetworkX graph
            G_nx = nx.Graph()

            # Add nodes with attributes to the NetworkX graph
            for v_idx, label in enumerate(temp_label):
                G_nx.add_node(v_idx, size=final_vertex_sizes[v_idx], label=label)

            # Add edges with weights to the NetworkX graph
            for e in G_igraph.es:
                source = e.source
                target = e.target
                weight = e['weight']
                G_nx.add_edge(source, target, weight=weight)
            

            # Adjust the node separation by specifying the k parameter (default is k=0.1)
            # Higher values of k will result in greater distances between nodes
            
            # Now G_nx is a NetworkX graph equivalent to G_igraph

            # Define a colormap
            #np.random.seed(42)
            colors = [(1, 1, 1), (1, 0, 0)]  # Red to white
            cm = LinearSegmentedColormap.from_list("my_colormap", colors, N=256)


            node_sizes = [G_nx.nodes[node]['size'] for node in G_nx.nodes()]
            node_labels = {node: G_nx.nodes[node]['label'] for node in G_nx.nodes()}

            # Get edge weights from the NetworkX graph
            edge_weights = [G_nx.edges[edge]['weight'] for edge in G_nx.edges()]

            # Create the figure and axis objects
            plt.figure(figsize=(8, 6))
            ax = plt.gca()
            prob=prob_terms_size_annotation.iloc[:,0].values
            
            
            # Draw the graph
            pos = nx.spring_layout(G_nx,k=0.9)  # You may need to adjust the layout algorithm
            nx.draw(G_nx, pos, ax=ax, with_labels=False, node_size=node_sizes, node_color=[cm(value) for value in prob], width=edge_weights)

            # Draw node labels above the nodes
            for node, (x, y) in pos.items():
                plt.text(x, y + 0.0, node_labels[node], horizontalalignment='center', color='blue', fontsize=11)

            # Draw color bar
            sm = plt.cm.ScalarMappable(cmap=cm)
            sm.set_array([])
            #cbar = plt.colorbar(sm, ticks=[0, 0.5, 1])
            plt.colorbar(sm, label='Probability of activity',
                         cax=ax.inset_axes([0.95, 0.1, 0.05, 0.8]),)
            #cbar.set_label('Node Values')

            plt.title('Top enriched gene-sets')

            # Save the plot to a PDF file
            plt.savefig(filename.replace('.pdf', '_%d_%s.pdf' % (n_terms, col)))
            
            # Close the plot to avoid displaying it
            plt.close()
            