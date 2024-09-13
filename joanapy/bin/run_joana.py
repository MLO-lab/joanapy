#!/usr/bin/env python
from pandas.core.dtypes import missing
from joanapy.enrichment_obj import ENRICHMENT_OBJ
from joanapy.joana import JOANA
from joanapy.joana_helper import transform_gmt_assignment_matrix
from joanapy.joana_helper import data_preproccessing
import os
import pandas as pd
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description='Run joana_single or joana_Cooperative based on the options provided')
    parser.add_argument('-o', '--omics1', help='File name for omics1')
    parser.add_argument('-o2', '--omics2', help='File name for omics2 (for cooperative mode)')
    parser.add_argument('-p', '--pathway', help='File name for pathway.gmt file')
    parser.add_argument('-d', '--dir', help='Path to input file directories and saving the results')
    parser.add_argument('-m', '--min_term_utilization', type=float,default=0.0, help='it removes pathways from gmt file if Minimum term utilization does not pass. It removes pathways which has atleast min_term_utilization measured genes inside')
    args = parser.parse_args()

    # Determine which function to call based on the options provided
    if args.omics2:
        run_joana_cooperative(args.omics1, args.omics2, args.pathway, args.dir, args.min_term_utilization)
    else:
        run_joana_single(args.omics1, args.pathway, args.dir, args.min_term_utilization)

def run_joana_single(filename_omics1,filename_pathway,dir,min_term_utilization=0.0):
    dir = os.path.normpath(dir)
    filename_output_single = os.path.join(dir,'results'+str(min_term_utilization), 'ResultSingle.csv')
    filename_assignment_matrix_gmt=os.path.join(dir,filename_pathway)
    filename_qvalues_first= os.path.join(dir,filename_omics1)
    #filename_qvalues_second= path+'omics2noHead.txt'
    
    min_term_size = 5
    max_term_size = 100000
    prior_pA = 1
    goodness_of_fit = True
    plot_components = True
    signif_threshold=0.1
    quantile=0.95
    if not os.path.exists(os.path.join(dir,'temp')):
        os.makedirs(os.path.join(dir,'temp'))
    
    res=data_preproccessing(filename_qvalues_first)
    
    #Cooperative Preproccessing
    if(res[0]=="coop"):
        res[1]['geneSymbol'].to_csv(os.path.join(dir,'temp','temp_geneSymbol.txt'), header=False, index=False)
        res[1]['qvalues1'].to_csv(os.path.join(dir,'temp','temp_qval1.txt'), header=False, index=False)
        res[1]['qvalues2'].to_csv(os.path.join(dir,'temp','temp_qval2.txt'), header=False, index=False)
        np.savetxt(os.path.join(dir,'temp','temp_missing1.txt'), res[2], fmt='%d')
        np.savetxt(os.path.join(dir,'temp','temp_missing2.txt'), res[3], fmt='%d')
        
    #Single Preproccessing
    else:
        res[1]['geneSymbol'].to_csv(os.path.join(dir,'temp','temp_geneSymbol.txt'), header=False, index=False)
        res[1]['qval'].to_csv(os.path.join(dir,'temp','temp_qval1.txt'), header=False, index=False)
        

    ###run Single

    filename_qvalues_first=os.path.join(dir,'temp','temp_qval1.txt')
    filename_gene_ids=os.path.join(dir,'temp','temp_geneSymbol.txt')
    filename_qvalues_first_filtered = os.path.join(dir,'temp','qvalues_first_species_filtered.txt')
    
    # read in first species
    gene_ids = pd.read_table(filename_gene_ids, header=None)
    print(min_term_size, max_term_size, min_term_utilization)
    filename_terms = os.path.join(dir,'temp','Human_GOALL_with_GO_iea_March_01_2021_symbol_terms_utilization' +str(min_term_size)+'_'+ str(max_term_size)+'_'+ str(min_term_utilization)+'.txt')
    filename_assignment_matrix = os.path.join(dir,'temp','Human_GOALL_with_GO_iea_March_01_2021_symbol_utilization'+ str(min_term_size)+'_'+ str(max_term_size)+'_'+ str(min_term_utilization)+'.txt')
    if not os.path.exists(os.path.dirname(filename_assignment_matrix)):
        os.makedirs(os.path.dirname(filename_assignment_matrix))

    ind_filtered_qvalues = transform_gmt_assignment_matrix(filename_assignment_matrix_gmt, gene_ids.values.squeeze().tolist(), filename_assignment_matrix, filename_terms, min_term_utilization=min_term_utilization)

    # filtere qvalues
    qvalues_first = pd.read_table(filename_qvalues_first, header=None)
    qvalues_first_filtered = qvalues_first.iloc[ind_filtered_qvalues]
    qvalues_first_filtered.to_csv(filename_qvalues_first_filtered, header=False, index=False)
    enrichment_obj_single = ENRICHMENT_OBJ(filename_assignment_matrix, filename_terms,
                                                filename_qvalues_first_filtered)
    joana_single = JOANA(enrichment_obj_single)
    joana_single.run(filename_output_single, goodness_of_fit=goodness_of_fit, plot_components=plot_components, prior_pA=prior_pA, min_term_size=min_term_size, max_term_size=max_term_size)
    joana_single.plot_hidden_significant_genes(os.path.dirname(filename_output_single), signif_threshold=signif_threshold, quantile=quantile) 
    joana_single.plot(filename_output_single.replace('.csv', '.pdf'))
    print(enrichment_obj_single.joana_output.sort_values(by='Single-Species', ascending=False).head())
    joana_single.plot_barplot(filename_output_single.replace('.csv', '_barplot.pdf'))



def run_joana_cooperative(filename_omics1,filename_omics2,filename_pathway,dir,min_term_utilization=0.0):
    path = os.path.normpath(dir)
    filename_output_cooperative = os.path.join(dir,'results'+str(min_term_utilization), 'ResultCooperative.csv')
    filename_assignment_matrix_gmt=os.path.join(dir,filename_pathway)
    filename_qvalues_first= os.path.join(dir,filename_omics1)
    filename_qvalues_second= os.path.join(dir,filename_omics2)
    
    
    min_term_size = 5
    max_term_size = 100000
    prior_pA = 1
    goodness_of_fit = True
    plot_components = True
    signif_threshold=0.1
    quantile=0.95
    if not os.path.exists(os.path.join(dir,'temp')):
        os.makedirs(os.path.join(dir,'temp'))
    
    res=data_preproccessing(filename_qvalues_first,filename_qvalues_second=filename_qvalues_second)
    #Cooperative Preproccessing
    if(res[0]=="coop"):
        res[1]['geneSymbol'].to_csv(os.path.join(dir,'temp','temp_geneSymbol.txt'), header=False, index=False)
        res[1]['qvalues1'].to_csv(os.path.join(dir,'temp','temp_qval1.txt'), header=False, index=False)
        res[1]['qvalues2'].to_csv(os.path.join(dir,'temp','temp_qval2.txt'), header=False, index=False)
        np.savetxt(os.path.join(dir,'temp','temp_missing1.txt'), res[2], fmt='%d')
        np.savetxt(os.path.join(dir,'temp','temp_missing2.txt'), res[3], fmt='%d')
    #Single Preproccessing
    else:
        res[1]['geneSymbol'].to_csv(os.path.join(dir,'temp','temp_geneSymbol.txt'), header=False, index=False)
        res[1]['qval'].to_csv(os.path.join(dir,'temp','temp_qval1.txt'), header=False, index=False)
        

    ###run Coop

    filename_qvalues_first=os.path.join(dir,'temp','temp_qval1.txt')
    filename_qvalues_second=os.path.join(dir,'temp','temp_qval2.txt')
    filename_gene_ids=os.path.join(dir,'temp','temp_geneSymbol.txt')
    filename_qvalues_first_filtered = os.path.join(dir,'temp','qvalues_first_species_filtered.txt')
    filename_qvalues_second_filtered = os.path.join(dir,'temp','qvalues_second_species_filtered.txt')
    filename_missing_first=os.path.join(dir,'temp','temp_missing1.txt')
    filename_missing_second=os.path.join(dir,'temp','temp_missing2.txt')
    filename_missing_first_filtered = os.path.join(dir,'temp','missing_first_species_filtered.txt')
    filename_missing_second_filtered = os.path.join(dir,'temp','missing_second_species_filtered.txt')

    # read in first species
    gene_ids = pd.read_table(filename_gene_ids, header=None)
    print(min_term_size, max_term_size, min_term_utilization)
    filename_terms = os.path.join(dir,'temp','Human_GOALL_with_GO_iea_March_01_2021_symbol_terms_utilization' +str(min_term_size)+'_'+ str(max_term_size)+'_'+ str(min_term_utilization)+'.txt')
    filename_assignment_matrix = os.path.join(dir,'temp','Human_GOALL_with_GO_iea_March_01_2021_symbol_utilization'+ str(min_term_size)+'_'+ str(max_term_size)+'_'+ str(min_term_utilization)+'.txt')
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


    enrichment_obj_cooperative = ENRICHMENT_OBJ(filename_assignment_matrix, filename_terms,
                                            filename_qvalues_first_filtered,
                                            filename_qvalues_second=filename_qvalues_second_filtered,
                                            filename_missing_first=filename_missing_first_filtered,
                                            filename_missing_second=filename_missing_second_filtered)
    
    joana_cooperative = JOANA(enrichment_obj_cooperative)
    joana_cooperative.run(filename_output_cooperative, goodness_of_fit=goodness_of_fit, plot_components=plot_components, prior_pA=prior_pA, min_term_size=min_term_size, max_term_size=max_term_size)
    joana_cooperative.plot_hidden_significant_genes(os.path.dirname(filename_output_cooperative), signif_threshold=signif_threshold, quantile=quantile) 
    joana_cooperative.plot(filename_output_cooperative.replace('.csv', '.pdf'))
    print(enrichment_obj_cooperative.joana_output.sort_values(by='Cooperative', ascending=False).head())
    joana_cooperative.plot_barplot(filename_output_cooperative.replace('.csv', '_barplot.pdf'))

#run_joana_single('omics1.txt','h.all.v6.2.symbols.gmt.txt','/Users/sareh/Desktop/testInput/hotTumorRNA/test/')
#run_joana_cooperative('omics1.txt','omics2.txt','h.all.v6.2.symbols.gmt.txt','/Users/sareh/Desktop/testInput/hotTumorMulti/')

if __name__=="__main__":
    main()