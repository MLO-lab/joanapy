import numpy as np 
import pandas as pd
from tqdm import trange



def transform_gmt_assignment_matrix(filename_gmt_assignment, gene_ids, filename_output, write_terms_to_file=True, write_term_summary=True):
    if not filename_gmt_assignment.endswith('.txt'):
        filename_gmt_assignment.replace(filename_gmt_assignment.split('.')[-1], 'txt')

    lines = tuple(open(filename_gmt_assignment, 'r'))
    
    term_id = list()
    term_description = list()
    gene_annotations = list()
    term_sizes = list()
    for i in range(0, len(lines)):
        if (lines[i].endswith('\n')):
            line_temp = lines[i][:-1]
        else:
            line_temp = lines[i]

        split = line_temp.split('\t')
        term_id.append(split[0].replace(' ', '_'))
        term_description.append(split[1].replace(' ', '_'))
        gene_annotations.append(list(split[2:]))
        if ' ' in gene_annotations[-1]:
            gene_annotations[-1].remove(' ')
        if '' in gene_annotations[-1]:
            gene_annotations[-1].remove('')
        # gene_annotations[-1] = np.asarray(gene_annotations[-1])
        term_sizes.append(len(gene_annotations[-1]))

    # create assignment matrix format
    list_assignment_matrix = list()
    for i in trange(len(gene_ids)):
        temp_assignments_gene = list()
        for j in range(len(gene_annotations)):
            if str(gene_ids[i]) in gene_annotations[j]:
                temp_assignments_gene.append(str(j))

        list_assignment_matrix.append(','.join(temp_assignments_gene))

    ind_gene_annotated = np.where(np.asarray(list_assignment_matrix) != '')[0]
    assignment_matrix = np.asarray(list_assignment_matrix)[ind_gene_annotated]

    # write to file
    with open(filename_output, "w") as myfile:
        for line in assignment_matrix:
            myfile.write(str(line) + '\n')

    if write_terms_to_file:
        pd.DataFrame(term_id).to_csv(filename_output.replace('.txt', '_terms.txt'), index=False, header=False)
    
    if write_term_summary:
        pd.DataFrame({'ID' : term_id, 'description' : term_description, 'size' : term_sizes, 'annotations' : ['/'.join(x) for x in gene_annotations]}).to_csv(filename_output.replace('.txt', '_term_summary.csv'), index=False, header=True)

    return(ind_gene_annotated)





    