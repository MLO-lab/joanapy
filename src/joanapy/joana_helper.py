import numpy as np 
import pandas as pd
from tqdm import trange



def transform_gmt_assignment_matrix(filename_gmt_assignment, gene_ids, filename_output_assignment, filename_output_terms, write_term_summary=True, min_term_utilization=0.):
    if not filename_gmt_assignment.endswith('.txt'):
        filename_gmt_assignment.replace(filename_gmt_assignment.split('.')[-1], 'txt')

    lines = tuple(open(filename_gmt_assignment, 'r'))
    
    term_id = list()
    term_description = list()
    gene_annotations = list()
    term_sizes = list()
    term_utilization = list()
    num_terms_dropped_utilization = 0
    num_terms_zero_size = 0
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

        # remove term if not genes are annotated and term size is zero
        if term_sizes[-1] == 0:
            term_id.pop()
            term_description.pop()
            gene_annotations.pop()
            term_sizes.pop()
            num_terms_zero_size += 1
            continue

        # filter for lowly utilized terms
        temp_number_annotated = len(np.intersect1d(gene_annotations[-1], gene_ids))
        term_utilization.append(temp_number_annotated / term_sizes[-1])
        if term_utilization[-1] < min_term_utilization:
            term_id.pop()
            term_description.pop()
            gene_annotations.pop()
            term_sizes.pop()
            term_utilization.pop()
            num_terms_dropped_utilization += 1
            continue

    print('Dropped %d terms due to non annotated genes in assignemnt.' %(num_terms_zero_size))
    print('Filtered %d terms, %d terms remain with a minimum term utilization of %s.' %(num_terms_dropped_utilization, len(term_id), "{:.2f}".format(min_term_utilization)))

    # create JOANA assignment matrix format
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
    with open(filename_output_assignment, "w") as myfile:
        for line in assignment_matrix:
            myfile.write(str(line) + '\n')

    pd.DataFrame(term_id).to_csv(filename_output_terms, index=False, header=False)
    
    if write_term_summary:
        pd.DataFrame({'ID' : term_id, 'description' : term_description, 'size' : term_sizes, 'utilization' : term_utilization, 'annotations' : ['/'.join(x) for x in gene_annotations]}).to_csv(filename_output_terms.replace('.txt', '_summary.csv'), index=False, header=True)

    return(ind_gene_annotated)



def load_assignmentMatrix(filename):
    lines = tuple(open(filename, 'r'))
    assignment = []

    for i in range(0, len(lines)):
        if (lines[i].endswith('\n')):
            line_temp = lines[i][:-1]
        else:
            line_temp = lines[i]

        split = line_temp.split(',')
        # assignment.append(map(int,split))
        assignment.append(np.asarray([int(x) for x in split]))

    return (assignment)


