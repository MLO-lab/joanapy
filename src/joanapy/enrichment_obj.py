import pandas as pd
import numpy as np
import dill
from joanapy.joana_helper import load_assignmentMatrix

class ENRICHMENT_OBJ:

    #TODO: will allow to add plotting functions usw and to save the object
    moment_fit_first = None
    moment_fit_second = None
    joana_output = None
    joana_output_upregulated = None
    joana_output_downregulated = None



    def __init__(self, filename_assignment_matrix, filename_terms, filename_qvalues_first, **kwargs):
        # single species
        self.filename_assignment_matrix = filename_assignment_matrix
        self.assignment_matrix = load_assignmentMatrix(filename_assignment_matrix)
        self.filename_terms = filename_terms
        self.terms = pd.read_table(self.filename_terms, header=None)
        self.filename_qvalues_first = filename_qvalues_first
        self.qvalues_first = pd.read_csv(self.filename_qvalues_first, header=None).values.squeeze()

        # cooperative
        self.filename_qvalues_second = kwargs.get("filename_qvalues_second", None)
        if not self.filename_qvalues_second is None:
            self.qvalues_second = pd.read_csv(self.filename_qvalues_second, header=None).values.squeeze()
        self.filename_missing_first = kwargs.get("filename_missing_first", None)
        if not self.filename_missing_first is None:
            self.missing_first = pd.read_table(self.filename_missing_first, header=None)
        self.filename_missing_second = kwargs.get("filename_missing_second", None)
        if not self.filename_missing_second is None:
            self.missing_second = pd.read_table(self.filename_missing_second, header=None)

        if self.filename_qvalues_second is None:
            self.type = 'single-species'
        else:
            self.type = 'cooperative'

    
    def save(self, filename):
        dill.dump(self, file = open(filename, "wb"))

    @staticmethod
    def load(filename):
        return(dill.load(open(filename, "rb")))



    