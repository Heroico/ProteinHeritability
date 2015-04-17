__author__ = 'heroico'
import csv

# look at gencode http://www.gencodegenes.org/data_format.html
COL_CHROMOSOME = 0
COL_FEATURE_TYPE = 1
COL_N_START_LOCATION = 2
COL_N_END_LOCATION = 3
COL_ENS_ID = 4
COL_GENE_NAME = 5
COL_GENE_TYPE = 6

KEY_ENSEMBLE_ID = "ensemble"
KEY_GENE_NAME = "name"

def ReadGeneCodeInput(file_name='Data/gencode.v12.V1.summary.protein'):
    gencodes = {}
    with open(file_name, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        read_first = False
        for row in reader:
            ensemble_version = row[COL_ENS_ID]
            ensemble = ensemble_version.split('.')[0]
            code = {KEY_ENSEMBLE_ID:ensemble, KEY_GENE_NAME: row[COL_GENE_NAME]}
            if not ensemble in gencodes:
                gencodes[ensemble] = code
            else:
                raise Exception('Duplicate ensemble id, check'+ensemble)
    return gencodes