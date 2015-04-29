__author__ = 'heroico'

import csv

def LoadPeople(sample_file_name):
    people_by_row = {}
    with open(sample_file_name, 'rb') as samples:
        reader = csv.reader(samples, delimiter='\t', quotechar='"')
        for row in reader:
            person_id = row[0]
            people_by_row[reader.line_num] = person_id
    return people_by_row


KEY_MRNA_GENE_NAME = "name"
KEY_MRNA_GENE_COL = "col"
KEY_MRNA_VALUES = "values"

def BuildGeneData(people_by_row, data_file_name):
    print data_file_name
    geneData = {}
    with open(data_file_name, 'rb') as file:
        reader = csv.reader(file, delimiter="\t", quotechar='"')
        read_first_row = False
        for row in reader:
            if not read_first_row:
                #print sorted(row)
                for col,gene in enumerate(row):
                    gene_item = {KEY_MRNA_GENE_NAME:gene, KEY_MRNA_GENE_COL:col}
                    geneData[gene] = gene_item
                read_first_row = True
            else:
                person_row = reader.line_num-1
                for gene, gene_item in geneData.iteritems():
                    values = None
                    if KEY_MRNA_VALUES in gene_item:
                        values = gene_item[KEY_MRNA_VALUES]
                    else:
                        values = {}
                        gene_item[KEY_MRNA_VALUES] = values
                    col = gene_item[KEY_MRNA_GENE_COL]
                    person = people_by_row[person_row]
                    values[person] = row[col]
    return geneData

import numpy

def Correlate(people_1, data_1, people_2, data_2):
    set_1 = []
    set_2 = []
    for gene_name, item_1 in data_1.iteritems():
        name_1 = item_1[KEY_MRNA_GENE_NAME]
        if name_1 in data_2:
            item_2 = data_2[name_1]

            values_1 = item_1[KEY_MRNA_VALUES]
            values_2 = item_2[KEY_MRNA_VALUES]
            for person, value in values_1.iteritems():
                if person in values_2:
                    set_1.append(float(value))
                    set_2.append(float(values_2[person]))
    pearson = numpy.corrcoef(set_1, set_2)[0,1]
    covariance = numpy.cov(set_1, set_2)[0,1]
    return pearson, covariance

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Correlation between predixcan and Affymetrix data.')
    parser.add_argument("--affy_results_file_name",
                        help="file with affyemtrix in predixcan format",
                        default="Data/affy-results.txt")

    parser.add_argument("--affy_sample_file_name",
                        help="Affymetrix people samples file",
                        default="Data/affy-samples.txt")

    parser.add_argument("--predixcan_results_file_name",
                        help="file with predixcan data",
                        default="Data/predixcan-results.csv")

    parser.add_argument("--predixcan_sample_file_name",
                        help="predixcan sample people file",
                        default="Data/predixcan-samples.txt")
    args = parser.parse_args()

    #-------------------------------------------------------------------------------------------------------------------
    affy_people_by_row = LoadPeople(args.affy_sample_file_name)
    affy_gene_data = BuildGeneData(affy_people_by_row, args.affy_results_file_name)

    #-------------------------------------------------------------------------------------------------------------------
    predy_people_by_row = LoadPeople(args.predixcan_sample_file_name)
    predy_gene_data = BuildGeneData(predy_people_by_row, args.predixcan_results_file_name)

    pearson, covariance = Correlate(affy_people_by_row, affy_gene_data, predy_people_by_row, predy_gene_data)
    print "Pearson:"+str(pearson)
    print "Cov:"+str(covariance)
