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

def PrintGene(gene_name,MRNA, people_by_predixcan_row, prefix):
    gene_item = MRNA[gene_name]
    with open(prefix+gene_name+".txt", "w+") as file:
        for num, person_id in people_by_predixcan_row.iteritems():
            values = gene_item[KEY_MRNA_VALUES]
            value = values[person_id]
            text = person_id + " " + value + "\n"
            file.write(text)

import numpy

def Correlate(people_1, data_1, people_2, data_2):
    set_1 = []
    set_2 = []
    gene_corr = []
    for gene_name, item_1 in data_1.iteritems():
        name_1 = item_1[KEY_MRNA_GENE_NAME]
        if name_1 in data_2:
            item_2 = data_2[name_1]

            values_1 = item_1[KEY_MRNA_VALUES]
            values_2 = item_2[KEY_MRNA_VALUES]

            flat_1 = []
            flat_2 = []
            for person, value in values_1.iteritems():
                if person in values_2:
                    value_1 = float(value)
                    value_2 = float(values_2[person])
                    set_1.append(value_1)
                    set_2.append(value_2)
                    flat_1.append(value_1)
                    flat_2.append(value_2)

            line = []
            pearson = 0
            covariance = 0
            samples = len(flat_1)
            if samples:
                pearson = numpy.corrcoef(flat_1, flat_2)[0,1]
                covariance = numpy.cov(flat_1, flat_2)[0,1]

            line.append(name_1)
            line.append(samples)
            line.append(pearson)
            line.append(covariance)
            gene_corr.append(line)

    pearson = numpy.corrcoef(set_1, set_2)[0,1]
    covariance = numpy.cov(set_1, set_2)[0,1]
    return gene_corr, pearson, covariance

def PrintCorrelation(gene_correlation, file_name):
    with open(file_name, "w+") as file:
        for row in gene_correlation:
            text = row[0] + " " + str(row[1]) + " " + str(row[2]) + " " + str(row[3]) + "\n"
            file.write(text)

def BuildIntersectionFiles(affy_people_by_row, affy_gene_data, predi_people_by_row, predi_gene_data, prefix):
    common_people = []
    affy_people = affy_people_by_row.values()
    predi_people = predy_people_by_row.values()
    for person_id in affy_people:
        if person_id in predi_people:
            common_people.append(person_id)
            predi_people.remove(person_id)

    common_genes = []
    affy_genes = affy_gene_data.keys()
    predi_genes = predi_gene_data.keys()
    for gene in affy_genes:
        if gene in predi_genes:
            common_genes.append(gene)
            predi_genes.remove(gene)

    affy_output = prefix+"affy.csv"
    with open(affy_output, "w+") as file:
        header = ",".join(common_genes)+"\n"
        file.write(header)

        for person_id in common_people:
            values=[]
            for gene in common_genes:
                gene_item = affy_gene_data[gene]
                values.append(gene_item[KEY_MRNA_VALUES][person_id])
            line = ",".join(values)+"\n"
            file.write(line)

    predy_output = prefix+"predi.csv"
    with open(predy_output, "w+") as file:
        header = ",".join(common_genes)+"\n"
        file.write(header)

        for person_id in common_people:
            values=[]
            for gene in common_genes:
                gene_item = predi_gene_data[gene]
                values.append(gene_item[KEY_MRNA_VALUES][person_id])
            line = ",".join(values)+"\n"
            file.write(line)

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

    parser.add_argument("--out",
                        help="file the gene ouput",
                        default="Out/predixcan-stats-results.txt")

    parser.add_argument("--affy_trace_gene",
                        help="affy gene to be output",
                        default=None)

    parser.add_argument("--intersection_output",
                       help="Will build files with intersection of affymetrix and predixcan data",
                       action='store_true',
                       default=False)

    parser.add_argument("--intersection_output_prefix",
                       help="Will build files with intersection of affymetrix and predixcan data",
                       default="IntermediateC/intersection_mrna_")

    args = parser.parse_args()

    #-------------------------------------------------------------------------------------------------------------------
    affy_people_by_row = LoadPeople(args.affy_sample_file_name)
    affy_gene_data = BuildGeneData(affy_people_by_row, args.affy_results_file_name)

    #-------------------------------------------------------------------------------------------------------------------
    predy_people_by_row = LoadPeople(args.predixcan_sample_file_name)
    predy_gene_data = BuildGeneData(predy_people_by_row, args.predixcan_results_file_name)

    if args.intersection_output:
        BuildIntersectionFiles(affy_people_by_row, affy_gene_data, predy_people_by_row, predy_gene_data, args.intersection_output_prefix)
    else:
        gene_correlation, pearson, covariance = Correlate(affy_people_by_row, affy_gene_data, predy_people_by_row, predy_gene_data)
        PrintCorrelation(gene_correlation, args.out)

        if args.affy_trace_gene is not None and len(args.affy_trace_gene):
            PrintGene(args.affy_trace_gene, affy_gene_data, affy_people_by_row, "Out/affy-")
