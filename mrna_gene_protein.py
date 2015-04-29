__author__ = 'heroico'

import csv

#-----------------------------------------------------------------------------------------------------------------------
KEY_GENE_NAME = "name"
KEY_GENE_ENSEMBLE = "ensemble"
KEY_GENE_PROTEINS = "proteins"

def BuildGeneToProteinRelationShip(file_name="Intermediate/HauseGeneToProtein.txt"):
    gene_to_protein = {}
    with open(file_name, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='"')
        for row in reader:
            cols = len(row)
            entry = {KEY_GENE_NAME: row[0], KEY_GENE_ENSEMBLE:row[1], KEY_GENE_PROTEINS:[]}
            for i in xrange(2,cols):
                entry[KEY_GENE_PROTEINS].append(row[i])
            gene_to_protein[entry[KEY_GENE_NAME]] = entry
    return gene_to_protein

#-----------------------------------------------------------------------------------------------------------------------
KEY_PERSON_ID = "id"
KEY_PERSON_PREDIXCAN_ROW = "predixcan_row"

def BuildPeopleInPhenoList(fam_file_name="Intermediate/hapmap_r23a.fam", predixcan_file_name="Data/predixcan-samples.txt"):
    people = []
    with open(fam_file_name, 'rb') as fam:
        reader = csv.reader(fam, delimiter=' ', quotechar='"')
        for row in reader:
            person = {KEY_PERSON_ID:row[1]}
            people.append(person)

    people_by_predixcan_row = {}
    with open(predixcan_file_name, 'rb') as samples:
        reader = csv.reader(samples, delimiter='\t', quotechar='"')
        for row in reader:
            person_id = row[0]
            candidates = [person for person in people if person[KEY_PERSON_ID] == person_id]
            candidates_num = len(candidates)

            if candidates_num > 1:
                raise "can't handle more than one candidate"
            elif candidates_num == 1:
                person = candidates[0]
                person[KEY_PERSON_PREDIXCAN_ROW] = reader.line_num
                people_by_predixcan_row[reader.line_num] = person
    return people, people_by_predixcan_row

#-----------------------------------------------------------------------------------------------------------------------
KEY_MRNA_GENE_NAME = "name"
KEY_MRNA_GENE_COL = "col"
KEY_MRNA_VALUES = "values"

def BuildMRNAData(selected_people_by_predixcan_row, gene_to_protein, predixcan_data_file_name):
    MRNA = {}
    with open(predixcan_data_file_name, 'rb') as file:
        reader = csv.reader(file, delimiter="\t", quotechar='"')
        read_first_row = False
        for row in reader:
            if not read_first_row:
                #print sorted(row)
                for col,gene in enumerate(row):
                    if gene in gene_to_protein:
                        mrna_item = {KEY_MRNA_GENE_NAME:gene, KEY_MRNA_GENE_COL:col}
                        MRNA[gene] = mrna_item
                read_first_row = True
            else:
                person_row = reader.line_num-1
                if person_row in selected_people_by_predixcan_row:
                    person = selected_people_by_predixcan_row[person_row]
                    for gene, mrna_item in MRNA.iteritems():
                        values = None
                        if KEY_MRNA_VALUES in mrna_item:
                            values = mrna_item[KEY_MRNA_VALUES]
                        else:
                            values = {}
                            mrna_item[KEY_MRNA_VALUES] = values
                        values[person[KEY_PERSON_ID]] = row[mrna_item[KEY_MRNA_GENE_COL]]
    return MRNA

#-----------------------------------------------------------------------------------------------------------------------
import process_pheno

KEY_PHENO_PERSON_ID = "person_id"
KEY_PHENO_VALUE = "value"

def LoadPheno(pheno_file_name):
    pheno = {}
    with open(pheno_file_name, 'rb') as file:
        reader = csv.reader(file, delimiter=" ", quotechar='"')
        for row in reader:
            pheno[row[1]] = row[2]
    return pheno

#-----------------------------------------------------------------------------------------------------------------------
import math
import numpy
def MRNAProteinAverageCorrelation(MRNA, gene_to_protein,  pheno_prefix):
    correlation = []
    for gene, mrna_item in MRNA.iteritems():
        proteins = gene_to_protein[gene][KEY_GENE_PROTEINS]
        mrna_values = []
        protein_values = []

        phenos = []
        for protein in proteins:
            name = process_pheno.PhenoFileName(pheno_prefix, protein)
            pheno = LoadPheno(name)
            phenos.append(pheno)

        for person_id, mrna_value in mrna_item[KEY_MRNA_VALUES].iteritems():
            if mrna_value == "NA":
                raise  "Can't handle invalid mrna value"

            pheno_value = 0
            pheno_count = 0

            for pheno in phenos:
                if person_id in pheno:
                    value = pheno[person_id]
                    if not "NA" in value:
                        float_value = float(value)
                        pheno_value += float_value
                        pheno_count += 1

            if pheno_count > 0:
                the_value = float(mrna_value)
                if not math.isnan(the_value):
                    mrna_values.append(the_value)
                    protein_values.append(pheno_value/pheno_count)

        pearson = numpy.corrcoef(mrna_values, protein_values)[0,1]
        line = []
        line.append(gene)
        line.append(pearson)
        line.append(len(mrna_values))
        correlation.append(line)
    return correlation

def MRNAProteinCorrelation(MRNA, gene_to_protein,  pheno_prefix):
    correlation = []
    for gene, mrna_item in MRNA.iteritems():
        proteins = gene_to_protein[gene][KEY_GENE_PROTEINS]

        for protein in proteins:
            name = process_pheno.PhenoFileName(pheno_prefix, protein)
            pheno = LoadPheno(name)

            protein_values = []
            mrna_values = []
            for person_id, pheno_value in pheno.iteritems():
                if "NA" in pheno_value:
                    continue

                if person_id in mrna_item[KEY_MRNA_VALUES]:
                    the_value = float(mrna_item[KEY_MRNA_VALUES][person_id])
                    the_pheno_value = float(pheno_value)
                    if not math.isnan(the_value) and not math.isnan(the_pheno_value):
                        mrna_values.append(the_value)
                        protein_values.append(the_pheno_value)

            pearson = "NA"
            if len(mrna_values) and len(protein_values):
                if len(mrna_values) == len(protein_values):
                    pearson = numpy.corrcoef(mrna_values, protein_values)[0,1]
                else:
                    raise Exception("Data does not comply")

            line = []
            line.append(gene+"-"+protein)
            line.append(pearson)
            line.append(len(mrna_values))
            correlation.append(line)

    return correlation

def PrintCorrelation(correlation,file_name):
    with open(file_name, "w+") as file:
        for line in correlation:
            text = line[0] + " " + str(line[1]) + " " + str(line[2]) + "\n"
            file.write(text)

def PrintGene(gene_name,MRNA, people_by_predixcan_row, prefix="Out/predixcandata-"):
    gene_item = MRNA[gene_name]
    with open(prefix+gene_name+".txt", "w+") as file:
        for num, person in people_by_predixcan_row.iteritems():
            person_id = person[KEY_PERSON_ID]
            values = gene_item[KEY_MRNA_VALUES]
            value = values[person_id]
            text = person_id + " " + value + "\n"
            file.write(text)

#-----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Correlations between proteins and mrna.')
    parser.add_argument("--gene_to_protein",
                        help="file with gene-to-protein relationships",
                        default="Intermediate/HauseGeneToProtein.txt")
    parser.add_argument("--predixcan_samples",
                        help="file with predixcan samples",
                        default = "Data/predixcan-samples.txt")
    parser.add_argument("--fam_file",
                        help="file with pheno fam information",
                        default="Intermediate/hapmap_r23a.fam")
    parser.add_argument("--predixcan_file",
                        help="file with predixcan data",
                        default="Data/predixcan-results.csv")
    parser.add_argument("--pheno_prefix",
                        help="prefix for phenotype input files",
                        default="./Intermediate/pheno_")
    parser.add_argument("--out",
                        help="name of file where the correlation (protein average) will be output",
                        default="./Out/hause-mrna-protein-correlation.txt")
    parser.add_argument("--out_ext",
                        help="name of file where the correlation will be output",
                        default="./Out/hause-mrna-protein-ext-correlation.txt")
    parser.add_argument("--trace_gene",
                        help="gene to be output",
                        default=None)
    args = parser.parse_args()

    #-------------------------------------------------------------------------------------------------------------------
    gene_to_protein_file = args.gene_to_protein
    gene_to_protein = BuildGeneToProteinRelationShip(gene_to_protein_file)

    #-------------------------------------------------------------------------------------------------------------------
    predixcan_sample_file = args.predixcan_samples
    fam_file = args.fam_file
    people, people_by_predixcan_row = BuildPeopleInPhenoList(fam_file, predixcan_sample_file)

    #-------------------------------------------------------------------------------------------------------------------
    predixcan_file = args.predixcan_file
    MRNA = BuildMRNAData(people_by_predixcan_row, gene_to_protein, predixcan_file)

    #-------------------------------------------------------------------------------------------------------------------
    pheno_prefix = args.pheno_prefix
    correlation = MRNAProteinAverageCorrelation(MRNA, gene_to_protein, pheno_prefix)

    output = args.out
    PrintCorrelation(correlation, output)

    #-------------------------------------------------------------------------------------------------------------------
    correlation_ext = MRNAProteinCorrelation(MRNA, gene_to_protein, pheno_prefix)
    output_ext = args.out_ext
    PrintCorrelation(correlation_ext, output_ext)

    #-------------------------------------------------------------------------------------------------------------------
    if args.trace_gene is not None and len(args.trace_gene):
        PrintGene(args.trace_gene, MRNA, people_by_predixcan_row)



