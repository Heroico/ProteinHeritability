__author__ = 'heroico'

import csv

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



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Correlations between proteins and mrna.')
    parser.add_argument("--gene_to_protein",
                        help="file with gene-to-protein relationships")
    parser.add_argument("--predixcan_samples",
                        help="file with predixcan samples")
    parser.add_argument("--fam_file",
                        help="file with pheno fam information")
    args = parser.parse_args()

    gene_to_protein_file = "Intermediate/HauseGeneToProtein.txt"
    if args.gene_to_protein and len(args.gene_to_protein):
        gene_to_protein_file = args.gene_to_protein
    gene_to_protein = BuildGeneToProteinRelationShip(gene_to_protein_file)

    predixcan_sample_file = "Data/predixcan-samples.txt"
    if args.predixcan_samples and len(args.predixcan_samples):
        predixcan_sample_file = args.predixcan_samples

    fam_file = "Intermediate/hapmap_r23a.fam"
    if args.fam_file and len(args.fam_file):
        fam_file = args.fam_file

    people, people_by_predixcan_row = BuildPeopleInPhenoList(fam_file, predixcan_sample_file)