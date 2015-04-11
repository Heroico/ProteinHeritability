__author__ = 'heroico'

#..............................................................................
# TODO:
# maybe discard NANS
#...........................................................................

import csv
import process_pheno

COL_INDIVIDUAL_LCL = 0
COL_SEX = 1
COL_POPULATION = 2
COL_FAMILY_ID = 3
COL_INDIVIDUAL_ID = 4
COL_DAD = 5
COL_MOM = 6

KEY_PHENO_NAME = "name"
KEY_PHENO_VALUES = "values"
KEY_PHENO_COLUMN = "column"

KEY_INDIVIDUAL_POPULATION = "population"
KEY_INDIVIDUAL_FAMILY_ID = "family_id"
KEY_INDIVIDUAL_ID = "individual_id"

def ReadWuPhenoInput(file_name='Data/pheno-protein-wu.csv',keep_only_europe=True):
    pheno_individuals = []
    phenos = {}
    with open(file_name, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        read_first = False
        for row in reader:
            #header row: some columns we'll skip, and then an individual name in each column
            if not read_first and len(row):
                read_first = True
                for index, item in enumerate(row):
                    if index <= COL_MOM:
                        pass
                    else:
                        pheno = {
                            KEY_PHENO_NAME:item,
                            KEY_PHENO_VALUES:{},
                            KEY_PHENO_COLUMN:index
                            }
                        phenos[item] = pheno

            # data row: get each value into its individual
            elif len(row):
                individual_info = {
                    KEY_INDIVIDUAL_POPULATION: row[COL_POPULATION],
                    KEY_INDIVIDUAL_FAMILY_ID: row[COL_FAMILY_ID],
                    KEY_INDIVIDUAL_ID: row[COL_INDIVIDUAL_ID]
                    }

                if keep_only_europe and not individual_info[KEY_INDIVIDUAL_POPULATION] == "CEU":
                    continue

                pheno_individuals.append(individual_info)
                for key,pheno in phenos.iteritems():
                    column = pheno[KEY_PHENO_COLUMN]
                    values = pheno[KEY_PHENO_VALUES]
                    values[individual_info[KEY_INDIVIDUAL_ID]] = row[column]

    return pheno_individuals, phenos

def FamPhenoIntersection(fam, pheno_individuals):
    intersection = []
    for item in fam:
        #print "'"+item[COL_FAM_IND]+"'"
        match = [x for x in pheno_individuals if  x[KEY_INDIVIDUAL_ID] == item[process_pheno.COL_FAM_IND]]
        if len(match) == 0:
            pass
        elif len(match) == 1:
            intersection.append(item)
        else:
            raise Exception('Cant handle duplicates')

    return intersection


from gene_code import ReadGeneCodeInput
from gene_code import KEY_GENE_NAME as GENCODE_KEY_GENE_NAME
from gene_code import KEY_ENSEMBLE_ID as GENCODE_KEY_ENSEMBLE_ID

def PrintWuGeneToProteinList(phenos, gene_codes, file_name= "IntermediateWu/WuGeneToProtein.txt"):
    output = {}
    for name, pheno in phenos.iteritems():
        ensemble = pheno[KEY_PHENO_NAME]
        if ensemble in gene_codes:
            gene_code = gene_codes[ensemble]
            line = None
            if not ensemble in output:
                line = []
                line.append(gene_code[GENCODE_KEY_GENE_NAME])
                line.append(ensemble)
                output[ensemble] = line
            else:
                line = output[ensemble]
            line.append(ensemble)
        else:
            raise Exception('Need gene code for every ensemble')
    output_lines = []
    for gene,line in output.iteritems():
        output_lines.append(line)
    output_lines =sorted(output_lines, key=lambda line: line[0])

    with open(file_name, 'w+') as out:
        for line in output_lines:
            out.write(" ".join(line)+"\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process phenotype from Wu for information needed further down in the chain.')
    parser.add_argument("--select_individuals",
                        help="build a selection of individuals to use with plink",
                        action='store_true')
    parser.add_argument("--build_phenotypes",
                        help="build phenotype files from intermediate",
                        action='store_true')
    parser.add_argument("--gene_to_protein_mode",
                       help="Will only build a table mapping gene to protein(must pass ensemble output)",
                       action='store_true')
    parser.add_argument("--gene_to_protein_output_name",
                        help="file output name",
                        default="IntermediateWu/WuGeneToProtein.txt")
    args = parser.parse_args()

    pheno_individuals, phenos = ReadWuPhenoInput()
    fam = process_pheno.ReadFamInput('Data/hapmap_r23a.fam')


    if args.gene_to_protein_mode:
        file_name = args.gene_to_protein_output_name
        gene_codes = ReadGeneCodeInput()
        PrintWuGeneToProteinList(phenos, gene_codes, file_name)

    if args.select_individuals:
        intersection = FamPhenoIntersection(fam, pheno_individuals)
        process_pheno.PrintIntersectionFamFile(intersection, file_name="IntermediateWu/selected.fam")

    if args.build_phenotypes:
        grm_ids = process_pheno.ReadGRMIDInput('IntermediateWu/hapmap_r23a.grm.id')
        process_pheno.PrintPhenotypeFiles(phenos,grm_ids,file_name_prefix='IntermediateWu/pheno_')