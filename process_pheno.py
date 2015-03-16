__author__ = 'heroico'

import csv

COL_ID = 0
COL_Protein = 1
COL_Gene = 2
COL_Ensembl_ID = 3
COL_Platform = 4
COL_Ab_Type = 5
COL_Dilution = 6
COL_Flag = 7
COL_Median_BCII = 8
COL_Median_SNR = 9
COL_Median_Signal_CV = 10
COL_Antibody_Quality = 11

KEY_PHENO_NAME = "name"
KEY_PHENO_VALUES = "values"

def ReadHausePhenoInput(file_name='Data/mmc5-i.csv'):
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
                    if index <= COL_Antibody_Quality:
                        pass
                    else:
                        pheno_individuals.append(item)
            # data row: get each value into its individual
            elif len(row):
                id = None
                protein = None
                for index, pheno_value in enumerate(row):
                    if index <= COL_Antibody_Quality:
                        if index == COL_ID:
                            id = pheno_value
                        if index == COL_Protein:
                            protein = pheno_value
                        pass
                    else:
                        name=id+"_"+protein
                        pheno = None
                        if not name in phenos:
                            pheno = {KEY_PHENO_NAME:name, KEY_PHENO_VALUES:{}}
                            phenos[name] = pheno
                        else:
                            pheno = phenos[name]
                        item = pheno_individuals[index-COL_Antibody_Quality-1]
                        values = pheno[KEY_PHENO_VALUES]
                        values[item] = pheno_value
    return pheno_individuals, phenos


COL_FAM_FAM = 0
COL_FAM_IND = 1
def ReadFamInput(file_name='Data/hapmap_r23a.fam'):
    fam = []
    with open(file_name, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='"')
        for row in reader:
            item = [row[COL_FAM_FAM], row[COL_FAM_IND]]
            fam.append(item)

    return fam

COL_GRM_FAM_ID = 0
COL_GRM_IND_ID = 1
def ReadGRMIDInput(file_name='Intermediate/hapmap_r23a.grm.id'):
    ids = []
    with open(file_name, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for row in reader:
            if len(row) == 2:
                item = [row[COL_GRM_FAM_ID], row[COL_GRM_IND_ID]]
                ids.append(item)
            else:
                raise Exception("Cant handle GRM Ids")
    return ids

def FamPhenoIntersection(fam, pheno_individuals):
    intersection = []
    for item in fam:
        #print "'"+item[COL_FAM_IND]+"'"
        match = [x for x in pheno_individuals if  x == item[COL_FAM_IND]]
        if len(match) == 0:
            pass
        elif len(match) == 1:
            intersection.append(item)
        else:
            raise Exception('Cant handle duplicates')

    return intersection

def PrintIntersectionFamFile(intersection,file_name='Intermediate/selected.fam'):
    with open(file_name, 'w') as out:
        for item in intersection:
            out.write(" ".join(item) + '\n')

def PrintPhenotypeFiles(phenos,grm_ids,file_name_prefix='Intermediate/pheno_'):
    for name, pheno in phenos.iteritems():
        pheno_name = name.replace('(','_').replace(')','_').replace('/','-')
        file_name = file_name_prefix+pheno_name+'.phen'
        with open(file_name, 'w+') as out:
            values = pheno[KEY_PHENO_VALUES]
            for ind in grm_ids:
                line = ind[COL_GRM_FAM_ID]+ " " + ind[COL_GRM_IND_ID] + " " + values[ind[COL_GRM_IND_ID]] + "\n"
                out.write(line)

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Process phenotype Hause for information needed further down in the chain.')
    parser.add_argument("--select_individuals",
                        help="build a selection of individuals to use with plink",
                        action='store_true')
    parser.add_argument("--build_phenotypes",
                        help="build phenotype files from intermediate",
                        action='store_true')
    args = parser.parse_args()

    pheno_individuals, phenos = ReadHausePhenoInput('Data/mmc5-i.csv')
    fam = ReadFamInput('Data/hapmap_r23a.fam')

    intersection = None
    if args.select_individuals:
        intersection = FamPhenoIntersection(fam, pheno_individuals)
        PrintIntersectionFamFile(intersection)

    if args.build_phenotypes:
        grm_ids = ReadGRMIDInput('Intermediate/hapmap_r23a.grm.id')
        PrintPhenotypeFiles(phenos,grm_ids)