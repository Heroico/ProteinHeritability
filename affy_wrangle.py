import csv

LogLevel = 0

def LoadIdToGene(file_name):
    id_to_gene_name = {}
    gene_name_to_id = {}
    repeated_ids = []
    repeated_genes = []
    NULL_genes = []
    with open(file_name, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for row in reader:
            if reader.line_num == 1:
                continue

            id = row[0]
            gene = row[1]

            valid = True
            if id in id_to_gene_name:
                valid = False
                repeated_ids.append(id)
                if LogLevel > 1: print "Repeated id:", id, " gene:", gene

            if gene in gene_name_to_id:
                valid = False
                repeated_genes.append(gene)
                if LogLevel > 1: print "Repeated gene:", gene, "id:", id

            if gene == "NULL":
                valid = False
                NULL_genes.append(gene)
                if LogLevel > 1: print "NULL gene, id:", id

            if valid:
                id_to_gene_name[id] = gene
                gene_name_to_id[gene] = id

    if LogLevel > 0: print "id to genes:",len(id_to_gene_name.keys())
    if LogLevel > 0: print "gene to ids:",len(gene_name_to_id.keys())
    if LogLevel > 0: print "Repeated ids:",len(repeated_ids)
    if LogLevel > 0: print "Repeated genes:",len(repeated_genes)
    if LogLevel > 0: print "NULL genes:",len(NULL_genes)

    return id_to_gene_name, gene_name_to_id

def LoadMatrix(id_to_gene_name, file_name):
    matrix_rows = []
    people = []
    header_row = []
    original_header = []
    original_first_row = []
    skipped = []
    valid_columns = []
    with open(file_name, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='"')
        for row in reader:
            if reader.line_num == 1:
                original_header = row
                cols = len(row)
                for i in xrange(1,cols-1):
                    id = row[i].translate(None, '"X')
                    if id in id_to_gene_name:
                        header_row.append(id_to_gene_name[id])
                        valid_columns.append(i+1) #remember that the file has the row name at the beggining, thewn the personid
                    else:
                        skipped.append(id)
                matrix_rows.append(header_row)
                if LogLevel > 1: print "Header:",len(header_row)
            else:
                if reader.line_num == 2:
                    original_first_row = row

                person_id = "NA"+row[1]
                people.append(person_id)

                person_row = []
                for i in valid_columns:
                    person_row.append(row[i])
                matrix_rows.append(person_row)
                if LogLevel > 1: print "Row length:", len(person_row)


    if LogLevel > 1: print "ValidColumns", len(valid_columns)
    if LogLevel > 1: print "Skipped:",len(skipped)
    if LogLevel > 2: print "\t".join(skipped)
    if LogLevel > 2: print "o header:", "\t".join(original_header[1:])
    if LogLevel > 2: print "o 1  row:", "\t".join(original_first_row[2:])
    if LogLevel > 2: print "  header:", "\t".join(header_row)
    if LogLevel > 2: print "  1  row:", "\t".join(matrix_rows[1])

    return matrix_rows, people

def PrintRows(rows, file_name):
    with open(file_name, 'w+') as out:
        for row in rows:
            out.write(" ".join(row)+"\n")

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Process Affymetrix data and tame it.')
    parser.add_argument("--matrix_file_name",
                        help="file with gene expresion matrix",
                        default="Data/affy-gen-protein-matrix.txt")
    parser.add_argument("--id_to_gene_file_name",
                        help="Affymetrix id to gene relationship file",
                        default="Data/affy-exon-array-trcid_to_gene.txt")
    parser.add_argument("--log_level",
                        help="debug log level, set a number",
                        default="0")
    parser.add_argument("--matrix_output",
                        help="output file with gene expresion matrix",
                        default="Data/affy-results.txt")
    parser.add_argument("--people_output",
                        help="output file specifying ids of each matrix row person",
                        default="Data/affy-samples.txt")
    args = parser.parse_args()

    LogLevel = int(args.log_level)

    id_to_gene_file_name = args.id_to_gene_file_name
    id_to_gene_name, gene_name_to_id = LoadIdToGene(id_to_gene_file_name)

    #-------------------------------------------------------------------------------------------------------------------
    matrix_file_name = args.matrix_file_name
    matrix_rows, people = LoadMatrix(id_to_gene_name, matrix_file_name)

    #-------------------------------------------------------------------------------------------------------------------
    output_matrix_file_name = args.matrix_output
    PrintRows(matrix_rows, output_matrix_file_name)

    output_people_file_name = args.people_output
    PrintRows(people, output_people_file_name)
