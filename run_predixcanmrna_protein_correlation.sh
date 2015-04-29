if [[ ! -d "./Intermediate" ]] ; then
    bash run.sh
fi

python process_pheno.py --gene_to_protein_mode

python mrna_gene_protein.py