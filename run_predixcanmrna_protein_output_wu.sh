if [[ ! -d "./IntermediateWu" ]] ; then
    bash run.sh
fi

python process_pheno_wu.py --gene_to_protein_mode

python mrna_gene_protein.py \
--gene_to_protein 'IntermediateWu/WuGeneToProtein.txt' \
--fam_file 'IntermediateWu/hapmap_r23a.fam' \
--pheno_prefix 'IntermediateWu/pheno_' \
--out './Out/wu-mrna-protein-correlation.txt' \
--out_ext './Out/wu-mrna-protein-ext-correlation.txt' \
--csv_output_mode \
--csv_output_prefix_name './Out/wu_'