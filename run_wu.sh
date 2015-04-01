#! /bin/bash

rm -rf IntermediateWu
mkdir IntermediateWu

python process_pheno_wu.py --select_individuals

plink --noweb --bfile Data/hapmap_r23a --maf 0.05 --keep IntermediateWu/selected.fam --make-bed --out IntermediateWu/hapmap_r23a

gcta64 --bfile IntermediateWu/hapmap_r23a --autosome --make-grm --out IntermediateWu/hapmap_r23a

python process_pheno_wu.py --build_phenotypes

for filename in ./IntermediateWu/pheno*.phen; do
  inname=${filename:2}
  outname="IntermediateWu/reml_${inname:21:-5}"
  gcta64 --grm IntermediateWu/hapmap_r23a --pheno $inname --reml --out $outname --thread-num 4
done

python process_reml.py --file_input_prefix './IntermediateWu/reml*' --reml_output Out/reml_results_wu.csv