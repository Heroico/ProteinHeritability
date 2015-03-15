#! /bin/bash

python process_pheno.py --select_individuals

plink --noweb --bfile Data/hapmap_r23a --maf 0.05 --keep Intermediate/selected.fam --make-bed --out Intermediate/hapmap_r23a

gcta64 --bfile Intermediate/hapmap_r23a --autosome --make-grm --out Intermediate/hapmap_r23a

python process_pheno.py --build_phenotypes

for filename in ./Intermediate/pheno*.phen; do
  inname=${filename:2}
  outname="Intermediate/reml_${inname:19:-5}"
  gcta64 --grm Intermediate/hapmap_r23a --pheno $inname --reml --out $outname --thread-num 4
done
