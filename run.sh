#! /bin/bash

python process_pheno.py

plink --noweb --bfile Data/hapmap_r23a --maf 0.05 --keep Intermediate/selected.fam --make-bed --out Intermediate/hapmap_r23a

gcta64 --bfile Intermediate/hapmap_r23a --autosome --make-grm --out Intermediate/hapmap_r23a
