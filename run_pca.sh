
plink --noweb --bfile Data/hapmap_r23a --maf 0.05  --make-bed --out IntermediatePCA/hapmap_r23a


gcta64 --bfile IntermediatePCA/hapmap_r23a --autosome --make-grm --out IntermediatePCA/hapmap_r23a