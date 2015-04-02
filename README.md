# ProteinHeritability

This repository contains several python scripts that can process protein phenotype data, and calculate its heritability.

The "run*.sh" are self contained bash scripts that command python scripts, gcta and plink to produce several calculations (most of them, protein heritability for the given data)

The python scripts themselves are general-purpose tools by themselves, although most of their code is input data wrangling.

# Assumptions

## Data

the user should create a subfolder within this repository, called:

```
data
```

containing the following files:

```
hapmap_r23a.bed
hapmap_r23a.bin
hapmap_r23a.fam
mmc5-i.csv #protein data from hauser
pheno-protein-wu.csv #protein data from Wu-Snyder
```

For correlating to mrna data, you need:

```
#generated with:
#python predict_gene_expression.py --dosages GEUVADIS/dosagefiles-hapmap2/ --weights DGN-WB_0.5.db --out out/results
predixcan-results.csv
#basically, samples files from GEUVADIS
predixcan-samples.txt
```

## GCTA

GCTA 1.24 should be on the shell path, and callable as gcta64
http://www.complextraitgenomics.com/software/gcta/

## PLINK

PLINK should be available on the shell path.
http://pngu.mgh.harvard.edu/~purcell/plink/

# Instructions

The protein heritability can be run as in:

```bash
$ bash run.sh
```

or:

```bash
$ bash run_wu.sh
```
