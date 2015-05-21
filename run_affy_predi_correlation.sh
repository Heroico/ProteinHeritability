if [[ ! -d "./IntermediateC" ]] ; then
    mkdir IntermediateC
fi

python predixcan_stats.py --intersection_output

cd Out

Rscript mrna_source_correlation.R