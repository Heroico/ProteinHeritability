if [[ ! -d "./IntermediateC" ]] ; then
    mkdir IntermediateC
    python predixcan_stats.py --intersection_output
fi

cd Out

Rscript mrna_source_correlation.R