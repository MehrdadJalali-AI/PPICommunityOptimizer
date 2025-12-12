#!/bin/bash
# Run pipeline on Gavin PPI + SGD GO data

python main.py \
  --mode gavin \
  --ppi gavin2006_socioaffinities_rescaled.txt \
  --go-file GO.txt \
  --go-use-symbol \
  --go-taxid 559292 \
  --alpha 0.5 \
  --overlap-tau 0.1 \
  --outdir outputs_gavin/ \
  --mcl-inflation 2.0 \
  --lea-population 30 \
  --lea-evaluations 500 \
  --random-seed 42

