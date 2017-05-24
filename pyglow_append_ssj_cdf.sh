#!/bin/bash
YEAR=2003
SAT=15
END=365
for DOY in $(seq 2 $END); do python pyglow_append_ssj_cdf.py $SAT $YEAR $DOY --clobber; done