#!/bin/bash
for DOY in $(seq 1 366)
do
	python pyglow_append_ssj_cdf.py 16 $DOY
	python pyglow_append_ssj_cdf.py 17 $DOY
	python pyglow_append_ssj_cdf.py 18 $DOY
done
