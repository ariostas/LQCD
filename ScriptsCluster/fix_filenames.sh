#!/bin/bash

cd /data/d02/arios/wil_16_60_aniso_cfgs

for file in *.lime*
do
	echo "Fixing filename for ${file}"
	mv ${file} ${file%%.*}_${file##*e}.lime
done

cd -