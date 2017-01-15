#!/bin/bash

cd /home/arios/Documents/LQCDConfigs/wil_16_60_aniso

for file in *.lime*
do
	echo "fixing name for ${file}"
	mv ${file} ${file%.*}_${file##*e}.lime
done

cd -