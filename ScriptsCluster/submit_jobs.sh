#!/bin/bash

LINES=100

cd /data/d02/arios/wil_16_60_aniso_cfgs/

for file in *.lime*
do
	echo "Making xml file for ${file%%.*}"
	sed "s/CONFIGNAME/${file%%.*}/g" /home/arios/scripts/xmlSEED.ini.xml > xml_${file%%.*}.ini.xml
done

ls *.ini.xml > xml_files.txt
split -l "${LINES}" xml_files.txt xmls
rm xml_files.txt

counter=0
countern=${LINES}
for file in xmls*
do
    echo "Making submit script for cfgs_$((counter+1))-${countern}"
    sed "s/SEED/${file}/g" /home/arios/scripts/invertSEED.sh > invert_${file}.sh
    chmod +x invert_${file}.sh
    echo "Submitting job for cfgs_$((counter+1))-${countern}"
    # qsub -N "cfgs_$((counter+1))-${countern}" invert_${file}.sh
    counter=${countern}
    countern=$((counter+LINES))
done

cd -