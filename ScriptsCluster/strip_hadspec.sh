#!/bin/bash

START_DIR=$PWD

cd /data/d02/arios/wil_16_60_aniso_cfgs/hadspec_nucleon_pol/
mkdir hadspec
cd hadspec

echo "Stripping hadspec data"

mkdir temp
cd temp

for x in {0..4}
do
	for y in {0..4}
	do
		/home/arios/chroma_utils/install/bin/strip_hadspec ../../hadspec_sh_${x}_sh_${y}_*.dat.xml

		for file in *
		do
			mv ${file} ../${file[@]%.SS}.sh_${x}_sh_${y}.SS
		done
	done
done

cd ..
rm -r temp

cd ${START_DIR}
