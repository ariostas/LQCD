#!/bin/bash

START_DIR=$PWD

cd /home/arios/Documents/LQCDConfigs/wil_16_60_aniso/hadspec_t21/
mkdir hadspec
cd hadspec

echo "Stripping hadspec data"

mkdir temp
cd temp

for x in {0..4}
do
	for y in {0..4}
	do
		/home/arios/ChromaCuda/chroma_utils/install/bin/strip_hadspec ../../hadspec_sh_${x}_sh_${y}_*.dat.xml

		for file in *
		do
			mv ${file} ../${file:0:-3}.sh_${x}_sh_${y}.SS
		done
	done
done

cd ..
rm -r temp

cd ${START_DIR}
