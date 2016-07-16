#!/bin/bash

START_DIR=$PWD

cd /home/arios/Documents/LQCDConfigs/wil_16_64_aniso
mkdir bar3ptfn
cd bar3ptfn

echo "Stripping bar3ptfn data"
mkdir temp
cd temp
for x in {0..2}
do
	for y in {0..2}
	do
		/home/arios/ChromaCuda/chroma_utils/install/bin/strip_bar3ptfn ../../bar3ptfn_sh_${x}_sh_${y}_*.dat
		for file in *
		do
			mv ${file} ${file}.DG$((${x}+1))_1.DG$((${y}+1))_1.SS
		done
		mv * ../
	done
done
cd ..
rm -rf temp

cd ${START_DIR}
