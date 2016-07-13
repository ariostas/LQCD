#!/bin/bash

START_DIR=$PWD

cd /home/arios/Documents/LQCDConfigs/cl3_16_48_b6p1_m0p2450
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
			mv ${file} ${file}.DG${x}_1.DG${y}_1.SS
		done
		mv * ../
	done
done
cd ..
rm -rf temp

cd ${START_DIR}
