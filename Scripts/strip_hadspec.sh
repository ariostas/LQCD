#!/bin/bash

START_DIR=$PWD

cd /home/arios/Documents/LQCDConfigs/cl3_16_48_b6p1_m0p2450
mkdir hadspec
cd hadspec

echo "Stripping hadspec data"

for x in {0..2}
do
	for y in {0..2}
	do
		/home/arios/ChromaCuda/chroma_utils/install/bin/strip_hadspec ../hadspec_sh_${x}_sh_${y}_*.dat.xml
	done
done

cd ${START_DIR}
