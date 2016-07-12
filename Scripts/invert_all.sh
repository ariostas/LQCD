#!/bin/bash

cd /home/arios/Documents/LQCDConfigs/cl3_16_48_b6p1_m0p2450

export LD_LIBRARY_PATH="/usr/local/cuda/lib64/:/usr/local/cuda/nvvm/lib64/"

for file in *.ini.xml
do
	echo "inverting config ${file:0:-8}"
    chroma -i ${file} -o ${file:0:-8}.out.xml
done
