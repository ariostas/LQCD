#!/bin/bash

cd /home/arios/Documents/LQCDConfigs/wil_16_64_aniso

export LD_LIBRARY_PATH="/usr/local/cuda/lib64/:/usr/local/cuda/nvvm/lib64/"

for file in *.ini.xml
do
	echo "inverting config ${file:0:-8}"
    chroma-jit -i ${file} -o ${file:0:-8}.out.xml
done
