#!/bin/bash

cd /home/arios/Documents/LQCDConfigs/wil_16_60_aniso

export LD_LIBRARY_PATH="/usr/local/cuda/lib64/:/usr/local/cuda/nvvm/lib64/"

echo "Producing configurations"
purgaug-jit -i /home/arios/Documents/temp_scripts/purgaug.ini.xml -o purgaug.out.xml
cd -
