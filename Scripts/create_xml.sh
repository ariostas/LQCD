#!/bin/bash

cd /home/arios/Documents/LQCDConfigs/wil_16_64_aniso

for file in *.lime
do
	echo "Making config for ${file:0:-5}"
	sed "s/CONFIGNAME/${file:0:-5}/g" /home/arios/Documents/temp_scripts/bar3ptfnSEED.ini.xml > bar3ptfn_${file:0:-5}.ini.xml	
done

cd -