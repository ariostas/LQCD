#!/bin/bash

cd /home/arios/Documents/LQCDConfigs/cl3_16_48_b6p1_m0p2450

for file in *.lime
do
	echo "Making config for ${file:0:-5}"
	sed "s/CONFIGNAME/${file:0:-5}/g" /home/arios/Documents/LQCD/Scripts/bar3ptfnSEED.ini.xml > bar3ptfn_${file:0:-5}.ini.xml	
done

cd -