#!/bin/bash

START_DIR=$PWD

cd /data/d02/arios/wil_16_60_aniso_cfgs/bar3ptfn_nucleon_pol/
mkdir bar3ptfn
cd bar3ptfn

echo "Stripping bar3ptfn data"
w_src=0
w_snk=0
mkdir temp
cd temp
for x in {0..4}
do
	for y in {0..4}
	do
		/home/arios/chroma_utils/install/bin/strip_bar3ptfn ../../bar3ptfn_sh_${x}_sh_${y}_*.dat
		for file in *
		do
			if [ ${x} -eq 0 ]
			then
				w_src=0
			elif [ ${x} -eq 1 ] || [ ${x} -eq 2 ]
			then
				w_src=1
			elif [ ${x} -eq 3 ] || [ ${x} -eq 4 ]
			then
				w_src=2
			fi

			if [ ${y} -eq 0 ]
			then
				w_snk=0
			elif [ ${y} -eq 1 ] || [ ${y} -eq 2 ]
			then
				w_snk=1
			elif [ ${y} -eq 3 ] || [ ${y} -eq 4 ]
			then
				w_snk=2
			fi
			mv ${file} "${file[@]%.SS}.DG${w_src}_1.DG${w_snk}_1.sh_${x}_sh_${y}.SS"
		done
		mv * ../
	done
done
cd ..
rm -rf temp

cd ${START_DIR}
