#! /bin/bash
#PBS -l nodes=4:ppn=24
#PBS -l walltime=3000:00:00
#PBS -l cput=3000:00:00
#PBS -j oe
#PBS -u arios

set -x

TOPDIR=/data/d02/arios/wil_16_60_aniso_cfgs/

cd ${TOPDIR}

. /opt/intel/composer_xe_2013.1.117/bin/compilervars.sh intel64
. /opt/intel/impi/4.1.0.024/intel64/bin/mpivars.sh

GEOM="2 2 2 12"
APRUN="mpirun -np 96 -perhost 24"
PROG=/home/arios/chroma/install/bin/chroma

[ -x $PROG ] || exit 1

linenum=1

while read line
do
    array=($line)
    if [[ "${array[0]}" != "#"* ]]; then 

        echo "" | ${APRUN} ${PROG} -i /data/d02/arios/wil_16_60_aniso_cfgs/${array[0]} -o /data/d02/arios/wil_16_60_aniso_cfgs/${array[0]%%.*}.out.xml -geom ${GEOM} 2>&1 | tee out_SEED.stdout
        sedstring="${linenum}s/^/#/"
        sed -i -e "${sedstring}" SEED

    fi

    linenum=$((linenum+1))

done < SEED

rm SEED
