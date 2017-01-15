#! /bin/bash
#PBS -l nodes=4:ppn=24
#PBS -l walltime=2000:00:00
#PBS -l cput=2000:00:00
#PBS -j oe

set -x

TOPDIR=/data/d02/arios/wil_16_60_aniso_cfgs/

cd ${TOPDIR}

# straight mpi runs - better but not great performance
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/wdetmold/opt/lib
.  /opt/intel/composer_xe_2013.1.117/bin/compilervars.sh intel64
. /opt/intel/impi/4.1.0.024/intel64/bin/mpivars.sh
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/software/lib/:$HOME/opt/lib/:$HOME/opt/usqcd/install/platypus/chroma/lib:$HOME/wdetmold/opt/lib/

GEOM="2 2 2 12"
APRUN="mpirun -np 96 -perhost 24"
PROG=/home/arios/chroma/install/bin/purgaug
INPUT=/home/arios/scripts/purgaug.ini.xml
# INPUT=su3_16_60_beta_5p8.ini.xml200

[ -x $PROG ] || exit 1

${APRUN} ${PROG} -i ${INPUT} -geom ${GEOM} 2>&1 | tee err.stdout 

