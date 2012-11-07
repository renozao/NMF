#!/bin/bash
#$ -cwd 
#$ -q opteron.q
#$ -pe mpich 5
#$ -M renaud@cbio.uct.ac.za
#$ -m bea

echo "Got $NSLOTS slots. $TMP/machines"

/opt/openmpi/bin/orterun -v -n $NSLOTS -hostfile $TMP/machines /usr/bin/R --slave -f mpi.R
