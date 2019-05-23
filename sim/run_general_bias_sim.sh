#!/bin/bash
#PBS -N hpy_bias
#PBS -o log/
#PBS -e log/
#PBS -l mem=3G
#PBS -l vmem=3G
#PBS -l walltime=00:20:00
#PBS -l nodes=1:ppn=10

cd $PBS_O_WORKDIR
mkdir -p log/

# First argument is filename for HPY parameters
FNAME=$1

echo Started: `date`
Rscript general_bias_sim.R $FNAME
echo Ended: `date`
