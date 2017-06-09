#!/bin/bash
###### Job name ######
#PBS -N planet
###### Output files ######
#PBS -o z.out
#PBS -e z.err
###### Number of nodes and cores ######
#PBS -l nodes=4:ppn=16
###### Queue name ######
#PBS -q medium
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh
module purge
module load openmpi/2.0.1_ic16.0

###### Run your jobs with parameters ######
if [ -n "$PBS_NODEFILE" ]; then
  if [ -f $PBS_NODEFILE ]; then
    NPROCS=`wc -l < $PBS_NODEFILE`
  fi
fi


###### Run your jobs with parameters ######
$OPENMPI_HOME/bin/mpirun -v -machinefile $PBS_NODEFILE -np $NPROCS ./pluto > this.out 

