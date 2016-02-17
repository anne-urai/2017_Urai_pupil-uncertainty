#!/bin/sh

# Use as: qsub sequentialJob 1, where 1 is the subject number. 
# Can easily be done in a loop: for i in {1..27}, do qsub -v sjnum=$i sequentialJob.sh; sleep 10; done

# on LISA, run this in parallel on nodes with 16 cores:
# qsub -v sjstart=1 sequentialJob.sh, where sjstart is the subject to begin on. Then check and resubmit at the next sj.

# embedded options to qsub - start with #PBS
# walltime: defines maximum lifetime of a job
# nodes/ppn: how many nodes? how many cores?

#PBS -S /bin/bash
#PBS -o ~/jobs
#PBS -e ~/jobs
#PBS -j oe
#PBS -q batch
#PBS -l walltime=100:00:00
#PBS -l nodes=1
#PBS -l pmem=1gbs

# -- run in the current working (submission) directory --
cd $PBS_O_WORKDIR
chmod g=wx $PBS_JOBNAME

for ((i=1; i<=16; i++)) ; do
(
                      filename=$(printf "$HOME/Dropbox/code/pupilUncertainty/serial-dependencies/data/2ifc_%s_sj%02d.txt" $whichFile $i); # get the name based on the input argument
                      python2.7 analysis.py -fr -n1000 $filename; # force rerun and dont plot
) &
done
wait