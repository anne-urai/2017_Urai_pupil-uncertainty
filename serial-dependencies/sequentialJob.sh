#!/bin/sh

# Use as: qsub sequentialJob 1, where 1 is the subject number. 
# Can easily be done in a loop: for i in {1..27}, do qsub -v sjnum=$i sequentialJob.sh; sleep 10; done

#PBS -q batch
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1 
#PBS -l pmem=1gb

# -- run in the current working (submission) directory --
cd $PBS_O_WORKDIR
chmod g=wx $PBS_JOBNAME

filename=$(printf "data/2ifc_rt_sj%02d.txt" $sjnum) # get the name based on the input argument
echo $filename

python2.7 analysis.py -fr -n1000 -p ~/Data/pupilUncertainty/Data/serialmodel/ \
$filename 1> ~/jobs/$PBS_JOBID.out 2> ~/jobs/$PBS_JOBID.err 

