#!/bin/sh

# Use as: qsub sequentialJob 1, where 1 is the subject number. 
# Can easily be done in a loop: for i in {1..27}, do qsub -v sjnum=$i sequentialJob.sh; sleep 10; done

# embedded options to qsub - start with #PBS
# walltime: defines maximum lifetime of a job
# nodes/ppn: how many nodes? how many cores?

#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=2
#PBS -l pmem=2gbs

# -- run in the current working (submission) directory --
cd $PBS_O_WORKDIR
chmod g=wx $PBS_JOBNAME

matlab -nodisplay -nodesktop -r "a1_PupilAnalysis($sjnum); exit" 1> ~/jobs/$PBS_JOBID.out 2> ~/jobs/$PBS_JOBID.err;
