#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -o /seq/mint/jishuxu/logs/
#$ -j y
#$ -l h_rt=1000:0:00
#$ -l h_vmem=8G
#$ -t 1-15
source /broad/software/scripts/useuse
use UGER
use .samtools-1.7
use .perl-5.20.3
export PATH=/broad/dejagerlab/cogdec/jishu/softwares/Miniconda3/bin:$PATH

TDHOME='/home/unix/jishuxu/mint/github/TransDecoder'
infile='/seq/mint/jishuxu/trinity/scRuns/run_list.txt'
line=`awk "NR==$SGE_TASK_ID" $infile`
fastafile=$(echo $line|cut -f 2 -d",")
cellname=$(echo $line|cut -f 1 -d",")
fastadir="/seq/mint/jishuxu/trinity/scRuns/fasta/"
echo ${fastadir}/${fastafile}
cd /home/unix/jishuxu/mint/trinity/scRuns/TransDecoder
${TDHOME}/TransDecoder.LongOrfs  -t ${fastadir}/${fastafile}
${TDHOME}/TransDecoder.Predict -t ${fastadir}/${fastafile}
