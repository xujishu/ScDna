#$ -S /bin/bash
#$ -o /seq/mint/jishuxu/logs/
#$ -j y
#$ -l h_rt=1000:0:00
#$ -l h_vmem=4G
#$ -t 1-15
source /broad/software/scripts/useuse
use UGER
use .samtools-1.7

export PATH=${PATH}:'/seq/mint/jishuxu/softwares/miniMap2/minimap2/'
infile='/seq/mint/jishuxu/trinity/scRuns/run_list.txt'
line=`awk "NR==$SGE_TASK_ID" $infile`
fastafile=$(echo $line|cut -f 2 -d",")
cellname=$(echo $line|cut -f 1 -d",")
fastadir="/seq/mint/jishuxu/trinity/scRuns/fasta/"
echo ${fastadir}/$fastafile
# trinity summary
TRINITYHOME='/home/unix/jishuxu/mint/github/trinityrnaseq/util'
output="/home/unix/jishuxu/mint/trinity/scRuns/trinity_summary/${cellname}.log"
${TRINITYHOME}/TrinityStats.pl ${fastadir}/$fastafile >$output
## samtools idxstats
bamfile="/seq/mint/jishuxu/trinity/scRuns/minimap2/${cellname}_T.bam"
samtools idxstats ${bamfile} >"/home/unix/jishuxu/mint/trinity/scRuns/trinity_summary/${cellname}_minimap2_T.txt" 
bamfile="/seq/mint/jishuxu/trinity/scRuns/minimap2/${cellname}_GT.bam"
samtools idxstats ${bamfile} >"/home/unix/jishuxu/mint/trinity/scRuns/trinity_summary/${cellname}_minimap2_GT.txt"
bamfile="/seq/mint/jishuxu/trinity/scRuns/minimap2/${cellname}_G.bam"
samtools idxstats ${bamfile} >"/home/unix/jishuxu/mint/trinity/scRuns/trinity_summary/${cellname}_minimap2_G.txt"

#samtools flagstat
bamfile="/seq/mint/jishuxu/trinity/scRuns/minimap2/${cellname}_T.bam"
samtools flagstat ${bamfile} >"/home/unix/jishuxu/mint/trinity/scRuns/trinity_summary/${cellname}_minimap2_flagstat_T.txt" 
bamfile="/seq/mint/jishuxu/trinity/scRuns/minimap2/${cellname}_GT.bam"
samtools flagstat ${bamfile} >"/home/unix/jishuxu/mint/trinity/scRuns/trinity_summary/${cellname}_minimap2_flagstat_GT.txt"
bamfile="/seq/mint/jishuxu/trinity/scRuns/minimap2/${cellname}_G.bam"
samtools flagstat ${bamfile} >"/home/unix/jishuxu/mint/trinity/scRuns/trinity_summary/${cellname}_minimap2_flagstat_G.txt"
