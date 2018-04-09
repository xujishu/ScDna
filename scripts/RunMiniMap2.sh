#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -o /seq/mint/jishuxu/logs/
#$ -j y
#$ -l h_rt=1000:0:00
#$ -l h_vmem=50G
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
rg="@RG\tID:${cellname}\tSM:${cellname}" 
#outdir="/seq/mint/jishuxu/trinity/scRuns/minimap2/"
outdir="/broad/hptmp/jishuxu/trinity/minimap2/"
minimap2 -ax splice -R ${rg}  /home/unix/jishuxu/mint/trinity/blastdb/GeTr/GRCh38.primary_assembly.genome.fa  ${fastadir}/${fastafile} |samtools view -Shu - |samtools sort -T /broad/hptmp/jishuxu/tmp/ - -o ${outdir}/${cellname}_G.bam 
samtools index ${outdir}/${cellname}_G.bam

minimap2 -ax splice  -R ${rg} /home/unix/jishuxu/mint/trinity/blastdb/GeTr/gencode.v27.transcripts.fa  ${fastadir}/${fastafile} |samtools view -Shu - |samtools sort -T /broad/hptmp/jishuxu/tmp/ - -o ${outdir}/${cellname}_T.bam 
samtools index ${outdir}/${cellname}_T.bam

minimap2 -ax splice -R ${rg} /home/unix/jishuxu/mint/trinity/blastdb/GeTr/GRCh38.Transcriptome.fa  ${fastadir}/${fastafile} |samtools view -Shu - |samtools sort -T /broad/hptmp/jishuxu/tmp/ - -o ${outdir}/${cellname}_GT.bam    
samtools index ${outdir}/${cellname}_GT.bam

