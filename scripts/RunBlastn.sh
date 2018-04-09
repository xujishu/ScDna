#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -o /broad/hptmp/jishuxu/log/
#$ -j y
#$ -l h_rt=1200:00:00
#$ -l h_vmem=100g
#$ -t 1-15
blastn='/seq/mint/jishuxu/softwares/ncbi-blast-2.7.1+/bin/blastn'
infile='/seq/mint/jishuxu/trinity/scRuns/run_list.txt'
line=`awk "NR==$SGE_TASK_ID" $infile`
fastafile=$(echo $line|cut -f 2 -d",")
cellname=$(echo $line|cut -f 1 -d",")
fastadir="/seq/mint/jishuxu/trinity/scRuns/fasta/"
echo ${fastadir}/$fastafile
outfile="/seq/mint/jishuxu/trinity/scRuns/blastn/${cellname}.GT.blastn.txt"

cmd="${blastn} -db  /seq/mint/jishuxu/trinity/blastdb/GeTr/GRCh38.Transcriptome.fa  -query ${fastadir}/$fastafile -outfmt 6 -num_alignments 10 -num_descriptions 1 -out ${outfile}"
echo ${cmd}
${cmd}
