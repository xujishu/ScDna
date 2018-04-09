#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -o /broad/hptmp/jishuxu/log/
#$ -j y
#$ -l h_rt=100:00:00
#$ -l h_vmem=8g
source /broad/software/scripts/useuse
use UGER

# download genome primary reference fasta and transcriptome fasta from gencode
cd /home/unix/jishuxu/mint/trinity/blastdb/GeTr/
wget  ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz ./
wget  ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz ./ 
zcat  *.fa.gz >GRCh38.Transcriptome.fa
/home/unix/jishuxu/mint/softwares/ncbi-blast-2.7.1+/bin/makeblastdb  -in GRCh38.Transcriptome.fa  -parse_seqids -dbtype nucl

