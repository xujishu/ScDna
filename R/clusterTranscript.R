#' Meta data 
source('~/Works/github/ScDna/R/functions.R')
meta<-read.csv('~/Documents/HCA/benchmarking/Sra_meta_groups_iPSC_730.csv')
#'  filter out bulk samples
meta.sc<-subset(meta,meta$cell !="population")
label<-paste(meta.sc$cell,meta.sc$lineage,sep='-')
meta.sc['label']<-label
#' load transcripts count
cnt<-read.csv('~/Documents/trinity/Docs/idxstat_T_all_protein_coding.csv',header=T)
#' parse data
cnt.dd<-cnt[,-c(1:2)]
cnt.hh<-colnames(cnt.dd)
mlist<-match(meta.sc$sra,cnt.hh)
cnt.dd<-cnt.dd[,mlist]
# convert greater than 1 value to 1
cnt.dd[cnt.dd>1] <- 1
# matrix
x<-as.matrix(cnt.dd)
hd<- hamming_binary(x)
nn<-nrow(x)
hd.r<-hd/nn
hd.log<-log(1-hd.r,base=2)
hd.tsne<-Rtsne(1-hd.r)
df<-data.frame("X1"=hd.tsne$Y[,1],"X2"=hd.tsne$Y[,2],meta.sc)
labels<-paste(df$cell,df$lineage)
df$labels<-labels
hd.grps <- buildSNNGraph(1-hd.r)
# extract memebership
hd.clusters <- cluster_fast_greedy(hd.grps)
hd.fc <-  membership(hd.clusters)
df$membership<-as.factor(hd.fc)
p.hd<-ggscatter(
  df,
  x = "X1",
  y = "X2",
  color = "membership",
  palette = "ucscgb",
  size = 3,
  alpha = 0.9,
  caption= "Tsne plot of gene matrix, color label represent cell type.",
  ggtheme = theme_minimal()
) + border()+
  xlab('Tsne 1')+
  ylab('Tsne 2')


# count matrix or TPM
tpm<-read.csv('~/Documents/HCA/pipeline_test/matrics/second_batch/HISAT2RSEM_merged_TPM.csv',header=T)
tpm.d<-tpm[,-c(1:2)]
tpm.d[tpm.d>0]<- 1
mlist<-match(meta.sc$sra,colnames(tpm.d))
tpm.d<- tpm.d[,mlist]
tpm.hd<-hamming_binary(as.matrix(tpm.d))
tpm.log<-log(tpm.hd+1,base=2)
tpm.tsne<-Rtsne(tpm.log)
df<-data.frame("X1"=tpm.tsne$Y[,1],"X2"=tpm.tsne$Y[,2],meta.sc)
labels<-paste(df$cell,df$lineage)
df$labels<-labels
p.tpm<-ggscatter(
  df,
  x = "X1",
  y = "X2",
  color = "labels",
  palette = "ucscgb",
  size = 3,
  alpha = 0.9,
  caption= "Tsne plot of gene matrix, color label represent cell type.",
  ggtheme = theme_minimal()
) + border()+
  xlab('Tsne 1')+
  ylab('Tsne 2')
tpm.grps <- buildSNNGraph(log(tpm.hd+1,base=2))
# extract memebership
tpm.clusters <- cluster_fast_greedy(tpm.grps)
tpm.fc <-  membership(tpm.clusters)