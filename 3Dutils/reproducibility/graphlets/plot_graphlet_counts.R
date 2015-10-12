
args=commandArgs(trailingOnly=TRUE)
gc_file=args[1] #gc_file="/srv/scratch/oursu/test/testgraph_code_orca.out"
out=args[2] #out="/srv/scratch/oursu/test/testgraph_code_orca.out.Rout"
annos=strsplit(args[3],',')[[1]] #should be files with in the same order as the nodes, stuff, one stuff per column, 0/1. For instance, to specify chromatin state, the file would have #nodes lines, and one column for each possible chromatin state, and 1/0 for each node, depending on whether it is in that chromatin state or not.

require(pheatmap)

graphlets=read.table(gc_file)
colnames(graphlets)=seq(from=0,to=dim(graphlets)[2]-1,by=1)
DEGREE=1

#If we have annotations, then plot a profile by	annotation
anno_profiles=list()
for (anno in annos){
 print(anno)
 anno_data=read.table(anno,header=TRUE)[,-c(1:3)]
 anno_profiles[[anno]]=data.frame(array(0,dim=c(ncol(anno_data),ncol(graphlets))))
 rownames(anno_profiles[[anno]])=colnames(anno_data)
 for (anno_type	in rownames(anno_profiles[[anno]])){
  print(anno_type)
  nodes=which(anno_data[,anno_type]==1)
  anno_profiles[[anno]][anno_type,]=colMeans(graphlets[nodes,])
 }
}


pdf(paste(out,'.pdf',sep=''))
#first, plot the degrees of the nodes
plot(ecdf(graphlets[,DEGREE]),ylab='ECDF',xlab='Node degree',main=paste('Node degree \n',basename(gc_file),sep=''))
#now, plot a heatmap of the graphlets
pheatmap(as.matrix(graphlets),cluster_cols=F,cluster_rows=F,main="Graphlet counts")
pheatmap(log(as.matrix(graphlets+1),base=10),cluster_cols=F,cluster_rows=F,main="Log10 (graphlet counts + 1) heatmap")
for (anno in annos){
    pheatmap(as.matrix(anno_profiles[[anno]]),cluster_cols=F,cluster_rows=F,main=paste("Mean graphlet counts\n",basename(anno),sep=''))
    pheatmap(as.matrix(log(anno_profiles[[anno]]+1,base=10)),cluster_cols=F,cluster_rows=F,main=paste("Mean log10(graphlet counts +1) counts\n",
    basename(anno),sep=''))
}
dev.off()
