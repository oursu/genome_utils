
args=commandArgs(trailingOnly=TRUE)
gc_file=args[1] #gc_file="/srv/scratch/oursu/test/testgraph_code_orca.out"
out=args[2] #out="/srv/scratch/oursu/test/testgraph_code_orca.out.Rout"

require(pheatmap)

graphlets=read.table(gc_file)
colnames(graphlets)=seq(from=0,to=dim(graphlets)[2]-1,by=1)
DEGREE=1


pdf(paste(out,'.pdf',sep=''))
#first, plot the degrees of the nodes
plot(ecdf(graphlets[,DEGREE]),ylab='ECDF',xlab='Node degree',main=paste('Node degree \n',basename(gc_file),sep=''))
#now, plot a heatmap of the graphlets
pheatmap(as.matrix(graphlets),cluster_cols=F,cluster_rows=F,main="Graphlet counts")
pheatmap(log(as.matrix(graphlets+1),base=10),cluster_cols=F,cluster_rows=F,main="Log10 (graphlet counts + 1) heatmap")
dev.off()
