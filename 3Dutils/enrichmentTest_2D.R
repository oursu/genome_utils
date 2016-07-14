
args=commandArgs(trailingOnly=TRUE)

f1=args[1]
f2=args[2]
outfile=args[3]

r1name=basename(f1)
r2name=basename(f2)

r1=read.table(f1)
r2=read.table(f2)
r2=data.frame(r2[,1],r2[,2],r2[,1],r2[,3],r2[,4],NA,NA)
colnames(r1)=colnames(r2)=c('chr1','frag1','chr2','frag2','count','pvalue','qvalue')

print('Read files')

#annotate distances and quantiles
#r1=data.frame(r1,dist=round(abs(as.numeric(r1[,2])-as.numeric(r1[,3]))/binSize))
#r2=data.frame(r2,dist=round(abs(as.numeric(r2[,2])-as.numeric(r2[,3]))/binSize))

r1names=paste(r1[,'chr1'],r1[,'frag1'],r1[,'chr2'],r1[,'frag2'])
r2names=paste(r2[,'chr1'],r2[,'frag1'],r2[,'chr2'],r2[,'frag2'])
rownames(r1)=r1names
rownames(r2)=r2names
r1_r2_common=unique(intersect(r1names,r2names))

r1_common=r1[r1_r2_common,]
r2_common=r2[r1_r2_common,]

print('Starting enrichment')

#enrichment function
get_enrichment=function(r1_common,r2_common,r1_col,r1_threshold,r1_larger,r2_col,r2_threshold,r2_larger){
	r1_common=data.frame(r1_common,quantile=rank(as.numeric(r1_common[,'count']),ties.method="random")/(dim(r1_common)[1]))
	r2_common=data.frame(r2_common,quantile=rank(as.numeric(r2_common[,'count']),ties.method="random")/(dim(r2_common)[1]))

	if (r1_larger){
	   r1_subset=which(as.numeric(r1_common[,r1_col])>=r1_threshold)
        }
	if (r1_larger==FALSE){
	   r1_subset=which(as.numeric(r1_common[,r1_col])<=r1_threshold)
	}
	if (r2_larger){
           r2_subset=which(as.numeric(r2_common[,r2_col])>=r2_threshold)
      	}
	if (r2_larger==FALSE){
	   r2_subset=which(as.numeric(r2_common[,r2_col])<=r2_threshold)
        }
	common=intersect(r1_subset,r2_subset)
	n=dim(r1_common)[1]
	m=data.frame(r1_stuff=c(length(common),length(setdiff(r1_subset,r2_subset))),
			not_r1_stuff=c(length(setdiff(r2_subset,r1_subset)),n-length(union(r1_subset,r2_subset))))
	out=fisher.test(m)
	return(out[['estimate']])
}

enrichments=c()
quantiles=seq(from=0,to=1,by=0.05)
for (q in quantiles){
    print(q)
    enrichments=c(enrichments,get_enrichment(r1_common,r2_common,'qvalue',0.05,FALSE,'quantile',q,TRUE))
}


outfilename=paste(outfile,r1name,'.VS.',r2name,sep='')
png(paste(outfilename,'.enrichment.png',sep=''))
par(mar=c(10,10,10,10))
plot(quantiles,enrichments,xlab=paste(r2name,'\nquantile'),ylab=paste('Fold enrichment in\n',r1name,'\ninteractions (0.05 FDR)'))
abline(h=1,col='gray')
dev.off()

png(paste(outfilename,'.scatterCounts_Log2.png',sep=''))
#ready for the scatterplot
#scatter counts
par(mar=c(10,10,10,10))
plot(log(r1_common[,'count'],base=2),log(r2_common[,'count'],base=2),xlab=paste(r1name,'\nlog2(counts)'),ylab=paste(r2name,'\nlog2(counts)'),pch='.')
dev.off()

print('===DONE')

