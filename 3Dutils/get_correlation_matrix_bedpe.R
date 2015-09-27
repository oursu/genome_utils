source('/srv/scratch/oursu/code/genome_utils/3Dutils/R_utils.R') #### change this
args=commandArgs(trailingOnly=TRUE)

mfile=args[1]
out=args[2]
#module load r/3.2.0

compute_corr=function(m,method_name,out){
	BED_ENTRIES=c(1:3)
	rownames(m)=c(1:(dim(m)[1]))
	bed=m[,BED_ENTRIES]
	m=m[,-BED_ENTRIES]
	chromosomes=as.character(unique(bed[,1]))
	for (chromo in chromosomes){
	    
	    #select chromosome items
	    print(chromo)
	    chr_items=which(as.character(bed[,1])==chromo)
	    chr_bed=bed[chr_items,]
	    chr_m=m[chr_items,]

	    #compute correlation
	    correlations=cor(t(chr_m),method=method_name)
	    
	    #reshape result
	    correlations_melted=matrix_wide_to_long_upperTri(correlations)
	    e1=chr_bed[as.character(correlations_melted[,'X1']),]
	    e2=chr_bed[as.character(correlations_melted[,'X2']),]

	    #convert to bedpe
	    cm_bedpe=data.frame(e1[,1],e1[,2],e1[,3],e2[,1],e2[,2],e2[,3],correlations_melted[,'value'])
	    
	    #write table
	    outfile=paste(out,'.',chromo,sep='')
	    write.table(cm_bedpe,file=outfile,sep='\t',row.names=F,col.names=F,quote=F)
	    gzip_file(outfile)
	}
}

print('Starting correlations')
m=read.table(mfile,header=TRUE)
print('Computing Pearson')
compute_corr(m,'pearson',paste(out,'.PearsonCorr',sep=''))
print('Computing Spearman')
compute_corr(m,'spearman',paste(out,'.SpearmanCorr',sep=''))
print('DONE computing correlations')
