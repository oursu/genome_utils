
args=commandArgs(trailingOnly=TRUE)

tf.f=args[1]
qtlpeak.f=args[2]
out=args[3]

testing=function(){
tf.f="/srv/gsfs0/projects/kundaje/users/oursu/crispr/CRISPR_screens_2016_08_15/results/tfs/TSS/TSS.BinnedDataTo1000.relativeToplasmid.atLeast5guides.E.overlapMatrix.gz"  
                                           
qtlpeak.f="/srv/gsfs0/projects/kundaje/users/oursu/crispr/CRISPR_screens_2016_08_15/results/processing/hit_files/TSS/TSS.BinnedDataTo1000.relativeToplasmid.atLeast5guides.CRKO.IDReffect.830.sig_elements.negHits.bed.gz"
out="/srv/gsfs0/projects/kundaje/users/oursu/crispr/CRISPR_screens_2016_08_15/results/tfs/TSS/TSS.BinnedDataTo1000.relativeToplasmid.atLeast5guides.E.IN.CRKO.IDReffect.830.sig_elements.negHits.enrichment.xls" 
}

compute_TF_enrichment=function(tf.f,qtlpeak.f,out){
tfdata=read.table(tf.f,header=TRUE)
tfdata[is.na(tfdata)]=0
tfs=setdiff(colnames(tfdata),c('chr','start','end','name','strand','score'))
qtlpeak=read.table(qtlpeak.f)
rownames(tfdata)=as.character(tfdata$name)
qtlpeaks=as.character(qtlpeak[,4])
#paste(as.character(qtlpeak[,1]),':',as.character(qtlpeak[,2]),'-',as.character(qtlpeak[,3]),sep='')
#restrict to only peaks with at least 1 TFBS
#peaksWithTFBS=which(rowSums(tfdata[,tfs])>0)
#tfdata=tfdata[peaksWithTFBS,]
#qtlpeaks=as.character(qtlpeak[which(as.character(qtlpeak[,1]) %in% as.character(tfdata$peak)),1])
#tfdata=data.frame(tfdata,QTLpeak=FALSE)
#tfdata[qtlpeaks,'QTLpeak']=TRUE

print(head(qtlpeaks))

n=dim(tfdata)[1]
result=data.frame(tf=tfs,
		A=NA,
		Overlap=NA,
		B=length(qtlpeaks),
		Total=n,
		enrichment=NA,
		enrichment_p=NA,
		confLow=NA,confHigh=NA)
rownames(result)=tfs
for (tf in tfs){
    tf_peaks=as.character(tfdata[which(as.numeric(as.character(tfdata[,tf]))>0),'name'])
    qtl_peaks_with_tf=intersect(tf_peaks,qtlpeaks)
    m=data.frame(TF=c(length(qtl_peaks_with_tf),length(tf_peaks)-length(qtl_peaks_with_tf)),
		notTF=c(length(qtlpeaks)-length(qtl_peaks_with_tf),n-length(tf_peaks)-length(qtlpeaks)+length(qtl_peaks_with_tf)))
    ftest=fisher.test(m)
    result[tf,'enrichment']=ftest$estimate
    result[tf,'enrichment_p']=ftest$p.value
    result[tf,'A']=length(tf_peaks)
    result[tf,'Overlap']=length(qtl_peaks_with_tf)
    result[tf,'confLow']=ftest$conf.int[1]
    result[tf,'confHigh']=ftest$conf.int[2]
}
result=result[order(result$enrichment_p),]
result=data.frame(result,BH=p.adjust(result$enrichment_p))
rownames(result)=NULL
write.table(result,
	file=out,
	sep='\t',quote=F,row.names=F,col.names=T)
outstuff=list()
outstuff[['result']]=result
outstuff[['N']]=n
return(outstuff)
}

### COMPUTE ENRICHMENTS
#######################
print('Computing enrichment')
outstuff=compute_TF_enrichment(tf.f,qtlpeak.f,out)
result=outstuff[['result']]
N=outstuff[['N']]
#replace complicated TF names by TF HGNCs
#hgnc.f='/home/oursu/ENCODE.hg19.TFBS.QC.metadata.jun2012 - TFs_SPP_pooled.tsv'
#hgnc=read.csv(hgnc.f,sep='\t')
#rownames(hgnc)=hgnc$FILENAME
result2=data.frame(result,TFname=gsub('idrOptimalBlackListFilt.','',gsub('.[a-zA-Z0-9]*.peaks','',result$tf)))
result2=data.frame(result2,name_percent=paste(result2$TFname,' | OR ',signif(result2$enrichment,2),' | ',
as.integer(100*result2$Overlap/result2$B),'% hits explained',sep=''))
sig_tfs=result2[which(as.numeric(as.character(result2$BH))<=0.05),'name_percent']

### PLOT ENRICHMENT
###################
print('Plotting')
require(ggplot2)
pdf(paste(out,'.pdf',sep=''),height=20,width=15)
res.cur=result2
res.cur$TFname=factor(res.cur$TFname,levels=res.cur$TFname[order(res.cur$enrichment)])
res.cur=cbind(res.cur,sig_after_BHcorrection=(as.numeric(as.character(res.cur$BH))<=0.05))
#print(ggplot(res.cur, aes(y=enrichment, x=TFname,col=sig_after_BHcorrection)) + coord_flip()+geom_point()+geom_errorbar(aes(ymax = res.cur$confLow, 
#ymin=res.cur$confHigh))+ylim(0,4)+
#ggtitle(basename(out)))
#print(res.cur)
print(ggplot(res.cur, aes(y=enrichment, x=TFname,col=sig_after_BHcorrection)) + coord_flip()+geom_point()+geom_errorbar(aes(ymax = res.cur$confLow, 
ymin=res.cur$confHigh))+ylim(0,20)+
ggtitle(basename(out)))

overlap_data=data.frame(Count=result2[,'Overlap'],TFname=result2[,'name_percent'],category='Hits in the set of interest')
b_minus_overlap_data=data.frame(Count=result2[,'B']-result2[,'Overlap'],TFname=result2$name_percent,category='Hits not in the set of interest')
a_minus_overlap_data=data.frame(Count=result2[,'A']-result2[,'Overlap'],TFname=result2$name_percent,category='Non-hits in the set of interest')
combined=rbind(overlap_data,b_minus_overlap_data,a_minus_overlap_data)
combined[,'TFname']=factor(combined[,'TFname'],levels=overlap_data[order(overlap_data[,'Count']),'TFname'])
sig_v=rep('Sig enrichment',times=dim(combined)[1])
sig_v[which(as.character(combined$TFname) %in% as.character(sig_tfs))]='Sig enrichment'
combined[,'category']=paste(sig_v,' | ',combined[,'category'],sep='')
#combined=data.frame(combined,sig=sig_v)
print(ggplot(combined, aes(y=Count, x=TFname,fill=category)) + coord_flip()+geom_bar(stat='identity')+ylim(0,N)+geom_hline(yintercept = result2[1,'B'])+
theme_bw()+scale_fill_manual(values=c('darkblue','lightblue','gray','darkred','pink','gray'))+
ggtitle(basename(out)))
print(result2)
dev.off()



