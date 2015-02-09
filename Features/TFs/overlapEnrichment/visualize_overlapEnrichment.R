read_enrichment_file=function(infile,SIG_THRESH,fillin){
  data=read.table(infile,header=TRUE)
  data.enrich=data.frame(enrichment=data[,'enrichment'])
  rownames(data.enrich)=data$tf
  to_remove=which(as.numeric(as.character(data$BH))>SIG_THRESH)
  #keep only the significant ones
  if (length(to_remove)>0){
    data.enrich[to_remove,'enrichment']=fillin
  }
  return(data.enrich)
}

optimal_ordering=function(m,meth){
  require(cba)
  d <- dist(as.matrix(m),method=meth)
  hc <- hclust(d)
  co <- order.optimal(d, hc$merge)
  m.optimalRows=as.matrix(m)[co$order,]
  return(m.optimalRows)
}

heatmap_enrichments=function(data,out,meth){
  data.optimal=optimal_ordering(data,meth)
  require(pheatmap)
  pdf(out)
  pheatmap(as.matrix(t(data.optimal)),cluster_rows=FALSE,fontsize=5,breaks=seq(from=-0.1,to=3,by=0.1),
        color=colorRampPalette(c("blue", "white", "red",'black'))(n = 31),cellwidth=5,cellheight=5)
  #pheatmap(as.matrix(t(data.optimal)),cluster_rows=FALSE,fontsize=5,
  #  cellwidth=5,cellheight=5,breaks=seq(from=0,to=1,by=0.01),
  #     color=colorRampPalette(c("white", "yellow","green","black"))(n = 100))
 
  dev.off()
}

overlapEnrichment_distalQTL=function(){
  enrichfiles='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/OverlapEnrichment/ENRICHDIR/ENRICHPREF'
  enrichments=c('TFBS_overlap_','Motif_overlap_','MotifCorrelatedLocal_overlap_') #add in disrupted motif overlaps, and hQTL overlaps

  hmarks=c('H3K4ME1','H3K4ME3','H3K27AC')
  for (enrich in enrichments){
    first=TRUE
    for (suffix in c('HMARK.QTLpeaks','LocalPeakIsHMARK.QTLpeaks_affectingDistalPeaks')){
      for (hmark in hmarks){
        f=gsub('ENRICHDIR',paste(enrich,'QTLpeaks0kb',sep=''),
                            gsub('ENRICHPREF',paste(enrich,hmark,'QTLpeaks0kb___.overlapEnrichIN',gsub('HMARK',hmark,suffix),sep=''),enrichfiles))
        cur_data=read_enrichment_file(f,0.05,-1)
        #cur_data=cur_data/max(cur_data[,1])
        if (suffix=='HMARK.QTLpeaks'){
          addon='Local'
        }
        if (suffix=='LocalPeakIsHMARK.QTLpeaks_affectingDistalPeaks'){
          addon='Distal'
        }
        rownames(cur_data)=gsub('bed.OverlapChIPseq','',
                            gsub('MotifMatch_','',
                              gsub('MergedPeaks_ChIPseq_','',
                                gsub('correlatedMotifs.motif.pouya.Motif.','',
                                  gsub('scanThresh0','',rownames(cur_data))))))
        if (first==FALSE){
          data=cbind(data,cur_data[rownames(data),])
          colnames(data)[ncol(data)]=paste(gsub('_',' ',enrich),hmark,' ',addon,sep='')
        }
        if (first==TRUE){
         data=cur_data
          first=FALSE
          colnames(data)[1]=paste(gsub('_',' ',enrich),hmark,' ',addon,sep='')
        }
      }
    }
    heatmap_enrichments(data,paste(dirname(f),'overlapEnrichmentHeatmap.pdf',sep=''),'euclidean')
  }
}

overlapEnrichment_distalQTL()

