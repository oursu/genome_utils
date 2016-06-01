require("reshape2",lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2/") ####### change this library


matrix_wide_to_long_upperTri=function(m){
	upperTriangle=upper.tri(m, diag=F)
        m.upperTriangle=m
        m.upperTriangle[!upperTriangle]=NA
        m_melted=na.omit(melt(m.upperTriangle, value.name ="value"))
        colnames(m_melted)<-c("X1", "X2", "value")
	return(m_melted)
}

gzip_file=function(f){
	system(paste('zcat -f ',f,' | gzip > ',f,'.gz',sep=''))
        system(paste('rm ',f,sep=''))
}