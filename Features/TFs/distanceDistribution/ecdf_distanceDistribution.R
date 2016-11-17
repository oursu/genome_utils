
args=commandArgs(trailingOnly=TRUE)

aggregated_values_f=args[1]
sig_f=args[2]
out=args[3]
mini=as.numeric(as.character(args[4]))
maxi=as.numeric(as.character(args[5]))

print(mini)
print(maxi)

d=read.table(aggregated_values_f)
colnames(d)=c('chr','start','end','name','value')
rownames(d)=as.character(d$name)

sig=read.table(sig_f)
colnames(sig)=c('chr','start','end','name')

common_sig=intersect(as.character(sig$name),as.character(d$name))
nonsig=setdiff(as.character(d$name),as.character(sig$name))

all_values=as.numeric(as.character(d[,'value']))
sig_values=as.numeric(as.character(d[common_sig,'value']))
nonsig_values=as.numeric(as.character(d[nonsig,'value']))

kstest_sig_vs_all=ks.test(all_values,sig_values)
kstest_sig_vs_nonsig=ks.test(nonsig_values,sig_values)

pdf(paste(out,'.pdf',sep=''))
plot(ecdf(all_values),xlim=c(mini,maxi),ylim=c(0,1),main='',col='black',xlab='',ylab='')
par(new=TRUE)
plot(ecdf(sig_values),xlim=c(mini,maxi),ylim=c(0,1),main=paste('KS test (sig vs all):','\n','D=',kstest_sig_vs_all[['statistic']],'\n','p=',kstest_sig_vs_all[['p.value']],sep=''),col='green',xlab=paste('Distance (bp)','\n',basename(out),sep=''),ylab='ECDF')

plot(ecdf(nonsig_values),xlim=c(mini,maxi),ylim=c(0,1),main='',col='gray',xlab='',ylab='')
par(new=TRUE)
plot(ecdf(sig_values),xlim=c(mini,maxi),ylim=c(0,1),main=paste('KS test (sig vs nonsig):','\n','D=',kstest_sig_vs_nonsig[['statistic']],'\n','p=',kstest_sig_vs_nonsig[['p.value']],sep=''),col='green',xlab=paste('Distance (bp)','\n',basename(out),sep=''),ylab='ECDF')

dev.off()

print(paste(out,'.pdf',sep=''))

stats=data.frame(test=c('sig_vs_all','sig_vs_nonsig'),size1=c(length(sig_values),length(sig_values)),size2=c(length(all_values),length(nonsig_values)),KStest.D=c(kstest_sig_vs_all[['statistic']],kstest_sig_vs_nonsig[['statistic']]),KStest.p=c(kstest_sig_vs_all[['p.value']],kstest_sig_vs_nonsig[['p.value']]))

write.table(stats,file=paste(out,'.txt',sep=''),row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')