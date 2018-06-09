
args=commandArgs(trailingOnly=TRUE)

aggregated_values_f=args[1]
sig_f=args[2]
out=args[3]
mini=as.numeric(as.character(args[4]))
maxi=as.numeric(as.character(args[5]))

print(mini)
print(maxi)

d=read.table(aggregated_values_f)
print(head(d))
colnames(d)=c('chr','start','end','name','value')
rownames(d)=as.character(d$name)
print(head(d))

sig=read.table(sig_f)
print(head(sig))
colnames(sig)=c('chr','start','end','name')

common_sig=intersect(as.character(sig$name),as.character(d$name))
nonsig=setdiff(as.character(d$name),as.character(sig$name))

all_values=as.numeric(as.character(d[,'value']))
sig_values=as.numeric(as.character(d[common_sig,'value']))
nonsig_values=as.numeric(as.character(d[nonsig,'value']))

print('sigs')
print(summary(sig_values))
print('nonsigs')
print(summary(nonsig_values))
kstest_sig_vs_all=ks.test(all_values,sig_values)
kstest_sig_vs_nonsig=ks.test(nonsig_values,sig_values)

pdf(paste(out,'.pdf',sep=''))
plot(ecdf(all_values),xlim=c(mini,maxi),ylim=c(0,1),main='',col='black',xlab='',ylab='')
par(new=TRUE)
plot(ecdf(sig_values),xlim=c(mini,maxi),ylim=c(0,1),main=paste('KS test (sig vs all):','\n','D=',kstest_sig_vs_all[['statistic']],'\n','p=',kstest_sig_vs_all[['p.value']],sep=''),col='green',xlab=paste('Distance (bp)','\n',basename(out),sep=''),ylab='ECDF')

plot(ecdf(nonsig_values),xlim=c(mini,maxi),ylim=c(0,1),main='',col='gray',xlab='',ylab='')
par(new=TRUE)
plot(ecdf(sig_values),xlim=c(mini,maxi),ylim=c(0,1),main=paste('KS test (sig vs nonsig):','\n','D=',kstest_sig_vs_nonsig[['statistic']],'\n','p=',kstest_sig_vs_nonsig[['p.value']],sep=''),col='green',xlab=paste('Distance (bp)','\n',basename(out),sep=''),ylab='ECDF')

require(ggplot2)
colors=c('lightblue','gray')
names(colors)=c('Hits','Non-hits')
dataset=data.frame(category=c(rep('Hits',times=length(sig_values)),rep('Non-hits',times=length(nonsig_values))),distance=c(sig_values,nonsig_values))
dataset=data.frame(dataset,distance_kb=as.numeric(as.character(dataset[,'distance']))/1000)
print(summary(dataset[which(as.character(dataset[,'category'])=='Hits'),]))
print(summary(dataset[which(as.character(dataset[,'category'])=='Non-hits'),]))

print(ggplot(dataset,aes(x=category,y=log(distance,10),fill=category)) +geom_boxplot())
print(ggplot(dataset,aes(x=category,y=distance_kb,fill=category)) +geom_boxplot(width=0.5)+theme(panel.border = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme_bw()+theme(axis.text=element_text(size=20),axis.title=element_text(size=20))+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("")+ylab("Distance to closest growth gene (kb)")+scale_fill_manual(values=colors))

#print(ggplot(dataset,aes(x=distance,fill=category)) +geom_histogram(alpha=0.2)+xlim(0,10000))
dev.off()

print(paste(out,'.pdf',sep=''))

stats=data.frame(test=c('sig_vs_all','sig_vs_nonsig'),size1=c(length(sig_values),length(sig_values)),size2=c(length(all_values),length(nonsig_values)),KStest.D=c(kstest_sig_vs_all[['statistic']],kstest_sig_vs_nonsig[['statistic']]),KStest.p=c(kstest_sig_vs_all[['p.value']],kstest_sig_vs_nonsig[['p.value']]))

write.table(stats,file=paste(out,'.txt',sep=''),row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')