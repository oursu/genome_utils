from optparse import OptionParser
import os
import math
from time import gmtime, strftime
import re
import numpy as np
'''
By Oana Ursu
oursu@stanford.edu
'''

def main():
    parser=OptionParser()
    parser.add_option('--bashrc',dest='bashrc',default='/srv/scratch/oursu/code/genome_utils/data_analysis/ChIA-PET/mango_bashrc')
    parser.add_option('--out_dir',dest='out_dir',default='/srv/scratch/oursu/3Dgenome/data/ChIA-PET/')
    parser.add_option('--metadata',dest='metadata',help='Columns are sra link, sample name, replicate group', default='/srv/scratch/oursu/3Dgenome/data/ChIA-PET/metadata/metadata_GSE59395')
    parser.add_option('--step',dest='step',help='Step: download, fastq, mango', default='mango')
    parser.add_option('--by_replicateGroup',dest='by_replGroup',action='store_true')
    parser.add_option('--blacklist',dest='blacklist',default='/srv/scratch/oursu/data/wgEncodeDacMapabilityConsensusExcludable.bed')
    opts,args=parser.parse_args()

    fastqdir=opts.out_dir+'/fastq/'
    os.system('mkdir '+fastqdir)
    mangodir=opts.out_dir+'/mango/'
    os.system('mkdir '+mangodir)

    d={}
    for line in open(opts.metadata,'r'):
        if line[0]=='#': #this is the header
            continue
        items=line.strip().split('\t')
        sra=items[0]
        samplename=items[1]
        replGroup=items[2]
        if replGroup not in d.keys():
            d[replGroup]={}
        d[replGroup][samplename]={}
        d[replGroup][samplename]['sra']=sra

    #===========================
    # SINGLE REPLICATES
    #===========================
    if not opts.by_replGroup:
        for line in open(opts.metadata,'r'):
            items=line.strip().split('\t')
            sra=items[0]
            samplename=items[1]
            replGroup=items[2]

            cmds=[]
            cmds.append('source '+opts.bashrc)

            #========================
            # Data download
            #========================
            if opts.step=='download':
                get_data(sra,fastqdir,'download')

            #========================
            # Conv to fastq
            #========================
            if opts.step=='fastq':
                get_data(sra,fastqdir,'fastq')

            #========================
            # Run mango on the data
            #========================
            if opts.step=='mango':
                run_mango(fastqdir+'/'+samplename+"_1.fastq.gz",fastqdir+'/'+samplename+"_2.fastq.gz",samplename,mangodir,cmds,opts.blacklist)

    #===============================
    # MERGED REPLICATES
    #===============================
    if opts.by_replGroup:
        for replGroup in d.keys():
            #if replGroup=='ChIAPET_GM12878_RAD21':
            #    continue
            cmds=[]
            cmds.append('source '+opts.bashrc)
            fq1s=''
            fq2s=''
            for samplename in d[replGroup].keys():
                fq1s=fq1s+' '+fastqdir+'/'+samplename+"_1.fastq.gz"
                fq2s=fq2s+' '+fastqdir+'/'+samplename+"_2.fastq.gz"

            #========================
            # Run mango on the data
            #========================
            if opts.step=='mango':
                run_mango(fq1s,fq2s,replGroup,mangodir,cmds,opts.blacklist)  
                


def run_mango(fq1,fq2,samplename,outdir,cmds,blacklist):
    outmango=outdir+'/'+samplename+'/'
    tmp=outmango+'/tmp'
    os.system("mkdir -p "+outmango)
    os.system("mkdir -p "+tmp)
    new_fq1=tmp+'/'+samplename+"_1.combined.fastq"
    new_fq2=tmp+'/'+samplename+"_2.combined.fastq"
    cmds.append("zcat -f "+fq1+" > "+new_fq1)
    cmds.append("zcat -f "+fq2+" > "+new_fq2)
    cmds.append("Rscript /software/mango/mango.R --fastq1 "+new_fq1+" --fastq2 "+new_fq2+" --outdir "+outmango+" --prefix "+samplename+".mango --argsfile ${ARGFILE} --chromexclude chrX,chrM,chrY --stages 1:5 --reportallpairs TRUE")
    #cmds.append('rm -r '+tmp)
    #cmds.append('rm '+outmango+'/*sam')
    #cmds.append('rm '+outmango+'/*fastq')
    #cmds.append('rm '+outmango+'/*bedpe')
    qsub_a_command(cmds,outmango+'/'+samplename+'.mango.script.sh','20G')

def get_data(sra,fastqdir,cmdtype):
    os.system('${SRASCRIPT} '+sra+' '+fastqdir+' '+cmdtype)

def qsub_a_command(cmds,shell_script_name,memory_number='20G'):
    #write a shell script (for reproducibility)
    f=open(shell_script_name,'w')
    print cmds
    for i in range(len(cmds)):
        f.write(cmds[i]+'\n')
    f.close()
    #make runnable
    os.system('chmod 711 '+shell_script_name)
    #Qsub the script
    os.system("qsub -l mem_free="+memory_number+" -l h_vmem="+memory_number+" -l h_rt=40:00:00 -o "+shell_script_name+'.o'+' -e '+shell_script_name+'.e'+' '+shell_script_name)



main()