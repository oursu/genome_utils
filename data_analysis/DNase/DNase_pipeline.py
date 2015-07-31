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
    parser.add_option('--smooth',dest='smooth',default='150')
    parser.add_option('--binSize',dest='binSize',default='1')
    parser.add_option('--normFactor',dest='normFactor',default='2150570000')
    parser.add_option('--chrSizes',dest='chrSizes',default='/mnt/lab_data/kundaje/projects/encodeenhancerpredict/mm10/mm10.chrom.sizes')
    parser.add_option('--genome',dest='genome',default='mm')
    parser.add_option('--bashrc',dest='bashrc',default='/srv/scratch/oursu/code/data_analysis/DNase/DNase_bashrc')
    parser.add_option('--out_dir',dest='out_dir',default='/mnt/lab_data/kundaje/projects/encodeenhancerpredict/DNaseI/results/')
    parser.add_option('--metadata',dest='metadata',help='Columns are sample name, replicate group, bam file', default='/mnt/lab_data/kundaje/projects/encodeenhancerpredict/histone_chipseq/doc/histone_chipseq_metadata.txt')
    parser.add_option('--step',dest='step',help='Step: tagAlign', default='tagAlign')
    parser.add_option('--by_replicateGroup',dest='by_replGroup',action='store_true')
    opts,args=parser.parse_args()

    
    os.system('mkdir '+opts.out_dir)
    tagAligndir=opts.out_dir+'/tagAlign/'
    os.system('mkdir '+tagAligndir)
    qcdir=opts.out_dir+'/QC'
    os.system('mkdir '+qcdir)
    peakdir=opts.out_dir+'/peaks/'
    os.system('mkdir '+peakdir)
    bigwigdir=opts.out_dir+'/bigwig/'
    os.system('mkdir '+bigwigdir)
    tagAligndirCenter=opts.out_dir+'/tagAlignCenteredHalfFragLen/'
    os.system('mkdir '+tagAligndirCenter)
    mergedRepldir=opts.out_dir+'/merged_replicates/'
    os.system('mkdir '+mergedRepldir)

    #map the data to the replicate groups so we can merge reads and call peaks together                                                                              
    sample2group={}
    for line in open(opts.metadata,'r').readlines():
        items=line.strip().split('\t')
        sampleid=items[0]
        replicateGroup=items[1]
        bam=items[2]
        if replicateGroup not in sample2group.keys():
            sample2group[replicateGroup]={}
            sample2group[replicateGroup]['bam']=set()
            sample2group[replicateGroup]['tagAlign']=set()
            sample2group[replicateGroup]['tagAlignCentered']=set()
        sample2group[replicateGroup]['bam'].add(bam)
        sample2group[replicateGroup]['tagAlign'].add(tagAligndir+'/'+sampleid+'.tagAlign.gz')
        sample2group[replicateGroup]['tagAlignCentered'].add(tagAligndirCenter+'/'+sampleid+'.centerByHalfFragLen.tagAlign.gz')

    #========================================================
    # SINGLE SAMPLES
    #========================================================
    if not opts.by_replGroup:
        for line in open(opts.metadata,'r').readlines():
            items=line.strip().split('\t')
            sampleid=items[0]
            replicateGroup=items[1]
            bam=items[2]
            tagAlign=tagAligndir+'/'+sampleid+'.tagAlign.gz'
            cmds=[]
            cmds.append('source '+opts.bashrc)

            #===================
            # make tagAlign file
            #===================
            if opts.step=='tagAlign':
                cmds.append('remove_chrM_qc30_2tagAlign.sh --inbam '+bam+' --outdir '+tagAligndir+' --sample '+sampleid+' --SE 1')
                qsub_a_command('qqqq'.join(cmds),tagAligndir+'/'+sampleid+'_2tagAlign.script.sh','qqqq','20G')

            #=============================
            # center tagAlign by fragLen/2 - this is for creating the bigwigs
            #=============================
            if opts.step=='centerByHalfFragLen':
                fragLenOver2=str(int(int(get_best_fragLen(qcdir+'/'+sampleid+'SPP_table.txt'))/2))
                cmds.append('slopBed -i '+tagAlign+' -g '+opts.chrSizes+' -l -'+fragLenOver2+' -r '+fragLenOver2+' -s | gzip > '+tagAligndirCenter+'/'+sampleid+'.centerByHalfFragLen.tagAlign.gz')
                qsub_a_command('qqqq'.join(cmds),tagAligndirCenter+'/'+sampleid+'_scenterByHalfFragLen.script.sh','qqqq','20G')

            #===================
            # run SPP
            #===================
            if opts.step=='SPP':
                sppcmd='Rscript $(which run_spp_nodups.R) -c='+tagAlign+' '+'-odir='+qcdir+' '+'-savp -out='+qcdir+'/'+sampleid+'SPP_table.txt'
                cmds.append(sppcmd)
                qsub_a_command('qqqq'.join(cmds),qcdir+'/'+sampleid+'_SPP.script.sh','qqqq','20G')
            
            #===================
            # MACS2
            #===================
            if opts.step=='MACS2':
                fragLenOver2=str(int(int(get_best_fragLen(qcdir+'/'+sampleid+'SPP_table.txt'))/2))
                peakfile=peakdir+'/'+sampleid

                if fragLen==0: #single-cut
                    #shift reads first
                    cmds.append('slopBed -i '+tagAlign+' -g '+opts.chrSizes+' -l 75 -r -75 -s > '+tagAlign+'.shift75bp.bed')
                    macscmd='macs2 callpeak -t '+tagAlign+'.shift75bp.bed'+' -f BED -n '+peakfile+' -g '+opts.genome+ ' -p 1e-2 --nomodel --shift 75 -B --SPMR --keep-dup all --call-summits'

                if fragLen!=0: #double cut, treat like chipseq
                    macscmd='macs2 callpeak -t '+tagAlign+' -f BED -n '+peakfile+' -g '+opts.genome+ ' -p 1e-2 --nomodel --shift '+fragLenOver2+' -B --SPMR --keep-dup all --call-summits'

                cmds.append(macscmd)
                qsub_a_command('qqqq'.join(cmds),peakdir+'/'+samplename+'_MACS.script.sh','qqqq','20G')
        

    
    #========================================================
    # MERGED REPLICATES
    #========================================================
    if opts.by_replGroup: #call peaks and bigwigs by replicate group
        for replGroup in sample2group.keys():
            cmds=[]
            cmds.append('source '+opts.bashrc)
            tagAligns=' '.join(list(sample2group[replGroup]['tagAlign']))
            centeredTagAligns=' '.join(list(sample2group[replGroup]['tagAlignCentered']))

            #figure out the fraglen
            #get median frag Len
            fragLenValues=[]
            for f in sample2group[replGroup]['tagAlign']:
                fragLenValues.append(get_best_fragLen(qcdir+'/'+re.sub('.tagAlign.gz','',os.path.basename(f))+'SPP_table.txt'))
            print fragLenValues
            print int(np.median(fragLenValues))

            if opts.step=='MACS2':
                combined_tagAlign=mergedRepldir+'/'+replGroup+'.merged.tagAlign.gz'
                cmds.append('zcat '+tagAligns+' | gzip > '+combined_tagAlign)
                cmds.append(macs2(combined_tagAlign,str(int(np.median(fragLenValues))/2),opts.genome,mergedRepldir+'/'+replGroup+'.peaks.narrowPeak'))
                qsub_a_command('qqqq'.join(cmds),mergedRepldir+'/'+replGroup+'_MACS2.script.sh','qqqq','20G')

            if opts.step=='bigwig':
                #these will be called on the centered data
                bigwig_from_tagAlign(replGroup,centeredTagAligns,opts.normFactor,opts.binSize,bigwigdir,opts.smooth)            
    


def bigwig_from_tagAlign(replGroup,tagAligns,normFactor,binSize,bigwigdir,smooth):
    suff=replGroup
    mtagAlign=bigwigdir+'/tmp/'+suff+'.tagAlign.gz'
    os.system('mkdir -p '+bigwigdir+'/tmp')
    cmds.append('cat '+tagAligns+' | gzip >  '+mtagAlign)
    cmds.append('bedToBam -i '+mtagAlign+' -g '+opts.chrSizes+' > '+mtagAlign+'.bam')
    cmds.append('samtools sort '+mtagAlign+'.bam '+mtagAlign+'.sorted')
    cmds.append('samtools index '+mtagAlign+'.sorted.bam')
    cmds.append('module load python_anaconda/2.2.0')
    cmds.append('bamCoverage --normalizeTo1x '+normFactor+' --smoothLength '+smooth+' --binSize '+binSize+' -b '+mtagAlign+'.sorted.bam'+' -o '+bigwigdir+'/'+suff+'.bw')
    print cmds
    qsub_a_command('qqqq'.join(cmds),bigwigdir+'/'+replGroup+'_bigwig.script.sh','qqqq','20G')

def get_best_fragLen(f):
    items=open(f,'r').readlines()[0].strip().split('\t')
    fs=items[2].split(',')
    ccs=items[3].split(',')
    ccs2=[]
    for item in ccs:
        ccs2.append(float(item))
    return int(fs[np.argmax(ccs2)])

def macs2(bed,shift,chrSizes,out):
    return 'macs2 callpeak -t '+bed+' -f BED -n '+out+' -g '+chrSizes+ ' -p 1e-2 --nomodel --shift '+shift+' -B --SPMR --keep-dup all --call-summits'

def qsub_a_command(cmd,shell_script_name,split_string=',',memory_number='20G'):
    #write a shell script (for reproducibility)
    f=open(shell_script_name,'w')
    #f.write("apo=\"'\"\n")
    cmds=cmd.split(split_string)
    print cmds
    for i in range(len(cmds)):
        f.write("cmd"+str(i)+"='"+cmds[i]+"'"+'\n')
        f.write('echo $cmd'+str(i)+'\n')
        f.write('eval $cmd'+str(i)+'\n')
    f.close()
    #make runnable
    os.system('chmod 711 '+shell_script_name)
    #Qsub the script
    os.system("qsub -l mem_free="+memory_number+" -l h_vmem="+memory_number+" -l h_rt=40:00:00 -o "+shell_script_name+'.o'+' -e '+shell_script_name+'.e'+' '+shell_script_name)



main()
