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
    parser.add_option('--normFactor',dest='normFactor',default='2150570000')
    parser.add_option('--chrSizes',dest='chrSizes',default='/mnt/lab_data/kundaje/projects/encodeenhancerpredict/mm10/mm10.chrom.sizes')
    parser.add_option('--genome',dest='genome',default='mm')
    parser.add_option('--bashrc',dest='bashrc',default='/srv/scratch/oursu/code/data_analysis/DNase/DNase_bashrc')
    parser.add_option('--out_dir',dest='out_dir',default='/mnt/lab_data/kundaje/projects/encodeenhancerpredict/DNaseI/results/')
    parser.add_option('--metadata',dest='metadata',help='Columns are fastq, bam, deduped bam', default='/mnt/lab_data/kundaje/projects/encodeenhancerpredict/DNaseI/data/mapped_metadata.txt')
    parser.add_option('--step',dest='step',help='Step: tagAlign', default='tagAlign')
    parser.add_option('--sample2group',dest='sample2group',default='/mnt/lab_data/kundaje/projects/encodeenhancerpredict/DNaseI/data/DNaseFastqMetadata_anshulDownloaded.namesMapped.replicateGroups.aligned')
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

    #map the data to the replicate groups so we can merge reads and call peaks together                                                                              
    sample2group={}
    for line in open(opts.sample2group,'r').readlines():
        items=line.strip().split('\t')
        sampleid=items[1]
        replicateGroup=items[2]
        #encode_group=items[0]
        #repl=items[2]
        if replicateGroup not in sample2group.keys():
            sample2group[replicateGroup]={}
            sample2group[replicateGroup]['bam']=set()
            sample2group[replicateGroup]['tagAlign']=set()
        sample2group[replicateGroup]['bam'].add('/srv/scratch/leepc12/run/ENCODE_enh_predict/'+sampleid.split('_Rep')[0]+'/out/'+sampleid+'.filt.nodup.srt.bam')
        sample2group[replicateGroup]['tagAlign'].add(tagAligndir+'/'+sampleid+'.filt.nodup.srt.tagAlign.gz')

    for line in open(opts.metadata,'r').readlines():
        items=line.strip().split()
        dedup=items[1]
        cmds=[]
        cmds.append('source '+opts.bashrc)
        samplename=os.path.basename(re.sub('.bam','',dedup))

        if opts.by_replGroup:
            continue

        #===================
        # make tagAlign file
        #===================
        if opts.step=='tagAlign':
            cmds.append('remove_chrM_qc30_2tagAlign.sh --indir '+os.path.dirname(dedup)+' --outdir '+tagAligndir+' --sample '+os.path.basename(re.sub('.bam','',dedup))+' --SE 1')
            qsub_a_command('qqqq'.join(cmds),tagAligndir+'/'+os.path.basename(re.sub('.bam','',dedup))+'_2tagAlign.script.sh','qqqq','20G')

        if opts.step=='shiftByFragLen':
            fragLenOver2=int(int(get_best_fragLen(qcdir+'/'+re.sub('.bam','',os.path.basename(dedup))+'SPP_table.txt'))/2)
            tagAlign=tagAligndir+'/'+os.path.basename(re.sub('.bam','',dedup))+'.tagAlign.gz'
            cmds.append('slopBed -i '+tagAlign+' -g '+opts.chrSizes+' -l '+str(fragLenOver2)+' -r -'+str(fragLenOver2)+' -s > '+tagAlign+'.shiftByHalfFragLen.bed')
            qsub_a_command('qqqq'.join(cmds),tagAligndir+'/'+os.path.basename(re.sub('.bam','',dedup))+'_2shiftedtagAlign.script.sh','qqqq','20G')
        
        #===================
        # make signal track
        #===================
        if opts.step=='bigwig':
            cpbam=opts.out_dir+'/tmp/'+samplename+'.bam'
            os.system('mkdir '+opts.out_dir+'/tmp')
            cmds.append('cp '+dedup+' '+cpbam)
            cmds.append('samtools index '+cpbam)
            cmds.append('module load python_anaconda/2.2.0')
            cmds.append('bamCoverage -b '+cpbam+' -o '+bigwigdir+'/'+samplename+'.bw')
            qsub_a_command('qqqq'.join(cmds),bigwigdir+'/'+samplename+'_bigwig.script.sh','qqqq','20G')

        #===================
        # run SPP
        #===================
        if opts.step=='SPP':
            sppcmd='Rscript $(which run_spp_nodups.R) -c='+dedup+' '+'-odir='+qcdir+' '+'-savp -out='+qcdir+'/'+os.path.basename(os.path.basename(re.sub('.bam','',dedup)))+'SPP_table.txt'
            cmds.append(sppcmd)
            qsub_a_command('qqqq'.join(cmds),qcdir+'/'+os.path.basename(re.sub('.bam','',dedup))+'_SPP.script.sh','qqqq','20G')
        
        #===================
        # MACS2
        #===================
        if opts.step=='MACS2':
            #macs2 callpeak -t <(${adjustedBed}) -f BED -n "${peakFile}" -g "${genomeSize}" -p 1e-2 --nomodel --shift "${fragLen}" -B --SPMR --keep-dup all --call-summits
            fragLen=int(open(qcdir+'/'+samplename+'SPP_table.txt','r').readlines()[0].strip().split('\t')[2].split(',')[0])
            peakfile=peakdir+'/'+samplename+'.peaks'
            bed=tagAligndir+'/'+samplename+'.tagAlign.gz'

            if fragLen==0: #single-cut
                #shift reads first
                cmds.append('slopBed -i '+bed+' -g '+opts.chrSizes+' -l 75 -r -75 -s > '+bed+'.shift75bp.bed')
                macscmd='macs2 callpeak -t '+bed+'.shift75bp.bed'+' -f BED -n '+peakfile+' -g '+opts.genome+ ' -p 1e-2 --nomodel --shift 75 -B --SPMR --keep-dup all --call-summits'

            if fragLen!=0: #double cut
                macscmd='macs2 callpeak -t '+bed+' -f BED -n '+peakfile+' -g '+opts.genome+ ' -p 1e-2 --nomodel --shift '+str(fragLen/2)+' -B --SPMR --keep-dup all --call-summits'

            cmds.append(macscmd)
            qsub_a_command('qqqq'.join(cmds),peakdir+'/'+samplename+'_MACS.script.sh','qqqq','20G')
        
    if opts.by_replGroup: #call peaks and bigwigs by replicate group
        for replGroup in sample2group.keys():
            cmds=[]
            cmds.append('source '+opts.bashrc)
            bams=' '.join(list(sample2group[replGroup]['bam']))
            tagAligns=' '.join(list(sample2group[replGroup]['tagAlign']))
            shiftedTagAligns=re.sub('tagAlign.gz','tagAlign.gz.shiftByHalfFragLen.bed',' '.join(list(sample2group[replGroup]['tagAlign'])))
            if opts.step=='MACS2':
                #get mean frag Len
                fragLenValues=[]
                #os.path.basename(os.path.basename(re.sub('.bam','',dedup)))
                for f in sample2group[replGroup]['bam']:
                    fragLenValues.append(get_best_fragLen(qcdir+'/'+re.sub('.bam','',os.path.basename(f))+'SPP_table.txt'))
                print fragLenValues
                print int(np.median(fragLenValues))
                combined_tagAlign=peakdir+'/'+replGroup+'.merged.tagAlign.gz'
                cmds.append('zcat '+tagAligns+' | gzip > '+combined_tagAlign)
                cmds.append(macs2(combined_tagAlign,str(int(np.median(fragLenValues))/2),opts.genome,peakdir+'/'+replGroup+'.peaks.narrowPeak'))
                qsub_a_command('qqqq'.join(cmds),peakdir+'/'+replGroup+'_MACS.script.sh','qqqq','20G')
            
            if opts.step=='bigwig':
                cpbam=opts.out_dir+'/tmp/'+replGroup+'.bam'
                os.system('mkdir '+opts.out_dir+'/tmp')
                cmds.append('samtools merge '+cpbam+' '+bams)
                cmds.append('samtools index '+cpbam)
                cmds.append('module load python_anaconda/2.2.0')
                cmds.append('bamCoverage -b '+cpbam+' -o '+bigwigdir+'/'+replGroup+'.bw')
                print cmds
                qsub_a_command('qqqq'.join(cmds),bigwigdir+'/'+replGroup+'_bigwig.script.sh','qqqq','20G')

            if opts.step=='shifted_bigwig':
                mtagAlign=opts.out_dir+'/tmp/'+replGroup+'.shifted.merged.tagAlign.gz'
                os.system('mkdir '+opts.out_dir+'/tmp')
                cmds.append('cat '+shiftedTagAligns+' | gzip >  '+mtagAlign)
                cmds.append('bedToBam -i '+mtagAlign+' -g '+opts.chrSizes+' > '+mtagAlign+'.bam')
                cmds.append('samtools sort '+mtagAlign+'.bam '+mtagAlign+'.sorted')
                cmds.append('samtools index '+mtagAlign+'.sorted.bam')
                cmds.append('module load python_anaconda/2.2.0')
                cmds.append('bamCoverage --normalizeTo1x '+opts.normFactor+' --binSize 10 -b '+mtagAlign+'.sorted.bam'+' -o '+bigwigdir+'/'+replGroup+'shifted.merged.bw')
                print cmds
                qsub_a_command('qqqq'.join(cmds),bigwigdir+'/'+replGroup+'_shiftedbigwig.script.sh','qqqq','20G')

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
