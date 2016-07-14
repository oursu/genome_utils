from optparse import OptionParser
import os
import math
from time import gmtime, strftime
import re
import gzip
import random
import numpy 
'''
Contact: Oana Ursu
oursu@stanford.edu
'''

def main():
    parser=OptionParser()
    
    parser.add_option('--out',dest='out',default='/srv/scratch/oursu/LongRange/NN/test/testMirrorPython')
    parser.add_option('--bedpe',dest='bedpe',default='/srv/scratch/oursu/LongRange/data/HiC/loops/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.bedpe')
    parser.add_option('--chr_sizes',dest='chrsizes',help='fa.fai',default='/mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa.fai')
    opts,args=parser.parse_args()

    print 'Building size dict'
    sizes={}
    CHR=0
    START=1
    END=2
    for line in open(opts.chrsizes,'r').readlines():
    	chromo,chromosize=line.strip().split('\t')[0:2]
    	sizes[re.sub('chr','',chromo)]=chromosize
    print sizes

    print 'Opening output file'
    out=open(opts.out,'w')

    print 'Mirroring'
    for line in open(opts.bedpe,'r').readlines():
    	e1_chr,e1_start,e1_end,e2_chr,e2_start,e2_end=line.strip().split('\t')[0:7]
    	smallest=0
    	if int(e2_end)<int(e1_start):
    		smallest=1
    	#now put the entries in order, starting with the smallest one. so 0 is the smallest, 1 is the largest
    	entries_init=[[e1_chr,int(e1_start),int(e1_end)],[e2_chr,int(e2_start),int(e2_end)]]
    	entries=[entries_init[smallest],entries_init[smallest-1]]
    	#get distance 
    	distance=int(abs(entries[0][END]-entries[1][START]))
    	#get entry sizes
    	entry_sizes=[abs(entries[0][START]-entries[0][END]),abs(entries[1][START]-entries[1][END])]
    	left_mirror=modify_position(entries[0][CHR],
    		int(entries[0][START]),
    		int(entries[0][END]),
    		distance,
    		entry_sizes[1],'left')
    	#only keep the example if it falls within the chromosome bounds
    	if left_mirror['start']<0:
    		continue
    	if left_mirror['end']>sizes[str(entries[0][CHR])]:
    		continue
    	right_mirror=modify_position(entries[1][CHR],
    		int(entries[1][START]),
    		int(entries[1][END]),
    		distance,
    		entry_sizes[0],'right')
    	if right_mirror['start']<0:
    		continue
    	if right_mirror['end']>sizes[str(entries[1][CHR])]:
    		continue
    	#now, write the mirrors to file
    	out.write(str(left_mirror['chromo'])+'\t'+str(left_mirror['start'])+'\t'+str(left_mirror['end'])+'\t'+str(entries[0][CHR])+'\t'+str(entries[0][START])+'\t'+str(entries[0][END])+'\n')
    	out.write(str(entries[1][CHR])+'\t'+str(entries[1][START])+'\t'+str(entries[1][END])+'\t'+str(right_mirror['chromo'])+'\t'+str(right_mirror['start'])+'\t'+str(right_mirror['end'])+'\n')


#modify (subtract/add position from the previous one)
def modify_position(chromo,start,end,distance,size_of_region,direction):
	new_pos={}
	new_pos['chromo']=chromo
	if direction=='left':
		new_pos['end']=start-distance
		new_pos['start']=new_pos['end']-size_of_region
	if direction=='right':
		new_pos['start']=end+distance 
		new_pos['end']=new_pos['start']+size_of_region
	return new_pos

main()