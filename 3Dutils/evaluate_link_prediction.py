
from optparse import OptionParser
import os
import math
from time import gmtime, strftime
import re
import numpy as np
import fnmatch
'''
By Oana Ursu
oursu@stanford.edu
'''

def main():
    parser=OptionParser()
    parser.add_option('--true_bedpe',dest='true_bedpe',default='/srv/scratch/oursu/test/test_class.bedpe.gz')
    parser.add_option('--predicted_bedpe',dest='predicted_bedpe',help='Must include a score column (8th column)',default='/srv/scratch/oursu/test/test_pred.bedpe.gz')
    opts,args=parser.parse_args()

    pred_bedpe=read_in_bedpe(opts.predicted_bedpe)
    true_bedpe=read_in_bedpe(opts.true_bedpe)

    
def synchronize_pred_true(pred_bedpe,true_bedpe):
    for 
    
def read_in_bedpe(bedpefile):
    cse1_cse2_v={}
    for line in gzip.open(bedpefile,'r'):
        items=line.strip().split('\t')
        if int(items[1])>int(items[4]):
            cse1=items[0]+':'+items[1]+'-'+items[2]
            cse2=items[3]+':'+items[4]+'-'+items[5]
        if int(items[1])<=int(items[4]):
            cse2=items[0]+':'+items[1]+'-'+items[2]
            cse1=items[3]+':'+items[4]+'-'+items[5]
        v=items[7]
        if cse1 not in cse1_cse2_v.keys():
            cse1_cse2_v[cse1]={}
        if cse2 in cse1_cse2_v[cse1].keys():
            print "pair "+cse1+' '+cse2+' appears multiple times in your dataset '+bedpefile+'. Each prediction should have only one entry. Exiting ...'
            exit
        if cse2 not in cse1_cse2_v[cse1].keys():
            cse1_cse2_v[cse1][cse2]=v
        return cse1_cse2_v
        
