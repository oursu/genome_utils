
import sys
sys.path.append('/srv/scratch/oursu/code/software/python/pyGPs-1.3.2/pyGPs/GraphExtensions/')
import nodeKernels
import numpy as np
from optparse import OptionParser
import os
import gzip

def main():
    parser=OptionParser('~/software/python/anaconda/bin/python /srv/scratch/oursu/code/genome_utils/3Dutils/diffuse_graph.py --beta 0.5 --matrix testmat.gz --matrix_format edgelist --make_symmetric')
    parser.add_option('--out',dest='out')
    parser.add_option('--beta',dest='beta')
    parser.add_option('--matrix',dest='m',default='')
    parser.add_option('--matrix_format',dest='matrix_format',help='for now, edgelist')
    parser.add_option('--make_symmetric',dest='sym',action='store_true')
    parser.add_option('--node_list',dest='node_list',help='These nodes will be added to the edge list, even if they have no edges',default='')
    opts,args=parser.parse_args()

    #Reading in matrix
    print 'Reading in matrix'
    if opts.matrix_format=='edgelist':
        m=edgelist2matrix(opts.m,opts.sym,opts.node_list)
        
    #Do diffusion kernel on it
    print 'Diffusing the matrix'
    k=nodeKernels.diffKernel(m, beta=float(opts.beta))
    print k.shape
    np.save(opts.out+'.beta'+opts.beta, m)
    print 'DONE!'
    

def edgelist2matrix(in_edgelist,make_sym,node_list):
    edge_dict={}
    if node_list!='':
        #pre-allocate the nodes now
        nodes={}
        node_c=0
        for line in gzip.open(node_list,'r'):
            node=line.strip().split('\t')[0]
            if node not in nodes.keys():
                nodes[node]=node_c
                node_c+=1
    #initialize a huge matrix
    m=np.zeros((node_c,node_c))
    c=1 ### remove soon
    for line in gzip.open(in_edgelist,'r'):
        print str(float(c/1000)) ## remove soon
        c+=1 ## remove soon
        items=line.strip().split('\t')
        n1=nodes[items[0]] 
        n2=nodes[items[1]]
        counts=float(items[2])
        m[n1,n2]=counts
        if make_sym:
            if m[n2,n1]!=0:
                if n1!=n2:
                    print "Trying to add an edge although the edge has already been added before. Exiting ..."
                    exit
            m[n2,n1]=counts
    return m   

main()
