
import sys
sys.path.append('/srv/scratch/oursu/code/software/python/pyGPs-1.3.2/pyGPs/GraphExtensions/')
sys.path.append('/srv/scratch/oursu/code/genome_utils/3Dutils/')
import nodeKernels
import numpy as np
from optparse import OptionParser
import os
import gzip
import python_3D_utils as utils3d

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
        m=utils3d.edgelist2matrix(opts.m,opts.sym,opts.node_list)
        
    #Do diffusion kernel on it
    print 'Diffusing the matrix'
    k=nodeKernels.diffKernel(m, beta=float(opts.beta))
    print k.shape
    np.save(opts.out+'.beta'+opts.beta, m)
    print 'DONE!'

main()
