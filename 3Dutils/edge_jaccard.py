from optparse import OptionParser
import numpy as np

def main():
    parser=OptionParser()
    parser.add_option('--out',dest='out',default='/srv/scratch/oursu/test/test_edge_jaccard')
    parser.add_option('--matrices',dest='ms',default='/srv/scratch/oursu/test/test_chr21.beta0.5.npy,/srv/scratch/oursu/test/test_chr21.beta0.5.npy',help='Comma delimited, .npy files, nodes should be aligned')
    parser.add_option('--thresh',dest='thresh',default='0.1')
    opts,args=parser.parse_args()

    thresh=float(opts.thresh)
    matrices=opts.ms.split(',')
    m1=np.load(matrices[0])
    m2=np.load(matrices[1])

    m1_bin=1*(m1>=thresh)
    m2_bin=1*(m1>=thresh)
    m1_e=np.count_nonzero(m1_bin)
    m2_e=np.count_nonzero(m2_bin)
    m1m2=np.multiply(m1_bin,m2_bin)
    m1m2_e=np.count_nonzero(m1m2)

    edges_union=m1_e+m2_e-m1m2_e
    edges_intersect=m1m2_e
    
    edge_jaccard=float(edges_intersect/edges_union)

    out=open(opts.out,'w')
    out.write('Edges_1\tEdges_2\tEdges_union\tEdges_intersect\tEdges_Jaccard\n')
    out.write(str(m1_e)+'\t'+str(m2_e)+'\t'+str(edges_union)+'\t'+str(edges_intersect)+'\t'+str(edge_jaccard)+'\n')
    out.close()

main()
