from optparse import OptionParser
import numpy as np
from scipy.stats import spearmanr

def main():
    parser=OptionParser()
    parser.add_option('--out',dest='out',default='/srv/scratch/oursu/test/test_spearman')
    parser.add_option('--matrices',dest='ms',default='/srv/scratch/oursu/test/test_chr21.beta0.5.npy,/srv/scratch/oursu/test/test_chr21.beta0.5.npy',help='Comma delimited, .npy files, nodes should be aligned')
    opts,args=parser.parse_args()

    matrices=opts.ms.split(',')
    m1=np.load(matrices[0])
    m2=np.load(matrices[1])

    print 'Computing Spearman correlation'
    sp=spearmanr(np.ndarray.flatten(m1),np.ndarray.flatten(m2))

    print 'Writing to file'
    out=open(opts.out,'w')
    out.write('Edge_Spearman_coeff\tEdge_Spearman_coeff_p\n')
    out.write(str(sp[0])+'\t'+str(sp[1])+'\n')
    out.close()

    print 'DONE!'

main()
