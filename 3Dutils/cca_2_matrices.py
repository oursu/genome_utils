import sys
sys.path.append("/srv/scratch/oursu/3Dgenome/src/kCCA/ChIAPET/pyrcca/")
import rcca
import numpy as np

def main():
    parser=OptionParser()
    parser.add_option('--out',dest='out')
    parser.add_option('--matrices',dest='ms',default='',help='Comma delimited, .npy files, nodes should be aligned')
    opts,args=parser.parse_args()

    matrices=opts.ms.split(',')
    m1=np.load(matrices[0])
    m2=np.load(matrices[1])

    # Set up Pyrcca
    cca = rcca.CCA(kernelcca=False, numCC=2, reg=0.5)
    # Find canonical components
    training=cca.train([m1,m2])


main()
