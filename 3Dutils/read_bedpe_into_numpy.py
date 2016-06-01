import numpy as np
import argparse
import gzip 

def main():
    parser = argparse.ArgumentParser(description='Store a bedpe file into a numpy array, indexed by entries in a bed file.')
    parser.add_argument('--bedpe',help='bedpe file')
    parser.add_argument('--bed',help='bed file')
    parser.add_argument('--out',help='out file, should end in .npy')
    parser.add_argument('--n1_col',type=int,default=1)
    parser.add_argument('--n2_col',type=int,default=4)
    parser.add_argument('--val_col',type=int,default=7)
    args = parser.parse_args()

    nodes,nodes_idx=read_nodes_from_bed(args.bed)
    nparray=edgelist2matrix(args.bedpe,args.n1_col,args.n2_col,args.val_col,nodes)
    np.save(args.out,nparray)
    
def read_nodes_from_bed(bedfile):
    nodes={}
    nodes_idx={}
    node_c=0
    for line in gzip.open(bedfile,'r'):
        items=line.strip().split('\t')
        chromo=items[0]
        start=items[1]
        end=items[2]
        node=items[3]
        if node in nodes.keys():
            print 'Error: node appears multiple times in your file!'
            sys.exit()
        if node not in nodes.keys():
            nodes[node]={}
            nodes[node]['idx']=node_c
            nodes[node]['chr']=chromo
            nodes[node]['start']=start
            nodes[node]['end']=end
            nodes_idx[node_c]=node
        node_c+=1
    return nodes,nodes_idx

#Assumes matrix is symmetric
def edgelist2matrix(in_edgelist,n1_col,n2_col,val_col,nodes):
    #returns a npy array of the counts.
    n=len(nodes.keys())
    c=0
    #keep records for row, col, value
    nparray=np.zeros((n,n))
    for line in gzip.open(in_edgelist,'r'):
        items=line.strip().split('\t')
        if items[n1_col] not in nodes.keys() or items[n2_col] not in nodes.keys():
            continue ###########
        n1=nodes[items[n1_col]]['idx']
        n2=nodes[items[n2_col]]['idx']
        counts=float(items[val_col])
        nparray[n1,n2]=counts
        nparray[n2,n1]=counts
    return nparray

main()

