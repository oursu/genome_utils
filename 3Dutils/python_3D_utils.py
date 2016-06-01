import numpy as np
import gzip

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
