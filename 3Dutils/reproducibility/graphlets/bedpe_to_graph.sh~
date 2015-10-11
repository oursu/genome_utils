
in_bedpe=$1
out_graph=$2

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: bedpe_to_graph_for_fascia.sh <in_bedpe>"
    echo "Assumes the graph is given as it is. So if you want to set a value threshold, you need to do it outside this code"
    exit
fi


#in_bedpe=/srv/scratch/oursu/3Dgenome/results/processed_data/HiC/counts/intra/GM12878_combined/5kb/chr21/chr21_5kb.RAWobserved.norm_SQRTVC.obsOverExp_SQRTVCexpected.bedpe.gz

#graph files are:
#number of nodes, number of edges, labelled nodes, edges (presumably edges are sorted by the first node, then the second one
zcat -f ${in_bedpe} | awk '{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' | sort -k1,1 -k2,2n | uniq | \
awk '{zeroindex=NR-1}{print $1":"$2"-"$3"\t"zeroindex}' | sort -k1b,1 > ${out_graph}_nodes
zcat -f ${in_bedpe} | awk '{print $1":"$2"-"$3"\t"$4":"$5"-"$6}' | sort -k1b,1 | \
join -1 1 -2 1 -o 1.2 2.2 ${out_graph}_nodes -  | sort -k2b,2 | \
join -1 1 -2 2 -o 2.1 1.2 ${out_graph}_nodes - > ${out_graph}_edges
numNodes=$(wc -l ${out_graph}_nodes | sed 's/ /\t/g' | cut -f1)
numEdges=$(wc -l ${out_graph}_edges | sed 's/ /\t/g' | cut -f1)
echo "${numNodes}" > ${out_graph}_notzip
echo "${numEdges}" >> ${out_graph}_notzip
echo "${numNodes} ${numEdges}" > ${out_graph}_notzip
zcat -f ${out_graph}_edges >> ${out_graph}_notzip
zcat -f ${out_graph}_notzip  > ${out_graph}
#zcat -f ${out_graph}_nodes | gzip > ${out_graph}_nodes.gz
rm ${out_graph}_notzip ${out_graph}_edges ${out_graph}_nodes
