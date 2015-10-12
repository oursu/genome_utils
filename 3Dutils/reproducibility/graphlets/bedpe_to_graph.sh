
in_bedpe=$1
out_graph=$2
output_type=$3

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: bedpe_to_graph.sh <in_bedpe> <out_graph> <output_type>"
    echo "Assumes the graph is given as it is. So if you want to set a value threshold, you need to do it outside of this code"
    echo "in_bedpe"
    echo "out_graph = prefix. Name of the file will be \${out_graph}.\${output_type}.graph"
    echo "output_type = can be fascia or orca, depening on which software will be used"
    exit
fi

out_graph=${out_graph}.${output_type}.graph
#graph files are:
#number of nodes, number of edges, labelled nodes, edges (presumably edges are sorted by the first node, then the second one
zcat -f ${in_bedpe} | awk '{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' | sort -k1,1 -k2,2n | uniq | \
awk '{zeroindex=NR-1}{print $1":"$2"-"$3"\t"zeroindex}' | sort -k1b,1 > ${out_graph}_nodes
zcat -f ${in_bedpe} | awk '{print $1":"$2"-"$3"\t"$4":"$5"-"$6}' | sort -k1b,1 | \
join -1 1 -2 1 -o 1.2 2.2 ${out_graph}_nodes -  | sort -k2b,2 | \
join -1 1 -2 2 -o 2.1 1.2 ${out_graph}_nodes - > ${out_graph}_edges
numNodes=$(wc -l ${out_graph}_nodes | sed 's/ /\t/g' | cut -f1)
numEdges=$(wc -l ${out_graph}_edges | sed 's/ /\t/g' | cut -f1)
if [[ ${output_type} == 'fascia' ]];
then
 echo "${numNodes}" > ${out_graph}
 echo "${numEdges}" >> ${out_graph}
fi
if [[ ${output_type} == 'orca' ]];
then
 echo "${numNodes} ${numEdges}" > ${out_graph}
fi
zcat -f ${out_graph}_edges >> ${out_graph}
zcat -f ${out_graph}_nodes | sed 's/:/\t/g' | sed 's/-/\t/g' | gzip > ${out_graph}_nodes.gz
rm ${out_graph}_edges ${out_graph}_nodes
