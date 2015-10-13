
in_bedpe=$1
out=$2
output_type=$3
GRAPHLET_NODES=$4
annos=$5
chrSizes=$6
make_bigwig=$7

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: bedpe_to_graphlet_counts.sh <in_bedpe> <out> <output_type> <annos> <chrSizes> <make_bigwig>"
    echo "Assumes the graph is given as it is. So if you want to set a value threshold, you need to do it outside of this code"
    echo "in_bedpe"
    echo "out = prefix."
    echo "output_type = can be fascia or orca, depening on which software will be used. For now, only orca is supported."
    echo "GRAPHLET_NODES = can be 4 or 5 for orca"
    echo "annos = additional annotation files, for plotting profiles of graphlets by annotation. should 1) have a header, 2) first 3 entries are bed coordinates, the next columns are annotations. e.g. for chromatin states, the next columns are each a chromatin state, and the values are 0/1 dpeending on whether the node contains that chromatin state or not. set as NA if you don't have this file"
    echo "chrSizes"
    echo "make_bigwig = Yes or No"
    exit
fi

out_graph=${out}.${output_type}.graph
out_counts=${out}.${output_type}.counts
out_plot=${out}.${output_type}.plot

#make the graph
echo "making graph"
bedpe_to_graph.sh ${in_bedpe} ${out} ${output_type} 

#count the graphlets
echo "counting graphlets"
orca.exe ${GRAPHLET_NODES} ${out_graph} ${out_counts}

#plot the graphlets
echo "plotting graphlet counts"
#if annos is something, then we need to make its feature file
outannos=""
if  [[ ${annos} != 'NA' ]];
then
 for annofile in $(echo ${annos} | sed 's/,/ /g');
 do
  outanno=${out_graph}_nodes.gz_FEATURE_$(basename ${annofile})
  make_1Dfeature.sh ${out_graph}_nodes.gz ${annofile} binary ${outanno}
 done
 outannos=$(echo "${outannos},${outanno}")
fi
outannos=$(echo ${outannos} | sed 's/^,//g')
echo ${outannos}

Rscript /srv/scratch/oursu/code/genome_utils/3Dutils/reproducibility/graphlets/plot_graphlet_counts.R ${out_counts} ${out_plot} ${outannos}

zcat -f ${out_counts} | gzip > ${out_counts}.gz
zcat -f ${out_graph} | gzip > ${out_graph}.gz
rm ${out_counts} ${out_graph} 

if [[ ${make_bigwig} == 'Yes' ]];
then
 echo "Making bigwig files"
 graphlet_counts_to_bigwig.sh ${out_graph}_nodes.gz ${out_counts}.gz ${out}.${output_type}_bigwig ${chrSizes} http://mitra.stanford.edu/kundaje/oursu/ /srv/www/kundaje/oursu/
fi

rm -r ${out}.${output_type}_bigwig*

echo "DONE graphlet analysis"
