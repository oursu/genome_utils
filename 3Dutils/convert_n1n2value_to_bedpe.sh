f=$1 #file with node1, node2, value
node_bed=$2 #bed file with coordinates of the nodes in $f, name column is the name of the node used in $f
out_bedpe=$3 #out file in bedpe format, shoud be gz 

if [[ "$#" -lt 2 ]]
then
    echo "USAGE: convert_n1n2value_to_bedpe.sh <f> <node_bed> <out_bedpe>"
    echo "f = file with node1, node2, value"
    echo "node_bed = bed file with coordinates of the nodes in $f, name column is the name of the node used in $f"
    echo "out_bedpe = out file in bedpe format, should be gz"                                                                                                      
    exit
fi


echo "====== Starting conversion from n1 n2 value to bedpe =========="
zcat -f ${node_bed} | sort -k4b,4 > ${out_bedpe}_node_bed_sorted
zcat -f ${f} | sort -k1b,1 | join -1 1 -2 4 -o 2.1 2.2 2.3 1.2 1.3 - ${out_bedpe}_node_bed_sorted | sort -k4b,4 | join -1 4 -2 4 -o 1.1 1.2 1.3 2.1 2.2 2.3 1.5 - ${out_bedpe}_node_bed_sorted | sed 's/ /\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNA\t"$7}' | gzip > ${out_bedpe}
rm ${out_bedpe}_node_bed_sorted
echo "====== DONE conversion from n1 n2 value to bedpe =========="
