
in_bedpe=$1
out_bedpe=$2
q=$3

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: bedpe_by_quantile_from_bedpe.sh <in_bedpe> <out_bedpe> <q>"
    echo "in_bedpe"
    echo "out_bedpe"
    echo "q = quantile as measured from the top (e.g. 0.1 is the top 10%)"
    exit
fi

numLines=$(zcat -f ${in_bedpe} | wc -l)
linesWanted=$(echo "${numLines} * ${q}" | bc -l | sed 's/[.]/\t/g' | cut -f1)
zcat -f ${in_bedpe} | sort -k8 -n | tail -n${linesWanted} | gzip > ${out_bedpe}
