
in_reads=$1
out_n1n2value=$2
mapq_threshold=$3
intra=$4
res=$5

if [[ "$#" -lt 4 ]]
then
    echo "USAGE: LA_reads_to_n1n2value_bins.sh <in_reads> <out_n1n2value> <mapq_threshold> <intra> <res>"
    echo "intra = can be intra or inter"
    exit
fi

R1=4
R2=8
MAPQ=10,11

if [[ ${intra} == 'intra' ]];
then
 zcat -f ${in_reads} | sed 's/ /\t/g' | awk '{if ($3==$7) print $0}'  > ${out_n1n2value}.tmp
fi
if [[ ${intra} == 'inter' ]];
then
zcat -f ${in_reads} | sed 's/ /\t/g' > ${out_n1n2value}.tmp
fi

zcat -f ${out_n1n2value}.tmp | cut -f3,${R1},7,${R2},${MAPQ} | \
awk -v mapq_t=${mapq_threshold} -v res=${res} '{r1=int($2/res)*res}{r2=int($4/res)*res}{if ($5>=mapq_t && $6>=mapq_t) print $1"\t"r1"\t"$3"\t"r2}' | \
sort | uniq -c | sed -e 's/ /\t/g' | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' | gzip > ${out_n1n2value}

rm ${out_n1n2value}.tmp 

