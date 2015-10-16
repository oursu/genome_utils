
in_reads=$1
out_n1n2value=$2
mapq_threshold=$3

F1=5
F2=9
MAPQ=10,11
zcat -f ${in_reads} | sed 's/ /\t/g' | cut -f${F1},${F2},${MAPQ} | awk -vmapq_t=${mapq_threshold} '{if ($3>=mapq_t && $4>=mapq_t) print $0}' | \
cut -f1,2 | sort | uniq -c | sed -e 's/ /\t/g' | awk '{print $2"\t"$3"\t"$1}' | gzip > ${out_n1n2value}


