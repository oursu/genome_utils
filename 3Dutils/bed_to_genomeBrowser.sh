
bed=$1
out=$2

zcat -f ${bed} | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"NR"\t."}' > ${out}.tmp
cat ${out}.tmp | sort -k1,1 -k2n,2n -k3n,3n | awk '{print $0"\t"NR"\t."}' > ${out}
bgzip -c ${out} > ${out}.gz
tabix -p bed ${out}.gz
rm ${out} ${out}.tmp
