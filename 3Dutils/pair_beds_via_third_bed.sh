
bed1=$1
bed2=$2
combination_bed=$3
out=$4

if [[ "$#" -lt 4 ]]
then
    echo "USAGE: pair_beds_via_third_bed <bed1> <bed2> <combination_bed> <out>"
    echo "bed1 = first bed (e.g. enhancer coordinates)"
    echo "bed2 = second bed (e.g. gene coordinates)"
    echo "combination_bed = bed for combining bed1 and bed2 (e.g. domains)"
    echo "out = bedpe file to write the interactions between items in bed1 and items in bed2"
    exit
fi

cols_bed1=$(zcat -f ${bed1} | awk '{print NF; exit}')
cols_bed2=$(zcat -f ${bed2} | awk '{print NF; exit}')
cols_combo=$(zcat -f ${combination_bed} | awk '{print NF; exit}')
zcat -f    ${bed1} | perl -lane '$" = "\t"; print "@F[0]:@F[1]-@F[2]\t@F[0..$#F]"' | sort   -k1b,1 > ${out}.bed1
zcat -f   ${bed2} | perl -lane '$" = "\t"; print "@F[0]:@F[1]-@F[2]\t@F[0..$#F]"' | sort   -k1b,1 > ${out}.bed2
bedtools window -w 0 -a   ${bed1} -b ${combination_bed} | awk -v cols=${cols_bed1} '{chr_comb=cols+1;start_comb=cols+2;end_comb=cols+3;print $1":"$2"-"$3"\t"$chr_comb":"$start_comb"-"$end_comb}' | sort -k2b,2 > ${out}.tmp1
bedtools window -w 0 -a   ${bed2} -b ${combination_bed} | awk -v cols=${cols_bed2} '{chr_comb=cols+1;start_comb=cols+2;end_comb=cols+3;print $1":"$2"-"$3"\t"$chr_comb":"$start_comb"-"$end_comb}' | sort -k2b,2 > ${out}.tmp2
join -1 2 -2 2 ${out}.tmp1 ${out}.tmp2 | sed 's/ /\t/g' | cut --complement -f1 | sort -k1b,1 | join -1 1 -2 1 - ${out}.bed1 | sed 's/ /\t/g' | cut --complement -f1 | sort -k1b,1 | join -1 1 -2 1 - ${out}.bed2 | sed 's/ /\t/g' | cut --complement -f1 | gzip > ${out}
rm ${out}.bed1 ${out}.bed2 ${out}.tmp1 ${out}.tmp2



