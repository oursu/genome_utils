
in_bed=$1
in_counts=$2
out_bigwig=$3
chrSizes=$4
WEBADDRESS=$5
WEBDIR=$6

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: graphlet_counts_to_bigwig.sh <in_bed> <in_counts> <out_bigwig> <chrSizes> <WEBADDRESS> <WEBDIR>"
    echo "in_bed = bed file with node coordinates"
    echo "in_counts = file with graphlet counts from orca"
    echo "out_bigwig = bigwig prefix"
    echo "chrSizes "
    echo "WEBADDRESS = http://mitra.stanford.edu/kundaje/oursu/"
    echo "WEBDIR = /srv/www/kundaje/oursu/"
    exit
fi

json=${out_bigwig}.json
zcat -f ${in_bed} | cut -f1-3 > ${out_bigwig}_tmp1
zcat -f ${in_counts} | sed 's/ /\t/g' > ${out_bigwig}_tmp2
countcols=$(zcat -f ${out_bigwig}_tmp2 | awk '{print NF; exit}')


echo "[{\"type\":\"native_track\",\"list\":[{\"name\":\"refGene\",\"mode\":\"full\"}]},{\"type\":\"coordinate_override\",\"coord\":\"chr9,36329955,chr9,37537411\"}," > ${json}.tmp
col=1
while [ "$col" -le "$countcols" ]; 
do
 echo $col
 colname=$(echo $col | awk '{name=$1-1}{print name}')
 pref=${out_bigwig}.${colname}
 bedgraph=${pref}.bedgraph
 bw=${pref}.bigwig
 zcat -f ${out_bigwig}_tmp2 | cut -f${col} | paste ${out_bigwig}_tmp1 - | sort -k1,1 -k2,2n > ${bedgraph}
 bedGraphToBigWig ${bedgraph} ${chrSizes} ${bw}
 #also make a json file, for easy loading
 echo "{\"type\":\"bigwig\",\"name\":\"$(basename ${pref})\",\"url\":\"${WEBADDRESS}/$(basename ${out_bigwig})/$(basename ${pref}).bigwig\",\"mode\":1,\"qtc\":{\"anglescale\":1,\"height\":40,\"summeth\":1,\"thtype\":0,\"pr\":0,\"pg\":0,\"pb\":255,\"smooth\":3}}," >> ${json}.tmp
 let "col+=1"
 rm ${bedgraph}
done
sed '$s/,$//' < ${json}.tmp > ${json}
echo "]" >> ${json}

mkdir -p ${WEBDIR}/$(basename ${out_bigwig})
mv ${out_bigwig}*bigwig ${WEBDIR}/$(basename ${out_bigwig})/
mv ${json} ${WEBDIR}/$(basename ${out_bigwig})/
rm ${json}.tmp ${out_bigwig}_tmp1 ${out_bigwig}_tmp2

