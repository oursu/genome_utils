in_bedpe=$1
translation_bed=$2
out_bedpe=$3
w=$4

if [[ "$#" -lt 2 ]]
then
    echo "USAGE: translate_bedpe.sh <in_bedpe> <translation_bed> <out_bedpe>"
    echo "in_bedpe = bedpe file input"
    echo "translation_bed = coordinates of entries to whch we want to translate"
    echo "out_bedpe = name of output file, should end in bedpe.gz"
    echo "w = window used for translating <in_bedpe> entries to the <translation_bed>"
    exit
fi

get_bedpe_entries(){
    local bedpe=$1
    local out=$2

    zcat -f ${bedpe} | cut -f1-3 | sort | uniq > ${out}.tmp
    zcat -f ${bedpe} | cut -f4-6 | sort| uniq >> ${out}.tmp
    zcat -f ${out}.tmp | sort | uniq | gzip > ${out}
    rm ${out}.tmp
}

echo "Starting translation"
get_bedpe_entries ${in_bedpe} ${out_bedpe}_entries.gz
bedtools window -w ${w} -a ${out_bedpe}_entries.gz -b ${translation_bed} | awk '{print $1":"$2"-"$3"\t"$4"\t"$5"\t"$6}' | sort -k1b,1 > ${out_bedpe}_translationFile
zcat -f ${in_bedpe} | awk '{print $1":"$2"-"$3"\t"$4":"$5"-"$6"\t"$0}' | \
sort -k1b,1 | join -e "NA" -t $'\t' -1 1 -2 1 ${out_bedpe}_translationFile - | cut --complement -f1 | \
sort -k4b,4 | join -e "NA" -t $'\t' -1 1 -2 4 ${out_bedpe}_translationFile - | cut --complement -f1,8-13 | gzip > ${out_bedpe}
rm ${out_bedpe}_entries.gz ${out_bedpe}_translationFile
echo "DONE translation, see/admire ${out_bedpe}"


