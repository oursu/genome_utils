in_bed=$1
out_bed=$2

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: bed_to_unique_bed_entries.sh <in_bedpe> <out_bed>"
    echo "in_bed = infile bed, can be multiple files, comma-delimited"
    echo "out_bed  = bed file with unique entries from the 2 entries from <in_bedpe>"
    exit
fi

beds=$(echo ${in_bed} | sed 's/,/ /g')
zcat -f ${beds} | cut -f1-3 | grep -v start | sort | uniq | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' | gzip > ${out_bed}
