in_bedpe=$1
out_bed=$2

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: bedpe_to_unique_bed_entries.sh <in_bedpe> <out_bed>"
    echo "in_bedpe = infile bedpe, can be multiple files, comma-delimited"
    echo "out_bed  = bed file with unique entries from the 2 entries from <in_bedpe>"
    exit
fi

bedpes=$(echo ${in_bedpe} | sed 's/,/ /g')
zcat -f ${bedpes} | awk '{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' | grep -v start | sort | uniq | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' | gzip > ${out_bed}
