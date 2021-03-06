in_n1n2value=$1
normfile=$2
out_n1n2value=$3
resolution=$4

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: LA_normalize_with_normFile.sh <in_n1n2value> <normfile> <out_n1n2value> <resolution>"
    echo "in_n1n2value = Lieberman-Aiden RAWobserved, coordinates represent beginning of the bin"
    echo "normfile = norm files, each row i is the value by which one should divide the value for the i-th fragment (coordinate/res+1)"
    echo "out_n1n2value = name of output file, also in n1n2value format. should be .gz"
    echo "resolution =  resolution used for HiC"
    exit
fi

normtmp=${out_n1n2value}_normtmp
zcat -f ${normfile} | awk '{print $0"\t"NR}' | sort -k2b,2 > ${normtmp}
zcat -f ${in_n1n2value} | awk -v res=${resolution} '{n1=$1/res+1}{n2=$2/res+1}{print $1"\t"$2"\t"$3"\t"n1"\t"n2}' | sort -k4b,4| join -1 2 -2 4 -o 2.1 2.2 2.3 2.4 2.5  1.1 -e "NA" -t $'\t' ${normtmp} - | sort -k5b,5 | join -1 2 -2 5 -o 2.1 2.2 2.3 2.4 2.5 2.6 1.1 -e "NA" -t $'\t' ${normtmp} - | awk '{norm=$3/($6*$7)}{print $1"\t"$2"\t"norm}' | gzip > ${out_n1n2value}
rm ${normtmp}

