in_n1n2value=$1
expectedfile=$2
out_n1n2value=$3
resolution=$4

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: LA_obsOverExp+with_expectedFile.sh <in_n1n2value> <expectedfile> <out_n1n2value> <resolution>"
    echo "in_n1n2value = Lieberman-Aiden RAWobserved, coordinates represent beginning of the bin"
    echo "expectedfile = expected file, each row i is the value by which one should divide the value for 2 fragments, separated by i-1 fragments"
    echo "out_n1n2value = name of output file, also in n1n2value format. should be .gz"
    echo "resolution =  resolution used for HiC"
    exit
fi

expectedtmp=${out_n1n2value}_expectedtmp
zcat -f ${expectedfile} | awk '{print $0"\t"NR}' | sort -k2b,2 > ${expectedtmp}
zcat -f ${in_n1n2value} | awk -v res=${resolution} '{dist=($2-$1)/res+1}{print $1"\t"$2"\t"$3"\t"dist}' | sort -k4b,4 | join -1 2 -2 4 -o 2.1 2.2 2.3 1.1 -e "NA" -t $'\t' ${expectedtmp} - | awk '{obsOverExp=$3/$4}{print $1"\t"$2"\t"obsOverExp}' |  gzip > ${out_n1n2value}


