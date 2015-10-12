
in_bedpe=$1
out_bedpe=$2
d_low=$3
d_high=$4

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: bedpe_by_distance_from_bedpe.sh <in_bedpe> <out_bedpe> <d_low> <d_high>"
    echo "in_bedpe"
    echo "out_bedpe"
    echo "d_low = Lower bound distance in bp"
    echo "d_high = Upper bound distance in bp"
    exit
fi

zcat -f ${in_bedpe} | awk -v low=${d_low} -v high=${d_high} '{dist=($5-$2)}{if (dist<0) dist=-dist}{if (dist>low && dist<high) print $0}' | gzip > ${out_bedpe}
