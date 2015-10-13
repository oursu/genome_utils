
in_bed=$1
featurefile=$2
featurescore=$3
out_features=$4


if [[ "$#" -lt 1 ]]
then
    echo "USAGE: make_1Dfeature.sh <in_bed> <featurefile> <featurescore> <out_features>"
    echo "in_bed = file of nodes"
    echo "featurefile = bed file with the feature as the name (4th) column" 
    echo "featurescore = can be count or binary, depending on how you want the feature to be reported"
    echo "out_features = name of out file (should end in .gz)"
    exit
fi

if [[ ${featurescore} == 'count' ]];
then
 binarize=''
fi

if [[ ${featurescore} == 'binary' ]];
then
 binarize="| awk '{binary=0}{if (\$1>0) binary=1}{print binary}'"
fi

options=$(zcat -f ${featurefile} | sed 's|/|\t|g' | cut -f4 | sort | uniq | sed 's/\n/ /g')
for o in ${options};
do 
 out=${out_features}.F${o}
 echo "${o}" > ${out}
 cmd="zcat -f ${featurefile} | grep -w ${o} | bedtools intersect -c -a ${in_bed} -b stdin | cut -f5 ${binarize} >> ${out}"
 eval $cmd
done

echo "chr_start_end" | sed 's/_/\t/g' > ${out_features}.tmp
zcat -f ${in_bed} | cut -f1-3 >> ${out_features}.tmp
paste ${out_features}.tmp ${out_features}.F* | gzip > ${out_features}
rm ${out_features}.tmp ${out_features}.F*

