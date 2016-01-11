bam=$1 #bam file
fragFile=$2
out_n1n2value=$3

if [[ "$#" -lt 3 ]]
then
    echo "USAGE: bam_to_n1n2value.sh <bam> <fragFile> <out_n1n2value>"
    echo "bam = input bam"
    echo "fragFile = bed file with the fragments"
    echo "out_n1n2value = name of outfile, should end in .gz"
    exit
fi

sbam=${bam}sortedbyReadName
reads=${out_n1n2value}.reads.bedpe.gz

samtools sort -n ${bam} ${sbam}
bedtools bamtobed -bedpe -i ${sbam}.bam | gzip > ${reads}
bedtools intersect -wo -a ${fragFile} -b ${reads} | awk '{print $8"\t"$9"\t"$10"\t"$5"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3}' | \
bedtools intersect -wo -a ${fragFile} -b stdin | awk '{print $11"\t"$12"\t"$13"\t"$1"\t"$2"\t"$3}' | sort | uniq -c | \
awk '{print $2"\t"$4"\t"$5"\t"$7"\t"$1}' | gzip > ${out_n1n2value}
rm ${reads}
rm ${sbam}

