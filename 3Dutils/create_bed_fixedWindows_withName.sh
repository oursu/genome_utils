chrSizes=$1 #chrSizes 
w=$2 #window size
out_bed=$3 #out bed file you get (will append to the name w${w}.gz)
naming=$4 #naming convention

if [[ "$#" -lt 2 ]]
then
    echo "Creates a bed file of genomic regions of fixed size, with each bed entry receiving a name (for instance, the fragment name)"
    echo "USAGE: create_bed_fixedWindows_withName.sh <chrSizes> <w> <out_bed> <naming>"
    echo "chrSizes = chrSizes file, for instance for human it's /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes "
    echo "w = window for the fixed size of the genomic regions. In bp"
    echo "out_bed = name of output file. Output file final name will be ${out_bed}.w${w}.gz"
    echo "naming = what to name each bed entry. Current options are endOfWindow and startOfWindow"
    exit
fi




echo "========= Starting to create window file =========="
naming_cmd=""
if [[ $naming == 'endOfWindow' ]];
then 
 naming_cmd="awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$3}' | "
fi

if [[ $naming == 'startOfWindow' ]];
then
 naming_cmd="awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$2}' | "
fi

cmd="bedtools makewindows -i winnum -w ${w} -g ${chrSizes} | ${naming_cmd} gzip > ${out_bed}.w${w}.gz"
eval $cmd

echo "========= DONE creating window file =========="
