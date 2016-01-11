
in_bed=$1
out_bed_pref=$2
w=$3
chrSizes=$4


if [[ "$#" -lt 4 ]]
then
    echo "USAGE: divide_bed_in_windows.sh <in_bed> <out_bed_pref> <w> <chrSizes>"
    echo "in_bed = in bed. should have 4 entries"
    echo "out_bed_pref = prefix of output bed files"
    echo "w = window size, in bp"
    echo "chrSizes = chrSizes file, for instance for human it's /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes"
    exit
fi

mkdir -p $(dirname ${out_bed_pref})

echo ${in_bed}
echo ${out_bed_pref}
echo "=="
#make windows
create_bed_fixedWindows_withName.sh ${chrSizes} ${w} ${out_bed_pref}_ws 'startOfWindow'
#intersect windows with the bed file
bedtools intersect -wo -a ${out_bed_pref}_ws.w${w}.gz -b ${in_bed} | awk -v out=${out_bed_pref} '{print $5"\t"$6"\t"$7"\t"$8>out".w"$2}'
cd $(dirname ${out_bed_pref})
for f in $(ls ${out_bed_pref}.w*);do gzip -f ${f};done
rm ${out_bed_pref}_ws.w${w}.gz
