
in_bedpe=$1
out_bedpe_pref=$2
w=$3
chrSizes=$4


if [[ "$#" -lt 4 ]]
then
    echo "USAGE: divide_bedpe_in_windows.sh <in_bedpe> <out_bedpe_pref> <w> <chrSizes>"
    echo "in_bedpe = in bedpe file. Must have 8 entries."
    echo "out_bedpe_pref = prefix of output bedpe files"
    echo "w = window size, in bp"
    echo "chrSizes = chrSizes file, for instance for human it's /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes"
    exit
fi

mkdir -p $(dirname ${out_bedpe_pref})

#make windows
create_bed_fixedWindows_withName.sh ${chrSizes} ${w} ${out_bedpe_pref}_ws 'startOfWindow'
#intersect windows with the bedpe file
bedpecols=$(zcat -f ${in_bedpe} | awk '{print NF; exit}')
bedtools pairtobed -type either -a ${in_bedpe} -b ${out_bedpe_pref}_ws.w${w}.gz | \
awk -v out=${out_bedpe_pref} '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8>out".w"$12}'
cd $(dirname ${out_bedpe_pref})
for f in $(ls ${out_bedpe_pref}.w*);do gzip -f ${f};done
rm ${out_bedpe_pref}_ws.w${w}.gz
