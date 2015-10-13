INFILE=$1
RES=$2
OUT=$3
CHR=$4
BL=$5
chrSizes=$6
norm=$7

MAP_CODE=/srv/scratch/oursu/code/genome_utils/3Dutils

if [[ "$#" -lt 7 ]]
then
    echo "USAGE: LA_to_FitHiC_make_fragsFile.sh <INFILE> <RES> <OUT> <CHR> <BL> <chrSizes> <norm>"
    echo "The frags file lists the fragments by their start, rather than by their midpoint"
    echo "INFILE = input file (should be the RAWobserved file of counts), e.g. /srv/scratch/oursu/3Dgenome/data/HiC/counts/intra/GM12878_combined/5kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_5kb.RAWobserved"
    echo "RES = resolution of bins (in bp)"
    echo "OUT = output name"
    echo "CHR = chromosome name"
    echo "BL = blacklist file for filtering regions out. Usually /srv/scratch/oursu/data/wgEncodeDacMapabilityConsensusExcludable.bed"
    echo "chrSizes"
    echo "norm = e.g. SQRTVC"
    exit
fi

echo "Starting to make frag file"
OUT=${OUT}.res${RES}.${CHR}
#make windows file 
create_bed_fixedWindows_withName.sh ${chrSizes} ${RES} ${OUT}_windows startOfWindow 
windowfile=${OUT}_windows.w${RES}.gz
zcat -f ${windowfile} | sort -k2b,2 | grep -w ${CHR} > ${windowfile}_unzip
#Now, get the marginalized count per fragment
zcat -f ${INFILE} | awk '{print $1"\t"$2"\t"int($3)}' > ${OUT}_interactions.tmp
zcat -f ${INFILE} | awk '{print $2"\t"$1"\t"int($3)}' >> ${OUT}_interactions.tmp
zcat -f ${OUT}_interactions.tmp | \
awk -v chromo=${CHR} '{a[$1]+=$3}END{for(i in a) print chromo"\tNA\t"i"\t"a[i]"\tNA"}' | gzip > ${OUT}.tmp.gz
#join windows file with the counts file, put a 0 where no count is available 
zcat -f ${OUT}.tmp.gz | sort -k3b,3 | join -a1 -1 2 -2 3 -e "0" -o 1.1 1.2 2.4 ${windowfile}_unzip - | sed 's/ /\t/g' | gzip > ${OUT}.tmp2.gz
#remove nodes from the blacklist
zcat -f ${OUT}.tmp2.gz | awk -v reso=${RES} '{end=$2+reso}{print $1"\t"$2"\t"end"\t"$3}' | sort | uniq | bedtools subtract -A -a stdin -b ${BL} > ${OUT}.tmp.rmBL
#get mappability of the bins
python ${MAP_CODE}/mappability_from_bed_bedout.py --regions ${OUT}.tmp.rmBL --read_length 101 --out ${OUT}.tmp.rmBL.mappability
zcat -f ${OUT}.tmp.rmBL.mappability | awk '{map=0}{if ($5>0.5) map=1}{print $1"\t0\t"$2"\t"$4"\t"map}' | gzip > ${OUT}.frags.gz
rm ${OUT}.tmp.rmBL ${OUT}.tmp.gz ${OUT}_interactions.tmp ${OUT}.tmp2.gz ${windowfile} ${windowfile}_unzip ${OUT}.tmp.rmBL.mappability
echo "DONE making frag file"
echo "making bias file"
normfile=$(echo ${INFILE} | sed 's/RAWobserved/'${norm}norm/)
zcat -f ${normfile} | awk -v res=${RES} -v chromo=${CHR} '{pos=(NR-1)*res}{print chromo"\t"pos"\t"$0}' | gzip > ${OUT}.bias${norm}.gz
echo "done bias file"



