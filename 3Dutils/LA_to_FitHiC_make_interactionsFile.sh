INFILE=$1
RES=$2
OUT=$3
CHR=$4
chrSizes=$5
BL=$6

if [[ "$#" -lt 6 ]]
then
    echo "USAGE: LA_to_FitHiC_make_interactionsFile.sh <INFILE> <RES> <OUT> <CHR> <chrSizes> <BL>"
    echo "INFILE = input file (should be the RAWobserved file of counts), e.g. /srv/scratch/oursu/3Dgenome/data/HiC/counts/intra/GM12878_combined/5kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_5kb.RAWobserved"
    echo "RES = resolution of bins (in bp)"
    echo "OUT = output name"
    echo "CHR = chromosome name"
    echo "chrSizes"
    echo "BL = e.g. /srv/scratch/oursu/data/wgEncodeDacMapabilityConsensusExcludable.bed"
    exit
fi

echo "Starting to make interactions file"
OUT=${OUT}.res${RES}.${CHR}.interactions
#make windows file                                                                                                                                                     
create_bed_fixedWindows_withName.sh ${chrSizes} ${RES} ${OUT}_windows startOfWindow
echo "create_bed_fixedWindows_withName.sh ${chrSizes} ${RES} ${OUT}_windows startOfWindow"
windowfile=${OUT}_windows.w${RES}.gz
zcat -f ${windowfile} | sort -k2b,2 | grep -w ${CHR} > ${windowfile}_unzip
#make bedpe file
convert_n1n2value_to_bedpe.sh ${INFILE} ${windowfile}_unzip ${OUT}_bedpe.gz
#remove blacklist from the bedpe file
bedtools pairtobed -type neither -a ${OUT}_bedpe.gz -b ${BL} |  awk '{if ($2!=$5) print $1"\t"$2"\t"$4"\t"$5"\t"int($8)}' | sed 's/chr//g' | gzip > ${OUT}.gz
rm ${windowfile}_unzip ${windowfile} ${OUT}_bedpe.gz
echo "DONE making interactions file"

