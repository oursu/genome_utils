infile=$1
out=$2
res=$3
chromo=$4
chrSizes=$5
BL=$6
l1l2norm=$7
U=$8
L=$9

if [[ "$#" -lt 7 ]]
then
    echo "USAGE: LA_ICE.sh <INFILE> <OUT> <RES> <CHR> <chrSizes> <BL> <l1l2norm> <U> <L>"
    echo "The frags file lists the fragments by their start, rather than by their midpoint"
    echo "INFILE = input file (should be the RAWobserved file of counts), e.g. /srv/scratch/oursu/3Dgenome/data/HiC/counts/intra/GM12878_combined/5kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_5kb.RAWobserved"
    echo "RES = resolution of bins (in bp)"
    echo "OUT = output name"
    echo "CHR = chromosome name"
    echo "BL = blacklist file for filtering regions out. Usually /srv/scratch/oursu/data/wgEncodeDacMapabilityConsensusExcludable.bed"
    echo "chrSizes"
    echo "l1l2norm = l1 or l2 for ice"
    exit
fi

predir=${out}.res${res}.${chromo}
mkdir -p ${predir}
pre=${predir}/$(basename ${out}).res${res}.${chromo}
icepref=${pre}.${l1l2norm}.ICE

#does not have to be symmetric here (ICE considers it symmetric)
echo "Making interactions"
LA_to_FitHiC_make_interactionsFile.sh ${infile} ${res} ${predir}/$(basename ${out}) ${chromo} ${chrSizes} ${BL}
echo "Making fragments"  
LA_to_FitHiC_make_fragsFile.sh ${infile} ${res} ${predir}/$(basename ${out}) ${chromo} ${BL} ${chrSizes} nothing

echo "ice"
python /srv/scratch/oursu/code/software/pipelines/bds/hic/ICE-with-sparseMatrix.py ${pre}.interactions.gz ${pre}.frags.forICE.gz ${l1l2norm} ${icepref} 0.5

echo "running fithic"
zcat -f ${icepref} | gzip > ${icepref}.gz
zcat -f ${icepref}.biases | gzip > ${icepref}.biases.gz
mkdir -p ${predir}/fithic
rm ${icepref} ${icepref}.biases ${pre}.frags.forICE.gz
python /srv/scratch/oursu/code/software/pipelines/bds/hic/fit-hic-LATEST_std.py -f ${pre}.frags.gz -i ${pre}.interactions.gz -t ${icepref}.biases.gz -o ${predir}/fithic -U ${U} -L ${L} -l $(basename ${out}).res${res}.${chromo}.${l1l2norm}.ICE
echo "DONE fithic"


