#!/bin/bash
usage()
{
cat <<EOF
usage: `basename $0` options
Filters out 1) reads not properly paired, 2) unmapped reads and their mates, 3) reads on the M chromosome. Then it converts to the end tagAlign format.
OPTIONS:
   -h           Show this message and exit
   --inbam FILE  [Required] Input bam
   --outdir DIR [Required] Output directory.
   --sample STR [Required] Sample name. Input files are read from <indir>/<sample>.bam and output is written in <outdir>/<sample>.tagAlign.bed.gzip.
   --SE Set to 1 if single end, else set to 0
EOF
}

ARGS=`getopt -o "hcns" -l "inbam:,outdir:,sample:,SE:" -- "$@"`
eval set -- "$ARGS"

INBAM=
OUTDIR=
SAMPLE=
PROPER_PAIR="-f3 "
while [ $# -gt 0 ]; do
    case $1 in
    -h) usage; exit;;
    --inbam) INBAM=$2; shift 2;;
    --outdir) OUTDIR=$2; shift 2;;
    --sample) SAMPLE=$2; shift 2;;
    --SE) SE=$2; shift 2;;
    --) shift; break;;
    *) usage; exit 1;;
    esac          
done

if [ $# -ne 0 ]; then
    usage; exit 1;
fi

if [[ -z $INBAM || -z $OUTDIR || -z $SAMPLE ]]; then
    usage; exit 1;
fi

if [[ $SE -eq 1 ]]; then
    PROPER_PAIR=""
fi

script_location=${OUTDIR}/${SAMPLE}_remove_chrM_qc30_2tagAlign_script.sh
outpref=${OUTDIR}/${SAMPLE}.tagAlign.gz
mkdir -p ${OUTDIR}

samtools view ${PROPER_PAIR} -q30 -F1548 -h ${INBAM} | grep -v chrM | samtools view -Sh - | awk '$4 !="*"{print $0}' | samtools view -S -b - | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; $5="1000" ; print $0}' | gzip > ${outpref}
#qsub -l h_vmem=3G -o ${script_location}.o -e ${script_location}.e ${script_location} 
