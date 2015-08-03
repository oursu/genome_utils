
fq1=$1 
fq2=$2
sample=$3
out=$4
argfile=$5

mkdir ${out}/tmp

s=${out}/${sample}_mango.script.sh
echo "#Script for running mango for sample ${sample}" > ${s}
echo "zcat -f ${fq1} > ${out}/tmp/$(basename ${fq1}).fastq" >> ${s}
echo "zcat -f ${fq2} > ${out}/tmp/$(basename ${fq2}).fastq" >> ${s}
echo "module add r/3.2.0" >> ${s}
echo "module load bowtie/1.1.1" >> ${s}
echo "module load bedtools/2.16.2" >> ${s}
echo "Rscript /software/mango/mango.R --fastq1 ${out}/tmp/$(basename ${fq1}).fastq --fastq2 ${out}/tmp/$(basename ${fq2}).fastq --outdir ${out}/ --prefix ${sample} --argsfile ${argfile} --chromexclude chrX,chrM --stages 1:5" >> ${s}
qsub -l h_vmem=20G -e ${s}.e -o ${s}.o ${s}

