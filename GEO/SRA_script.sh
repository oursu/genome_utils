
#general command for downloading is
#/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/GEO/SRA_script.sh ${sra} /srv/gsfs0/projects/kundaje/users/oursu/Predict3D/HiCData/Fastq/ download
#/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/GEO/SRA_script.sh ${sra} /srv/gsfs0/projects/kundaje/users/oursu/Predict3D/HiCData/Fastq/ fastq


sra_path=$1
outdir=$2
step=$3
mkdir -p ${outdir}
sra_name=$(basename ${sra_path})
outname=${outdir}/${sra_name}

#here come the commands
if [[ ${step} == 'download' ]]; then
 s=${outname}_downloadScript.sh
 echo "module load sratoolkit/2.4.2" > ${s}
 echo "wget --directory-prefix ${outdir} ${sra_path}" >> ${s}
 chmod 755 ${s}
 ${s}
fi
if [[ ${step} == 'fastq' ]]; then
 s=${outname}_FastqScript.sh
 echo "module load sratoolkit/2.4.2" > ${s}
 echo "fastq-dump --split-3 --outdir ${outdir} ${outname}" >> ${s}
 echo "gzip $(echo ${outname} | sed 's/.sra//g')_1.fastq " >> ${s}
 echo "gzip $(echo ${outname} | sed 's/.sra//g')_2.fastq" >> ${s}
 chmod 755 ${s}
 qsub -l h_vmem=100G -o ${s}.o -e ${s}.e ${s}
fi
