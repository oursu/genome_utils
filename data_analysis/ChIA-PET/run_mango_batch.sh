
ARGFILE=/srv/scratch/oursu/code/genome_utils/data_analysis/ChIA-PET/mango_argfile
SRA_SCRIPT=/srv/scratch/oursu/code/genome_utils/GEO/SRA_script.sh

metadata=$1
datadir=$2
outdir=$3
step=$4 #download, fastq, mango

#first, one round of deleting the old scripts
while read line
do
    #skip lines starting with #
    if [[ $line == [#]* ]];then
       echo "skipping line $line because it starts with #"
    else
    read -a items <<< "$line"
    replicateGroup=${items[2]}
    #get rid of any script that existed before for this sample
    cat ${outdir}/${replicateGroup}_mango.batchScript.*sh >> ${outdir}/${replicateGroup}_mango.batchScript.sh.old
    fi
done < ${metadata}

#===============================
#now, actually process the data
while read line
do
    #skip lines starting with #
    if [[ $line == [#]* ]];then
       echo "skipping line $line because it starts with #"
    else
    read -a items <<< "$line"
    sra=${items[0]}
    sample=${items[1]}
    replicateGroup=${items[2]}
    s=
    #get rid of any script that existed before for this sample
    #==================
    #download sra files
    #==================
    if [[ ${step} == 'download' ]]; then
     ${SRA_SCRIPT} ${sra} ${outdir} download

    fi
done < ${metadata}