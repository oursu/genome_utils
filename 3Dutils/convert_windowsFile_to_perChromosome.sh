wfile_bed=$1 #bed file with windows

if [[ "$#" -lt 2 ]]
then
    echo "Splits a file of genomic regions into separate files, one per chromosome"
    echo "USAGE: convert_windowsFile_to_perChromosome.sh <wfile_bed>"
    echo "wfile_bed = bed file of genomic windows"
    exit
fi



outpref=$(ls ${wfile_bed} | sed 's/.gz//g')_"byChromosome"

zcat -f ${wfile_bed} | awk -v out=${outpref} '{print $1"\t"$2"\t"$3"\t"$4>out$1}'
for f in $(ls ${outpref}*);do gzip $f;done
