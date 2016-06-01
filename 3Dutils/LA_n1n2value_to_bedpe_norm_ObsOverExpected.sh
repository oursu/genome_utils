in_n1n2value=$1
norm=$2
obsOverExpected=$3
out_bedpe=$4
resolution=$5
chrSizes=$6
chromo=$7

if [[ "$#" -lt 1 ]]
then
    echo "USAGE: LA_n1n2value_to_bedpe_norm_obsOverExpected.sh <in_n1n2value> <norm> <obsOverExpected> <out_bedpe> <resolution> <chrSizes> <chromo>"
    echo "in_n1n2value = RAWobserved file"
    echo "norm = normalization scheme. Can be: 'nothing', 'KR', 'SQRTVC', 'VC'"
    echo "obsOverExpected = can be 'nothing' (which means keep entries as they are), RAWexpected, KRexpected, SQRTVCexpected, VCexpected"
    echo "out_bedpe = the name of the output bedpe file, without .gz"
    echo "resolution = HiC resolution, in bp"
    echo "chrSizes = chrSizes file, for instance for human it's /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes"
    echo "chromo = chromosome for the analysis"
    exit
fi

out_bedpe=${out_bedpe}.norm_${norm}.obsOverExp_${obsOverExpected}.bedpe
normf=$(echo ${in_n1n2value} | sed 's/RAWobserved/'${norm}norm/)
normalized_file=${out_bedpe}_${norm}norm.n1n2value.gz
obsOverExp_file=${out_bedpe}_${norm}norm.n1n2value.obsOverExp.gz
expected=$(echo ${in_n1n2value} | sed 's/RAWobserved/'${obsOverExpected}/)


#make files for the windows used
create_bed_fixedWindows_withName.sh ${chrSizes} ${resolution} ${out_bedpe}_windowsFile startOfWindow
convert_windowsFile_to_perChromosome.sh ${out_bedpe}_windowsFile.w${resolution}.gz

#operations on the values
#normalization
if [[ ${norm} == 'nothing' ]];
then
 cp ${in_n1n2value} ${normalized_file}
fi
if [[ ${norm} != 'nothing' ]];
then
 LA_normalize_with_normFile.sh ${in_n1n2value} ${normf} ${normalized_file} ${resolution}
fi
#obsOverExpected
if [[ ${obsOverExpected} != 'nothing' ]];
then
 LA_obsOverExp_with_expectedFile.sh ${normalized_file} ${expected} ${obsOverExp_file} ${resolution}
fi
if [[ ${obsOverExpected} == 'nothing' ]];
then
 mv ${normalized_file} ${obsOverExp_file}
fi 
echo "done values operations"

#conversion to bedpe
convert_n1n2value_to_bedpe.sh ${obsOverExp_file} ${out_bedpe}_windowsFile.w${resolution}_${chromo}.gz ${out_bedpe}.gz

#removing unnecessary files
#rm ${normalized_file} ${obsOverExp_file} ${out_bedpe}_windowsFile.w${resolution}*.gz
 
