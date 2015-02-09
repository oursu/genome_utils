## FUNCTIONS
#############

OverlapRegionsWithTFBS () {
	regionfile=$1
	TFData=$2
	TFData_fileRegex=$3
	outdir=$4
	overlapWindow=$5
	aggregatedOver=$6
	out_prefix=$7
	bashrc=$8
	#make outdir if it doesn't exist
	mkdir -p ${outdir}
	#go through each TF dataset and compute overlap using bedtools
	chips=$(ls ${TFData}/${TFData_fileRegex})
	for chip in ${chips};
   	do
    	chip_base=$(basename ${chip} | sed 's/.gz//g')
    	outname=${outdir}/${out_prefix}.${chip_base}_OVERLAP_$(basename ${aggregatedOver}).tfbsCount

    	s=${outname}.sh
    	echo "source ${bashrc}" > ${s}
    	echo "echo ${chip_base} > ${outname}" >> ${s}
    	#windowbed between our regions and the TF: bedtools window -w ${overlapWindow} -c -a ${regionfile} -b -
    	#keep regions for which we found >=1 TFBS: awk '{if (\$5>0) print \$0}' 
    	#aggregate these over the aggregatedOver bed file: bedtools intersect -c -a ${peakfile} -b - 
    	#keep just the counts: cut -f5 
    	echo "zcat -f ${chip} | sed 's/chr//g' | \
    	bedtools window -w ${overlapWindow} -c -a ${regionfile} -b - | \
    	awk '{if (\$5>0) print \$0}' | \
    	bedtools intersect -c -a ${aggregatedOver} -b - | \
    	cut -f5 >> ${outname}" >> ${s} 
    	chmod 777 ${s}
    	echo $s
    	qsub -o ${s}.o -e ${s}.e ${s}	
    done
}


#Get arguments
regionfile=$1 #where to look for overlaps
TFData=$2 #directory with Tfs for which to compute the overlap
TFData_fileRegex=$3 #TF data name regex
outdir=$4 #directory where to write output
overlapWindow=$5 #window of overlap relative to the $regionfile
aggregatedOver=$6 #aggregate overlaps over this bed file
bashrc=$7 #bashrc that loads bedtools
out_prefix=$8
step=$9 #step: computeOverlaps, overlapsToMatrix, overlapEnrichment
sig_aggregatedOver=${10} #the significant set with which to overlap

if [[ $step == "computeOverlaps" ]];then
	echo "======computing overlaps"
	OverlapRegionsWithTFBS ${regionfile} ${TFData} ${TFData_fileRegex} ${outdir} ${overlapWindow} ${aggregatedOver} ${out_prefix} ${bashrc}
fi

overlapMatrix=${outdir}/${out_prefix}.overlapMatrix
if [[ $step == "overlapMatrix" ]];then
	echo "Making the overlap matrix"
	#combine all TFBS overlaps into 1 file
	echo chr_start_end_peak | sed 's/_/\t/g' > ${overlapMatrix}.tmp
	zcat -f ${aggregatedOver} >> ${overlapMatrix}.tmp
	paste ${overlapMatrix}.tmp ${outdir}/${out_prefix}*tfbsCount | gzip > ${overlapMatrix}.gz
	#remove scripts used for the previous step
	rm ${outdir}/${out_prefix}*tfbsCount 
	cat ${outdir}/${out_prefix}*sh > ${outdir}/${out_prefix}.totalScripts 
	rm ${outdir}/${out_prefix}*sh ${outdir}/${out_prefix}*sh.o ${outdir}/${out_prefix}*sh.e ${overlapMatrix}.tmp
fi

if [[ $step == "enrichmentAnalysis" ]];then
	source ${bashrc}
	#perform enrichment analysis
	echo "computing significance"
	${R_WITH_PACKAGES} ${enrichRcode} ${overlapMatrix}.gz ${sig_aggregatedOver} ${outdir}/${out_prefix}.overlapEnrichIN$(basename ${sig_aggregatedOver})
fi


