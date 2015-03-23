#1=output script
#2=interaction
#3=bed of old coordinates
#4=bed of new coordinates
#5=outfile
#interaction has only the interacting fragments
TranslateInteraction () {
	s=$1
    intfile=$2
	inbed=$3
	outbed=$4
	outfile=$5
	sourcefile=$6
	echo "source ${sourcefile}" >> ${s}
	echo "bedtools intersect -wo -b ${outbed} -a ${inbed} | awk '{print \"chr\"\$1\"_\"\$4\"\t\"\$5\"_\"\$8}' | gzip > ${outfile}_transl" >> ${s}
	echo "zcat ${intfile} | sort -k1b,1 > ${intfile}_i" >> ${s}
	echo "zcat ${outfile}_transl | sort -k1b,1 > ${outfile}_transl_t" >> ${s}
	echo "join -1 1 -2 1 ${intfile}_i ${outfile}_transl_t | sort -k2b,2 | join -1 2 -2 1 - ${outfile}_transl_t | awk '{print \$3\"\t\"\$4\"\t\"\$2\"\t\"\$1}' | gzip > ${outfile}" >> ${s}
}

#1=bed file
#2=bedpe file
#3=outfile
#4=window around the bed file to look for intersections
#5=columns chosen;can be first,second,all (all is default). bedpe order of columns will stay unchanged
#produces a file of overlaps with the following format: 
#<bedpe entries> <bed file entries> for any windowed overlap
window_bedpe_with_bed (){
 bed=$1
 bedpe=$2
 out=$3
 w=$4
 columns_chosen=$5 #can be first,second,all (all is default). bedpe order of columns should stay unchanged

 echo $bed
 echo $bedpe 
 echo $out
 echo $w
 echo $columns_chosen
 module load bedtools/2.21.0
 numBedCols=$(awk '{print NF; exit}' ${bed})
 if [ ${columns_chosen} == "first" ] || [ ${columns_chosen} == "all" ];then 
  bedtools window -w ${w} -a ${bedpe} -b ${bed} > ${out}.overlapsCol1.bedpe
 fi
 if [ ${columns_chosen} == 'second' ] || [ ${columns_chosen} == 'all' ];then 
  zcat -f ${bedpe} | awk '{print $4"\t"$5"\t"$6"\t"$0}' | \
  bedtools window -w ${w} -a - -b ${bed} | cut --complement -f1-3 > ${out}.overlapsCol2.bedpe
 fi
 if [ ${columns_chosen} == 'first' ];then
  zcat -f ${out}.overlapsCol1.bedpe > ${out}
  rm ${out}.overlapsCol1.bedpe
 fi
 if [ ${columns_chosen} == 'second' ];then
  zcat -f ${out}.overlapsCol2.bedpe > ${out}
  rm ${out}.overlapsCol2.bedpe
 fi
 if [ ${columns_chosen} == 'all' ];then
  zcat -f ${out}.overlapsCol1.bedpe ${out}.overlapsCol2.bedpe | sort | uniq > ${out}
  rm ${out}.overlapsCol1.bedpe ${out}.overlapsCol2.bedpe
 fi
}
#example
#gwas=/srv/gsfs0/projects/kundaje/users/pgreens/projects/enh_gene_link_gwas/results/roadmap_EUR.AD_Igap_stage1/EUR.AD_Igap_stage1_sig_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed
#bedpe_file=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/QTLs/chr8.DistalQTLs_expandYRILD.bedpe
#outtest=/home/oursu/testW
#window_bedpe_with_bed ${gwas} ${bedpe_file} ${outtest} 0 first


translate_bedpe (){
 bedpe_file=$1
 output_bed=$2
 w=$3
 outname=$4
 keepUntranslated=$5 #keep_untranslated

 #make short bedpe - for now just keep the name, then bed file of the interacting regions
 #TODO: keep all entries from the bedpe
 zcat -f ${bedpe_file} | awk '{print $1":"$2"-"$3"\t"$4":"$5"-"$6"\t"$7}' | sed 's/chr//g' > ${outname}.shortbedpe
 zcat -f ${outname}.shortbedpe | cut -f1 > ${outname}.bedpe.bed.tmp
 zcat -f ${outname}.shortbedpe | cut -f2 >> ${outname}.bedpe.bed.tmp
 zcat -f ${outname}.bedpe.bed.tmp | sort | uniq | sed 's/:/\t/g' | sed 's/-/\t/g' | awk '{print $0"\t"$1":"$2"-"$3}' > ${outname}.bedpe.bed
 rm ${outname}.bedpe.bed.tmp 
 #overlap this bed file with the desired bedfile
 module load bedtools/2.21.0
 zcat -f ${output_bed} | \
 awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' | \
 bedtools window -w ${w} -a ${outname}.bedpe.bed -b - | cut -f4,8 | sort -k1b,1 > ${outname}.nameMapping
 #now, a simple join will replace our interactions with the new coordinates
 zcat -f ${outname}.shortbedpe | sort -k1b,1 | join -1 1 -2 1 ${outname}.nameMapping - | \
 sed 's/ /\t/g' | \
 cut --complement -f1 | \
 sort -k2b,2 | \
 join -1 1 -2 2 ${outname}.nameMapping - | \
 sed 's/ /\t/g' | \
 cut --complement -f1 | \
 awk '{gsub(":","\t",$2)"\t"$1"\t"$3}' | \
 sort | uniq > ${outname}
}

bedpe_file=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/QTLs/chr8.DistalQTLs_expandYRILD.bedpe
outname=/home/oursu/testTranslate
output_bed=/srv/gsfs0/projects/snyder/oursu/histoneQTL/networks/Data/GeneAnnoJudith.bed
w=0





