#TODO: write some tests for these

#=====================================================
#take entries from bedpe and mirror them
#e.g. for entry (DNA1,DNA2), will get (DNA1,DNA1-|DNA1-DNA2|) and (DNA2,DNA2+|DNA1-DNA2|)
#TODO: make sure the 2 entries are on the same chromosome
#TODO: combine files in a nicer way
#TODO: add slop. the problem is that now i get already <0 values when mirroring, 
#and I would need to not make the mirror in the first place when this happens, similar when chromosome ends are surpassed
mirror_bedpe(){
	local bedpe=$1
	local out=$2
	local combineMirrors=$3
	local CHRSIZES=$4

	combineMirrors=1
	out=/srv/scratch/oursu/LongRange/NN/test/testMirror
	bedpe=/srv/scratch/oursu/LongRange/data/HiC/loops/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.bedpe
	#split into DNA1>DNA2 and DNA1<DNA2 based on start
	zcat -f ${bedpe} | \
	awk -v outfile=${out} '{diff=$3-$5}{if (diff>=0) print $0>outfile"firstEntryAfterSecond.bedpe"}{if (diff<0) print $0>outfile"firstEntryBeforeSecond.bedpe"}'
	#add appropriate mirror left and right (slop the mirror as it gets created)
	for f in ${out}firstEntryAfterSecond.bedpe ${out}firstEntryBeforeSecond.bedpe;
	do
	 if [ -e ${f} ];
	 then
	  if [ ${f} == "${out}firstEntryBeforeSecond.bedpe" ];
	  then
	  entry1_chr=1
	  entry1_start=2
	  entry1_end=3
	  entry2_chr=4
	  entry2_start=5
	  entry2_end=6
	  fi
	  if [ ${f} == "${out}firstEntryAfterSecond.bedpe" ];
	  then
	   entry2_chr=1
	   entry2_start=2
	   entry2_end=3
	   entry1_chr=4
	   entry1_start=5
	   entry1_end=6
	  fi
	  #mirror left
	  zcat -f ${f} | 
	  awk -v e1_chr=${entry1_chr} -v e1_start=${entry1_start} -v e1_end=${entry1_end} \
	  -v e2_chr=${entry2_chr} -v e2_start=${entry2_start} -v e2_end=${entry2_end} \
	  'function abs(x){return ((x < 0.0) ? -x : x)} {diff=abs($e1_end-$e2_start)}{leftEnd=$e1_start-diff}{leftStart=leftEnd-abs($e2_start-$e2_end)}{print $e1_chr"\t"leftStart"\t"leftEnd"\t"$e1_chr"\t"$e1_start"\t"$e1_end}' > ${f}.mirrorLeft.bedpe
	  #cat ${f}.mirrorLeft.bedpe | cut -f1-3 | bedtools slop -b 0 -i stdin -g ${CHRSIZES} > ${f}.mirrorLeftSlop
	  #mirror right
	  zcat -f ${f} | \
	  awk -v e1_chr=${entry1_chr} -v e1_start=${entry1_start} -v e1_end=${entry1_end} \
	  -v e2_chr=${entry2_chr} -v e2_start=${entry2_start} -v e2_end=${entry2_end} \
	  'function abs(x){return ((x < 0.0) ? -x : x)} {diff=abs($e1_end-$e2_start)}{rightStart=$e2_end+diff}{rightEnd=rightStart+abs($e1_start-$e1_end)}{print $e2_chr"\t"$e2_start"\t"$e2_end"\t"$e2_chr"\t"rightStart"\t"rightEnd}' > ${f}.mirrorRight.bedpe
	 fi
	done
	if [ ${combineMirrors} == 1 ];
	 then
	 zcat -f ${out}firstEntry*Second.mirror*.bedpe > ${out}
	 rm ${out}firstEntry*Second.mirrorLeft.bedpe ${out}firstEntry*Second.mirrorRight.bedpe
	fi
}

#=====================================================
#1=bed file
#2=bedpe file
#3=outfile
#4=window around the bed file to look for intersections
#5=columns chosen;can be first,second,all (all is default). bedpe order of columns will stay unchanged
#produces a file of overlaps with the following format: 
#<bedpe entries> <bed file entries> for any windowed overlap
window_bedpe_with_bed (){
 local bed=$1
 local bedpe=$2
 local out=$3
 local w=$4
 local columns_chosen=$5 #can be first,second,all (all is default). bedpe order of columns should stay unchanged

 echo $bed
 echo $bedpe 
 echo $out
 echo $w
 echo $columns_chosen
 #module load bedtools/2.21.0
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

#============================================================
translate_bedpe (){
 local bedpe_file=$1
 local output_bed=$2
 local w=$3
 local outname=$4
 local keepUntranslated=$5 #keep_untranslated

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






