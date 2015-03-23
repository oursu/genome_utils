source_file=/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/Bedtools_like/Bedtools_for_pairs.sh
source ${source_file}
# datasets
#=========
DATA=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/
gwasdir=${DATA}/gwasOverlaps
qtldir=${DATA}/QTLs/
mkdir -p ${qtldir}
mkdir -p ${gwasdir}
#distal_qtls=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Distal/DistalQTLs.2015-02-02
distal_qtls=/srv/gsfs0/projects/kundaje/users/oursu/test/Distals_test
local_qtls=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Local/LocalQTLs.2015-02-25
#==========================
#Processing of QTL datasets
#==========================
#create file of local QTLs as follows: COMMENT THIS ONCE DONE!!!!!!
peakanno=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Data/CombinedPeakAnno.gz
peakanno_tss=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Data/CombinedPeakAnno.pointTSS
#making the merged regulatory elements, and lists mapping
regElts=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/Annotations/peaks.hmarks.dhs.bed.gz
annos=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/Annotations/
#- snps to regulatory elements
#- peaks to regulatory elements
genes=/srv/gsfs0/projects/snyder/oursu/histoneQTL/networks/Data/GeneAnnoJudith.bed
zcat -f ${peakanno} | grep RNA | awk '{start=$2+1500}{end=$2+1501}{print $1"\t"start"\t"end"\t"$4}' > ${peakanno_tss}
zcat -f ${peakanno} | grep -v RNA > ${peakanno}.notss 
zcat -f  ${peakanno}.notss ${peakanno_tss} > ${peakanno}.tssPlusPeaks
module load bedtools/2.21.0
zcat -f ${regElts} | sed 's/chr//g' | bedtools window -w 0 -a ${peakanno}.tssPlusPeaks -b - | cut -f4,8 | sort -k1b,1 > ${regElts}.2peaks.mapping
echo "regElt_Symbol" | sed 's/_/\t/g' > ${regElts}.toSymbol
zcat -f ${regElts} | bedtools window -w 0 -a - -b ${genes}.TSS.bed | \
cut -f4,8 | sort -k2b,2 | sed 's/chr//g' | \
join -1 2 -2 1 - /srv/gsfs0/projects/snyder/oursu/histoneQTL/networks/results/ensgenenames2.sorted | \
sed 's/ /\t/g' | awk '{print $2"\t"$3}' >> ${regElts}.toSymbol
zcat -f ${qtldir}/ALL.LocalQTLs.bedpe | sed 's/,/\t/g' | \
cut -f7,10 | sort -k2b,2 | join -1 2 -2 1 - ${regElts}.2peaks.mapping | \
sed 's/ /\t/g' | cut -f2-3 |  sort -k2b,2 | uniq | \
bedtools groupby -i - -g 2 -c 1 -o first > ${regElts}.toSNPs
echo "regElt_SNPandSymbol" | sed 's/_/\t/g' > ${regElts}.toSNPs_and_Symbol
zcat -f ${regElts}.toSymbol ${regElts}.toSNPs | sort -k2b,2 | uniq | \
bedtools groupby -i - -g 1 -c 2 -o collapse >> ${regElts}.toSNPs_and_Symbol
#label merged peaks with chromatin states
cs=/srv/gs1/projects/snyder/jzaugg/histoneQTL/chromatinStates/C*Core_mnemonics.bed
zcat -f /srv/gs1/projects/snyder/jzaugg/histoneQTL/chromatinStates/C*Core_mnemonics.bed | \
bedtools window -w 0 -a ${regElts}  -b - | sed 's/_/\t/g' | cut -f4,8 | \
sed 's/chr//g' | sort -k1b,1 | sort -k2b,2n | awk '!a[$1]++' > ${regElts}.cs 
# merge all enhancers
zcat -f /srv/gs1/projects/snyder/jzaugg/histoneQTL/chromatinStates/C*Core_mnemonics.bed | \
grep -v '13_\|14_\|15_' | sort -k1,1 -k2,2n | \
bedtools merge -i - > ${annos}/enhancerStates_mergedAcross7YRI.bed
zcat -f ${annos}/enhancerStates_mergedAcross7YRI.bed | \
bedtools window -w 0 -a ${regElts}  -b - | sed 's/_/\t/g' | cut -f4,8 | \
sed 's/chr//g' | awk '{print $1"\tENH"}' > ${regElts}.cs.enhancerStates_mergedAcross7YRI.nodeAnno
# get TSSs from Judith's file of genes
cat ${peakanno_tss} | awk '{print "chr"$0}' | \
bedtools window -w 0 -u -a ${regElts}  -b - | \
sed 's/_/\t/g' | cut -f4,8 | \
sed 's/chr//g' | awk '{print $1"\tTSS"}' > ${regElts}.cs.3kbaroundTSS_gencode.nodeAnno
# annotate as enh, tss and enh-tss by windowing these 2 files with the regElts
echo "node_TSSvsENH" | sed 's/_/\t/g' > ${regElts}.cs.mergedEnhancers_and_3kbaroundTSSgencode.nodeAnno
cat ${regElts}.cs.enhancerStates_mergedAcross7YRI.nodeAnno ${regElts}.cs.3kbaroundTSS_gencode.nodeAnno |
sort | uniq | \
sort -k1b,1 | bedtools groupby -g 1 -c 2 -o collapse >> ${regElts}.cs.mergedEnhancers_and_3kbaroundTSSgencode.nodeAnno


zcat -f /srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Local/SNPQTLmatrix/SNPQTLmatrix.*.gz | \
awk '{endsnp=$1+1}{if ($9=="pass") print $2"\t"$1"\t"endsnp"\t"$10"_"$4"\t"$3}' | sort -k4b,4 > ${local_qtls}.tmp
zcat -f ${peakanno}.tssPlusPeaks | sort -k4b,4 > ${peakanno}.sort
join -1 4 -2 4 ${local_qtls}.tmp ${peakanno}.sort | sed 's/ /\t/g' | \
awk '{print "chr"$2"\t"$3"\t"$4"\tchr"$6"\t"$7"\t"$8"\t"$5",chr"$2":"$3"-"$4",Local,"$1",chr"$6":"$7"-"$8}' > ${local_qtls}
rm ${peakanno}.sort ${local_qtls}.tmp
#now, make bedpe files for the local and distal QTLs
zcat -f ${distal_qtls} | awk -F "\t" -v qtldirec=${qtldir} '{plusone=$7+1}{if ($21=="QTL") print "chr"$18"\t"$7"\t"plusone"\tchr"$18"\t"$8"\t"$9"\t"$1",chr"$18":"$7"-"plusone",Distal,"$2",chr"$18":"$8"-"$9>qtldirec"/chr"$18".DistalQTLs.bedpe"}'
zcat -f ${local_qtls} | awk -F "\t" -v qtldirec=${qtldir} '{print $0>qtldirec"/"$1".LocalQTLs.bedpe"}' 

#======================
# LD expansion of QTLs
#====================== 
#Make a bedpe file with both types of QTLs, and expand with LD
for chromo in {1..22};
do
 #Distals
 out=${qtldir}/chr${chromo}.DistalQTLs_expandYRILD.bedpe
 s=${out}.sh
 echo "python /srv/gsfs0/projects/kundaje/users/oursu/code/ChromatinVariationCode/gwas/expandLD_by_column.py --out ${out}.tmp --chromo ${chromo} --snpcol 7 --input ${qtldir}/chr${chromo}.DistalQTLs.bedpe" > ${s}
 echo "echo \"chrLDSNP_startLDSNP_endLDSNP_chrPeak_startPeak_endPeak_QTLchr:start-end,LocalDistal,Peak,chrPeak:startPeak-endPeak\" | sed 's/_/\t/g' > ${out}.columns" >> ${s}
 echo "zcat -f ${out}.tmp | cut --complement -f4-6 > ${out}" >> ${s}
 echo "rm ${out}.tmp" >> ${s}
 chmod 711 ${s}
 #qsub -o ${s}.o -e ${s}.e ${s}
 #Locals
 out=${qtldir}/chr${chromo}.LocalQTLs_expandYRILD.bedpe
 s=${out}.sh
 echo "python /srv/gsfs0/projects/kundaje/users/oursu/code/ChromatinVariationCode/gwas/expandLD_by_column.py --out ${out}.tmp --chromo ${chromo} --snpcol 7 --input ${qtldir}/chr${chromo}.LocalQTLs.bedpe" > ${s}
 echo "echo \"chrLDSNP_startLDSNP_endLDSNP_chrPeak_startPeak_endPeak_QTLchr:start-end,LocalDistal,Peak,chrPeak:startPeak-endPeak\" | sed 's/_/\t/g' > ${out}.columns" >> ${s}
 echo "zcat -f ${out}.tmp | cut --complement -f4-6 > ${out}" >> ${s}
 echo "rm ${out}.tmp" >> ${s}
 chmod 711 ${s}
 qsub -o ${s}.o -e ${s}.e ${s}
done
#put them all together into 1 file
zcat -f ${qtldir}/chr*.DistalQTLs_expandYRILD.bedpe | grep rs* > ${qtldir}/ALL.DistalQTLs_expandYRILD.bedpe
zcat -f ${qtldir}/chr*.LocalQTLs_expandYRILD.bedpe | grep rs* > ${qtldir}/ALL.LocalQTLs_expandYRILD.bedpe
zcat -f ${qtldir}/ALL.DistalQTLs_expandYRILD.bedpe ${qtldir}/ALL.LocalQTLs_expandYRILD.bedpe > ${qtldir}/ALL.LocalandDistalQTLs_expandYRILD.bedpe
# ===========================
# Overlap with GWAS
#============================
gwases=$(echo $(ls /srv/gsfs0/projects/kundaje/users/pgreens/projects/enh_gene_link_gwas/results/roadmap_EUR*/*_sig_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed))
#gwases="/srv/gsfs0/projects/kundaje/users/pgreens/projects/enh_gene_link_gwas/results/roadmap_EUR.AD_Igap_stage1/EUR.AD_Igap_stage1_sig_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed /srv/gsfs0/projects/kundaje/users/pgreens/projects/enh_gene_link_gwas/results/roadmap_EUR.AdvancedAMD/EUR.AdvancedAMD_sig_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed /srv/gsfs0/projects/kundaje/users/pgreens/projects/enh_gene_link_gwas/results/roadmap_EUR.c4d_cad/EUR.c4d_cad_sig_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed"
#bedpe_interest=${qtldir}/ALL.DistalQTLs_expandYRILD.bedpe
bedpe_interest=${qtldir}/ALL.LocalandDistalQTLs_expandYRILD.bedpe
gwas=/srv/gsfs0/projects/kundaje/users/pgreens/projects/enh_gene_link_gwas/results/roadmap_EUR.AD_Igap_stage1/EUR.AD_Igap_stage1_sig_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed
for gwas in ${gwases};
do
 out=${gwasdir}/overlapQTL_$(basename ${gwas} | sed 's/.bed//g')/overlapQTL_$(basename ${gwas}).bedpe
 mkdir -p $(dirname ${out})
 s=${out}.sh
 #overlap snps with the gwas
 echo "module load bedtools/2.21.0" > ${s}
 echo "source ${source_file}" >> ${s}
 echo "window_bedpe_with_bed ${gwas} ${bedpe_interest} ${out}.overlappers.bedpe 0 first" >> ${s} #this gives us the QTL SNPs to keep (we'll keep them + all their LD)
 #go back to the LD file and pick out all LDed SNPs with this one
 echo "zcat -f ${out}.overlappers.bedpe | cut -f7 | sort | uniq > ${out}.qtlsnpsKeep" >> ${s}
 echo "zcat -f ${bedpe_interest} | sort -k7b,7 | \
 join -1 1 -2 7 ${out}.qtlsnpsKeep - | sed 's/ /\t/g' | \
 awk '{print \$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$7\"\t\"\$1}' > ${out}" >> ${s}
 #traslate affected peak to the merged peaks
 echo "zcat -f ${out} | sed 's/,/\t/g' | \
 sort -k10b,10 | \
 join -1 10 -2 1 - ${regElts}.2peaks.mapping | sed 's/ /\t/g' | sed 's/:/\t/g' | sed 's/-/\t/g' | \
 awk '{print \$2\"\t\"\$3\"\t\"\$4\"\t\"\$16\"\t\"\$17\"\t\"\$18\"\t\"\$8\",\"\$9\":\"\$10\"-\"\$11\",\"\$12\",\"\$16\":\"\$17\"-\"\$18\",\"\$1}' | sort | uniq > ${out}.affected2RegElements.bedpe" >> ${s}
 #translate SNP to the local peaks, then to the merged peaks
 echo "zcat -f ${out}.affected2RegElements.bedpe | \
 bedtools window -w 0 -a - -b ${qtldir}/ALL.LocalQTLs.bedpe | \
 sed 's/,/\t/g' | \
 awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$9\"\t\"\$10\"\t\"\$11\"\t\"\$21}' | 
 sort -k12b,12 | \
 join -1 12 -2 1 - ${regElts}.2peaks.mapping | sed 's/ /\t/g' | sort | uniq > ${out}.Reg2RegElements.long.bedpe" >> ${s}
 echo "zcat -f ${out}.Reg2RegElements.long.bedpe | \
 awk '{print \$13\"\t\"\$5\":\"\$6\"-\"\$7\"\t\"\$12\"\t\"\$10}' | \
 sed 's/:/\t/g' | sed 's/-/\t/g' | sed 's/_/\t/g' | \
 awk '{dist=int(((\$2+\$3)/2-(\$5+\$6)/2)/1000)}{posdist=dist;if (dist<0) posdist=-dist}{if (posdist!=0) print \$1\":\"\$2\"-\"\$3\"\t\"\$4\":\"\$5\"-\"\$6\"\t\"posdist\"\t\"\$7\"\t\"\$9}' | \
 sed 's/H3K4ME1/chromatin/g' | sed 's/H3K4ME3/chromatin/g' | sed 's/H3K27AC/chromatin/g' | sed 's/dhs/chromatin/g' | \
 sort | uniq > ${out}.Reg2RegElements" >> ${s}
 echo "zcat -f ${out}.Reg2RegElements | \
 awk '{print \"chr\"\$1\"\tchr\"\$2}' |  perl -lane 'printf qq[%s\n], join q[ ], sort @F' | sed 's/ /\t/g' | \
 sort | uniq | awk '{print \$0\"\t1\"}' > ${out}.Reg2RegElements.pairwiseTrack" >> ${s}
 chmod 755 ${s}
 qsub -o ${s}.o -e ${s}.e ${s}
done



















#=======================
# Overlap with each GWAS 
#=======================
overlap_bedpe_with_gwas () {
 #peyton's GWAS files have as name the tag SNP, so that's why you get the same SNP name mapped to multiple genomic coordinates
 gwases_input=$1
 outpref=$2
 bedpe_file=$3
 w=$4
 delimiter=$5
 gwases=$(echo ${gwases_input} | sed 's/'${delimiter}'/ /g')
 keepOnlyColumn=$6
 colOfName=$7
 outbed=$8
 for gwas in ${gwases};
 do
  echo ${gwas}
  out=${outpref}_$(basename ${gwas})
  s=${out}.sh
  echo "module load bedtools/2.21.0" > ${s}
  echo "echo \"#overlapType_#chrGWAS_#startGWAS_#endGWAS_#tagSNPGWAS_#chrQTL-LDSNP_#startQTL-LDSNP_#endQTL-LDSNP_#chrPeak_#startPeak_#endPeak_#QTL2Peak_#R2QTL-LDSNP_#LDSNP_#overlappedBases\" | sed 's/_/\t/g' > ${out}.columns" > ${s}
  for chromo in {1..22};
  do
   bedpe_current=$(echo ${bedpe_file} | sed 's/CHROMO/'${chromo}'/g')
   overlap_bed_with_bedpe ${s} ${gwas} ${bedpe_current} ${out} ${w}
  done
  #in out, we have a set of snp-peaks, where we know that the snp overlaps a gwas hit - make a bedpe from that
  bedpe_out=${out}.QTLoverlapsGWAS.SNP2Peak.bedpe
  echo "zcat -f ${out} | ${keepOnlyColumn} cut -f6-11,13 > ${bedpe_out}" >> ${s}
  #(also make edges from gwas tag snp to our snps[translated to merged peaks])

  #now, we keep GWAS-SNP with position, LD-QTLSNP with position, SNP-affectedPeak
  #looking for only column 1 overlaps, since we were overlapping SNPs
  echo "echo \"#chrGWAS_#startGWAS_#endGWAS_#tagSNPGWAS_#chrPeak_#startPeak_#endPeak_#QTL2Peak_#GWAS\" | sed 's/_/\t/g' > ${out}.summary.columns" >> ${s}
  echo "zcat -f ${out} | ${keepOnlyColumn}cut -f2-5,9-11,${colOfName} | sort | uniq | awk -v gwasname=$(basename ${gwas}) '{print \$0\"\t\"gwasname}' >> ${out}.summary" >> ${s}
  #make a signal track
  echo "zcat -f ${out}.summary | awk '{print \$1\":\"\$2\"-\"\$3\"\t\"\$5\":\"\$6\"-\"\$7\"\t1\"}' | sort | uniq > ${out}.pairwiseTrack" >> ${s}
  source_file=/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/Bedtools_like/Bedtools_for_pairs.sh
  echo "source ${source_file}" >> ${s}
  translated_bedpe=${bedpe_out}.2mergedPeaks.bedpe
  echo "translate_bedpe ${bedpe_out} ${outbed} 0 ${translated_bedpe}" >> ${s}
  genes=/srv/gsfs0/projects/snyder/oursu/histoneQTL/networks/Data/GeneAnnoJudith.bed
  #${genes}.TSS.bed
  #echo "translate_bedpe ${translated_bedpe} ${genes}.TSS.bed 0 ${translated_bedpe}.toGenes.bedpe" >> ${s}
  
  echo "zcat -f ${translated_bedpe} | \
  awk '{print \$1\":\"\$2\"-\"\$3\"\t\"\$4\":\"\$5\"-\"\$6\"\t1\"}' > ${translated_bedpe}.pairwiseTrack" >> ${s}
  chmod 755 ${s}
  qsub -o ${s}.o -e ${s}.e ${s}
 done
}

#gwases_input=$(echo $(ls /srv/gsfs0/projects/kundaje/users/pgreens/projects/enh_gene_link_gwas/results/roadmap_EUR*/*_sig_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed) | sed 's/ /DELIMIT/g')
outname_input=${gwasdir}/overlapGWAS_distalQTLs_
qtls_withLD_input=${qtldir}/chrCHROMO.DistalQTLs_expandYRILD.bedpe 
overlap_bedpe_with_gwas ${gwases_input} ${outname_input} ${qtls_withLD_input} 0 DELIMIT "grep overlapsCol1 | " 13 /srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/Annotations/peaks.hmarks.dhs.bed.gz























#====

# annotations
peakanno=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Data/CombinedPeakAnno.gz
genes=/srv/gsfs0/projects/snyder/oursu/histoneQTL/networks/Data/GeneAnnoJudith.bed
zcat -f ${genes} | awk '{start=$2+2000}{end=$2+2001}{print "chr"$1"\t"start"\t"end"\t"$4}' > ${genes}.TSS.bed




to_network (){
 source_file=/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/Bedtools_like/Bedtools_for_pairs.sh
 script=$1
 bedpe_current=$2
 outbed=$3
 outputname=$4
 echo "source ${source_file}" > ${script}
 echo "translate_bedpe ${bedpe_current} ${outbed} 0 ${outputname}" >> ${script}
 chmod 755 ${script}
 qsub -o ${script}.o ${script}.e ${script}	
}


