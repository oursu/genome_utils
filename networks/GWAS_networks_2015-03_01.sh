source_file=/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/Bedtools_like/Bedtools_for_pairs.sh
source ${source_file}
# datasets
#=========
DATA=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/
gwasdir=${DATA}/gwasOverlaps
qtldir=${DATA}/QTLs/
mkdir -p ${qtldir}
mkdir -p ${gwasdir}
distal_qtls=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Distal/DistalQTLs.2015-02-02
local_qtls=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Local/LocalQTLs.2015-02-25
# all marks and RNA annotations
annos=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/Annotations/
peakanno=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Data/CombinedPeakAnno.gz
# ENS annotation of TSSs (with gene symbols added if possible)
peakanno_RNA=/srv/gs1/projects/snyder/jzaugg/histoneQTL/peakAnnotation/peakLocations.ALL.RNA.txt
tss=${annos}/pointTSS.bed
ens2symbols=/srv/gsfs0/projects/snyder/oursu/histoneQTL/networks/results/ensgenenames2.sorted
regElements=${annos}/regElts.mergedHistonesDhs.bed
peaks2regElements=${annos}/peaks2regElts.mergedHistonesDhs.txt
peaks2TSS=${annos}/peaks2TSS.txt
TSSRegElt=${annos}/TSSRegEltNodeAnno.txt


#==========================
#Annotations
#==========================
#make TSS annotations
zcat -f ${peakanno_RNA} | sed 's/ /\t/g' | awk '{start=$3+1500}{end=$3+1501}{print $2"\t"start"\t"end"\t"$5}' | \
sort -k4b,4 | join -1 4 -2 1 -a1 -o 1.1 1.2 1.3 0 2.2 - ${ens2symbols} | grep -v gene_id | sed 's/ /\t/g' > ${tss}
#make the merged elements enhancers
module load bedtools/2.21.0
zcat -f ${peakanno} | grep -v RNA | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - | \
awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > ${regElements}
#map individual elements to the regulatory elements
zcat -f ${peakanno} | grep -v RNA | cut -f1-4 | \
bedtools window -w 0 -a - -b ${regElements} | cut -f4,8 | sort -k1b,1 > ${peaks2regElements}
# map histonepeaks/dhs to the genes
zcat -f ${peakanno} | cut -f1-4 | \
bedtools window -w 0 -a - -b ${tss} | cut -f4,8,9 | sort -k1b,1 > ${peaks2TSS}
#color nodes by TSS or reg elts
zcat -f ${peaks2TSS} | awk '{print $3"\tTSS"}' > ${TSSRegElt}.tmp
zcat -f ${regElements} | awk '{print $4"\tRegElt"}' >> ${TSSRegElt}.tmp
echo "node_TSSRegElt" | sed 's/_/\t/g' > ${TSSRegElt}
zcat -f ${TSSRegElt}.tmp | awk '!_[$1]++' >> ${TSSRegElt}
rm ${TSSRegElt}.tmp
#===========================

#==============
#make qtl data
#==============
zcat -f /srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Local/SNPQTLmatrix/SNPQTLmatrix.*.gz | \
awk -F "\t" -v qtldirec=${qtldir} '{endsnp=$1+1}{if ($9=="pass") print $2"\t"$1"\t"endsnp"\t"$3","$2":"$1"-"endsnp",Local,"$10"_"$4",NA">qtldirec"/chr"$2".LocalQTLs.bed"}'
zcat -f /srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Distal/DistalQTLs.2015-02-02 | \
awk -F "\t" -v qtldirec=${qtldir} '{plusone=$7+1}{print "chr"$18"\t"$7"\t"plusone"\t"$1","$18":"$7"-"plusone",Distal,"$15","$2>qtldirec"/chr"$18".DistalQTLs.bed"}'
rm ${qtldir}/chrchr.DistalQTLs.bed
#===============

#======================
# LD expansion of QTLs
#====================== 
#Make a bed file with both types of QTLs, and expand with LD
columns=chrLDSNP_startLDSNP_endLDSNP_QTL,chr:start-end,LocalDistal,LocalPeak,DistalPeak
for chromo in {1..22};
do
 #Distals
 out=${qtldir}/chr${chromo}.DistalQTLs_expandYRILD.bed
 s=${out}.sh
 echo "python /srv/gsfs0/projects/kundaje/users/oursu/code/ChromatinVariationCode/gwas/expandLD_by_column.py --out ${out}.tmp --chromo ${chromo} --snpcol 4 --input ${qtldir}/chr${chromo}.DistalQTLs.bed" > ${s}
 echo "echo \"${columns}_R2_LDSNP\" | sed 's/_/\t/g' > ${out}.columns" >> ${s}
 echo "zcat -f ${out}.tmp | cut --complement -f4-6 > ${out}" >> ${s}
 echo "rm ${out}.tmp" >> ${s}
 chmod 711 ${s}
 #qsub -o ${s}.o -e ${s}.e ${s}
 #Locals
 out=${qtldir}/chr${chromo}.LocalQTLs_expandYRILD.bed
 s=${out}.sh
 echo "python /srv/gsfs0/projects/kundaje/users/oursu/code/ChromatinVariationCode/gwas/expandLD_by_column.py --out ${out}.tmp --chromo ${chromo} --snpcol 4 --input ${qtldir}/chr${chromo}.LocalQTLs.bed" > ${s}
 echo "echo \"${columns}_R2_LDSNP\" | sed 's/_/\t/g' > ${out}.columns" >> ${s}
 echo "zcat -f ${out}.tmp | cut --complement -f4-6 > ${out}" >> ${s}
 echo "rm ${out}.tmp" >> ${s}
 chmod 711 ${s}
 #qsub -o ${s}.o -e ${s}.e ${s}
done
#put them all together into 1 file
zcat -f ${qtldir}/chr*.DistalQTLs_expandYRILD.bed | grep rs*  > ${qtldir}/ALL.DistalQTLs_expandYRILD.bed
zcat -f ${qtldir}/chr*.LocalQTLs_expandYRILD.bed | grep rs* > ${qtldir}/ALL.LocalQTLs_expandYRILD.bed
zcat -f ${qtldir}/ALL.DistalQTLs_expandYRILD.bed ${qtldir}/ALL.LocalQTLs_expandYRILD.bed > ${qtldir}/ALL.LocalandDistalQTLs_expandYRILD.bed
echo "${columns}_R2_LDSNP" | \
sed 's/_/\t/g' > ${qtldir}/ALL.LocalandDistalQTLs_expandYRILD.bed.columns
#===========================

#==============
# Convert to our regulatory elements (merged histone + dhs), TSSs
#==========================================================
#convert to regulatory elements the name of the local peak and the name of the distal peak
LOCAL_PEAKCOL=7
DISTAL_PEAKCOL=8
columns=chrLDSNP_startLDSNP_endLDSNP_QTL,chr:start-end,LocalDistal,LocalPeak,DistalPeak
pref=${qtldir}/ALL.LocalandDistalQTLs_expandYRILD
zcat -f ${qtldir}/ALL.LocalandDistalQTLs_expandYRILD.bed | \
sed 's/,/\t/g' | sort -k${LOCAL_PEAKCOL}b,${LOCAL_PEAKCOL} | \
join -e "NA" -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 2.2 -a1  -1 ${LOCAL_PEAKCOL} -2 1 - ${peaks2regElements} | \
sed 's/ /\t/g' | \
sort -k${DISTAL_PEAKCOL}b,${DISTAL_PEAKCOL} | \
join -e "NA" -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 2.2 -a1  -1 ${DISTAL_PEAKCOL} -2 1 - ${peaks2regElements} | \
sed 's/ /\t/g' > ${pref}.2regElts.bed 
echo "${columns},R2,LDSNP,RegEltLocal,RegEltDistal" | \
sed 's/_/\t/g' | sed 's/,/\t/g' > ${pref}.2regElts.bed.columns
#convert to TSS the name of the local peak and the name of the distal peak
zcat -f ${pref}.2regElts.bed | \
sed 's/,/\t/g' | sort -k${LOCAL_PEAKCOL}b,${LOCAL_PEAKCOL} | \
join -e "NA" -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 2.2 -a1 -1 ${LOCAL_PEAKCOL} -2 1 - ${peaks2TSS} | \
sed 's/ /\t/g' | \
sort -k${DISTAL_PEAKCOL}b,${DISTAL_PEAKCOL} | \
join -e "NA" -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 2.2 -a1 -1 ${DISTAL_PEAKCOL} -2 1 - ${peaks2TSS} | \
sed 's/ /\t/g' > ${pref}.2regElts.2TSS.bed 
echo "${columns},R2,LDSNP,RegEltLocal,RegEltDistal,TSSLocal,TSSDistal" | \
sed 's/_/\t/g' | sed 's/,/\t/g' > ${pref}.2regElts.bed.2TSS.columns

#========================
# GWAS
#=========================
make_net (){
  infile=$1
  outfile=$2
  pval_thresh=$3
  LDthresh=$4
  script=$5

  annos=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/Annotations/
  peaks2TSS=${annos}/peaks2TSS.txt
  tss=${annos}/pointTSS.bed
  regElements=${annos}/regElts.mergedHistonesDhs.bed

  echo "zcat -f ${infile} | grep Local | \
  awk '{if ((\$19<=${pval_thresh}) && (\$9>=${LDthresh})) print \$0}' | \
  awk '{n1=\$18;n2=\$11}{if (\$13!=\"NA\") n2=\$13}{print n1\"\t\"n2\"\t\"\$6\"\t\"\$7}' | \
  sed 's/_/\t/g' | \
  sed 's/H3K4ME1/chromatin/g' | sed 's/H3K4ME3/chromatin/g' | sed 's/H3K27AC/chromatin/g' | sed 's/dhs/chromatin/g' | \
  cut -f1-4 | sort | uniq > ${outfile}" >> ${script}

  echo "zcat -f ${infile} | grep Distal | \
  awk '{if ((\$19<=${pval_thresh}) && (\$9>=${LDthresh})) print \$0}' | \
  awk '{n1=\$11;n2=\$12}{if (\$13!=\"NA\") n1=\$13}{if (\$14!=\"NA\") n2=\$14}{print n1\"\t\"n2\"\t\"\$6\"\t\"\$8}' | \
  sed 's/_/\t/g' | \
  sed 's/H3K4ME1/chromatin/g' | sed 's/H3K4ME3/chromatin/g' | sed 's/H3K27AC/chromatin/g' | sed 's/dhs/chromatin/g' | \
  cut -f1-4 | sort | uniq >> ${outfile}" >> ${script}

  #annotate the genes with their symbol
  echo "zcat -f ${peaks2TSS} | sort -k2b,2 > ${outfile}.TSSsorted" >> ${script}
  echo "zcat -f ${outfile} | \
  sort -k1b,1 | join -1 2 -2 1 -a2 -o 2.2 2.3 2.4 1.3 2.1 ${outfile}.TSSsorted - | sed 's/ /\t/g' | awk '{print \$4\"\t\"\$1\"\t\"\$2\"\t\"\$3}' | \
  sort -k2b,2 | join -1 2 -2 2 -a2 -o 2.1 2.3 2.4 1.3 2.2 ${outfile}.TSSsorted - | sed 's/ /\t/g' | awk '{print \$1\"\t\"\$4\"\t\"\$2\"\t\"\$3}' | \
  sort | uniq > ${outfile}.withGeneNames.tmp" >> ${script}
  #annotate the genes with coordinates, so that we can compute distances between elements in our networks
  echo "zcat -f ${tss} | cut -f1-3,5 | zcat -f - ${regElements} | cut -f1-4 | sort -k4b,4 > ${outfile}.allElementsBed" >> ${script}
  echo "zcat -f ${outfile}.withGeneNames.tmp | grep Distal | sort -k1b,1 | \
  join -1 1 -2 4 -a1 -o 1.1 1.2 1.3 1.4 2.1 2.2 2.3 - ${outfile}.allElementsBed | sed 's/ /\t/g' | sort -k2b,2 | \
  join -1 2 -2 4 -a1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 2.1 2.2 2.3 - ${outfile}.allElementsBed | sed 's/ /\t/g' | \
  awk '{dist=int(((\$6+\$7)/2-(\$9+\$10)/2)/1000)}{posdist=dist;if (dist<0) posdist=-dist}{if (posdist!=0) print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"posdist\" kb\"}' > ${outfile}.withGeneNames" >> ${script}
  echo " zcat -f ${outfile}.withGeneNames.tmp | grep Local >> ${outfile}.withGeneNames" >> ${script}

  #and make nicer names for the nodes
  echo "zcat -f ${outfile}.withGeneNames | cut -f1 > ${outfile}.nodes" >> ${script}
  echo "zcat -f ${outfile}.withGeneNames | cut -f2 >> ${outfile}.nodes" >> ${script}
  echo "zcat -f ${outfile}.nodes | sort | uniq | grep -v ":" | awk '{print \$1\"\t\"\$1}' > ${outfile}.nodeLabels" >> ${script}
  echo "zcat -f ${outfile}.nodes | sort | uniq | grep ":" | \
  sed 's/:/\t/g' | sed 's/-/\t/g' | \
  awk '{n=((int(\$2+\$3)/2000)/1000)}{printf \"%s:%s-%s\tchr%s: %.3f Mb\n\", \$1,\$2,\$3,\$1,n}' >> ${outfile}.nodeLabels" >> ${script}
  echo "rm ${outfile}.withGeneNames.tmp ${outfile}.TSSsorted"
}

columns=chrLDSNP_startLDSNP_endLDSNP_QTL,chr:start-end,LocalDistal,LocalPeak,DistalPeak
gwases=$(echo $(ls /srv/gsfs0/projects/kundaje/users/pgreens/projects/hqtl_gwas/enrich_results/roadmap_EUR*/EUR*_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed))
pref=${qtldir}/ALL.LocalandDistalQTLs_expandYRILD
bed_interest=${pref}.2regElts.2TSS.bed
#gwas=/srv/gsfs0/projects/kundaje/users/pgreens/projects/hqtl_gwas/enrich_results/roadmap_EUR.CD/EUR.CD_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed
for gwas in ${gwases};
do
 out=${gwasdir}/overlapQTL_$(basename ${gwas} | sed 's/.bed//g')/overlapQTL_$(basename ${gwas}).LD08_p5.bed
 mkdir -p $(dirname ${out})
 s=${out}.sh
 GWAS_SNP_COL=18
 #overlap snps with the gwas
 echo "module load bedtools/2.21.0" > ${s}
 echo "source ${source_file}" >> ${s}
 echo "zcat -f $(echo ${gwas} | sed 's/no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed/no_TAG_SNPs_in_LD.bed/g') | sort -k4b,4 > ${out}.pvals" >> ${s}
 echo "bedtools window -w 0 -a ${bed_interest} -b ${gwas} | \
 sort -k${GWAS_SNP_COL}b,${GWAS_SNP_COL} | \
 join -1 4 -2 ${GWAS_SNP_COL} -o 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10 2.11 2.12 2.13 2.14 2.15 2.16 2.17 2.18 1.5  ${out}.pvals - | \
 sed 's/ /\t/g' > ${out}.overlapGWAS.bed" >> ${s} #this gives us the QTL SNPs to keep (we'll keep them + all their LD)
 echo "echo ${columns},R2,LDSNP,RegEltLocal,RegEltDistal,TSSLocal,TSSDistal_GWASchr_GWASstart_GWASend_GWAStagSNPID_GWASassociationPval | \
 sed 's/_/\t/g' | sed 's/_/\t/g' > ${out}.overlapGWAS.bed.columns" >> ${s}
 make_net ${out}.overlapGWAS.bed ${out}.overlapGWAS.network 0.00001 0.8 ${s}
 chmod 755 ${s}
 qsub -o ${s}.o -e ${s}.e ${s}
done

#make 1 big GWAS table
for gwas in ${gwases};
do
  out=${gwasdir}/overlapQTL_$(basename ${gwas} | sed 's/.bed//g')/overlapQTL_$(basename ${gwas}).bed
  zcat -f ${out}.overlapGWAS.bed | \
  awk -v gwasname=$(echo $(basename ${gwas}) | sed 's/_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed//g') '{print $0"\t"gwasname}' > ${out}.GWAStable.txt
done
columns=chrLDSNP_startLDSNP_endLDSNP_QTLSNP,QTLchr:start-end,QTLType,LocalPeak,DistalPeak
echo "${columns},R2,LDSNP,RegEltLocal,RegEltDistal,TSSLocal,TSSDistal_GWASchr_GWASstart_GWASend_GWAStagSNPID_GWASassociationPval_Disease" | \
sed 's/_/\t/g' | sed 's/,/\t/g' > ${gwasdir}/GWAStable.total.txt
zcat -f ${gwasdir}/overlapQTL_*/*GWAStable.txt >> ${gwasdir}/GWAStable.total.txt
zcat -f ${gwasdir}/GWAStable.total.txt | head -n1 > ${gwasdir}/GWAStable.total.GWASsig.txt
zcat -f ${gwasdir}/GWAStable.total.txt | awk '{if ($19<=0.00000001) print $0}' >> ${gwasdir}/GWAStable.total.GWASsig.txt




# make motif table
cols_we_want=c('SNP','affected.gene','snp.pos','localpeak','spearmanCorr_MotifScoreVSpeak',
  'pearsonCorr_MotifScoreVSpeak','maxObsMotif','motifName',
  'genomeWide_match.p.sp','crossTF_pBH.spearman','affected.peak.chr',
  'affected.peak.start','affected.peak.end')

motif_data=rbind(total_distal[,cols_we_want],
  total_local_all[,cols_we_want])


















awk -v reone=${re1} -v retwo=${re2} -v gone=${g1} -v gtwo=${g2} '{n1=reone;n2=retwo}{if ()}
# ===========================
# Overlap with GWAS
#============================
make_net (){
 infile=$1
 s=$2
 pval_thresh=$3
 LDthresh=$4

 PVAL_COL=17
 LD_COL=9
 #make the network files
 # - add all edges of SNP-local peak
 echo "zcat -f ${infile} | sed 's/,/\t/g' | \
 awk '{if ((\$${PVAL_COL}<=${pval_thresh}) && (\$${LD_COL}>=${LDthresh})) print \$16\"\t\"\$7\"\t\"\$11}' | sed 's/_/\t/g' |\
 cut -f1-3 | sed 's/H3K4ME1/chromatin/g' | sed 's/H3K4ME3/chromatin/g' | sed 's/H3K27AC/chromatin/g' | sed 's/dhs/chromatin/g' | \
 sort | uniq | awk '{print \$0\"Local\"}' > ${infile}.network" >> ${s}
 # - for distal, also add edges of localpeak ->distalpeak, if distalpeak is RNA
  echo "zcat -f ${infile} | sed 's/,/\t/g' | \
  awk '{if ((\$${PVAL_COL}<=${pval_thresh}) && (\$${LD_COL}>=${LDthresh})) print \$0}' | \
  grep Distal | \
  awk '{print \$7\"\t\"\$8\"\t\"\$12}' | sed 's/_/\t/g' |\
  cut -f1-3 | sed 's/H3K4ME1/chromatin/g' | sed 's/H3K4ME3/chromatin/g' | sed 's/H3K27AC/chromatin/g' | sed 's/dhs/chromatin/g' | \
  sort | uniq | \
  awk '{print \$0\"Distal\"}' >> ${infile}.network" >> ${s}
  echo "zcat -f ${infile}.network.tmp | grep Distal | sed 's/:/\t/g' | sed 's/-/\t/g' | \
  awk '{dist=int(((\$2+\$3)/2-(\$5+\$6)/2)/1000)}{posdist=dist;if (dist<0) posdist=-dist}{if (posdist!=0) print \$1\":\"\$2\"-\"\$3\"\t\"\$4\":\"\$5\"-\"\$6\"\t\"\$7\"\t\"posdist\" kb\"}' > ${infile}.network.dist" >> ${s}
  echo "zcat -f ${infile}.network.tmp | grep Local | sed 's/:/\t/g' | sed 's/-/\t/g' |\
  awk '{print \$1\"\t\"\$2\":\"\$3\"-\"\$4\"\t\"\$5\"\t\"}' >> ${infile}.network.dist" >> ${s}
  echo "zcat -f ${regElts}.toSymbol.sorted1 > ${infile}.network.nodeLabels.tmp" >> ${s}
  echo "zcat -f ${infile}.network | cut -f1 | grep rs | \
  awk '{print \$1\"\t\"\$1}' >> ${infile}.network.nodeLabels.tmp" >> ${s}
  #add in normal names of nodes, if they can't be converted to genes
  echo "zcat -f ${infile}.network | cut -f1,2 | \
  awk '{for(i=1;i<=NF;i++){print \$i}}' - | sort | uniq | \
  awk '{print \$1\"\t\"\$1}' >> ${infile}.network.nodeLabels.tmp" >> ${s}
  echo "zcat -f ${infile}.network.nodeLabels.tmp |\
  awk '!_[\$1]++' > ${infile}.network.nodeLabels" >> ${s}
  echo "rm ${infile}*.tmp"
}

make_simple_net (){
 infile=$1
 script=$2
 outf=$3

 echo "zcat -f ${infile} | grep Distal | cut -f4 | \
 sed 's/H3K4ME1/chromatin/g' | sed 's/H3K4ME3/chromatin/g' | sed 's/H3K27AC/chromatin/g' | sed 's/dhs/chromatin/g' | \
 sed 's/,/\t/g' | cut -f4,5,9 | sed 's/_/\t/g' | \
 cut -f1-3 | sort | uniq | sed 's/:/\t/g' | sed 's/-/\t/g' | \
 awk -v outfile=${outf} '{dist=int(((\$2+\$3)/2-(\$5+\$6)/2)/1000)}{posdist=dist;if (dist<0) posdist=-dist}{if (posdist!=0) print \$1\":\"\$2\"-\"\$3\"\t\"\$4\":\"\$5\"-\"\$6\"\t\"\$7\"\t\"posdist\" kb\">outfile\".chr\"\$1}'" > ${script} 
 chmod 777 ${script}
 qsub -o ${script}.o -e ${script}.e ${script}
}

make_simple_net ${bed_interest} ${bed_interest}.makeNet.sh ${bed_interest}.network














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


