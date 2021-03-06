usage()
{
cat "usage: `basename $0` options
Code by Oana Ursu: oursu@stanford.edu
Distance distribution comparison
OPTIONS:
   -h     Show this message and exit
   -a all regions to be considered
   -f comma-delimited bed files for filtering regions in -a. set to NA by default.
   -g aggregatedOver a bed file that tells you how to combine the regions in -a into elements
   -s sigfile significant (a subset of -g) elements to contrast with the whole set of aggregatedOver
   -d datset to which you want to compare distances
   -o out
   -b bashrc
   -w window for intersecting -a with -g
"
}

ALL=''
FILTER='NA'
AGGREGATED_OVER=''
SIGFILE=''
DISTANCE_TO=''
OUT=''
BASHRC=''
w='0'

while getopts "ha:f:g:s:d:o:b:w:" opt
do
    case $opt in
	h)
            usage; exit;;
	a)
            ALL=$OPTARG;;
	f)
            FILTER=$OPTARG;;
	g)
            AGGREGATED_OVER=$OPTARG;;
	s)
            SIGFILE=$OPTARG;;
	d) 
            DISTANCE_TO=$OPTARG;;
	o)
            OUT=$OPTARG;;
	b)
            BASHRC=$OPTARG;;
	w)
	    w=$OPTARG;;
	?)
            usage
	    exit 1;;
    esac    
done

testing(){
    BASHRC=/srv/gsfs0/projects/snyder/oursu/software/git/crispr_screens/bashrc_crispr_screens_scg4
    data=/srv/gsfs0/projects/kundaje/users/oursu/crispr/CRISPR_screens_2016_08_15/results/
    pref=${data}processing/hit_files/TSS/TSS.BinnedDataTo1000.relativeToplasmid.atLeast5guides
    ALL=${pref}.CRISPRi.castle_p_BH.0.05.allGuides.Concordant.CRKOcutsite.bed.gz
    FILTER=${pref}.allElements.bed.gz
    AGGREGATED_OVER=${FILTER}
    SIGFILE=${pref}.CRISPRi.castle_p_BH.0.05.sig_elements.allHits.bed.gz
    DISTANCE_TO=/srv/gsfs0/projects/kundaje/users/oursu/crispr/data/K562_DNase_model_june8_binaryPeaks_numPeaks_400000_pvalThresh_5_corePeakSize_500_totalPeakSize_1000_jitter_15.txt.gz
    OUT=/srv/gsfs0/projects/kundaje/users/oursu/test/testdistance
    w=0
}

source ${BASHRC}
if [[ ${FILTER} != 'NA' ]];
then
    filters_combined=${OUT}.filters.bed
    zcat -f $(echo ${FILTER} | sed 's/,/ /g') > ${filters_combined}
else
    filters_combined=${ALL}
fi

#filtering step
bedtools window -w 0 -a ${ALL} -b ${filters_combined} | cut -f1-4 | sort | uniq | sort -k1,1 -k2,2n | gzip > ${OUT}.all.filtered.bed.gz

#compute distances to dataset
distances=${OUT}.all.filtered.D2.$(basename ${DISTANCE_TO}).gz
zcat -f ${DISTANCE_TO} | sort -k1,1 -k2,2n | grep -v start | cut -f1-3 | bedtools closest -d -a ${OUT}.all.filtered.bed.gz -b stdin | cut -f1-4,8 | sort | uniq | gzip > ${distances}

#accumulate with aggregatedOver
aggregatedDistances=${OUT}.all.filtered.D2.$(basename ${DISTANCE_TO}).agg_$(basename ${AGGREGATED_OVER}).gz
bedtools window -w ${w} -a ${distances} -b ${AGGREGATED_OVER} | awk '{print $6"\t"$7"\t"$8"\t"$9"\t"$5}' | sort -k5 -n | awk '!seen[$4]++' | gzip > ${aggregatedDistances}

#ECDF test
Rscript ${genome_utils_path}/Features/TFs/distanceDistribution/ecdf_distanceDistribution.R ${aggregatedDistances} ${SIGFILE} ${OUT} 0 2000000

#rm ${filters_combined} ${OUT}.all.filtered.bed.gz ${distances} ${aggregatedDistances} 
