usage()
{
cat "usage: `basename $0` options
Code by Oana Ursu: oursu@stanford.edu
Distance distribution comparison
OPTIONS:
   -h     Show this message and exit
   -s metadata for samples all guides, aggregated over, sigfile
   -d metadata distances (just a set of bed files one per line)
   -b bashrc
   -f metadata filters (just a set of bed files one per line)
   -o OUT
   -w window for intersecting samples with aggregator
"
}

META_SAMPLES=''
META_DISTANCES=''
BASHRC=''
META_FILTER=''
OUT=''
w='0'


while getopts "hs:d:b:f:o:w:" opt
do
    case $opt in
	h)
            usage; exit;;
	s)
            META_SAMPLES=$OPTARG;;
	d)
            META_DISTANCES=$OPTARG;;
	b)
            BASHRC=$OPTARG;;
	f)
            META_FILTER=$OPTARG;;
	o) 
            OUT=$OPTARG;;
	w)
	    w=$OPTARG;;
	?)
            usage
	    exit 1;;
    esac    
done

echo "running"

zcat -f ${META_SAMPLES}

while read line
do
    echo "reading"
    #skip lines starting with #                                                                                                                             
    if [[ $line == [#]* ]];then
        echo "skipping line $line because it starts with #"
    else
        read -a items <<< "$line"
	naming_sample=${items[0]}
        all=${items[1]}
	agg=${items[2]}
	sig=${items[3]}
	echo "read line"
	while read line
	do
	    if [[ $line == [#]* ]];then
		echo "skipping line $line because it starts with #"
	    else
		read -a items <<< "$line"
		naming_distance=${items[0]}
		distance_file=${items[1]}
		while read line
		do
		    if [[ $line == [#]* ]];then
			echo "skipping line $line because it starts with #"
		    else
			read -a items <<< "$line"
			naming_filter=${items[0]}
			filterfiles=${items[1]}
			#outpref=${OUT}/agg_$(basename ${agg})/TOTAL_$(basename ${all})/SIG_$(basename ${sig})/FILTER_to_$(basename ${filterfiles})/DISTANCE_to_$(basename ${distance_file})/Distance
			outpref=${OUT}/FilterTo_${naming_filter}/DistanceTo_${naming_distance}/FilterTo_${naming_filter}.DistanceTo_${naming_distance}.${naming_sample}
			mkdir -p $(dirname ${outpref})
			#now, run the distance code
			s=${outpref}.sh
			echo "source ${BASHRC}" > ${s}
			echo "distance_distribution.sh -a ${all} -f ${filterfiles} -g ${agg} -s ${sig} -d ${distance_file} -o ${outpref} -b ${BASHRC} -w ${w}" >> ${s}
			qsub -o ${s}.o -e ${s}.e ${s}
		    fi
		done < ${META_FILTER}
	    fi
	done < ${META_DISTANCES}
    fi
done < ${META_SAMPLES}
