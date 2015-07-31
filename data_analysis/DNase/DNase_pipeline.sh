usage()
{
usage: `basename $0` options                                                                                                                                                               
DNase pipelime
OPTIONS:                                                                                                                                                                                   
   -h     Show this message and exit                                                                                                                                                       
   --aligned_bam [Required] The bam of aligned data                                                                                                                                             
   --counts [Required] Column with the counts                                                                                                                                              
   --f1 [Required] Column with the name of fragment 1                                                                                                                                      
   --f2 [Required] Column with the name of fragment 2                                                                                                                                      
   -c     Overwrite [0]                                                                                                                                                                    
}

INFILE=
COUNTS_COL=
F1_COL=
F2_COL=
CLEAN=''

ARGS=`getopt -o "hc" -l "count_data:,counts:,f1:,f2:" -- "$@"`
eval set -- "$ARGS"

while [ $# -gt 0 ]; do
    case $1 in
        -h) usage; exit;;
 --count_data) INFILE=$2;shift 2;;
        --counts) COUNTS_COL=$2;shift 2;;
        --f1) F1_COL=$2; shift 2;;
        --f2) F2_COL=$2; shift 2;;
        -c) CLEAN=1; shift;;
        --) shift; break;;
    esac
done

echo "Preparing input===================="
echo "Count data: ${INFILE}"
echo "COUNTS_COL: ${COUNTS_COL}"
echo "F1: ${F1_COL}"
echo "F2: ${F2_COL}"

