#!/bin/sh 

## For submitting a sample query to the hail database 

set -euo pipefail

help() {

echo "For generating a sample query from the hail database.

      Usage:
      sh submit_sample_query.sh -s sample1,sample2 ... -p pedigree in *.fam format 
      
      Options:
         -s comma delimited list of samples to query, i.e./ a trio
         -p pedigree file (optional)
         -l in the form of chr:start-end (optional)
         -h help
      "

}

## Initalizing variables 

NUM_NON_REF=1
SAMPLES=""
PEDIGREE=""
LOCUS=""

while getopts ":s:p:l:h" opt; do
  case ${opt} in
    s)
      SAMPLES="${OPTARG}"
      ;;
    p)
      PEDIGREE="${OPTARG}"
      ;;
    l)
      LOCUS="${OPTARG}"
      ;;
    h)
      help
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

if [[ -z "${SAMPLES}" ]]; then 
   printf "ERROR: no samples provided to -s.\n"
   help
   exit 0
fi

IFS=',' read -r -a sarray <<< "${SAMPLES}"
SAMPLE_ARGS=""
OUTDIR="hail-query-"
for sample in "${sarray[@]}"; do 
   SAMPLE_ARGS+="${sample} "
   OUTDIR+="${sample}-"
done  
OUTDIR+="out"
echo "Results will be written to ${OUTDIR}"

module load hail

rundate=`date +'%m%d%y'`
mkdir -p ${OUTDIR} ${OUTDIR}/sample_query_${rundate}_logs 
printf "\nSamples to query: ${SAMPLE_ARGS}\n\n" 
printf "Using pedigree file: ${PEDIGREE}\n\n"
printf "Initializing sample query swarm file: sample-query-ped-${rundate}.swarm\n\n"
echo "#SWARM -g 25 -t 8 --gres=lscratch:100 --module hail --logdir sample_query_${rundate}_logs" > sample-query-ped-${rundate}.swarm
echo "#SWARM -g 25 -t 8 --gres=lscratch:100 --module hail --logdir sample_query_${rundate}_logs" > combinepairs-${rundate}.swarm

DB_PATHS=(
  "/data/Kastner_PFS/WGS/cohort_db/version_091725/db1.gatk.091725.mt"
  "/data/Kastner_PFS/WGS/cohort_db/version_091725/db2.gatk.091725.mt"
  "/data/Kastner_PFS/WGS/cohort_db/version_091725/db3.gatk.091725.mt"
  "/data/Kastner_PFS/WGS/cohort_db/version_091725/db4.gatk.091725.mt"
  "/data/Kastner_PFS/WGS/cohort_db/bcftools_091725/db1.bcftools.091725.mt"
  "/data/Kastner_PFS/WGS/cohort_db/bcftools_091725/db2.bcftools.091725.mt"
  "/data/Kastner_PFS/WGS/cohort_db/bcftools_091725/db3.bcftools.091725.mt"
  "/data/Kastner_PFS/WGS/cohort_db/bcftools_091725/db4.bcftools.091725.mt"
)

if [[ -z "${LOCUS}" ]]; then 
  printf "No locus provided, querying all variants for samples.\n\n"

  for db in "${DB_PATHS[@]}"; do
    echo "python3-hail hail-sample-ped-query-fxn.py --db ${db} \
      --samples ${SAMPLE_ARGS} \
      --pedigree ${PEDIGREE} \
      --outdir /lscratch/\${SLURM_JOB_ID}" >> sample-query-ped-${rundate}.swarm
  done

  for num in $(seq 1 4); do 
    echo "python3-hail hail-sample-ped-query-combinepairs.py --g_ht_path ${OUTDIR}/db${num}.gatk.091725.query.ht \
      --b_ht_path ${OUTDIR}/db${num}.bcftools.091725.query.ht \
      --outdir ${OUTDIR}" >> combinepairs-${rundate}.swarm
  done 

else
  printf "Querying variants at locus: ${LOCUS}\n\n"
  chr=$(echo "${LOCUS}" | cut -d':' -f1)

  DB=""

  if [[ "${chr}" == @(chr1|chr2|chr3|chr4) ]]; then 
    echo "Searching db1 for locus ${LOCUS}..."
    DB=1
  elif [[ "${chr}" == @(chr5|chr6|chr7|chr8|chr9|chr10) ]]; then 
    echo "Searching db2 for locus ${LOCUS}..."
    DB=2
  elif [[ "${chr}" == @(chr11|chr12|chr13|chr14|chr15|chr16|chr17) ]]; then 
    echo "Searching db3 for locus ${LOCUS}..."
    DB=3
  elif [[ "${chr}" == @(chr18|chr19|chr20|chr21|chr22|chrX|chrY) ]]; then 
    echo "Searching db4 for locus ${LOCUS}..."
    DB=4
  else 
    printf "ERROR: locus ${LOCUS} not recognized. Please use format chr:start-end with valid chromosome.\n"
    exit 0
  fi

  echo "python3-hail hail-sample-ped-query-fxn.py --db /data/Kastner_PFS/WGS/cohort_db/version_091725/db${DB}.gatk.091725.mt \
    --samples ${SAMPLE_ARGS} \
    --pedigree ${PEDIGREE} \
    --locus ${LOCUS} \
    --outdir /lscratch/\${SLURM_JOB_ID}" >> sample-query-ped-${rundate}.swarm
  echo "python3-hail hail-sample-ped-query-fxn.py --db /data/Kastner_PFS/WGS/cohort_db/bcftools_091725/db${DB}.bcftools.091725.mt \
    --samples ${SAMPLE_ARGS} \
    --pedigree ${PEDIGREE} \
    --locus ${LOCUS} \
    --outdir /lscratch/\${SLURM_JOB_ID}" >> sample-query-ped-${rundate}.swarm
  echo "python3-hail hail-sample-ped-query-combinepairs.py --g_ht_path ${OUTDIR}/db${DB}.gatk.091725.query.ht \
    --b_ht_path ${OUTDIR}/db${DB}.bcftools.091725.query.ht \
    --outdir ${OUTDIR}" >> combinepairs-${rundate}.swarm

fi


query_jid=$(swarm sample-query-ped-${rundate}.swarm)
swarm --dependency=afterok:${query_jid} combinepairs-${rundate}.swarm