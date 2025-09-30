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
         -h help
      "

}

## Initalizing variables 

NUM_NON_REF=1
SAMPLES=""
PEDIGREE=""

while getopts ":s:p:h" opt; do
  case ${opt} in
    s)
      SAMPLES="${OPTARG}"
      ;;
    p)
      PEDIGREE="${OPTARG}"
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

query_jid=$(swarm sample-query-ped-${rundate}.swarm)
swarm --dependency=afterok:${query_jid} combinepairs-${rundate}.swarm