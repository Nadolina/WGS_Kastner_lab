#!/bin/sh 

## For parsing input arguments to alignment_cleanup.sh to allow for single or batch submissions 

Help()
{
	# Display help 
	echo "To automate the processes of converting BAMs to uBAMs, marking adapters and aligning to the reference genome (currently GRCh38)."

	echo "To run this pipeline use the following command:"
	echo "sh pre-process-pipe.sh -b [original bam] OR -f [text file of bam paths]"

	echo "-b Full path to original bam."
	echo "-f Text file containing batch of BAMs. Full path required, one per line"

}

while getopts "b:f:h" option; do
   case $option in
	b) mbam=$OPTARG ;;
	f) batch=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

rundate=$(date +%Y%m%d)
scriptpth="$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')"
sourcedir="$(dirname $scriptpth)"

## Parsing inputs -----
if [ -z "${mbam}" ] && [ -z "${batch}" ]
then
	echo "No bam or batch file provided, please provide one using -b or -f"
	exit 1
elif [ -n "${batch}" ] && [ -z "${mbam}" ]
then
	echo "Batch file provided: ${batch}. BAMs to be processed:"
	cat $batch

	## If a batch file is provided, create a swarm file and submit to swarm
    echo "#SWARM -g 80 -t 16 --gres=lscratch:400 --time=1-00:00:00" > alignment_cleanup_${rundate}.swarm
    mkdir -p alignment_cleanup_${rundate}_logs

    for bam in $(cat $batch); do
        echo "Submitting job for BAM: $bam"
        echo "sh ${sourcedir}/alignment_cleanup.sh -b ${bam}" >> alignment_cleanup_${rundate}.swarm
    done

    swarm --module GATK/4.6.0.0 --logdir alignment_cleanup_${rundate}_logs alignment_cleanup_${rundate}.swarm

elif [ -n "${mbam}" ] && [ -z "${batch}" ]
then 

	## If a single bam is provided, submit alignment_cleanup.sh directly to slurm with the bam as input
	echo "Single bam provided: ${mbam}"
    sbatch --mem=80g --cpus-per-task=24 --gres=lscratch:400 --time=1-00:00:00 ${sourcedir}/alignment_cleanup.sh -b ${mbam}
else 
	echo "Please provide either a single bam (-b) or a batch file (-f), not both."
	exit 1
fi
