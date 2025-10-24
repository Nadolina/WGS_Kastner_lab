#!/bin/sh 

## For parsing input arguments to pre-process-pipe.sh and submitting to slurm

Help()
{
	# Display help 
	echo "To automate the processes of converting BAMs to uBAMs, marking adapters and aligning to the reference genome (currently GRCh38)."

	echo "To run this pipeline use the following command:"
	echo "sh pre-process-pipe.sh -b [original bam] OR -f [text file of bam paths]"

	echo "-b Full path to original BAM for realignment."
	echo "-f Text file containing batch of BAMs. Full path required, one per line."

}

while getopts "b:f:h" option; do
   case $option in
	b) alnbam=$OPTARG ;;
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
if [ -z "$alnbam" ] && [ -z "$batch" ]
then
	echo "No bam or batch file provided, please provide one using -b or -f"
	exit 1
elif [ -n "$batch" ] && [ -z "$alnbam" ]
then 
	echo "Batch file provided: ${batch}. BAMs to be processed:"
	cat $batch

	## If a batch file is provided, create a swarm file and submit to swarm
    echo "#SWARM -g 48 -t 16 --gres=lscratch:400 --time=1-00:00:00" > swarm_pre_process_${rundate}.txt
    mkdir -p pre_process_${rundate}_logs

    for bam in $(cat $batch); do
        echo "Submitting job for BAM: $bam"
        echo "sh ${sourcedir}/pre-process-pipe.sh -b ${bam}" >> swarm_pre_process_${rundate}.txt
    done

    swarm --module picard/3.2.0 --module bwa-mem2/2.2.1 --module fastqc/0.12.1 --module samtools/1.19 --logdir pre_process_${rundate}_logs swarm_pre_process_${rundate}.txt

elif [ -n "$alnbam" ] && [ -z "$batch" ]
then 

	## If a single bam is provided, submit pre-process-pipe.sh directly to slurm with the bam as input
	echo "Single bam provided: ${alnbam}"
    sbatch --mem=48g --cpus-per-task=16 --gres=lscratch:400 --time=1-00:00:00 --output=pre_process_%j.out --error=pre_process_%j.err ${sourcedir}/pre-process-pipe.sh -b ${alnbam}
else 
	echo "Please provide either a single bam (-b) or a batch file (-f), not both."
	exit 1
fi

