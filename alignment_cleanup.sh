#!/bin/sh

# The purpose of this script is to take the merged bam alignment output from pre-process-pipe.sh and prepare it for variant calling. 
# References for this include: 
#   https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery

Help()
{
	# Display help 
	echo "To automate the processes cleaning up a bwa-mem2 alignment to HG38 prior to performing variant calling. This script assumes the bam has already undergone MergeBamAlignments."
	echo ""
	echo "To run this pipeline use the following command:"
	echo ""
	echo "		sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] --time=[] alignment_cleanup.sh [arguments]"
	echo ""
	echo "	You must supply just one of -b, -o or -l. The multiple input options are designed to provide some flexiblity to the location of a bam file."
	echo "	-b merged bam alignment"
	echo " 	-l locations; instead of directly providing the merged bam, you can provide the batch list containing one sample ID on each line, corresponding to that sample's folder name."
	echo "	   This option will then be looking for "/[sample_ID]/[sample ID]_out/[sample ID].mergedaln.bam", and will move into [sample ID]_out, so if your directories do not comply with that structure, use -b instead."
	echo "	-o path to orig bam; this is based on the assumption that in the first step (pre-process-pipe.sh) the user looped through a list of paths to original bams located in a directory other than the working directory."
	echo "	   This option will then parse that path name to find the presumed working directory name for a given sample and identify the merged bam; i.e. /[sample_ID]/[sample ID]_out/[sample ID].mergedaln.bam." 
	echo "	-r start script after markduplicates spark, no input just pass the flag; ie./ if your previous run fails but generated the *markdups_sort.bam correctly (OPTIONAL)"
	echo "	-h help"
	echo "" 
	echo "You can loop through a batch to pass locations or BAMs to -l and -o in the same way as pre-process-pipe.sh. The -b is mostly for passing a single BAM."
	echo ""
	echo "For WGS, I recommend providing at least --mem=96g and --cpus-per-task=24, gres=lscratch:400 and --time=1-06:00:00, assuming a bam of around 60GB."

}

restart='true'
while getopts "b:o:l:rh" option; do
   case $option in
	b) mbam=$OPTARG ;;
	o) origbam=$OPTARG ;;
	l) location=$OPTARG ;;
	r) restart='false' ;;
	h) # display Help
         Help
         exit;;
   esac
done

if [ -z "$location" ] && [ -z "$origbam" ]
then
	sample=`basename $mbam | cut -f1 -d'.'`
	printf "Working directory: ${PWD}\n"
	printf "Merged alignment bam: $mbam\n"
elif [ -z "$mbam" ] && [ -z "$location" ]
then
	sample=`basename $origbam | cut -f1 -d'.'`
	printf "\nDirectory for sample: ${sample}\n"
	printf "No location or merged bam provided, original bam provided: $origbam\n"
	cd ${sample}/${sample}_out/
	printf "\nWorking directory: ${PWD}\n"
	mbam=${sample}.mergedaln.bam
	printf "Merged alignment bam under the name: ${mbam}\n\n"
else
	sample=$(echo $location)
	printf "The location of the working directory was provided as: $location. Moving into ${location}/${location}_out to perform alignment post-processing.\n"
	cd $location/"${location}_out"
	mbam=`ls $location.mergedaln.bam`
	printf "Merged alignment bam under the name: ${mbam}\n"
fi 

module load GATK/4.6.0.0 || fail "Could not load GATK"

if [[ ! -d /lscratch/${SLURM_JOBID} ]]; then mkdir -p /lscratch/${SLURM_JOBID}; fi

markdupssam=`echo $mbam | sed -E 's/mergedaln.sam|mergedaln.bam|merged.bam|merged.sam/markdups_sort.bam/g'`


##The sets of known variants are required for the BQSR modeling. 
knownindels=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.known_indels.vcf.gz
knownsnps=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.dbsnp138.vcf.gz
millsgold=/data/Kastner_PFS/references/HG38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta

set -x

if ${restart}  ##Mark duplicates takes a long time, this modularizes it a bit in the event something happens with the model but the mark duplicates bam was created just fine. 
then
	gatk MarkDuplicatesSpark --java-options "-Xms60G -Xmx60G" \
	-I ${mbam} \
	-O $markdupssam \
	-M markdup_metrics.txt \
	--create-output-bam-index true \
	--remove-all-duplicates true \
	--spark-master local[12] 2> err.${SLURM_JOB_ID}.txt
else
	continue 
fi 

gatk BaseRecalibrator --java-options "-Xms12G -Xmx12G -XX:ParallelGCThreads=4" \
	-R $ref \
	-I $markdupssam \
	-O ${sample}_bqsr_1.table \
	--known-sites $knownindels \
	--known-sites $knownsnps \
	--known-sites $millsgold \
	--tmp-dir /lscratch/${SLURM_JOB_ID} > log-bqsr-${SLURM_JOB_ID}.txt 2>&1

gatk ApplyBQSR --java-options "-Xms12G -Xmx12G -XX:ParallelGCThreads=4" \
	-bqsr ${sample}_bqsr_1.table \
	-I $markdupssam \
	-O ${sample}_bqsr_1.bam \
	--tmp-dir /lscratch/${SLURM_JOB_ID} > log-ApplyBQSR-${SLURM_JOB_ID}.txt

gatk BaseRecalibrator --java-options "-Xms12G -Xmx12G -XX:ParallelGCThreads=4" \
	-R $ref \
	-I ${sample}_bqsr_1.bam \
	-O ${sample}_bqsr_2.table \
	--known-sites $knownindels \
	--known-sites $knownsnps \
	--known-sites $millsgold \
	--tmp-dir /lscratch/${SLURM_JOB_ID} > log-bqsr-2-${SLURM_JOB_ID}.txt 2>&1 

# Generating some QC plots so we can see how BQSR improves the quality scoring. 
gatk AnalyzeCovariates --java-options "-Xms12G -Xmx12G -XX:ParallelGCThreads=4" \
	-before ${sample}_bqsr_1.table \
	-after ${sample}_bqsr_2.table \
	-plots AnalyzeCovariates.pdf 2>&1 

rm "${markdupssam}" ##removing the intermediate mark duplicates bam 
set +x


