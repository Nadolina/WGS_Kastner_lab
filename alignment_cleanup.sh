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
	echo "		sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] --time=days-hours:minutes:seconds alignment_cleanup.sh [arguments]"
	echo ""
	echo "	You must supply -b or -l, but do not supply both."
	echo "	-b merged bam alignment"
	echo " 	-l batch list of locations; instead of directly providing the merged bam, you can provide the batch list containing one sample ID on each line, corresponding to that sample's folder name."
	echo "	   This option will then be looking for "/[sample_ID]_out/[sample ID]_out/[sample ID].mergedaln.bam", and will move into [sample ID]_out, so if your directories do not comply with that structure, use -b instead."
	echo "	   Note that the script is still designed to just take one bam at a time, so the user will have to loop through the batch list themselves, this -l just permits for that. I will internalize looping capabilities in the future."
	echo "	   For now use something like:" 
	echo "		while read sample; do sbatch --mem=128g --cpus-per-task=32 --gres=lscratch:500 --time=1-05:30:00 $SCRIPTS/alignment_cleanup.sh -l $sample ; done < HC-batch-091324.txt"
	echo "	-r start script after markduplicates spark, no input just pass the flag; ie./ if your previous run fails but generated the *markdups_sort.bam correctly (OPTIONAL)"
	echo "	-h help"
	echo ""
	echo "For WGS, I recommend providing at least --mem=96g and --cpus-per-task=16 and --time=1-00:00:00, assuming a bam of around 60GB, scale appropriately."

}

restart='true'
while getopts "b:l:rh" option; do
   case $option in
	b) mbam=$OPTARG ;;
	l) location=$OPTARG ;;
	r) restart='false' ;;
	h) # display Help
         Help
         exit;;
   esac
done

if [ -z "$location" ]
then
	printf "no location provided, bam was provided: $alnbam\n"
else
	cd $location/"${location}_out"
	mbam=`ls $location.mergedaln.bam`
fi 

module load GATK/4.6.0.0 || fail "Could not load GATK"

if [[ ! -d /lscratch/${SLURM_JOBID} ]]; then mkdir -p /lscratch/${SLURM_JOBID}; fi

markdupssam=`echo $mbam | sed -E 's/mergedaln.sam|mergedaln.bam|merged.bam|merged.sam/markdups_sort.bam/g'`

printf "${PWD}/$markdupssam\n"

##The sets of known variants are required for the BQSR modeling. 
knownindels=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.known_indels.vcf.gz
knownsnps=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.dbsnp138.vcf.gz
millsgold=/data/Kastner_PFS/references/HG38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta

set -x

if ${restart}  ##Mark duplicates takes a long time, this modularizes it a bit in the event something happens with the model but the mark duplicates bam was created just fine. 
then
	gatk MarkDuplicatesSpark -I ${mbam} -O $markdupssam -M markdup_metrics.txt --create-output-bam-index true --remove-all-duplicates true --tmp-dir /lscratch/${SLURM_JOB_ID} --spark-master local[12] 2> err.${SLURM_JOB_ID}.txt
else
	continue 
fi 

gatk BaseRecalibrator -R $ref -I $markdupssam -O bqsr_1.table --known-sites $knownindels --known-sites $knownsnps --known-sites $millsgold --tmp-dir /lscratch/${SLURM_JOB_ID} > log-bqsr-${SLURM_JOB_ID}.txt 2>&1
gatk ApplyBQSR -bqsr bqsr_1.table -I $markdupssam -O bqsr_1.bam --tmp-dir /lscratch/${SLURM_JOB_ID} > log-ApplyBQSR-${SLURM_JOB_ID}.txt

gatk BaseRecalibrator -R $ref -I bqsr_1.bam -O bqsr_2.table --known-sites $knownindels --known-sites $knownsnps --known-sites $millsgold --tmp-dir /lscratch/${SLURM_JOB_ID} > log-bqsr-2-${SLURM_JOB_ID}.txt 2>&1 

# Generating some QC plots so we can see how BQSR improves the quality scoring. 
gatk AnalyzeCovariates -before bqsr_1.table -after bqsr_2.table -plots AnalyzeCovariates.pdf 2>&1 

rm $markdupssam ##removing the intermediate mark duplicates bam 
set +x


