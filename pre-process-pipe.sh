#!/bin/sh

# The purpose of this script is to revert aligned BAMs (i.e./ by NISC or elsewhere) to unaligned BAMs, such that we can customize alignment and variant calling for the Kastner lab. 
# References for this include:
#	https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM
#	https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
#	https://s3.amazonaws.com/broad-references/broad-references-readme.html

set -e

Help()
{
	# Display help 
	echo "Submit pre-process-pipe.sh to SLURM using submit-pre-process-pipe.sh"

}

while getopts "b:h" option; do
   case $option in
	b) alnbam=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

module load picard/3.2.0 || fail "Could not load picard module"
module load bwa-mem2/2.2.1 || fail "Could not load bwa-mem2"
module load fastqc/0.12.1 || fail "Could not load fastqc"
module load samtools/1.19 || fail "Could not load samtools"

orig_dir=$(dirname ${alnbam})
printf "The source directory of the bam is $orig_dir.\n"

prefix=$(basename $alnbam | cut -f1 -d'.')
printf "Creating sub-directory $prefix in $PWD for outputs.\n"

ref=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta ##the bwa mem2 index files are present in this directory as well 

## Create output directories
mkdir -p $PWD/${prefix}/${prefix}_out
mkdir -p $PWD/${prefix}/qc_out
mkdir -p $PWD/${prefix}/logs
OUTDIR=$(echo $PWD/${prefix}/${prefix}_out)
LOGS=$(echo ${PWD}/${prefix}/logs)

## Confirming index exists for bam and creating otherwise 
if [ -f "${alnbam}.bai" ]; then 	
	echo "index found: ${alnbam}.bai"
else 
	echo "creating index"
	samtools index ${alnbam}
fi

if [[ ! -d /lscratch/${SLURM_JOB_ID} ]]; then mkdir -p /lscratch/${SLURM_JOB_ID}; fi

## Generating list of IDs from read group lines in a bam 
RGlist=$(samtools view -H ${alnbam} | grep '^@RG' | awk -F'\t' '{for (i=1; i<=NF; i++) if ($i ~ /^ID:/) print substr($i, 4)}')


## FUNCTIONS -----

revert() { ## Converts the original aligned bam to an unaligned bam

	##ATTRIBUTES_TO_CLEAR are annotations in an aligned bam that are not compatible with an unaligned bam/future alignment and must be removed 
	##Outputting a per-read-group bam to enable parallelization of the next steps.
	java -Xmx16G -jar $PICARDJAR RevertSam \
     		I=${alnbam} \
			OUTPUT_BY_READGROUP=true \
			OUTPUT=${OUTDIR} \
     		ATTRIBUTE_TO_CLEAR=BD \
			ATTRIBUTE_TO_CLEAR=BI \
     		ATTRIBUTE_TO_CLEAR=XS \
     		ATTRIBUTE_TO_CLEAR=XT \
     		ATTRIBUTE_TO_CLEAR=XN \
     		ATTRIBUTE_TO_CLEAR=OC \
     		ATTRIBUTE_TO_CLEAR=OP \
     		ATTRIBUTE_TO_CLEAR=sd \
     		TMP_DIR=/lscratch/${SLURM_JOB_ID} \
     		> ${LOGS}/log-RevertSam-${prefix}-${SLURM_JOB_ID}.txt 2>&1

}

mkAdapters_aln() { ## Piping MarkIlluminaAdapters through SamToFastq and BWA-MEM2 in one swarm file.

	echo "#SWARM -g 42 -t 16 --gres=lscratch:300 --module GATK/4.6.0.0 --module bwa-mem2 --time=16:00:00 --logdir ${LOGS}" > ${OUTDIR}/swarm-mkAdapter-aln.txt

	## Looping through each read group ID to create the swarm commands for adapter marking and alignment
	for ID in ${RGlist}; do
	echo "java -Xmx20G -jar $PICARDJAR MarkIlluminaAdapters \
				I=$OUTDIR/${ID}.bam \
				O=/dev/stdout \
				M=${prefix}/qc_out/${ID}.mkAdapter.metrics.txt \
				TMP_DIR=/lscratch/\${SLURM_JOB_ID} \
				2> ${LOGS}/log-mkadapters-${ID}-${SLURM_JOB_ID}.txt |\
			java -Xmx20G -jar $PICARDJAR SamToFastq \
				I=/dev/stdin \
				FASTQ=/dev/stdout \
				CLIPPING_ACTION=2 CLIPPING_ATTRIBUTE=XT INTERLEAVE=true NON_PF=true TMP_DIR=/lscratch/\${SLURM_JOB_ID} 2> ${LOGS}/log-samtofastq-${ID}-\${SLURM_JOB_ID}.txt | \
			bwa-mem2 mem -t $SLURM_CPUS_PER_TASK -M -p $ref -o ${OUTDIR}/${ID}.aln.bam /dev/stdin 2> ${LOGS}/log-bwamem2-${ID}-\${SLURM_JOB_ID}.txt" >> ${OUTDIR}/swarm-mkAdapter-aln.txt

	done

}


## PRE-PROCESSING AND ALIGNMENT ------ 

set -x ##Using set -x and +x to print the function calls. 
revert
rbams=""
for ID in $(RGlist); do 
	rbams+="${OUTDIR}/${ID}.bam "
done ## Listing reverted per-read-group bams. 

mkAdapters_aln
mkAdapters_jid=$(swarm ${OUTDIR}/swarm-mkAdapter-aln.txt)

echo "Creating swarm file for sorting aligned and reverted bams..."
echo "Creating swarm file for merging BAM alignments..."
echo "#SWARM -g 10 --gres=lscratch:200 --module GATK/4.6.0.0 --time=8:00:00 --logdir ${LOGS}" > ${OUTDIR}/sort-${prefix}-${SLURM_JOB_ID}.swarm
echo "#SWARM -g 24 -t 4 --gres=lscratch:300 --module GATK/4.6.0.0 --time=12:00:00 --logdir ${LOGS}" > ${OUTDIR}/mergebamaln-${prefix}-${SLURM_JOB_ID}.swarm

for ID in ${RGlist}; do

    echo "Adding ${ID}.aln.bam and ${ID}.bam to sort swarm..."
    echo "gatk --java-options '-Xmx8G' SortSam \
            -I ${OUTDIR}/${ID}.aln.bam \
            -O ${OUTDIR}/${ID}.aln.sort.bam \
            --SORT_ORDER queryname \
            --TMP_DIR /lscratch/\${SLURM_JOB_ID}" \
			2> ${LOGS}/log-sort-alnbam-${ID}-\${SLURM_JOB_ID}.txt >> ${OUTDIR}/sort-${prefix}-${SLURM_JOB_ID}.swarm
    echo "gatk --java-options '-Xmx8G' SortSam \
            -I ${OUTDIR}/${ID}.bam \
            -O ${OUTDIR}/${ID}.sort.bam \
            --SORT_ORDER queryname \
            --TMP_DIR /lscratch/\${SLURM_JOB_ID}" \
			2> ${LOGS}/log-sort-bam-${ID}-\${SLURM_JOB_ID}.txt >> ${OUTDIR}/sort-${prefix}-${SLURM_JOB_ID}.swarm

	echo "gatk --java-options '-Xmx20G' MergeBamAlignment \
		R=$ref \
		UNMAPPED_BAM=${OUTDIR}/${ID}.bam \
		O=${OUTDIR}/${ID}.mergedaln.bam \
		ALIGNED_BAM=${OUTDIR}/${ID}.aln.bam \
		CREATE_INDEX=true \
		ADD_MATE_CIGAR=true \
		CLIP_ADAPTERS=false \
		CLIP_OVERLAPPING_READS=true \
		INCLUDE_SECONDARY_ALIGNMENTS=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		ATTRIBUTES_TO_RETAIN=XS \
		ATTRIBUTES_TO_RETAIN=XA \
		ATTRIBUTES_TO_RETAIN=pa \
		TMP_DIR=/lscratch/\${SLURM_JOB_ID} \
		2> ${LOGS}/log-mergebaln-${ID}-\${SLURM_JOB_ID}.txt" >> ${OUTDIR}/mergebamaln-${prefix}-${SLURM_JOB_ID}.swarm
done

sort_jid=$(swarm --dependency afterok:${mkAdapters_jid} ${OUTDIR}/sort-${prefix}-${SLURM_JOB_ID}.swarm)
merge_jid=$(swarm --dependency afterok:${sort_jid} ${OUTDIR}/mergebamaln-${prefix}-${SLURM_JOB_ID}.swarm)

revertbams=""
alnbams=""
sortbams=""
sortalnbams=""
mbams=""
for ID in ${RGlist}; do 
	mbams+="${OUTDIR}/${ID}.mergedaln.bam "
	revertbams+="${OUTDIR}/${ID}.bam "
	alnbams+="${OUTDIR}/${ID}.aln.bam "
	sortbams+="${OUTDIR}/${ID}.sort.bam "
	sortalnbams+="${OUTDIR}/${ID}.aln.sort.bam "
done ## Listing reverted per-read-group bams. 

fullmerged_jid=$(sbatch --mem=16g \
	-c 4 \
	--gres=lscratch:300 \
	--time=4:00:00 \
	--dependency afterok:${merge_jid} \
	--wrap="module load samtools ; samtools merge -f -@4 -o ${OUTDIR}/${prefix}.fullmerged.bam --write-index ${mbams}")

## removing intermediate files if the samtools merge is successful
sbatch --dependency afterok:${fullmerged_jid} --wrap="rm ${revertbams} ${sortbams} ${sortalnbams} ${alnbams}"


set +x
