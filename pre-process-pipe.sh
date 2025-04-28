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
	echo "To automate the processes of converting BAMs to uBAMs, marking adapters and aligning to the reference genome (currently GRCh38)."
	echo "For WGS, I recommend providing at least --mem=48g and --cpus-per-task=16 and --time=1-00:00:00, assuming a bam of around 60GB, scale appropriately."
	echo ""
	echo "To run this pipeline use the following command:"
	echo "sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] pre-process-pipe.sh -b [original bam]"
	echo "OR"
	echo "sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] pre-process-pipe.sh -l [location]"

	echo "You need to pass EITHER -b or -l, but NOT BOTH."
	echo "The -l is retained as a legacy option, and was used for looping through original BAMs that had been copied to the working directory but had variable naming schemes."
	echo "The -l [locations] option will look in the folder for anything that matches *.bam, so ensure the bam of interest is the only *.bam in the folder provided."
	echo "The -b is preferable because you can point to an original BAM anywhere on biowulf, i.e./ Kastner_PFS, (pending permissions) without having to copy large BAMs."
	echo ""
	echo "You can pass a location or BAM directly to -l and -b, but you can also loop through a file of either, with one location/path per line."
	echo "while read path; do sbatch [OPTIONS] pre-process-pipe.sh -b $\{path} ; done < batch.txt"
	echo ""

}

while getopts "b:l:h" option; do
   case $option in
	b) alnbam=$OPTARG ;;
	l) location=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

if [ -z "$location" ]
then
	printf "no location provided, bam was provided: $alnbam\n"
else
	cd $location 
	alnbam=`ls *bam`
	printf "The bam is $alnbam\n"
fi 

module load picard/3.2.0 || fail "Could not load picard module"
module load bwa-mem2/2.2.1 || fail "Could not load bwa-mem2"
module load fastqc/0.12.1 || fail "Could not load fastqc"
module load samtools/1.19 || fail "Could not load samtools"

orig_dir=`dirname $alnbam`
printf "The source directory of the bam is $orig_dir.\n"

prefix=`basename $alnbam | cut -f1 -d'.'`
printf "Creating sub-directory $prefix in $PWD for outputs.\n"

ref=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta ##the bwa mem2 index files are present in this directory as well 

mkdir -p $PWD/${prefix}/${prefix}_out
OUTDIR=`echo $PWD/${prefix}/${prefix}_out`

mkdir -p $PWD/${prefix}/qc_out

index=`ls ${orig_dir}/${prefix}*.bai`
if [ -f $index ]; then 	
	echo "index found: $index"
else 
	echo "creating index"
	samtools index $alnbam
fi

if [[ ! -d /lscratch/${SLURM_JOB_ID} ]]; then mkdir -p /lscratch/${SLURM_JOB_ID}; fi

## FUNCTIONS -----

revert() { ## Converts the original aligned bam to an unaligned bam

	##ATTRIBUTES_TO_CLEAR are annotations in an aligned bam that are not compatible with an unaligned bam/future alignment and must be removed 

	java -Xmx16G -jar $PICARDJAR RevertSam \
     		I=$alnbam \
     		O=$OUTDIR/$prefix.reverted.bam \
     		SANITIZE=true \
     		ATTRIBUTE_TO_CLEAR=BD \
			ATTRIBUTE_TO_CLEAR=BI \
     		ATTRIBUTE_TO_CLEAR=XS \
     		ATTRIBUTE_TO_CLEAR=XT \
     		ATTRIBUTE_TO_CLEAR=XN \
     		ATTRIBUTE_TO_CLEAR=OC \
     		ATTRIBUTE_TO_CLEAR=OP \
     		ATTRIBUTE_TO_CLEAR=sd \
     		# ATTRIBUTE_TO_CLEAR=xq \
			# ATTRIBUTE_TO_CLEAR=XQ \
     		TMP_DIR=/lscratch/${SLURM_JOB_ID} \
     		> $OUTDIR/log-RevertSam-${prefix}-${SLURM_JOB_ID}.txt 2>&1

}

mkAdapters() { # Using MarkIlluminaAdapters to mark the illumina adapters in our bam file before performing alignment. 

	java -Xmx8G -jar $PICARDJAR MarkIlluminaAdapters \
     		I=$OUTDIR/${prefix}.reverted.bam \
     		O=$OUTDIR/${prefix}.mkAdapter.bam\
     		M=${prefix}/qc_out/${prefix}.mkAdapter.metrics.txt \
     		TMP_DIR=/lscratch/${SLURM_JOB_ID} \
    		 > $OUTDIR/log-mkAdapter-${prefix}-${SLURM_JOB_ID}.txt 2>&1
}

alignment() { ## Converting the marked adapters bam into a sam and piping that through BWA-MEM2 against HG38 and merging the unaligned reads to the end of the bam with mergebamalignments. 

	echo "Converting sam to an interleaved fastq and piping the output to bwa-mem2." 

	java -Xmx16G -jar $PICARDJAR SamToFastq \
                I=$OUTDIR/${prefix}.mkAdapter.bam \
                FASTQ=/dev/stdout \
                CLIPPING_ACTION=2 CLIPPING_ATTRIBUTE=XT INTERLEAVE=true NON_PF=true TMP_DIR=/lscratch/${SLURM_JOB_ID} 2> ${OUTDIR}/log-samtofastq-${SLURM_JOB_ID}.txt | \
	bwa-mem2 mem -t $SLURM_CPUS_PER_TASK -M -p $ref /dev/stdin 2> ${OUTDIR}/log-bwamem2-${SLURM_JOB_ID}.txt  | \
	java -Xmx16G -jar $PICARDJAR MergeBamAlignment R=$ref UNMAPPED_BAM=${OUTDIR}/${prefix}.reverted.bam ALIGNED_BAM=/dev/stdin O=$OUTDIR/${prefix}.mergedaln.bam \
		CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA ATTRIBUTES_TO_RETAIN=pa TMP_DIR=/lscratch/${SLURM_JOB_ID}

}

## PRE-PROCESSING AND ALIGNMENT ------ 

set -x ##Using set -x and +x to print the function calls. 
revert
fastqc -o ${prefix}/qc_out ${OUTDIR}/${prefix}.reverted.bam &
mkAdapters 
alignment
samtools stats -@${SLURM_CPUS_PER_TASK} ${OUTDIR}/${prefix}.mergedaln.bam > ${prefix}/qc_out/${prefix}.mergedaln.stats ## generating stats of the original and new alignments. 
samtools stats -@${SLURM_CPUS_PER_TASK} ${alnbam} > ${prefix}/qc_out/${prefix}.originalaln.stats
rm $OUTDIR/${prefix}.mkAdapter.bam $OUTDIR/${prefix}.reverted.bam 
set +x 
