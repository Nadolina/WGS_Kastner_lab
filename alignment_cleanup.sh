#!/bin/sh

# The purpose of this script is to take the merged bam alignment output from pre-process-pipe.sh and prepare it for variant calling. 
# References for this include: 
#   https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery

Help()
{
	# Display help 
	echo "To automate the processes cleaning up a bwa-mem2 alignment to HG38 prior to performing variant calling. This script assumes the bam has already undergone MergeBamAlignments."
	echo ""
	echo "Run cleanup via the submissions script submit-alignment-cleanup.sh."

}

restart='true'
while getopts "b:h" option; do
   case $option in
	b) bam=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

if [[ -z "${bam}" ]]; then
	echo "No bam provided, please provide one using -b"
	exit 1
else
	sample=$(basename ${bam} | sed 's/.bam//g')
	OUTDIR="$PWD/${sample}/${sample}_out"
	LOGS="$PWD/${sample}/logs"
	mbam="${OUTDIR}/${sample}.fullmerged.bam"

	echo "Merged alignment bam: ${mbam}"
	echo "Output directory: ${OUTDIR}"
fi

module load GATK/4.6.0.0 || fail "Could not load GATK"
module load samtools/1.19 || fail "Could not load samtools"

if [[ ! -d /lscratch/${SLURM_JOBID} ]]; then mkdir -p /lscratch/${SLURM_JOBID}; fi
rundate=$(date +%Y%m%d)

##The sets of known variants are required for the BQSR modeling. 
knownindels=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.known_indels.vcf.gz
knownsnps=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.dbsnp138.vcf.gz
millsgold=/data/Kastner_PFS/references/HG38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta

set -x

mkdir -p /lscratch/${SLURM_JOB_ID}
gatk MarkDuplicatesSpark --java-options "-Xms60G -Xmx60G" \
	-I ${mbam} \
	-O ${OUTDIR}/${sample}.marked_dups.bam \
	-M ${OUTDIR}/markdup_metrics.txt \
	--create-output-bam-index true \
	--remove-all-duplicates true \
	--tmp-dir /lscratch/${SLURM_JOB_ID} \
	--spark-master local[12] 2> ${LOGS}/markdups_${rundate}.err

>${OUTDIR}/${sample}.alt_contigs.list
samtools view -H ${mbam} |\
	grep -e chrM -e alt -e random -e decoy -e HLA -e Un -e EBV |\
	awk '{print $2}' |\
	sed 's/SN://g' >> ${OUTDIR}/${sample}.alt_contigs.list

echo "#SWARM -g 16 -t 4 --gres=lscratch:100 --time=4:00:00 --module GATK/4.6.0.0 --logdir ${LOGS}" > ${OUTDIR}/bqsr_1.${rundate}.swarm
echo "#SWARM -g 6 --gres=lscratch:50 --module GATK/4.6.0.0 --logdir ${LOGS}" > ${OUTDIR}/applyBQSR.${rundate}.swarm
echo "#SWARM -g 16 --gres=lscratch:100 --time=4:00:00 --module GATK/4.6.0.0 --logdir ${LOGS}" > ${OUTDIR}/bqsr_2.${rundate}.swarm

## For looping through chromosomes for parallel BQSR.
chrs=$(echo $(seq 1 1 22) X Y)

for chr in ${chrs}; do 

	## Generating a swarm file to run BQSR in parallel per chromosome. 
	echo "gatk BaseRecalibrator --java-options '-Xms12G -Xmx12G' \
		-R ${ref} \
		-I ${OUTDIR}/${sample}.marked_dups.bam \
		-O ${OUTDIR}/chr${chr}.${sample}_bqsr_1.table \
		--intervals chr${chr} \
		--known-sites ${knownindels} \
		--known-sites ${knownsnps} \
		--known-sites ${millsgold} \
		--tmp-dir /lscratch/\${SLURM_JOB_ID}" >> ${OUTDIR}/bqsr_1.${rundate}.swarm

	echo "gatk ApplyBQSR --java-options '-Xms4G -Xmx4G' \
		-bqsr ${OUTDIR}/${sample}_combined_bqsr_1.table \
		-I ${OUTDIR}/${sample}.marked_dups.bam \
		-O ${OUTDIR}/chr${chr}.${sample}_bqsr_1.bam \
		--intervals chr${chr} \
		--tmp-dir /lscratch/\${SLURM_JOB_ID}" >> ${OUTDIR}/applyBQSR.${rundate}.swarm

	echo "gatk BaseRecalibrator --java-options '-Xms12G -Xmx12G' \
		-R $ref \
		-I ${OUTDIR}/chr${chr}.${sample}_bqsr_1.bam \
		-O ${OUTDIR}/chr${chr}.${sample}_bqsr_2.table \
		-L "chr${chr}" \
		--known-sites $knownindels \
		--known-sites $knownsnps \
		--known-sites $millsgold \
		--tmp-dir /lscratch/\${SLURM_JOB_ID} ;\
		gatk AnalyzeCovariates --java-options '-Xmx4G' \
		-before ${OUTDIR}/chr${chr}.${sample}_bqsr_1.table \
		-after ${OUTDIR}/chr${chr}.${sample}_bqsr_2.table \
		-plots ${OUTDIR}/chr${chr}_AnalyzeCovariates.pdf" >> ${OUTDIR}/bqsr_2.${rundate}.swarm

done

echo "gatk BaseRecalibrator --java-options '-Xms12G -Xmx12G' \
		-R ${ref} \
		-I ${OUTDIR}/${sample}.marked_dups.bam \
		-O ${OUTDIR}/${sample}.alt_contigs.bqsr_1.table \
		-L ${OUTDIR}/${sample}.alt_contigs.list \
		--known-sites ${knownindels} \
		--known-sites ${knownsnps} \
		--known-sites ${millsgold} \
		--tmp-dir /lscratch/\${SLURM_JOB_ID}" >> ${OUTDIR}/bqsr_1.${rundate}.swarm

echo "gatk ApplyBQSR --java-options '-Xms4G -Xmx4G' \
		-bqsr ${OUTDIR}/${sample}_combined_bqsr_1.table \
		-I ${OUTDIR}/${sample}.marked_dups.bam \
		-L ${OUTDIR}/${sample}.alt_contigs.list \
		-O ${OUTDIR}/${sample}.alt_contigs.bqsr_1.bam \
		--tmp-dir /lscratch/\${SLURM_JOB_ID}" >> ${OUTDIR}/applyBQSR.${rundate}.swarm

## Initating swarm jobs with dependencies -------

bqsr_1_jid=$(swarm ${OUTDIR}/bqsr_1.${rundate}.swarm)

## Generating combined BQSR table from contig-scattered tables for use in ApplyBQSR step.
BQSRTABLES=""
for chr in ${chrs}; do 
	BQSRTABLES+="-I ${OUTDIR}/chr${chr}.${sample}_bqsr_1.table "
done
BQSRTABLES+="-I ${OUTDIR}/${sample}.alt_contigs.bqsr_1.table"
gather_jid=$(sbatch --mem=5g --dependency afterok:${bqsr_1_jid} \
	--wrap="module load GATK/4.6.0.0 ; gatk GatherBQSRReports ${BQSRTABLES} \
		-O ${OUTDIR}/${sample}_combined_bqsr_1.table \
		2> ${LOGS}/gatherBQSRReports_${rundate}.err")

applyBQSR_jid=$(swarm --dependency afterok:${gather_jid} ${OUTDIR}/applyBQSR.${rundate}.swarm)
bqsr_2_jid=$(swarm --dependency afterok:${applyBQSR_jid} ${OUTDIR}/bqsr_2.${rundate}.swarm)

BAMs=""
for chr in ${chrs}; do 
	BAMs+="${OUTDIR}/chr${chr}.${sample}_bqsr_1.bam "
done
BAMs+="${OUTDIR}/${sample}.alt_contigs.bqsr_1.bam"

echo "Merging chrosome-scattered BQSR BAMs into a single BAM for downstream variant calling."
echo "samtools merge -f -@4 -o ${OUTDIR}/${sample}_bqsr_1.bam --write-index ${BAMs}"
merge_jid=$(sbatch --mem=16g -c 4 --gres=lscratch:100 --time=4:00:00 --dependency afterok:${bqsr_2_jid} \
	--wrap="module load samtools ; samtools merge -f -@4 -o ${OUTDIR}/${sample}_bqsr_1.bam --write-index ${BAMs}")

## Removing intermediaries -------

# intermediates1=$(ls ${OUTDIR}/*.mergedaln.bam ${OUTDIR}/*mergedaln.bai ${OUTDIR}/*.marked_dups.bam)
# intermediates2=$(for CHR in $(echo $(seq 1 1 22) X Y); do ls ${OUTDIR}/chr${CHR}.${sample}_bqsr_1.ba* ; done)

intermediates=""
for f in ${OUTDIR}/*.mergedaln.ba* ; do 
	intermediates+="${f} "
done
intermediates+="${OUTDIR}/*.marked_dups.ba* "
for CHR in $(echo $(seq 1 1 22) X Y); do 
	intermediates+="${OUTDIR}/chr${CHR}.${sample}_bqsr_1.ba* "
done

sbatch --dependency afterok:${merge_jid} \
	--wrap="echo 'Final BQSR BAM merged finished successfully; deleting intermediate files' && rm -f ${intermediates}"

set +x


