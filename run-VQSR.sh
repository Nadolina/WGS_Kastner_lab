#!/bin/sh

## References:
##    https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
##    https://hpc.nih.gov/training/gatk_tutorial/vqsr.html#optimized-script-for-variantrecalibrator

set -e

Help()
{
	# Display help 
	echo "Running variant score recalibration on the combined VCF."
	echo ""
	echo "To run this pipeline use the following command:"
	echo "sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] run-VQSR.sh -v [directory of genotyped VCFs]"

}

while getopts "v:h" option; do
   case $option in
	v) vcfdir=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

module load GATK/4.6.0.0
module load bcftools/1.19
module load picard/3.2.0

rundate=`date +'%m%d%y'`

ref=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta

mkdir -p ${PWD}/VQSR_logs_${rundate}

## We at this point have one VCF representing each chromosome across our batch. We need to combine them all before VQSR. 
## In this loop I am creating the portion of the gatherVCFs command line defining all the VCF locations 
vcflist=`ls ${vcfdir} | grep '.vcf.gz$' | sort -V`  
gatherinput=""
for vcf in ${vcflist};
do
    gatherinput+="I=${vcfdir}/${vcf} " 
done 

echo "java -jar $PICARDJAR GatherVcfs ${gatherinput} 0=${PWD}/merged-$rundate.vcf.gz" 2> ${PWD}/VQSR_logs_${rundate}/gathervcfs-${rundate}-$SLURM_JOB_ID.log
java -jar $PICARDJAR GatherVcfs ${gatherinput} O=${PWD}/merged-$rundate.vcf.gz 2> ${PWD}/VQSR_logs_${rundate}/gathervcfs-${rundate}-${SLURM_JOB_ID}.log

set -x
bcftools index -t merged-${rundate}.vcf.gz 
bcftools stats merged-${rundate}.vcf.gz > merged-${rundate}.vcf.stats

mkdir -p filtered_VCFs_${rundate}

## Running a round of VQSR on the SNPs and indels separately, as recommended by Biowulf and GATK documentation. 
gatk --java-options "-XX:ParallelGCThreads=2 -Xmx8G" VariantRecalibrator \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 \
  -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
  -tranche 95.0 -tranche 94.0 \
  -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  -AS \
  -R ${ref} \
  -V merged-$rundate.vcf.gz \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
  /fdb/GATK_resource_bundle/hg38/hapmap_3.3.hg38.vcf.gz  \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 \
  /fdb/GATK_resource_bundle/hg38/1000G_omni2.5.hg38.vcf.gz \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 \
  /fdb/GATK_resource_bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
  /fdb/GATK_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP\
  -mode SNP -O filtered_VCFs_${rundate}/merged_SNP1.recal --tranches-file filtered_VCFs_${rundate}/output_SNP1.tranches \
  --rscript-file filtered_VCFs_${rundate}/output_SNP1.plots.R \
  --tmp-dir /lscratch/${SLURM_JOB_ID} 2> ${PWD}/VQSR_logs_${rundate}/slurm-vqsrsnps-${SLURM_JOB_ID}.log 

gatk --java-options "-Xmx8G -XX:ParallelGCThreads=2" VariantRecalibrator \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 \
  -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
  -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
  -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  -R ${ref} \
  -V merged-$rundate.vcf.gz \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 \
  /fdb/GATK_resource_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
  /fdb/GATK_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
  -mode INDEL -O filtered_VCFs_${rundate}/merged_indel1.recal --tranches-file filtered_VCFs_${rundate}/output_indel1.tranches \
  --rscript-file filtered_VCFs_${rundate}/output_indel1.plots.R \
  --tmp-dir /lscratch/${SLURM_JOB_ID} 2> ${PWD}/VQSR_logs_${rundate}/slurm-vqsr-indels-${SLURM_JOB_ID}.log 

## Running apply VQSR with the SNP recal file, then passing the vcf generated from this command to the applyVQSR for indels, so we have a combined snp-indels output vcf. 
gatk --java-options "-Xmx8G -XX:ParallelGCThreads=2" ApplyVQSR \
  -V merged-$rundate.vcf.gz \
  --recal-file filtered_VCFs_${rundate}/merged_SNP1.recal \
  -mode SNP \
  -AS \
  --tranches-file filtered_VCFs_${rundate}/output_SNP1.tranches \
  --truth-sensitivity-filter-level 99.9 \
  --create-output-variant-index true \
  -O filtered_VCFs_${rundate}/SNP.recalibrated_99.9.${rundate}.vcf.gz 2> ${PWD}/VQSR_logs_${rundate}/slurm-applyvqsr-snps-${SLURM_JOB_ID}.log

gatk --java-options "-Xmx8g -XX:ParallelGCThreads=2" ApplyVQSR \
  -V filtered_VCFs_${rundate}/SNP.recalibrated_99.9.${rundate}.vcf.gz \
  -mode INDEL \
  --recal-file filtered_VCFs_${rundate}/merged_indel1.recal \
  --tranches-file filtered_VCFs_${rundate}/output_indel1.tranches \
  --truth-sensitivity-filter-level 99.9 \
  --create-output-variant-index true \
  -O filtered_VCFs_${rundate}/indel.SNP.recalibrated_99.9.${rundate}.vcf.gz 2> ${PWD}/VQSR_logs_${rundate}/slurm-applyvqsr-indels-${SLURM_JOB_ID}.log
 
set +x 
