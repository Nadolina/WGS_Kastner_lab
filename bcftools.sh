#!/bin/sh

## attempting to run bcftools mpileup, parallelized by chromosome, across our five samples

set -e

Help()
{
	# Display help 
	echo "Running bcftools mpileup across chromosomes across our test batch."
	echo ""
	echo "To run this pipeline use the following command:"
	echo "sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] bcftools.sh -b [batch files with sample IDs]"

}

while getopts "b:h" option; do
   case $option in
	b) batchfile=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

module load bcftools/1.19
module load samtools/1.19

rundate=`date +'%m%d%y'`

ref=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta
thousandgAF=/data/Kastner_PFS/references/1000genomes/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.AFs.tab.gz
thousandgHDR=/data/Kastner_PFS/references/1000genomes/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.AFs.hdr

echo "#SWARM -g 32 -t 16 --time=12:00:00" > mpileup-${rundate}.swarm
echo "#SWARM -g 32 -t 8 --time=06:00:00" > annotate-${rundate}.swarm
echo "#SWARM -g 32 -t 8 --time=06:00:00" > call-${rundate}.swarm
mkdir -p ${PWD}/bcftools_logs
mkdir -p ${PWD}/bcftools_out

rm ${PWD}/bams-bcftools-call-${rundate}.txt && touch ${PWD}/bams-bcftools-call-${rundate}.txt
while IFS="" read -r id || [ -n "$id" ]
do 
    printf "${PWD}/${id}/${id}_out/bqsr_1.bam\n" >> ${PWD}/bams-bcftools-call-${rundate}.txt
	if [ ! -f ${PWD}/${id}/${id}_out/bqsr_1.csi ]; then
		samtools index -c -@4 -o ${PWD}/${id}/${id}_out/bqsr_1.csi ${PWD}/${id}/${id}_out/bqsr_1.bam
	fi
done < ${batchfile}

chrlist=($(seq 1 1 22) "X" "Y")

for value in "${chrlist[@]}" 
do 
    chrnum="chr${value}"

    printf "bcftools mpileup -b ${PWD}/bams-bcftools-call-${rundate}.txt\
		-f ${ref} \
		-r ${chrnum} \
		-a AD,DP,SP \
		-o bcftools_out/mpileup-${chrnum}-${rundate}.vcf.gz \
		-Oz --threads 8 \
		2> ${PWD}/bcftools_logs/mpileup-${chrnum}-${rundate}-\${SLURM_JOB_ID}.log; bcftools index -t --threads 8 -o bcftools_out/mpileup-${chrnum}-${rundate}.vcf.gz.tbi bcftools_out/mpileup-${chrnum}-${rundate}.vcf.gz\n" >> mpileup-${rundate}.swarm
	printf "bcftools annotate -a $thousandgAF \
		-h $thousandgHDR \
		-c CHROM,POS,REF,ALT,REF_AN,REF_AC \
		-r ${chrnum} \
		-o bcftools_out/annotate-${chrnum}-${rundate}.vcf.gz \
		-Oz --threads 8 \
		bcftools_out/mpileup-${chrnum}-${rundate}.vcf.gz \
		2> ${PWD}/bcftools_logs/annotate-${chrnum}-${rundate}-\${SLURM_JOB_ID}.log\n" >> annotate-${rundate}.swarm
	printf "bcftools call -o bcftools_out/call-${chrnum}-${rundate}.vcf.gz \
		-mv -a GQ,GP \
		--prior-freqs REF_AN,REF_AC \
		-Oz --threads 8 \
		-v bcftools_out/annotate-${chrnum}-${rundate}.vcf.gz \
		2> ${PWD}/bcftools_logs/bcftools-${chrnum}-${rundate}-\${SLURM_JOB_ID}.log\n" >> call-${rundate}.swarm

done 

mpileup_jid=$(swarm --module bcftools --gres=lscratch:200 -g 32 -t 16 --logdir ${PWD}/bcftools_logs mpileup-${rundate}.swarm)
echo "BCFtools annotate will run after the mpileup swarm jobs with the JOBID ${mpileup_jid} are complete."

annotate_jid=$(swarm --module bcftools --dependency afterok:$mpileup_jid --gres=lscratch:100 -g 32 -t 8 --logdir ${PWD}/bcftools_logs annotate-${rundate}.swarm)
echo "BCFtools call will run after the annotate swarm jobs with the JOBID ${annotate_jid} are complete."

swarm --module bcftools --dependency afterok:$annotate_jid --gres=lscratch:100 -g 32 -t 8 --logdir ${PWD}/bcftools_logs call-${rundate}.swarm