#!/bin/sh

##References
# 1. https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels

set -e 

Help()
{
	# Display help 
	echo "To automate the processes of calling variants from WGS that has undergone standard GATK pre-processing recommendations, alignment (currently GRCh38), and standard alignment cleanup (mark duplicates, BQSR)."
	echo ""
	echo "To run this pipeline use the following command:"
	echo "sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] variant_calling_GATK.sh -b [batchfile] "
    echo "-b    This is a textfile of the IDs (assuming you are following the prescribed directory structure, and you have working directories named with their IDs only), with one ID on each line."  
    echo "-o    This is a textfile of the original bam paths, assuming you have bams in locations other than the working directory and created this file for use in pre-process-pipe.sh. One path per line."
    echo "      The goal of this was to just reduce the need to create intermediate files and for continuity with the same batch file." 

}

while getopts "b:o:h" option; do
   case $option in
	b) batchfile=$OPTARG ;;
    o) origbams=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

module load GATK/4.6.0.0 || fail "Could not load GATK module"
module load bcftools/1.19 || fail "Could not load bcftools module"
module load picard/3.2.0 || fail "Could not load picard module"

rundate=`date +'%m%d%y'`  ## for use in various file and log names 
chrlist=($(seq 1 1 22) "X" "Y") 
ref=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta
dbsnp=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.dbsnp138.vcf.gz

scriptpth="$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')"
sourcedir="$(dirname $scriptpth)"

if [ -z $batchfile ]
then
	batch=""
	while read path
	do
		prefix=`basename $path | cut -f1 -d'.'`
		batch+="${prefix} "
	done < ${origbams}
else
	batch=`cat ${batchfile}`
fi 

## FUNCTIONS -----

## Haplotype caller function generates a swarm command for each sample passed
generate_HC_swarm() { ## pass the input bam and ID number as positional arguments 

    inbam=$1 
    id=$2 

    echo "#SWARM -t 8 -g 32 --time=12:00:00" > HC-${rundate}-${id}.swarm ## Initating a swarm file 

    for value in "${chrlist[@]}" ## Looping through a list of chromosome names to generate lines for a swarm script that will run haplotype caller on each chromosome in parallel.
    do 
        chrnum="chr${value}"

        echo "TWD=/lscratch/\${SLURM_JOB_ID}; \
                gatk --java-options "-Xmx4g" HaplotypeCaller  \
                -R ${ref} \
                -I ${inbam} \
                -L ${chrnum} \
                -O ${PWD}/${id}/${id}_out/gvcfs_${rundate}/${chrnum}.g.vcf.gz \
                -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation --pcr-indel-model NONE -ERC GVCF \
                --native-pair-hmm-threads 4 \
                --tmp-dir /lscratch/\${SLURM_JOB_ID} 2> ${PWD}/${id}/${id}_out/logs_${rundate}/log-HC-${chrnum}-\${SLURM_JOB_ID}.txt" >> HC-${rundate}-${id}.swarm
    done

}

## This function creates a swarm to combine the GVCFs for each chromosome across the samples in our batch 
generate_combineGVCF_swarm() {

    echo "#SWARM -t 8 -g 32 --time=12:00:00" > ${PWD}/combinedGVCFs-${rundate}.swarm

    mkdir -p ${PWD}/combinedGVCFs_${rundate}
    mkdir -p ${PWD}/conbinedGVCFs_logs

    for value in "${chrlist[@]}"
    do 
        chrnum="chr${value}"
        
        gvcflist="" ## Creating a list of the paths to GVCFs to be combined 
        for id in ${batch}
        do 
            chrgvcf="${PWD}/${id}/${id}_out/gvcfs_${rundate}/${chrnum}.g.vcf.gz"
            gvcflist+="-V ${chrgvcf} " 
        done
        ## using the batchfile of sample IDs to combine chromosome gvcfs across samples 



        ## one swarm line per chromosome 
        echo "TWD=/lscratch/\${SLURM_JOB_ID}; \
                gatk --java-options "-Xmx4g" CombineGVCFs \
                -R $ref ${gvcflist} \
                -G StandardAnnotation -G AS_StandardAnnotation \
                -L ${chrnum} \
                -O ${PWD}/combinedGVCFs_${rundate}/${chrnum}.combined.g.vcf.gz \
                --tmp-dir /lscratch/\${SLURM_JOB_ID} 2> ${PWD}/combinedGVCFs_logs/${chrnum}-${rundate}-\${SLURM_JOB_ID}.log" >> ${PWD}/combinedGVCFs-${rundate}.swarm 
    done 
}

generate_genotypeGVCFs_swarm() {  ##generating another swarm to genotype the batch chromosomes in parallel

    echo "#SWARM -t 8 -g 32 --time=12:00:00" > ${PWD}/genotypeGVCFs-${rundate}.swarm

    mkdir -p ${PWD}/genotypedVCFs_${rundate}
    mkdir -p ${PWD}/genotypeGVCFs_logs

    for value in "${chrlist[@]}" 
    do 
        chrnum="chr${value}"
        chrgvcf="${PWD}/combinedGVCFs_${rundate}/chr${value}.combined.g.vcf.gz" 

        echo "TWD=/lscratch/\${SLURM_JOB_ID}; \
            gatk --java-options '-Xmx4g' GenotypeGVCFs \
            --tmp-dir /lscratch/\${SLURM_JOB_ID} \
            -R $ref \
            -V $chrgvcf \
            -G StandardAnnotation \
            -G AS_StandardAnnotation \
            --dbsnp $dbsnp \
            -L ${chrnum} \
            -O ${PWD}/genotypedVCFs_${rundate}/${chrnum}.genotypevcf.${rundate}.vcf.gz \
            2> ${PWD}/genotypeGVCFs_logs/${chrnum}-${rundate}-${SLURM_JOB_ID}.log" >> ${PWD}/genotypeGVCFs-${rundate}.swarm
    done 

}

# HAPLOTYPE CALLER ------
# Running haplotype caller on each chromosome of each sample
# A separate swarm script is generated for each sample by this loop

set -x 

jids=""
for id in ${batch}
do 
    bam="${PWD}/${id}/${id}_out/bqsr_1.bam"

    mkdir -p ${PWD}/${id}/${id}_out/logs_${rundate} ##Creating the necessary subdirectories for each sample 
    mkdir -p ${PWD}/${id}/${id}_out/gvcfs_${rundate} 

    generate_HC_swarm ${bam} ${id}

    jid=$(swarm --module GATK --gres=lscratch:200 -g 32 -t 8 --logdir ${PWD}/${id}/${id}_out/logs_${rundate} ${PWD}/HC-${rundate}-${id}.swarm)
    echo $jid
    jids+=`echo "$jid "`
done 

jidlist=`echo $jids | sed 's/ /:/g'` ## Creating a list of jobids to use as dependencies for the next job, so that combineGVCFs does not start running until all HC jobs are complete 
printf "Combining GVCFs dependent on following jobs completing: $jids\n" 

## COMBINING GVCFs ------
## GVCFs combined by chromosome, i.e./ we combined the GVCFs for chrN from each sample in the batch, resulting in 24 GVCFs representing all samples 
## One swarm job per chromosome

generate_combineGVCF_swarm
combine_jid=$(swarm --module GATK --dependency afterok:$jidlist --gres=lscratch:200 -g 32 -t 8 --logdir ${PWD}/combinedGVCFs_logs/ ${PWD}/combinedGVCFs-${rundate}.swarm) ## a second jid to hold genotyping until all combineGVCF jobs are done 
echo $combine_jid
printf "Haplotype caller swarms completed, combining GVCFs across chromosomes. Genotyping will began upon completion of $combine_jid\n"

## GENOTYPING GVCFs ------
## Genotyping the chromosome-combined GVCFs 
## One swarm job per chromosome 

generate_genotypeGVCFs_swarm
genotype_jid=$(swarm --module GATK --dependency afterok:$combine_jid --gres=lscratch:200 -g 32 -t 8 --logdir ${PWD}/genotypeGVCFs_logs ${PWD}/genotypeGVCFs-${rundate}.swarm)
echo $genotypejid

sbatch --mem=96g --cpus-per-task=24 --gres=lscratch:400 --time=1-06:00:00 --dependency=afterok:$genotype_jid ${sourcedir}/run-VQSR.sh -v ${PWD}/genotypedVCFs_${rundate}

set +x 



