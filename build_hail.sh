#!/bin/sh 

## For building a hail database from multiple batches after they are annotated with VEP. 

set -e

Help()
{
	# Display help 
    echo ""
	echo "This script is intended to prepare VEP-annotated VCFs and build a hail database from a number of batches (~100 WGS samples in total is ideal). Build_hail.sh expects the output as generated by run-VEP.sh, which is per-batch folders of uncompressed VCFs in the format chr1.filtered.vcf, chr2.filtered.vcf ... chrN.filtered.vcf. Assuming the inputs are coming from the Kastner pipeline, do not combine GATK and bcftools into the same hail database. Build one database for each GATK and bcftools,for a given group of batches.
    "

	echo "To run this pipeline use the following command:"
	echo "sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] build_hail.sh [folder names of batches]
    "

}

while getopts "h" option; do
   case $option in
	h) # display Help
         Help
         exit;;
   esac
done

module load bcftools/1.19

batches=()
printf "
The batches provided include:
"
for folder in "$@"; do 
    batches+="${folder} "
    echo "  $folder"
done

rundate=`date +'%m%d%y'`
chrlist=($(seq 1 1 22) "X" "Y")

touch index-${rundate}.swarm
touch concat-${rundate}.swarm
touch split-${rundate}.swarm
echo "#SWARM -g 32 -t 8 --time=4:00:00" > index-${rundate}.swarm
echo "#SWARM -g 32 -t 8 --time=4:00:00" > concat-${rundate}.swarm
echo "#SWARM -g 32 -t 8 --time=4:00:00" > split-${rundate}.swarm

## List of fields to exclude from INFO header and variant record 
## There are fields listed in INFO that are present in the CSQ, so I remove them from INFO before CSQ split, which will add them back
exclude=$(echo 'INFO/AC,INFO/AN,INFO/AS_FilterStatus,INFO/AS_RAW_BaseQRankSum,INFO/AS_RAW_MQ,INFO/AS_RAW_MQRankSum,INFO/AS_RAW_ReadPosRankSum,INFO/AS_ReadPosRankSum,INFO/AS_SB_TABLE,INFO/AS_SOR,INFO/AS_VQSLOD,INFO/AS_culprit,INFO/BaseQRankSum,
INFO/DB,INFO/ExcessHet,INFO/FS,INFO/InbreedingCoeff,INFO/MLEAC,INFO/MLEAF,INFO/MQ,INFO/MQRankSum,INFO/NEGATIVE_TRAIN_SITE,INFO/POSITIVE_TRAIN_SITE,INFO/QD,INFO/RAW_MQandDP,INFO/ReadPosRankSum,INFO/SOR,INFO/VQSLOD,INFO/culprit,
INFO/gnomADg,INFO/gnomADg_AC,INFO/gnomADg_AN,INFO/gnomADg_nhomalt_joint,INFO/gnomADg_nhomalt,INFO/gnomADg_AC_joint,INFO/gnomADg_AF_joint,INFO/gnomADg_AF_grpmax,INFO/gnomADg_AC_XY,INFO/gnomADg_non_par,INFO/ClinVar,INFO/ClinVar_CLNSIG,
INFO/ClinVar_CLNREVSTAT,INFO/ClinVar_CLNDN,INFO/GERP')

mkdir -p ${PWD}/buildhail_logs

csqsplit_vcfs=()

for folder in $batches; 

    do echo "
    $folder:
    "
    
    prefix=$(basename ${folder} | sed 's/-VEP-VCFs//g')
    printf "
    Concatenated chr-VCFs for ${folder} to be found in ${folder}/${prefix}.concat.vcf.gz.
    "

    vcflist=()

    ## swarm for preparing batch-chromosome VCFs with compression and indexing 
    for chr in "${chrlist[@]}";

        do if test -f "${folder}/chr${chr}.filtered.vcf.gz" ; then 
            printf "
            ${folder}/chr${chr}.filtered.vcf.gz found to be compressed. Creating index now.
            "
            printf "bcftools index -f $folder/chr${chr}.filtered.vcf.gz\n" >> index-${rundate}.swarm
            vcflist+="$folder/chr${chr}.filtered.vcf.gz " 
        elif test -f "${folder}/chr${chr}.filtered.vcf"; then
            printf "
            ${folder}/chr${chr}.filtered.vcf.gz found uncompressed. Compressing and creating index.
            "
            printf "bgzip ${folder}/chr${chr}.filtered.vcf; bcftools index -f ${folder}/chr${chr}.filtered.vcf.gz\n" >> index-${rundate}.swarm
            vcflist+="$folder/chr${chr}.filtered.vcf.gz " 
        else 
            printf "
            Expected VCF file ${folder}/chr${chr}.filtered.vcf[.gz] not found.
            "
        fi
    done
    
    ## generating swarm for concatenating VCFs in a batch 
    echo "
    Concatenating chr-VCFs with following command: 
    
    bcftools concat -Oz -o ${folder}/${prefix}.concat.vcf.gz -W ${vcflist}
    "

    printf "bcftools concat -Oz -o ${folder}/${prefix}.concat.vcf.gz -W ${vcflist}\n" >> concat-${rundate}.swarm

    ## Generating swarm for  splitting CSQ lines and removing unnecessary fields in batch VCFs 
    printf "bcftools +split-vep -c- -d ${folder}/${prefix}.concat.vcf.gz | \
        bcftools annotate -x "${exclude}" -Oz -o ${folder}/${prefix}.csqsplit.vcf.gz -W\n" >> split-${rundate}.swarm
    
    csqsplit_vcfs+="${folder}/${prefix}.csqsplit.vcf.gz "

done

## Preparing script to merge csq-split batch VCFs after swarms are complete 
printf "#!/bin/sh\n\nmodule load bcftools\n" > merge-${rundate}.sh 
printf "bcftools merge -m none -Oz -o multibatch.merged.${rundate}.vcf.gz --write-index ${csqsplit_vcfs}\n" >> merge-${rundate}.sh 

index_jid=$(swarm --module bcftools --gres=lscratch:200 -g 32 -t 8 --logdir buildhail_logs index-${rundate}.swarm)
echo "Concatenation will run after the compression and indexing swarm jobs with the JOBID ${index_jid} are complete."

concat_jid=$(swarm --module bcftools --gres=lscratch:300 --dependency afterok:$index_jid -g 32 -t 8 --logdir buildhail_logs concat-${rundate}.swarm)
echo "Concatenating chr-VCFs into whole batch VCFs. Will begin splitting consequence lines after swarm jobs with the JOBID ${concat_jid}."

split_jid=$(swarm --module bcftools --gres=lscratch:300 --dependency afterok:${concat_jid} -g 32 -t 8 --logdir buildhail_logs split-${rundate}.swarm)
echo "Performing bcftools split-CSQ on the swarm jobs with the JOBID ${split_jid}."

sbatch --mem=24g -c 8 --gres=lscratch:50 --dependency=afterok:${split_jid} merge-${rundate}.sh

