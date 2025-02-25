#!/bin/sh

## REFERENCES 
## 1. https://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html)
## 2. 

set -e 

Help()
{
	# Display help 
	echo "
    To automate annotation of batch VCFs." 
    echo "
    Use the following command:

        sbatch [parameters] run-VEP.sh -v [batch VCF]

        -h help

        "

}

while getopts "v:h" option; do
   case $option in
	v) VCF=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

module load VEP/112
module load samtools/1.19
module load bcftools/1.19
module load GATK/4.6.0.0

chrlist=($(seq 1 1 22) "X" "Y") 
ref=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta


## identifying the file type and deciding whether to compress and/or index 
if [[ $(file -b --mime-type $VCF) == 'application/x-gzip' ]]
then
    printf "BGZIP file recognized, performing tabix indexing.\n"
    bcftools index -f -t $VCF 
elif [[ $(file -b --mime-type $VCF)  == 'text/plain' ]]
then
    printf "Uncompressed file recognized, performing bgzip and indexing.\n"
    bgzip $VCF ; bcftools index -f -t "${VCF}.gz"
    VCF="${VCF}.gz"
else 
    printf "Unrecognized file type, pass *.vcf or *.vcf.gz\n"
fi

prefix=$(echo ${VCF} | sed 's/.vcf.gz//g')

## ALT allele deletions in GATK are signified with an asterisk (*), but VEP employs a dash (-)
## ALT deletions left as an * will not be annotated 
## Here, I normalize the VCF (as per recommendations in 1,2,3) and substitute asterisks in the ALT column for dashes.

bcftools norm \
    --threads ${SLURM_CPUS_PER_TASK} \
    -m-any \
    --check-ref x \
    -f ${ref} \
    ${VCF} \
    --force | \
    sed 's/*\t/-\t/g' | \
    bgzip > ${prefix}.norm.vcf.gz

bcftools index ${prefix}.norm.vcf.gz

deletions=$(bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -H ${VCF} | \
    cut -f5 | \
    grep -w '*' | \
    wc -l)

echo "
Number of deletions (*) observed in ${VCF}: ${deletions}"

echo "#SWARM -t 16 -g 64 --time=12:00:00 --logdir ${SLURM_JOB_ID}.tmpdir" > "${prefix}-VEP.swarm"
echo "#SWARM -t 4 -g 16 --time=4:00:00 --logdir ${SLURM_JOB_ID}.tmpdir" > "${prefix}-VEP-filter.swarm"
mkdir -p ${prefix}-VEP-logs 
mkdir -p ${prefix}-VEP-VCFs 

for value in "${chrlist[@]}"
do 
    chrnum="chr${value}"

    echo "TWD=/lscratch/\${SLURM_JOB_ID}; \
        tabix -h ${prefix}.norm.vcf.gz ${chrnum} | vep \
            --fasta  ${ref} \
            --format vcf \
            --species human \
            --assembly GRCh38 \
            --vcf \
            -o ${prefix}-VEP-VCFs/${chrnum}.vcf \
            --force_overwrite \
            --no_stats \
            --offline --cache --dir ${VEP_CACHEDIR} \
            --plugin CADD,${VEP_CACHEDIR}/CADD_1.4_GRCh38_whole_genome_SNVs.tsv.gz,${VEP_CACHEDIR}/CADD_1.4_GRCh38_InDels.tsv.gz \ 
            --custom file=${VEP_CACHEDIR}/gnomad.genomes.v4.0.sites.noVEP.vcf.gz,short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AC%AN%AF_afr%AF_amr%AF_asj%AF_eas%AF_fin%AF_nfe%AF_oth \
            --plugin REVEL,${VEP_CACHEDIR}/revel_GRCh38_1.3.tsv.gz \
            --plugin SpliceAI,snv=${VEP_CACHEDIR}/spliceai_scores.masked.snv.hg38.vcf.gz,indel=${VEP_CACHEDIR}/spliceai_scores.masked.indel.hg38.vcf.gz \
            --plugin MaxEntScan,${VEP_CACHEDIR}/MaxEntScan \
            --verbose 2> ${prefix}-VEP-logs/${chrnum}.log" >> ${prefix}-VEP.swarm
    
    echo "TWD=/lscratch/\${SLURM_JOB_ID}; \
        filter_vep -i ${prefix}-VEP-VCFs/${chrnum}.vcf \
        -f 'gnomADg_AF_afr < 0.1 and gnomADg_AF_amr < 0.1 and gnomADg_AF_asj < 0.1 and gnomADg_AF_eas < 0.1 and gnomADg_AF_fin < 0.1 and gnomADg_AF_nfe < 0.1' \
        --only_matched \
        --format vcf \
        --output_file ${prefix}-VEP-VCFs/${chrnum}.filtered.vcf \
        --force_overwrite \
        2> ${prefix}-VEP-logs/${chrnum}.filter.log" >> ${prefix}-VEP-filter.swarm 

done 

vep_jid=$(swarm --module VEP/112,samtools/1.19 -t 16 -g 64 --gres=lscratch:200 --time=1-06:00:00 ${prefix}-VEP.swarm)

swarm --dependency afterok:${vep_jid} --module VEP/112 -t 4 -g 16 --gres=lscratch:50 --time=04:00:00 ${prefix}-VEP-filter.swarm



