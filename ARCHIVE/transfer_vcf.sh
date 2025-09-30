#!/bin/sh

### INPUTS -----------------------

Help()
{
	# Display help 
    echo ""
	echo "For transferring unannotated per-sample VCFs to their respective folders in /data/Kastner_PFS/WGS/.
    "

	echo "To run this pipeline use the following command:"
	echo "sbatch [OPTIONS] -g [GATK VCF] -b [BCFTOOLS VCF] -o [ORIGINAL BATCH TEXT FILE]

    The original batch text file is required to identify the YEAR directory a given sample belogns to. 
    "

}

while getopts "g:b:o:h" option; do
   case $option in
    g) g=$OPTARG ;;
    b) b=$OPTARG ;;
    o) orig=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

## Ensuring the batch file was passed, and exiting otherwise. 
if [ -z "${orig}" ]; then 
    printf "ERROR Original batch text file not passed, exiting script.\n"
    exit 1
elif [ ! -f "${g}" ]; then 
    printf "ERROR ${g} does not exist, exiting script.\n"
    exit 1 
elif [ ! -f "${b}" ]; then 
    printf "ERROR ${b} does not exist, exiting script.\n"
    exit 1 
fi

module load bcftools/1.19 

# IDENTIFYING SAMPLE CONGRUENCE AND STORAGE DESTINATION --------------------

## Printing lists of samples represented in the batch and VCFs to confirm they are shared 
batch_samples=$(cat ${orig} | awk -F '/' '{print $NF}' | sed 's/.bam//g' | sort)
g_samples=$(bcftools query -l ${g} | sort)
b_samples=$(bcftools query -l ${b} | sort)

printf "\nSamples identified:\n"
echo "batch file: $(echo ${batch_samples[@]})"
echo "gatk VCF: $(echo ${g_samples[@]})" 
echo "bcftools VCF: $(echo ${b_samples[@]})"

## If the lists of batch and VCF samples are not the same, the user is prompted to indicate whether they want to continue with preparing the transfer 
if [[ $batch_samples != $g_samples ]]; then
    printf "\nDifferent samples found between batch and ${g}, do you want to continue running the transfer?.\n"
    read -p "(y/n)?" choice
    case "$choice" in 
        y|Y ) printf "User acknowledges samples may be different, continuing with transfer.\n" ;;
        n|N ) printf "no, exiting script\n" ; exit 1 ;;
        *) printf "reply y or no, exiting script\n" ; exit 1 ;;
    esac 
fi

if [[ $batch_samples != $b_samples ]]; then
    printf "Different samples found between batch and ${b}, do you want to continue running the transfer?\n"
    read -p "(y/n)?" choice
    case "$choice" in 
        y|Y ) printf "User acknowledges samples may be different, continuing with transfer.\n" ;;
        n|N ) printf "no, exiting script\n" ; exit 1 ;;
        *) printf "reply y or no\n" ; exit 1 ;;
    esac
fi

batch1=$(grep '/data/Kastner_PFS/WGS/20*' ${orig})
batch2=$(grep -v '/data/Kastner_PFS/WGS/20*' ${orig})

declare -A date_dict
no_date=()

for line in ${batch}
do

    set +e 

    if echo ${batch1} | grep -q ${line}
        then 
            ## Isolating the sample name from the original BAM file 
            sample=$(echo ${line} | awk -F '/' '{print $NF}' | sed 's/.bam//g' )

            ## Parsing path to get original date of receipt 
            mod_date=$(echo ${line} | cut -f 5 -d'/')

            echo $sample $mod_date
    elif echo ${batch2} | grep -q ${line} 
        then
            sample=$(echo ${line} | awk -F '/' '{print $NF}' | sed 's/.bam//g' )

            ## Identifying the year the original BAM was modified or created - for sorting in ketu 
            mod_date=$(date -r ${line} | awk -F ' ' '{print $NF}')
            echo $sample $mod_date 
    fi

    if [ -z ${mod_date} ]; then
        printf "No origin date identified for ${sample}.\n"
        no_date+="${sample} "
        continue
    else
        printf "Origin date of ${sample} identified as ${mod_date}.\n"
        date_dict[${sample}]=${mod_date}
    fi

    set -e

done 

# rundate=`date +'%m%d%y'`

# prefix=$( echo $orig | sed 's/.txt//g')
# mkdir -p vcf-transfer-${prefix}-logs
# echo "#SWARM -g 5 --time=4:00:00 --logdir vcf-transfer-${prefix}-logs" > vcf-transfer-${prefix}.swarm

# for key in "${!date_dict[@]}"; do
#     if echo "${g_samples}" | grep -q -w "${key}"; then 
#         echo "mkdir -p /data/Kastner_PFS/WGS/${date_dict[$key]}/${key}; bcftools view -s ${key} -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${date_dict[$key]}/${key}/${key}.gatk.vcf.gz ${g}" >> vcf-transfer-${prefix}.swarm
#         echo "mkdir -p /data/Kastner_PFS/WGS/${date_dict[$key]}/${key}; bcftools view -s ${key} -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${date_dict[$key]}/${key}/${key}.bcftools.vcf.gz ${b}" >> vcf-transfer-${prefix}.swarm
#     fi 
# done

# for s in ${no_date}; do 
#     if echo "${g_samples}" | grep -q -w "${sample}"; then 
#         echo "bcftools view -s ${s} -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${s}.gatk.vcf.gz ${g}" >> vcf-transfer-${prefix}.swarm
#         echo "bcftools view -s ${s} -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${s}.bcftools.vcf.gz ${b}" >> vcf-transfer-${prefix}.swarm
#     fi
# done 

# # STORING ANNOTATED VCF --------------

# ## Need to identify unfiltered annotated VCFs 
# ## Finding the directory for VEP output - confirming it exists or exiting 
# gdir=$(ls $(dirname ${g}) | grep -E 'indel.SNP.recalibrated_99.9.[0-9]{6}+\-VEP\-VCFs')
# bdir=$(ls $(dirname ${b}) | grep -E 'filter-snpsindels-[0-9]{6}+\-VEP\-VCFs')

# if [ -z  ${gdir} ]; then
#     printf "\nGATK VEP output directory not found, exiting script.\n"
#     exit 1 
# elif [ -z ${bdir} ]; then
#     printf "\nBcftools VEP output directory not found, exiting script.\n"
#     exit 1  
# else 
#     printf "\nVEP output directories found: 
#     GATK: ${gdir}
#     bcftools: ${bdir}\n"
# fi

# ## Listing chromosome vcfs 
# chrlist=$(echo $(seq 1 1 22) "X" "Y")
# g_vcfs=()
# b_vcfs=()

# for chr in ${chrlist}; do 
#     g_vcf=$(ls $(dirname ${g})/${gdir}/chr${chr}.vcf*)
#     b_vcf=$(ls $(dirname ${b})/${bdir}/chr${chr}.vcf*)

#     g_vcfs+="${g_vcf} "
#     b_vcfs+="${b_vcf} "
# done 

# ## preparing output directories in Kastner_PFS
# batchdate=$(echo ${orig} | sed 's/batch-//g ; s/.txt//g')
# outdir=/data/Kastner_PFS/WGS/cohort_db/batch_VCFs/${batchdate}
# g_prefix=$( echo ${g} | awk -F '/' '{print $NF}' | sed 's/.vcf.gz//g' )
# b_prefix=$( echo ${b} | awk -F '/' '{print $NF}' | sed 's/.vcf.gz//g')
# mkdir -p ${outdir}

# ## Concatenating chromosome VCFs and outputting to Kastner_PFS 
# bcftools concat -a -Oz -o ${outdir}/${g_prefix}.concat.vcf.gz ${g_vcfs}
# bcftools concat -a -Oz -o ${outdir}/${b_prefix}.concat.vcf.gz ${b_vcfs}







