#!/bin/sh 

# INPUTS -----------------------

Help()
{
	# Display help 
    echo ""
	echo "After processing, unannotated batch VCFs are transferred to Kastner_PFS for permanent storage. Individual sample VCFs are similarly placed in their designated WGS/[YEAR]/[SAMPLE] directories.
    "

	echo "To run this pipeline use the following command:"
	echo "sbatch [OPTIONS] -g [GATK VCF] -b [BCFTOOLS VCF] -o [ORIGINAL BATCH TEXT FILE]

    The original batch text file is required to identify the YEAR directory a given sample belongs to. 
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

## Confirming samples in batch file and VCFs are the same 

batch_samples=$( cat ${orig} | awk -F '/' '{print $NF}' | sed 's/.bam//g' | sort )
g_samples=$(bcftools query -l ${g} | sort)
b_samples=$(bcftools query -l ${b} | sort)

printf "\nSamples identified:\n"
echo "batch file: $(echo ${batch_samples[@]})"
echo "gatk VCF: $(echo ${g_samples[@]})" 
echo "bcftools VCF: $(echo ${b_samples[@]})"

## User input to move forward with script if some samples do not match 

printf "
WARNING any discongruent samples between the VCF(s) and batch filewill need to be tranferred manually. 
"
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

## Discriminating between paths in batch file pointing to the prescribed WGS storage location (/data/Kastner_PFS/WGS/[YEAR]) or to older storage locations 
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
    elif echo ${batch2} | grep -q ${line} 
        then
            sample=$(echo ${line} | awk -F '/' '{print $NF}' | sed 's/.bam//g' )

            ## Identifying the year the original BAM was modified or created - for sorting in ketu 
            mod_date=$(date -r ${line} | awk -F ' ' '{print $NF}')
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

rundate=`date +'%m%d%y'`

prefix=$( echo $orig | sed 's/.txt//g')
mkdir -p vcf-transfer-${prefix}-logs
echo "#SWARM -g 5 --time=4:00:00 --logdir vcf-transfer-${prefix}-logs" > vcf-transfer-${prefix}.swarm

for key in "${!date_dict[@]}"; do
    if echo "${g_samples}" | grep -q -w "${key}"; then 
        echo "mkdir -p /data/Kastner_PFS/WGS/${date_dict[$key]}/${key}; bcftools view -s ${key} -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${date_dict[$key]}/${key}/${key}.gatk.vcf.gz ${g}" >> vcf-transfer-${prefix}.swarm
        echo "mkdir -p /data/Kastner_PFS/WGS/${date_dict[$key]}/${key}; bcftools view -s ${key} -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${date_dict[$key]}/${key}/${key}.bcftools.vcf.gz ${b}" >> vcf-transfer-${prefix}.swarm
    fi 
done

for s in ${no_date}; do 
    if echo "${g_samples}" | grep -q -w "${sample}"; then 
        echo "bcftools view -s ${s} -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${s}.gatk.vcf.gz ${g}" >> vcf-transfer-${prefix}.swarm
        echo "bcftools view -s ${s} -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${s}.bcftools.vcf.gz ${b}" >> vcf-transfer-${prefix}.swarm
    fi
done 