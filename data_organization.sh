#!/bin/sh 

## GOAL is to take a batch file and from it:
## 1. transfer the BQSR BAM to ketu
## 2. transfer the original BAMs to ketu  
## 2. transfer the new sample VCFs to Kastner_PFS/WGS 
## 3. transfer the original sample VCFs to Kastner_PFS/WGS

Help()
{
    # Display help 
    echo "
    To automate offloading files of samples that have been re-processed fully." 
    echo "
    Use the following command:

        sbatch [parameters] globus_transfer.sh -o [batch] -g [GATK VCF (unannotated)] -b [bcftools VCF (unannotated)]

        -o batch file of paths to original bams 
        -g gatk batch VCF
        -b bcftools batch VCF
        -h help

        "

}

while getopts "o:g:b:h" option; do
   case $option in
        o) batch=${OPTARG} ;;
        g) gatk=${OPTARG} ;;
        b) bcftools=${OPTARG} ;;
        h) # display Help
         Help
         exit;;
   esac
done

module load bcftools/1.19

## FUNCTIONS -----------------------------------------

## Given the path to an original BAM used in re-alignment (per pre-process-pipe.sh and alignment_cleanup.sh), run_batch() wil:
## 1. identify the corresponding BQSR BAM and BAI
## 2. create the appropriate directory in ketu 
## 3. add the transfer commands for the original and BQSR files to the globus batch. 
function run_batch() {

    bam=$1 
    bai=$( ls $(dirname $bam)/*.bai )

    ## Confirming bam passed exists or exiting the function
    if [[ -f ${bam} && -f ${bai} ]]; then 
        printf "${bam} and ${bai} found.\n"
    else 
        printf "ERROR: ${bam} NOT found.\n"
        continue
    fi 

    ## Identifying sample ID and BQSR bams
    sample=$(echo ${bam} | awk -F '/' '{print $NF}' | sed 's/.bam//g' )
    newbam=${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bam 
    newbai=${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bai

    if [[ "$bam" == *"/data/Kastner_PFS/WGS/"* ]]; then 
        mod_date=$(echo ${bam} | cut -f 5 -d'/')
    else
        mod_date=$(date -r ${bam} | awk -F ' ' '{print $NF}')
    fi 

    ## Determining whether a directory has been created for this sample in ketu; creating if not.
    status=$(globus stat c1cb6b0e-4b64-11ef-b90c-c30b78766494:/mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample})
    if [[ ${status} == "Nothing found at"* ]]; then 
        printf "${sample} not found.\n"
        globus mkdir c1cb6b0e-4b64-11ef-b90c-c30b78766494:/mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample}
    else
        printf "Folder found for ${sample} in ketu, see status below:\n"
        printf "${status}\n"
    fi

    if [[ -f ${bam} && -f ${bai} ]]; then 
        echo ""
        echo "Adding ${bam} and ${bai} to globus transfer: globus-transfer-${rundate}-${batch_prefix}.txt"
        echo "--no-recursive ${bam} /mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample}/${sample}.bam" >> globus-transfer-${rundate}-${batch_prefix}.txt
        echo "--no-recursive ${bai} /mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample}/${sample}.bai" >> globus-transfer-${rundate}-${batch_prefix}.txt
    else 
        printf "ERROR: ${bam} or ${bai} does not exist for ${sample}. Check path is correct or if BAM has already been transfered to ketu.\n"
    fi

    ## Transferring a batch of original bams assumes they have all completed re-alignment so we check for the completion of BQSR BAMs.
    if [[ -f ${newbam} && -f ${newbai} ]]; then 
        echo ""
        echo "Adding ${newbam} and ${newbai} to globus transfer: globus-transfer-${rundate}-${batch_prefix}.txt"
        echo "--no-recursive ${newbam} /mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample}/${sample}_bqsr_1.bam" >> globus-transfer-${rundate}-${batch_prefix}.txt
        echo "--no-recursive ${newbai} /mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample}/${sample}_bqsr_1.bai" >> globus-transfer-${rundate}-${batch_prefix}.txt
    else
        printf "ERROR: ${newbam} or ${newbai} does not exist for ${sample}. Confirm sample re-alignment is complete. Skipping ${sample}.\n"
    fi

}

function get_orig_vcfs() { 

    bam=$1
    dir=$(dirname ${bam})
    sample=$(echo ${bam} | awk -F '/' '{print $NF}' | sed 's/.bam//g')

    ## filtering out BAM paths already pointing to Kastner_PFS/WGS
    ## original VCFs, if there can stay, or won't be found there 
    if [[ "${bam}" == *"/data/Kastner_PFS/WGS/"* ]]; then
        continue 
    else
        ## Gathering lists of standard VCFs from NISC 
        filtvcf=$(ls ${dir}/*filtered.vcf.gz*)
        gvcf=$(ls ${dir}/*g.vcf.gz*)
        ## Finding origin date 
        mod_date=$(date -r  ${dir}/${sample}.filtered.vcf.gz | awk -F ' ' '{print $NF}')    

        if [ -n "${filtvcf}" ]; then 
            for file in ${filtvcf}
                do 
                    echo "mv ${file} /data/Kastner_PFS/WGS/${mod_date}/${sample}/" >> orig_vcf_transfer_${rundate}-${batch_prefix}.txt
                done
        else 
            echo "Original filtered VCFs not found for ${sample}."
        fi

        if [ -n "$gvcf" ]; then 
            for file in ${gvcf}
                do 
                    echo "mv ${file} /data/Kastner_PFS/WGS/${mod_date}/${sample}/" >> orig_vcf_transfer_${rundate}-${batch_prefix}.txt
                done
        else 
            echo "Original GVCFs not found for ${sample}."
        fi

    fi 

}

## EXECUTING FUNCTIONS -------------------------------

batch_prefix=$(echo ${batch} | sed 's/.txt//g')

if [[ -n ${batch} ]]; then

    ## Initializing globus transfer batch file
    rundate=`date +'%m%d%y'`
    >globus-transfer-${rundate}-${batch_prefix}.txt
    >orig_vcf_transfer_${rundate}-${batch_prefix}.txt

    ## Looping through BAM paths in batch file to:
    ##  1. generate the globus transfer batch file 
    ##  2. generate a text file to move original VCFs to their appropriate sub-directory in /data/Kastner_PFS/WGS/[YEAR]/[SAMPLE]
    while read bam; 
        do 
            run_batch ${bam}
            get_orig_vcfs ${bam}
    done < ${batch}

else
    echo "Batch file was not passed, checking other arguments."
fi

## If the batch text file and the GATK and/or the bcftools VCF were passed generate VCFs to distribute to sample directories in Kastner_PFS/WGS/[YEAR]
if [[ -f ${batch} ]] && ([[ -f ${gatk} ]] || [[ -f ${bcftools} ]]); then 

    echo "${batch} and VCFs provided - starting sample VCF preparation." 

    ## Initating a swarm file for generating the sample VCFs 
    echo "#SWARM -g 12 -t 4 --gres=lscratch:50 --time 04:00:00 --logdir vcf_transfer_${batch_prefix}_logs" > vcf_transfer_${batch_prefix}.swarm

    while read bam; do 

        ## Variables needed to make appropriate destination directories
        dir=$(dirname ${bam})
        sample=$(echo ${bam} | awk -F '/' '{print $NF}' | sed 's/.bam//g' )

        ## Parsing path of original BAM to find correct date of origin 
        ## BAMs may have been moved to a new directory /data/Kastner_PFS/[YEAR] so their modification date doesn't reflect their true year of origin
        ## New directories sorted by date - so the path can be used to find a date 
        if [[ "$bam" == *"/data/Kastner_PFS/WGS/"* ]]; then 
            mod_date=$(echo ${bam} | cut -f 5 -d'/')
        else
            mod_date=$(date -r ${bam} | awk -F ' ' '{print $NF}')
        fi 

        ## Creating the sample's directory if it does not exist already
        if [ -d /data/Kastner_PFS/WGS/${mod_date}/${sample} ]; then 
            printf "Directory exists: /data/Kastner_PFS/WGS/${mod_date}/${sample}\n"
        else 
            printf "Directory DOES NOT exist: /data/Kastner_PFS/WGS/${mod_date}/${sample}\n"
            printf "Making directory."
            mkdir -p /data/Kastner_PFS/WGS/${mod_date}/${sample}
        fi

        ## Adding commands to swarm file if gatk and/or bcftools batch VCFs were passed
        if [ -n ${gatk} ]; then 
            ## Creating an individual VCF for the sample 
            echo ""
            echo "Adding command to vcf_transfer_${batch_prefix}.swarm to generate VCF for sample ${sample} from GATK variant calls: 
                bcftools view -s ${sample} --threads 4 -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${mod_date}/${sample}/gatk.${sample}.vcf.gz ${PWD}/{gatk}"
            echo "bcftools view -s ${sample} --threads 4 -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${mod_date}/${sample}/gatk.${sample}.vcf.gz ${PWD}/${gatk}" >> vcf_transfer_${batch_prefix}.swarm
        fi 

        if [ -n ${bcftools} ]; then
            echo ""
            echo "Adding command to vcf_transfer_${batch_prefix}.swarm to generate VCF for sample ${sample} from bcftools variant calls:
                bcftools view -s ${sample} --threads 4 -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${mod_date}/${sample}/bcftools.${sample}.vcf.gz ${PWD}/${bcftools}"
            echo "bcftools view -s ${sample} --threads 4 -U -c 1 -Oz --write-index -o /data/Kastner_PFS/WGS/${mod_date}/${sample}/bcftools.${sample}.vcf.gz ${PWD}/${bcftools}" >> vcf_transfer_${batch_prefix}.swarm
        fi

    done < ${batch}
else 
    echo "No batch VCFs provided, exiting program without preparing sample VCFs."
fi

## Instructions for output batch and swarm files. 
printf "\nThe file globus-transfer-${rundate}-${batch_prefix}.txt will contain any BQSR *bam/*bai and original *bam/*bai files available to transfer to ketu for the samples in ${batch}.\n"
printf "Run the command 

        globus transfer --no-verify-checksum --batch globus-transfer-${rundate}-${batch_prefix}.txt {origin endpoint globus ID} {destination endpoint globus ID}

    The origin endpoint should be e2620047-6d04-11e5-ba46-22000b92c6ec. Use this command to verify, look for 'NIH HPC Data Transfer (Biowulf)': 

        globus endpoint search NIH 

    Your ketu globus endpoint will be different. Use the below to get your ID:

        globus endpoint search nhgrishell02 

    Add --dry-run to the end of the globus transfer command to test the validity of the transfer. Remove the flag (as shown above) to submit the transfer.

"
printf "\nThe file vcf_transfer_${batch_prefix}.swarm contains a series of bcftools commands to generate per-sample VCFs from the batch VCFs of GATK and bcftools variant calls.\n"
printf "Run the following command to run the swarm

        swarm --module bcftools vcf_transfer_${batch_prefix}.swarm

"
