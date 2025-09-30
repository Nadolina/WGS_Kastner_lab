#!/bin/sh

set -e 

Help()
{
        # Display help 
        echo "
    To automate offloading files to ketu by globus commandline." 
    echo "
    Use the following command:

        sbatch [parameters] globus_transfer.sh -b [original batch]

        -b batch file of paths to original bams 
        -h help

        "

}

while getopts "b:h" option; do
   case $option in
        b) batchfile=${OPTARG} ;;
        h) # display Help
         Help
         exit;;
   esac
done

## Initiating the batch file that will be passed to the globus transfer command 
rundate=`date +'%m%d%y'`
>globus-transfer-${rundate}-${batchfile}

set +e

batch1=$(grep '/data/Kastner_PFS/WGS/20*' ${batchfile})
batch2=$(grep -v '/data/Kastner_PFS/WGS/20*' ${batchfile})

set -e

fullbatch=$(cat ${batchfile})


for line in $fullbatch
do 

    set +e 

    ## if-fi checks for formatting of file name
    ## some files have been moved in Kastner_PFS and so their "modification date" does not reflect their true date of receipt
    if echo ${batch1} | grep -q ${line}
        then 
            ## Isolating the sample name from the original BAM file 
            sample=$(echo ${line} | awk -F '/' '{print $NF}' | sed 's/.bam//g' )

            ## Parsing path to get original date of receipt 
            mod_date=$(echo ${line} | cut -f 5 -d'/')
            printf "The bam file ${line} was created in ${mod_date}.\n"
    elif echo ${batch2} | grep -q ${line} 
        then
            sample=$(echo ${line} | awk -F '/' '{print $NF}' | sed 's/.bam//g' )

            ## Identifying the year the original BAM was modified or created - for sorting in ketu 
            mod_date=$(date -r ${line} | awk -F ' ' '{print $NF}')
            printf "The bam file ${line} was created in ${mod_date}.\n"
    fi

    set -e

echo $mod_date $sample
        
    ## globus mkdir does not have the -p function like linux to check for the existence of a file or folder 
    ## using the "stat" output to see whether a folder already exists for a given sample, and creating it if not 
    status=$(globus stat c1cb6b0e-4b64-11ef-b90c-c30b78766494:/mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample})
    if [[ ${status} == "Nothing found at"* ]]; then 
        printf "${sample} not found.\n"
        globus mkdir c1cb6b0e-4b64-11ef-b90c-c30b78766494:/mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample}
    else
        printf "Folder found for ${sample} in ketu, see status below:\n"
        printf "${status}\n"
    fi

    ## Conditionals to identify the existence of BQSR BAM and *bai files for a sample, and appending to the globus transfer batch file  
    if [[ -f ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bam && -f ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bai ]]; then 
        printf "Found bam for sample ${sample}: ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bam\n"
        printf "Found index: ${sample}: ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bai\n"
        ls ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bam
        ls ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bai

        echo "--no-recursive ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bam /mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample}/${sample}_bqsr_1.bam" >> globus-transfer-${rundate}-${batchfile}
        echo "--no-recursive ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bai /mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample}/${sample}_bqsr_1.bai" >> globus-transfer-${rundate}-${batchfile}
    
    elif [ -f ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bam ]; then 
        printf "Only found bam for sample ${sample}, but no index: ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bam\n"
        ls ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bam

        echo "--no-recursive ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bam /mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample}/${sample}_bqsr_1.bam" >> globus-transfer-${rundate}-${batchfile}

    elif [ -f ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bai ]; then 
        printf "Only found an index for BQSR bam: ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bai\n"
        printf "Confirm BAM exists or was already transferred.\n"
        ls ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bai

        echo "--no-recursive ${PWD}/${sample}/${sample}_out/${sample}_bqsr_1.bai /mnt/brajukanm/ketu_labs/Kastner/WGS/${mod_date}/${sample}/${sample}_bqsr_1.bai" >> globus-transfer-${rundate}-${batchfile}

    else
        printf "No bam or index found for sample ${sample}.\n"
    fi

done 

## Instructions for output batch file. 
printf "\nThe file globus-transfer-${rundate}-${batchfile} will contain any BQSR *bam and *bai files available to transfer to ketu.\n"
printf "Run the command 

        globus transfer --no-verify-checksum --batch globus-transfer-${rundate}-${batchfile} {origin endpoint globus ID} {destination endpoint globus ID}

    The origin endpoint should be e2620047-6d04-11e5-ba46-22000b92c6ec. Use this command to verify, look for 'NIH HPC Data Transfer (Biowulf)': 

        globus endpoint search NIH 

    Your ketu globus endpoint will be different. Use the below to get your ID:

        globus endpoint search nhgrishell02 

"