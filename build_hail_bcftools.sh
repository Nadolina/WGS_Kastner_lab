#!/bin/sh

set -e

### INPUTS -----------------------

Help()
{
	# Display help 
    echo ""
	echo "This script is very much like build_hail_2.sh, but is designed for the bcftools calls, not GATK. 
    "

	echo "To run this pipeline use the following command:"
	echo "sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] build_hail_bcftools.sh [folder names of batches]
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

## Parsing the input for the list of folders containing annotated-filtered VCFs 
batches=()
printf "
The batches provided include:
"
for folder in "$@"; do 
    batches+="${folder} "
    echo "  $folder"
done

### FUNCTIONS -----------------------

prep_vcfs() {

    ## list of chromosome numbers included in each hail database 
    db=$1 

    ## Tracking which db we are on 
    db_num=$(( $db_num + 1 )) 
    printf "Preparing VCFs for db${db_bum}.\n"

    ## Initializing list of VCFs to concatenate (to form one hail matrix table) for a given db
    vcflist_db=()

    printf "\nIdentifying compression status of chromosomes in $batches and preparing for concatenation/CSQ line splitting:\n"
    for chr in $db
    do

        ## Initializing list of VCFs to merge for a chromosome with a db (i.e./ chr3 from batchA, chr3 from batchB, chr3 from batchC)
        vcflist_chr=()

        ## Looping through each chromosome and batch to perform compression, indexing and prepare lists of VCFs for merging and concatenation 
        ## Generating swarm scripts 
        for folder in ${batches}
        do if test -f "${folder}/chr${chr}.filtered.vcf.gz" ; then 
            printf "${folder}/chr${chr}.filtered.vcf.gz found to be compressed. Creating index now.\n"
            printf "bcftools index -f $folder/chr${chr}.filtered.vcf.gz\n" >> index-bcftools-${rundate}.swarm

            vcflist_chr+="$folder/chr${chr}.filtered.vcf.gz " 

        elif test -f "${folder}/chr${chr}.filtered.vcf"; then
            printf "${folder}/chr${chr}.filtered.vcf.gz found uncompressed. Compressing and creating index.\n"
            bgzip -c ${folder}/chr${chr}.filtered.vcf > ${folder}/chr${chr}.filtered.vcf.gz 
            # printf "bgzip ${folder}/chr${chr}.filtered.vcf && bcftools index -f ${folder}/chr${chr}.filtered.vcf.gz\n" >> index-${rundate}.swarm
            printf "bcftools index -f ${folder}/chr${chr}.filtered.vcf.gz\n" >> index-bcftools-${rundate}.swarm

            vcflist_chr+="$folder/chr${chr}.filtered.vcf.gz " 

        else 
            printf "Expected VCF file ${folder}/chr${chr}.filtered.vcf[.gz] not found.\n"
        fi
        done

        ## Merging VCFs for a given chromosome across the batches, the splitting the CSQ line and removing an unnecessary files
        ## Writes output to a folder of the name buildhail_bcftools_[run date]
        printf "Merging VCFs ${vcflist_chr}\n"
        printf "bcftools merge -m none ${vcflist_chr} | bcftools +split-vep -c- -d | bcftools annotate --threads 4 -x "${exclude}" -Oz -o buildhail_bcftools_${rundate}/chr${chr}.split.vcf.gz --write-index\n"  >> split-bcftools-${rundate}.swarm 

        vcflist_db+="buildhail_bcftools_${rundate}/chr${chr}.split.vcf.gz "

    done 

    printf "CSQ-split VCFs to concatenate include ${vcflist_db}"
    printf "bcftools concat --threads 4 -Oz -o buildhail_bcftools_${rundate}/db${db_num}.concat.${rundate}.vcf.gz ${vcflist_db}\n" >> concat-bcftools-${rundate}.swarm

}

## VARIABLES AND SETUP -----------------------

rundate=`date +'%m%d%y'`
mkdir -p buildhail_bcftools_${rundate}_logs 
mkdir -p buildhail_bcftools_${rundate}

# Defining the divison of chromosomes to pass to function
db1=$(seq 1 4)
db2=$(echo $(seq 5 10))
db3=$(echo $(seq 11 17))
db4=$(echo $(seq 18 22) "X" "Y")

printf "\nThe database will be divided by chromosome as such:\n
    Database 1: chr1-4
    Database 2: chr5-10
    Database 3: chr11-17
    Database 4: chr18-22,X,Y\n"

## Initializing swarm files to produce from prep_vcf
echo "#SWARM -g 5 --time=4:00:00" > index-bcftools-${rundate}.swarm
echo "#SWARM -g 12 -t 4 --time=4:00:00" > split-bcftools-${rundate}.swarm
echo "#SWARM -g 12 -t 4 --time=4:00:00" > concat-bcftools-${rundate}.swarm

## Fields in the VCF that are not useful and get removed
exclude=$(echo 'INFO/IDV,INFO/VDB,INFO/RPBZ,INFO/MQBZ,INFO/BQBZ,INFO/MQSBZ,INFO/SCBZ,INFO/MQ0F,INFO/CSQ')
echo $exclude

### RUNNING FUNCTIONS AND SWARMS -----------------------

## prep_vcfs loop increases db_num value to track which db we are on and name outputs accordingly 
db_num=0

prep_vcfs "${db1[@]}"
prep_vcfs "${db2[@]}"
prep_vcfs "${db3[@]}"
prep_vcfs "${db4[@]}"

## Running swarm with job IDs dependencies to ensure the swarms run successively 
index_jid=$(swarm --module bcftools/1.19 -b 3 --logdir buildhail_bcftools_${rundate}_logs index-bcftools-${rundate}.swarm)
split_jid=$(swarm --module bcftools/1.19 --dependency afterok:${index_jid} --logdir buildhail_bcftools_${rundate}_logs split-bcftools-${rundate}.swarm)
swarm --module bcftools/1.19 --dependency afterok:${split_jid} --logdir buildhail_bcftools_${rundate}_logs concat-bcftools-${rundate}.swarm


