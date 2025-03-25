#!/bin/sh

## REFERENCES 
## 1. https://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html)
## 2. https://bioinformatics.stackexchange.com/questions/22124/variants-from-multiple-tools-normalization-before-or-after-annotation-with-vep
## 3. https://www.ensembl.info/2020/05/26/normalising-variants-to-standardise-ensembl-vep-output/
## 4. https://samtools.github.io/bcftools/bcftools.html#terminology
## 5. https://gatk.broadinstitute.org/hc/en-us/articles/360035531912-Spanning-or-overlapping-deletions-allele#:~:text=We%20use%20the%20term%20spanning,heterozygous%20or%20homozygous%20variant%20form
## 6. https://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#input


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


# identifying the file type and deciding whether to compress and/or index 
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

# ALT allele deletions in GATK are signified with an asterisk (*), but VEP employs a dash (-)
# ALT deletions left as an * will not be annotated 
# Here, I normalize the VCF (as per recommendations in 1,2,3) and substitute asterisks in the ALT column for dashes (see 4,5,6 for more info on the * and - notation in bcftools and vep).

bcftools sort ${VCF} | \
    bcftools norm --threads 4 -m-any --check-ref x -f ${ref} --force | \
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

if test -f ${prefix}-VEP.swarm 
then
    printf "${prefix}-VEP.swarm found - overwritting.\n"
    rm ${prefix}-VEP.swarm
fi

echo "#SWARM -g 10 --time=1-00:00:00 --logdir ${prefix}-VEP-logs" > "${prefix}-VEP.swarm"
echo "#SWARM -g 5 --time=4:00:00 --logdir ${prefix}-VEP-logs" > ${prefix}-VEP-filter.swarm
mkdir -p ${prefix}-VEP-logs 
mkdir -p ${prefix}-VEP-VCFs 

##Generating a swarm file where there is each line is the annotation of one chromosome 
for value in "${chrlist[@]}"
do 
    chrnum="chr${value}"

    echo "tabix -h ${prefix}.norm.vcf.gz ${chrnum} | vep \
        --fasta /data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta \
        --format vcf \
        --species human \
        --assembly GRCh38 \
        --vcf \
        -o ${prefix}-VEP-VCFs/${chrnum}.vcf \
        --force_overwrite \
        --offline --cache --dir ${VEP_CACHEDIR} \
        --hgvs \
        --check_existing \
        --plugin HGVSReferenceBase \
        --plugin CADD,${VEP_CACHEDIR}/CADD_1.4_GRCh38_whole_genome_SNVs.tsv.gz,${VEP_CACHEDIR}/CADD_1.4_GRCh38_InDels.tsv.gz \
        --custom file=/data/Kastner_PFS/references/HG38/gnomad.genomes.v4.0.sites.noVEP.subset_annots.vcf.gz,short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AC%AN%nhomalt_joint%nhomalt%AC_joint%AF_joint%AF_grpmax%AC_XY%non_par%sift_max%polyphen_max \
        --plugin SpliceAI,snv=${VEP_CACHEDIR}/spliceai_scores.masked.snv.hg38.vcf.gz,indel=${VEP_CACHEDIR}/spliceai_scores.masked.indel.hg38.vcf.gz \
        --plugin MaxEntScan,${VEP_CACHEDIR}/MaxEntScan \
        --plugin REVEL,${VEP_CACHEDIR}/revel_GRCh38_1.3.tsv.gz \
        --plugin AlphaMissense,file=${VEP_CACHEDIR}/AlphaMissense_hg38.tsv.gz,cols=all \
        --custom file=/fdb/clinvar/vcf_GRCh38/clinvar_20220517.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
        --custom file=/fdb/VEP/loftee/2022-06/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw,short_name=GERP,format=bigwig \
        --flag_pick \
        --pick_order rank \
        --canonical \
        --no_stats --verbose 2> ${prefix}-VEP-logs/${chrnum}.log" >> ${prefix}-VEP.swarm

        # --plugin GeneBe \

    # a second swarm for filtering post annotation
    echo "filter_vep -i ${prefix}-VEP-VCFs/${chrnum}.vcf \
        -f 'gnomADg_AF_grpmax < 0.1 or not gnomADg_AF_grpmax' \
        --only_matched \
        --format vcf \
        --output_file ${prefix}-VEP-VCFs/${chrnum}.filtered.vcf \
        --force_overwrite \
        2> ${prefix}-VEP-logs/${chrnum}.filter.log" >> ${prefix}-VEP-filter.swarm 

done 

vep_jid=$(swarm --module VEP/112,samtools/1.19 ${prefix}-VEP.swarm)

swarm --dependency afterok:${vep_jid} --module VEP/112 --time=04:00:00 ${prefix}-VEP-filter.swarm



