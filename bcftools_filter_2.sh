#!/bin/sh

set -e

Help()
{
	# Display help 
	echo "Combining the genotyped VCFs from each chromosome."
	echo ""
	echo "To run this pipeline use the following command:"
	echo "sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] bcftools_filter.sh -d [directory of bcftools called VCFs]"

}

while getopts "d:h" option; do
   case $option in
	d) vcfdir=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

module load R/4.4.1
module load bcftools/1.19
module load picard/3.2.0
module load vcflib
# module load vcftools

rundate=`date +'%m%d%y'`
cd ${vcfdir}
printf "Working in ${vcfdir}.\n" 

## merging the VCF and getting the number of variants identified across all chromosomes and samples 
vcflist=`ls | grep '.vcf.gz$' | grep '^call'| sort -V`

mergevcf=merged-bcftools-$rundate.vcf.gz

gatherinput=""
for vcf in ${vcflist};
do
    gatherinput+="I=${vcf} " 
done 

echo "java -jar $PICARDJAR GatherVcfs ${gatherinput} 0=$mergevcf" 2> gathervcfs-bcftools-${rundate}.log
java -jar $PICARDJAR GatherVcfs ${gatherinput} O=$mergevcf 2> gathervcfs-bcftools-${rundate}.log

bcftools index -t $mergevcf

samples=`bcftools query -l $mergevcf`
samplestr=`echo $samples | sed -E 's/\n/,/g'`

#printf "The number of variants in %s is: " "$samplestr" >> gathervcfs-bcftools-$rundate.log ; bcftools view -H $mergevcf | wc -l >> gathervcfs-bcftools-$rundate.log 

# Random sampling the VCF to create a subset we can identify filtering thresholds from 

bcftools view --types 'snps' $mergevcf | vcfrandomsample -r 0.01 | bgzip > subset-bcftools-snps-$rundate.vcf.gz 
bcftools index subset-bcftools-snps-$rundate.vcf.gz 
subsetsnps=subset-bcftools-snps-$rundate.vcf.gz

bcftools view --types 'indels' $mergevcf | vcfrandomsample -r 0.01 | bgzip > subset-bcftools-indels-$rundate.vcf.gz 
bcftools index subset-bcftools-indels-$rundate.vcf.gz 
subsetindels=subset-bcftools-indels-$rundate.vcf.gz

## Generating the statistics to plot in R 
mkdir -p stats 
out="stats/snp-out"

bcftools query -f '%CHROM\t%POS\t%QUAL\t%DP\t%MQBZ' subset-bcftools-snps-$rundate.vcf.gz -o snp-metrics-$rundate.tsv 
bcftools query -f '%CHROM\t%POS\t%QUAL\t%DP\t%IMF' subset-bcftools-indels-$rundate.vcf.gz -o indels-metrics-$rundate.tsv 

echo 'Rscript -e "rmarkdown::render('/data/Kastner_PFS/scripts/pipelines/WGS_Kastner_lab/bcftools_filter_thresholds_2.rmd',params=list(myargs = '${PWD}'),output_dir='${PWD}')"'
Rscript -e "rmarkdown::render('/data/Kastner_PFS/scripts/pipelines/WGS_Kastner_lab/bcftools_filter_thresholds_2.rmd',params=list(myargs = '${PWD}'),output_dir='${PWD}')"

qual=`cut -f2 -d ',' parameters_snps.csv| awk 'NR==2'`
mindepth=`cut -f2 -d ',' parameters_snps.csv| awk 'NR==3'`
maxdepth=`cut -f2 -d ',' parameters_snps.csv| awk 'NR==4'`
minMQBZ=`cut -f2 -d ',' parameters_snps.csv| awk 'NR==5'`
maxMQBZ=`cut -f2 -d ',' parameters_snps.csv | awk 'NR==6'`

filters=`echo "QUAL < $qual || INFO/DP < $mindepth || INFO/DP > ${maxdepth} || MQBZ < ${minMQBZ}"`

echo "bcftools view -e "${filters}" --min-af 0.05 --types snps -O z -o filter-maf-snps-$rundate.vcf.gz $mergevcf "
bcftools view -e "${filters}" --min-af 0.05 --types snps -O z -o filter-maf-snps-$rundate.vcf.gz $mergevcf 
echo "bcftools view -e "${filters}" --types snps -O z -o filter-nomaf-snps-$rundate.vcf.gz $mergevcf" 
bcftools view -e "${filters}" --types snps -O z -o filter-nomaf-snps-$rundate.vcf.gz $mergevcf

iqual=`cut -f2 -d ',' parameters_indels.csv| awk 'NR==2'`
imindepth=`cut -f2 -d ',' parameters_indels.csv| awk 'NR==3'`
imaxdepth=`cut -f2 -d ',' parameters_indels.csv| awk 'NR==4'`
imf=`cut -f2 -d ',' parameters_indels.csv| awk 'NR==5'`

ifilters=`echo "QUAL < $iqual || INFO/DP < $imindepth || INFO/DP > ${imaxdepth} || IMF < 0.1"`

echo "bcftools view -e "${ifilters}" --types indels -O z -o filter-indels-$rundate.vcf.gz $mergevcf"
bcftools view -e "${ifilters}" --types indels -O z -o filter-indels-$rundate.vcf.gz $mergevcf 