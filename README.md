# Whole genome variant calling pipeline 

This README details the four modules of the WGS pipeline developed for rare variant calling, the choices made in development, and how to use the various scripts. The modules are as such:
1. Revert the original BAM and re-align to GRCh38.
2. Cleanup the alignment to prepare it for variant calling.
3. Variant calling and filtering with GATK.
4. Variant calling and filtering with bcftools.

Modules 3 and 4 can happen in parallel. 

<br />

<div align="center">
  <img src="https://github.com/user-attachments/assets/4fa4b648-50c5-44e8-8d72-dd2d4d546bca"> 
</div>

## Dependencies

```
picard/3.2.0
bwa-mem/2 2.2.1
fastqc/0.12.1
samtools/1.19
```

## Where do we start?

```
[sample]/
  [sample].bam
  [sample].bam.bai
  [sample].filtered.vcf.gz
  [sample].filtered.vcf.gz
  [sample].g.vcf.gz
  [sample].g.vcf.gz.tbi 
```

Most of our files come from NISC. We typically get three files of importance (the bam, and the two vcf.gz files), along with their indexes, which are usually necessary for performing downstreams analysis but can be recreated easily. Some people in our group use the gVCF or VCF for annotation and variant interpretation. Others perform their own variant calling (generating their own gVCF and VCF) from the BAM. The problem with the majority of our BAM files, is that they are aligned to GRCh37. This pipeline is an effort to re-align everything to GRCh38, which is a newer and more accurate reference genome. Additionally, the goal of this pipeline is to internally standardize variant calling workflows. _So, we need to start by extracting the sequence reads from the aligned BAM providied by NISC_.


## Inputs 

1. The [sample].bam and its index, [sample].bam.bai
2. GRCh38 with decoys and other accompanying files (/data/Kastner_PFS/references/HG38/) (1)

With the current pipeline configuration, I recommend setting up your working directory as such: 

```
[sample 1]/
  *bam
  *bam.bai
[sample 2]/
  *bam 
  *bam.bai 
[sample n]/
  *bam 
  *bam.bai
batch-[rundate].txt
```
Where the batch text file contains one sample ID per line, which corresponds to their directory names


## The pipeline 

### Module 1: Re-alignment 

The process of extract reads from an aligned BAM mostly adheres to the instructions in (2). This extracts the sequences, and removes alignment features, producing an unaligned BAM (uBAM). Then, we mark any adapters which may be artifacts from sequencing and can interfere with alignment. MarkIlluminaAdapters produces a new BAM with these demarcations, as well as a metrics file describing the adapters found, but for space, the demarcated BAM is deleted, and the metrics are retained. In the last step, I pipeline the adapter-marked bam through SamToFastq, bwa-mem2 alignment to GRCh38 with decoys, and MergeBamAlignment to concatenate the unmapped reads back to the final bam. 

Samtools stats is run on both the original and final bam. Fastqc is run on the reverted bam to look at the quality of the sequence reads. 

<div align="center">
  <img src="https://github.com/user-attachments/assets/5df34f9f-9167-4ac3-8ec9-7742378c3a49">
</div>

```
batch --mem=[] --cpus-per-task=[] --gres=lscratch:[] pre-process-pipe.sh -b [original bam]
  OR
sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] pre-process-pipe.sh -l [location]

  You need to pass EITHER -b or -l, but NOT BOTH.
  The -l [locations] option will look in the folder for anything that matches *.bam, so ensure the bam of interest is the only *.bam in the folder provided.

  -l pass the path to the directory containing the bam; allows user to loop through a text file containing locations (like a batch)
  -b pass the bam file 
```
<br />
The -l option was implemented here as a temporary solution to passing batches or pointing to BAMs in another directory besides the working directory. For the -l to work on a batch, you currently need to loop through a textfile with one location per line. It will work best if you use the prescribed working directory structure in "Inputs", where sample directories are named for sample IDs. This is because the -l looks for anything named *bam in the directory provided. For example: 

```
while read loc; do sbatch --mem=129g --cpus-per-task=16 --gres=lscratch:400  $SCRIPTS/pre-process-pipe.sh -l $loc ; done < HC-batch-091324.txt
```

It is on my to-do list to wrap these scripts such as to accept a batch list as input. 

### Module 2: Alignment cleanup 

This module follows instructions from GATK's pre-variant calling recommendations (4). Alignment cleanup involves marking duplicate read pairs and base quality score recalibration. Duplicate read pairs can arise from the same DNA fragments in sequencing, and unless marked variant callers may consider them as independent sequences, which would be inaccurate. Here we fully remove duplicate read pairs. 

Base recalibration models the base quality scores using prior knowledge of variant sites, to identify patterns of sequencing bias. The model is applied to the bam to correct these biases. GATK documentation then recommends running a second round of BQSR for quality assurance. The second model isn't applied, but instead we compare the models with CovariateAnalysis, which produces a pdf summary of the results. The documentation does not explain what the QC value is, but I assume it is to confirm that a the first model was adequate for correction and that a second would be of little benefit, or in fact introduce over-correction. 

<div align="center">
  <img src="https://github.com/user-attachments/assets/0cd58f81-0d37-46cc-a711-7d775abab496">
</div>

<br />

```
sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] --time=days-hours:minutes:seconds alignment_cleanup.sh [arguments]

  -l pass the path to the directory containing the bam; allows user to loop through a text file containing locations (like a batch)
  -b merged bam alignment
  -r start script after markduplicates spark, no input just pass the flag; ie./ if your previous run fails but generated the *markdups_sort.bam correctly (OPTIONAL)
  -h help
```
The -l here works in the same way as in module 1, where you have to loop through the textfile. 

### Module 3: Variant calling with GATK 

GATK variant calling mostly subscribes to the recommendations in GATK's documentation and tutorials (5). At this point, the pipeline becomes more parallelized. I generate Biowulf swarms to run HaplotypeCaller on each chromosome of each sample in parallel. Once the haplotype called gVCFs are done, I combine chromosome gVCFs across samples with GATKs CombineGVCF, resulting in 24 gVCFs. Documentation recommends using GenomicsDB for this gVCF gathering step, but for ease of use I chose CombineGVCFs. The chromosome-combined gVCFs are then genotyped. 

After genotyping, we then need to filter the called variants. GATK has a program for this called Variant Quality Score Recalibration, which is akin to BQSR. Again using prior knowledge from curated datasets like dnsnp and 1000Genomes, and annotations in our VCF, VariantRecalibrator tries to model variant scores that are likely to be true variants. This model is then applied back to the VCFs to divide variants into likelihood "tranches". 

<div align="center">
  <img width=900 src="https://github.com/user-attachments/assets/8545bc01-8023-477e-a8dc-f161ff941500">
</div>

## Cumulative outputs

```
[sample]/
  [sample]_out/
    [sample].reverted.bam
    [sample].mergedaln.bam
    [sample].mkAdapter.metrics.txt
    [sample].originalaln.stats
    [sample].mergedaln.stats
  qc_out/
    [sample].reverted_fastqc.html
    [sample].reverted_fastqc.zip
    
```
##  Batch size and memory allocation 

## Laundry list 

1. Implement a batching solution for modules one and two. 

## References 

1. https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
2. https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM
3. https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
4. https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
5. https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels


