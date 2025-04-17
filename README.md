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


<details> 

<summary>
Getting started with the piepline - click on the arrow dropdown to review dependencies and inputs for this pipeline. 
</summary>

## Dependencies

```
picard/3.2.0
bwa-mem/2 2.2.1
fastqc/0.12.1
samtools/1.19
GATK/4.6.0.0
bcftools/1.19
R
```

**You will also need the dependencies managed by a renv that can be found here: https://github.com/Nadolina/WGS-renv.git.** If you plan to copy the pipeline to another location, please follow the recommendations in (10) to clone the renv. You would then also need to modify the path to the renv in the Rmd. If you are in the Kastner group and using the shared installation of the pipeline, you do not need to clone the renv as there is already a shared clone. 

If you so choose, you can install the R packages manually instead, but this does not guarantee versioning. You will need: 

```
ggplot2
vcfR
readr
tidyverse
data.table
dplyr
kableExtra
stringr
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

</details>

## The pipeline 

Click on the drop down arrows to view a module. 

<details>
<summary>Module 1: Re-alignment to GRCh38</summary>

### Module 1: Re-alignment 

The process of extracting reads from an aligned BAM mostly adheres to the instructions in (2). This extracts the sequences, and removes alignment features, producing an unaligned BAM (uBAM). Then, we mark any adapters which may be artifacts from sequencing and can interfere with alignment. MarkIlluminaAdapters produces a new BAM with these demarcations, as well as a metrics file describing the adapters found, but for space, the demarcated BAM is deleted, and the metrics are retained. In the last step, I pipeline the adapter-marked bam through SamToFastq, bwa-mem2 alignment to GRCh38 with decoys, and MergeBamAlignment to concatenate the unmapped reads back to the final bam. 

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
The -l option was implemented here as a temporary solution to passing batches, but the bam still needs to be copied into your working directory (the location) because the script moves into that directory. For the -l to work on a batch, you currently need to loop through a textfile with one location per line. It will work best if you use the prescribed working directory structure in "Inputs", where sample directories are named for sample IDs. This is because the -l looks for anything named *bam in the directory provided. For example: 

```
while read loc; do sbatch --mem=129g --cpus-per-task=16 --gres=lscratch:400  $SCRIPTS/pre-process-pipe.sh -l $loc ; done < HC-batch-091324.txt
```

It is on my to-do list to wrap these scripts such as to accept a batch list as input. 

</details>

<details>
<summary>Module 2: Alignment cleanup</summary>

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

</details>

<details>
<summary>Module 3: Variant calling with GATK</summary>

### Module 3: Variant calling with GATK 

GATK variant calling mostly subscribes to the recommendations in GATK's documentation and tutorials (5). At this point, the pipeline becomes more parallelized. I generate Biowulf swarms to run HaplotypeCaller on each chromosome of each sample in parallel (A). I combine chromosome gVCFs produced by HaplotypeCaller across samples with GATKs CombineGVCFs, resulting in 24 gVCFs (B). Documentation recommends using GenomicsDB for this gVCF gathering step, but for ease of use I chose CombineGVCFs. The chromosome-combined gVCFs are then genotyped (C). 

After genotyping, we then need to filter the called variants. GATK has a program for this called Variant Quality Score Recalibration, which is akin to BQSR. Again using prior knowledge from curated datasets like dnsnp and 1000Genomes, and annotations in our VCF, VariantRecalibrator tries to model variant scores that are likely to be true variants. This model works is data greedy, so we combined all the gentotyped CHRn-VCFs into a single VCF (D). Then we generate one VariantRecalibrator model for SNPs and another for indels (E,F), as recommended by GATK documentation (6). The models are applied back to the VCF to organize variants into tranches, effectively filtering them. The SNP model is applied first (G), and then the indel model (H), resulting in a fully variant quality score recalibrated VCF. 

<div align="center">
  <img src="https://github.com/user-attachments/assets/c7eea530-4b99-40a7-9576-b9506fbb4042">
</div>

<br />

```
sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] variant_calling_GATK.sh -b [batchfile]

  -b textfile of locations of directories containing original BAM; one per line 
  -h help
```

As with previous modules, the -b mostly anticipates the structure prescribed in "Inputs", where each directory is named for it's sample ID. Unlike previous modules, -b can take the whole textfile as an input, rather than having to loop through it on the command line.

</details>

<details>
<summary>Module 4: Variant calling with bcftools </summary>

### Module 4: Variant calling with bcftools

<img align="right" src="https://github.com/user-attachments/assets/d4e7ec8d-7904-4bd1-b688-49e314005de4">

Lierature suggests using two or more variant callers, because there are several highly accurate SNP and indel callers available, and concordance between multiple callers lends confidence to  calls. In addition to GATK, this workflow performs variant calling with bcftools (7). This module also starts with the base recalibrated BAMs from the batch. Using the swarm functionality again, we run bcftools mpileup per chromosome, across samples, producing 24 VCFs. These VCFs are annotated with known alleles frequencies from the 1000 Genomes project, because the --prior-freqs flag is used in the bcftools call command to improve calling performance. 

<br />

Unlike GATK, bcftools does not have a model to perform filtering. There are various approaches to filtering, but most tutorials and documentation recommend a hard-filtering approach. My approach was largely informed by the recommendations in (8,9). In nearly all the VCFs we have received from NISC, you can see the same set of hard filters have been applied to SNPs and indels. These are very common filtering parameters and values. 

  ![image](https://github.com/user-attachments/assets/52d42b91-a062-4e3f-a992-4e35288c556b)

I wanted to make our filtering workflow a little more adaptable. Instead of applying the same filtering thresholds to all of our batches, I aim to identify thresholds that are tailored to the batch, to try to mitigate changes and biases in the batch. I take a subset (1%) of SNPs and indels from each batch, and extract their variant metrics. I identify the extreme outlier values in the subset using just quantiles, and apply those back to the whole batch of variants as filtering thresholds. For SNPs, we use the following metrics: variant quality, depth and mapping quality Z-score (MQBZ). MQBZ looks at the differences in mapping quality at heterozygous sites, because significant differences at the ALT allele could indicate a false positive. For indels, quality and depth are also used, as well as IMF, which is the fraction of reads supporting an indel.  

In the example below, you can see the range of quality values varies between Batch 1 and 2, with Batch 1 exhibiting on average higher quality values. If we were apply the same quality filter to both batches, we would remove more variants in Batch 2 than 1, and there is an increased chance we are removing true variants. The motivation is to try to accomodate those metric variations between batches, to retain as many true variants as possible. While this does mean the filtering thresholds won't be standardize and this can complicate methods reporting a bit, I expect that metrics should hold relatively steady between batches and so thresholds should be similar. This approach is more of a contingency than anything. 

<div align="center">
  <img src="https://github.com/user-attachments/assets/62bd4413-f982-483e-96db-c898d5c36ae1">
</div>

<br />

Because this approach varies with each batch, I wanted to ensure we had QC outputs describing the thresholds identified. So, the filtering program produces a report on the subset statistics, in the form of an HTML file so it can be viewed in any browser. Below is an example of one of the figures in the HTML report. 

<br />

<div align="center">
  <img src="https://github.com/user-attachments/assets/5c433fa2-e04e-44f8-9749-027a0f671d5d">
</div>

<br />

```
sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] bcftools.sh -b [batch files with sample IDs]

-b  the same batch file as passed to GATK, one sample ID/location per line 
```

<br />


</details>

<details>
<summary>Module 5: Annotation with VEP </summary>

### Module 5: Annotation with VEP

This VEP annotation (11) script simply accepts a VCF as input, which means it can be run on any VCF that you would like. The script will:
1. check for compression and indexing and perform both if needed
2. sort and normalize the VCF according to recommendations from ensembl and other users (13,14)
3. generates a swarm for to annotate chromosomes in parallel
4. performs VEP filtering on each chromosome

Annotations employed include:
* HGVS
* CADD
* gnomADg and gnomAD joint
* SpliceAI
* MaxEntScan
* REVEL
* AlphaMissense
* ClinVar
* GERP

I also flag the "picked" variant according to their impact ranking. Please refer to (15) for more information on this --pick_flag. In addition to this I have also added a "CANONICAL" flag, which will indicate with just a "YES" if a given variant falls on a canonical transcript (11). 

After annotation, I apply a very simple filter of 'gnomADg_AF_grpmax < 0.1 or not gnomADg_AF_grpmax'. This just means I only retain variants that have a maximum population allele frequency (grpmax) of 0.1, or if they do have a grpmax AF at all, because absense from gnomAD may indicate the variant is rare. I also apply --only_matched, which will only retain consequence blocks that pass the prescribed filters, if a variant has more than one block. This is more useful if a variant has multiple recorded consequences, i.e./ "intron_variant&non_coding_transcript_variant" and "upstream_gene_variant", and one of your filters is to retain only intron variants. Then, the --only_matched would remove the upstream_gene_variant block. Allelic frequenices for a given variant do not change between blocks, so this does not really apply but I have included it nonetheless. 

After annotation and filtering you will have a *vcf and *filtered.vcf for each chromosome for a batch of samples. 

This VEP script just takes a VCF as input, so it can be applied to one or more samples from any variant caller. I pass batch VCFs from both GATK and bcftools to run-VEP.sh. VEP does not have sophisticated threading, parallelization or lscratch space, so I recommend only allocating 10-15g of memory across the biowulf default of 2 threads. 

```
sbatch --mem=[] --cpus-per-task=[] run-VEP.sh -v [VCF]

  -v VCF
  -h help
```

</details>

<details>
<summary>Module 6: Building the Hail database </summary>

### Module 6: Building the Hail database 

Hail is a database software designed specifically for managing genomic data, because standard dataframe options like SQL are not as sophisticated for this work (16). 1000 genomes, gnomAD and the UK Biobank among other large consortia employ Hail. Users can import a multi-sample VCF into Hail and construct a matrixtable, which can represent varying numbers of records and identifiers. Hail has extensive functionality for querying VCF data, aggregating statistics and performing annotation. 

Even though the MatrixTable can accomodate much larger quanities of data than we store in the Kastner lab, I chose to split the database into four groups of chromosomes: chr1-4, chr5-10, chr11-17 and chr18-22,X,Y. 


</details>

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

## Laundry list 

1. Implement a batching solution for modules one and two. The -l in modules 1 and 2 serves as a bit of a work-around to batching. The problem with the current set up is you need to copy the bam and pre-make working directories for each sample. I think I will wrap this such that it a) allows the user to point to a bam not in the working directory and b) has a second mandatory flag for the sample ID, so that creation of the sample ID-named working directory is standardized. 
2. Parallelize modules 1 and 2. I think we could break alignment by chromosome or read group and carry that break through base recalibration.
3. I created a renv for the rmarkdown script at the end, which is just in my own directories and is still linked for now, but it should be moved to the Kastner scripts directory. 

## References 

1. https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
2. https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM
3. https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
4. https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
5. https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels
6. https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR
7. https://samtools.github.io/bcftools/bcftools.html
8. https://www.htslib.org/workflow/filter.html
9. https://speciationgenomics.github.io/filtering_vcfs/
10. https://rstudio.github.io/renv/articles/renv.html
11. https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic
12. https://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html
13. https://bioinformatics.stackexchange.com/questions/22124/variants-from-multiple-tools-normalization-before-or-after-annotation-with-vep
14. https://www.ensembl.info/2020/05/26/normalising-variants-to-standardise-ensembl-vep-output/
15. https://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options
16. https://blog.hail.is/introtohail/


