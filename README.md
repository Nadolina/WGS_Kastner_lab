# Whole genome variant calling pipeline 

This README details the modules of the WGS pipeline developed for rare variant calling, the choices made in development, and how to use the various scripts. The modules are as such:
1. Revert the original BAM and re-align to GRCh38.
2. Cleanup the alignment to prepare it for variant calling.
3. Variant calling and filtering with GATK.
4. Variant calling and filtering with bcftools.
5. Separate annotation of GATK and bcftools calls with VEP.
6. Construction of hail databases for each GATK and bcftools.
7. Querying the database for samples, genes or other cohort characteristics.
8. Transferring data from your working directory to long-term storage 

Modules 3 and 4 can happen in parallel. 

<br />

<div align="center">
  <img src="https://github.com/user-attachments/assets/afd59b84-626c-4ffd-bb64-7e564739ec88"
">

The coloured boxes in this figure correspond to additional flowcharts describing subprocesses in the section 'The Pipeline'. 

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
VEP/112
python3

python packages:
argparse
hail
pandas
time
datetime

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

This pipeline is primarily run using a text file that contains the path to one original BAM file on each line. This list will be referred to as the batch. You can run individual BAMs through pre-process-pipeline.sh and alignment_cleanup.sh but a batch will be required from variant calling onwards. Here is an example snippet of one batch:

```
/data/Kastner_PFS/seq_data_storage/wgs_data/NISC_wgs_072023_dir/5915_dir/5915.bam
/data/Kastner_PFS/seq_data_storage/wgs_data/592/592.bam
/data/Kastner_PFS/seq_data_storage/wgs_data/NISC_wgs_072023_dir/5922/5922.bam
/data/Kastner_PFS/seq_data_storage/wgs_data/5924/5924.bam
```
You will also need to make sure the script is pointing to a directory containing the GRCh38 reference with decoys and alternate sequences (1). The necessary files are currently stored in /data/Kastner_PFS/references/HG38/. 

The original input to the pipeline was based on the user copying the BAMs of interest to their working directory in the prescribe format shown below. For legacy purpose, and more options for the user, I have decided to leave this functionality, even though I __strongly recommend__ just passing a textfile of paths pointing to the original BAMs. 

```
[sample]
  [sample].bam
  [sample].bai 
```


</details>

## The Pipeline 

Click on the drop down arrows to view a module. 

<details>
<summary>Module 1: Re-alignment to GRCh38</summary>

### Module 1: Re-alignment 

The process of extracting reads from an aligned BAM mostly adheres to the instructions in (2). This extracts the sequences, and removes alignment features, producing an unaligned BAM (uBAM). Then, we mark any adapters which may be artifacts from sequencing and can interfere with alignment. MarkIlluminaAdapters produces a new BAM with these demarcations, as well as a metrics file describing the adapters found, but for space, the demarcated BAM is deleted, and the metrics are retained. In the last step, I pipe the adapter-marked bam through SamToFastq, bwa-mem2 alignment to GRCh38 with decoys, and MergeBamAlignment to concatenate the unmapped reads back to the final bam. 

Samtools stats is run on both the original and final bam. Fastqc is run on the reverted bam to look at the quality of the sequence reads. 

<div align="center">
  <img src="https://github.com/user-attachments/assets/5df34f9f-9167-4ac3-8ec9-7742378c3a49">
</div>

```
batch --mem=[] --cpus-per-task=[] --gres=lscratch:[] pre-process-pipe.sh -b [original bam]
  OR
sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] pre-process-pipe.sh -l [location]

  The -l [locations] option will look in the folder for anything that matches *.bam, so ensure the bam of interest is the only *.bam in the folder provided.

  -l pass the path to the directory containing the bam
  -b pass the original bam file 
```
<br />
The -l option was implemented as an earlier solution to batching, when the program still required that you copy BAM files to your work space in the format shown below. I left this option available for anyone who may prefer to use this format. Basically you would structure your working directory with one directory per sample, containing the original bam. Then, your batch file would contain a list of paths pointing to these sample folders in your working directory, one per line. 

```
[sample]
  [sample].bam
  [sample].bai
```

If you want to run a batch with either -l or -b, you will need to loop through the batch textfile. Alternatively, you can run the program without looping, with just one BAM at a time. 
```
while read location; do sbatch --mem=48g --cpus-per-task=8 --gres=lscratch:400  $SCRIPTS/pre-process-pipe.sh -l $location ; done < HC-batch-091324.txt
while read sample; do sbatch --mem=48g --cpus-per-task=8 --gres=lscratch:400  $SCRIPTS/pre-process-pipe.sh -b $sample ; done < HC-batch-091324.txt
```

It is on my to-do list to wrap these scripts such as to accept a batch list as input without requiring looping. 

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

  -l location of sample working directory 
  -b merged bam alignment
  -o path to the original BAM 
  -r start script after markduplicates spark, no input just pass the flag; ie./ if your previous run fails but generated the *markdups_sort.bam correctly (OPTIONAL)
  -h help
```
Like with module 1, you can loop through a batch file or just run per usual with one bam or location. 

</details>

<details>
<summary>Module 3: Variant calling with GATK</summary>

### Module 3: Variant calling with GATK 

GATK variant calling mostly subscribes to the recommendations in GATK's documentation and tutorials (5). At this point, the pipeline becomes more parallelized. I generate Biowulf swarms to run HaplotypeCaller on each chromosome of each sample in parallel (A). I combine chromosome gVCFs produced by HaplotypeCaller across samples with GATK's CombineGVCFs, resulting in 24 gVCFs (B). Documentation recommends using GenomicsDB for this gVCF gathering step, but for ease of use I chose CombineGVCFs. The chromosome-combined gVCFs are then genotyped (C). 

After genotyping, we then need to filter the called variants. GATK has a program for this called Variant Quality Score Recalibration, which is akin to BQSR. Again using prior knowledge from curated datasets like dnsnp and 1000Genomes, and annotations in our VCF, VariantRecalibrator tries to model variant scores that are likely to be true variants. This model is data greedy, so we combine all the gentotyped CHRn-VCFs into a single VCF (D). Then we generate one VariantRecalibrator model for SNPs and another for indels (E,F), as recommended by GATK documentation (6). The models are applied back to the VCF to organize variants into tranches, effectively filtering them. The SNP model is applied first (G), and then the indel model (H), resulting in a fully variant quality score recalibrated VCF. 

<div align="center">
  <img src="https://github.com/user-attachments/assets/c7eea530-4b99-40a7-9576-b9506fbb4042">
</div>

<br />

```
sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] variant_calling_GATK.sh -b [batchfile]

  -b This is a textfile of the IDs (assuming you are following the prescribed directory structure, and you have working directories named with their IDs only), with one ID on each line.
  -o This is a textfile of the original bam paths, assuming you have bams in locations other than the working directory and created this file for use in pre-process-pipe.sh. One path per line.
  -h help

```

The -b is comparable to -l in previous modules. The -o requires a textfile of paths pointing to the original BAMs for a given sample, which the program will parse to identify the files it requires. You do not need to loop through the textfile as with previous modules. 

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

I also flag the "picked" variant according to their impact ranking. Please refer to (15) for more information on the --pick_flag. In addition to this I have also added a "CANONICAL" flag, which will indicate with just a "YES" if a given variant falls within a canonical transcript (11). 

After annotation, I apply a very simple filter of 'gnomADg_AF_joint < 0.2 or not gnomADg_AF_joint'. So, we only retain variants that have a maximum joint (gnomAD exome and genome) allele frequency of 0.2, or if they do not have a joint AF at all, because absense from gnomAD may indicate the variant is rare. I also apply --only_matched, which will only retain consequence blocks that pass the prescribed filters, if a variant has more than one block. This is more useful if a variant has multiple recorded consequences, i.e./ "intron_variant&non_coding_transcript_variant" and "upstream_gene_variant", and one of your filters is to retain only intron variants. Then, the --only_matched would remove the upstream_gene_variant block. Allelic frequenices for a given variant do not change between blocks, so this does not really apply but I have included it nonetheless. 

After annotation and filtering you will have a *vcf and *filtered.vcf for each chromosome for a batch of samples. 

This VEP script just takes a VCF as input, so it can be applied to one or more samples from any variant caller. I pass batch VCFs from both GATK and bcftools to run-VEP.sh. VEP does not have sophisticated threading, parallelization or lscratch space usage, so I recommend only allocating 10-15g of memory across the biowulf default of 2 threads. 

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

Even though the MatrixTable can accomodate much larger quanities of data than we store in the Kastner lab, I chose to split the database into four groups of chromosomes: chr1-4, chr5-10, chr11-17 and chr18-22,X,Y. This was mostly to improve the speed for querying the database for a gene, although it does not improve the speed of a sample query. In querying for either genes or samples, the isolation of that specific gene or sample(s) is fairly quick. The resulting query is temporarily stored as a new matrix table. However, I add more annotations (pLI, LOEUF, mis-Z-score, bcftools calls, etc.) to the query, and also perform some reformating for easier filtering when the output TSV is view in excel. 

Batch VCFs need to be modified before building these four hail databases. To perform these modifications, I have written build_hail_2.sh and build_hail_bcftools.sh. These are a little redundant, and will be consolidated at some point, per my laundry list. I only wrote them separately because the naming scheme for bcftools files and bcftools annotations are different from GATK. The scripts:
1. accept a list of folders containing batch annotated VCF files
2. merge chromosomes across batches annotated VCFs (i.e./ 1 with 1, 2 with 2, etc.)
3. convert consequence blocks into a Hail compatible format (bcftools csq-split pluggin)
4. remove unnecessary annotations (mostly from variant calling)
5. concatenate formatted chromosome VCFs into the four prescribed sub-database VCFs

The script will submit three series of swarm commands, one to perform any necessary compression and indexing, then to perform the CSQ-splitting and the last for the concatenation. The result will be a folder, named buildhail_[rundate] for build_hail_2.sh (for annotated GATK calls) or buildhail_bcftools_[rundate]. Both will contain:

    chr1.split.vcf.gz
    chr1.split.vcf.gz.csi
    chr2.split.vcf.gz
    chr2.split.vcf.gz.csi
    ...
    chrN.split.vcf.gz
    chrN.split.vcf.gz.csi 
    db1.concat.[rundate].vcf.gz
    db2.concat.[rundate].vcf.gz
    db3.concat.[rundate].vcf.gz
    db4.concat.[rundate].vcf.gz 

```
sbatch [OPTIONS] build_hail_2.sh [batch annotation folder 1] [batch annotation folder 2] ... [batch annotation folder n]

Example using prescribed pipeline folders:
  * running from my wgs working directory

sbatch --mem=10g -c 6 --gres=lscratch:200 --time=04:00:00 /data/Kastner_PFS/scripts/pipelines/WGS_Kastner_lab/build_hail_2.sh \
  filtered_VCFs_011425_2/indel.SNP.recalibrated_99.9.011425-VEP-VCFs \
  filtered_VCFs_011425/indel.SNP.recalibrated_99.9.011425-VEP-VCFs \
  filtered_VCFs_012525/indel.SNP.recalibrated_99.9.012525-VEP-VCFs \
  filtered_VCFs_110824/indel.SNP.recalibrated_99.9.110824-VEP-VCFs

You can pass as many folders of annotated VCFs as you'd like. 
```

__NOTE__ that these "db" VCFs are not the final database. But, they are the final VCF from which the four hail databases are constructed. To construct the Hail databases, you must: 
1. Make a new directory in /data/Kastner_PFS/WGS/cohort_db of the name version_[rundate]. You can see version_031925 as an example. In that example you will see the 4 VCFs as well as 4 folders with the same name but a '.mt' suffix instead of '.vcf.gz'. These are the matrix tables (sub-databases) constructed from the four VCFs.
2. Move the VCFs to your new directory.
3. Assuming the batches you passed to build_hail_2.sh and build_hail_bcftools.sh have not been included in previous versions of the database, merge the previous versions' VCFs with yours. Once merged, you will have your four cohort representive VCFs from which you can finally build your matrix tables. For example:

```
bcftools merge -m none -Oz -o [merged VCF] db1.concat.031725.vcf.gz db1.concat.[rundate].vcf.gz
bcftools merge -m none -Oz -o [newname] db2.concat.031725.vcf.gz db2.concat.[rundate].vcf.gz
etc.
```
  
4. Start an sinteractive session and allocation at least 40g of memory and 16 cores.
5. Load the hail module and start an interactive python session.

```
module load hail
ipython
```
6. Import and initiate hail with memory.

```
import hail as hl
hl.init(spark_conf={'spark.driver.memory': '25g'})
```
7. For each of the four VCFs, perform the following command to import the VCF and save it into a matrix table. 

```
hl.import_vcf('[merged VCF]',force_bgz=True,reference_genome='GRCh38',array_elements_required=False).write('db1.concat.[rundate].mt',overwrite=True)
```
8. repeat for bcftools VCFs in a different folder in cohort_db, called something like 'bcftools_rundate'.
9. I unfortunately have the python hail query scripts (hail-gene-query.py and hail-sample-query.py) hard-coded with these database names, so edit them with the new matrix table file paths here:

   ![image](https://github.com/user-attachments/assets/2d0ff666-00df-4b08-9762-3fda759888ff)

Use 'quit()' to exit the interactive python session. 

At this point you should have eight matrix tables (between GATK and bcftools) ready for use independently or with the query scripts provided! 

</details>

<details>
<summary>Module 7: Querying the hail database</summary>

### Module 7: Querying the Hail database 

Hail has extensive functionality for querying matrix tables containing genomic data. I have coded two python programs for sample and gene querying, but I strongly encourage you to review their interactive python documentation (17). The programs are hail-gene-query.py and hail-sample-query.py, which are located in /data/Kastner_PFS/scripts/pipelines/WGS_Kastner_Lab. I also recommend consulting their forum for any questions you may have; there is a very active community of users and Hail engineers troubleshooting a variety of problems on the forum (18). 

Hail-gene-query.py will collect any variants in the gene of interest, and aggregate a list of heterozygous and alt-homozygous samples as called by GATK and/or bcftools. The script adds some further annotation, and reformats the Hail table into a pandas dataframe that can be easily exported to a TSV with a name 'hail-[GENE SYMBOL]-[DATE].tsv'. This should only take a few minutes to collect. 

The second script is hail-sample-query.py. __NOTE__ that this sample query only outputs 'MODERATE' and 'HIGH' impact variants, per VEP's classification. This output can be quite substantial, and so I chose to apply this preliminary filtering. This adds the same annotations, and will again aggregate across samples at each variant site to list the heterozygous and alt-homozygous samples called by GATK and/or bcftools. This will of course be limited to just the samples you pass. You can pass multiple samples to the program as a space-delimited list. Since the script will go through each sub-database, as opposed to just one for a specific gene, this program does take a bit more time. 

I wanted to implement a way to further filter the hail-sample-query.py outputs. So, I coded an additional parameter -n (--n_non_ref), which allows the user to tell the program to only retain variants where the sum of non-reference genoetypes between the prescribed samples is equal to or greater than a given number. For example, if you have a trio, with a known effected parent and patient, you could pass -n 2, so that you only look at sites where at least two of the three samples contain alternate genotypes. You need to pass a minimum of -n 1. 

```
module load hail

python3-hail hail-gene-query.py -g [GENE SYMBOL, i.e/MEFV]

python3-hail hail-samply-query.py -s [SAMPLE ID 1] [SAMPLE ID 2] ... [SAMPLE ID N] -n [number of non-ref genotypes]
ex.  python3-hail hail-samply-query.py -s 1122 1123 1124 -n 2

```

</details>

<details>
<summary>Module 8: Transferring completed VCFs and BAMs to long-term storage</summary>
  
### Module 8: Transferring completed VCFs and BAMs to long-term storage 

Storing our re-aligned samples and variants is an essential final step in this workflow. All data from modules 1-5 should be in your working directory (i.e: /data/[USERNAME]/working_dir/). Some of the files you produced need to be retained. See the table below for a list of retained files and their storage destination. 

| File | Destination | Path |
|------|-------------|------|
|[SAMPLE]_bqsr_1.bam | ketu | /mnt/[USERNAME]/ketu_labs/Kastner/WGS/[YEAR ORIGINAL DATA PRODUCED]/[SAMPLE] |




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
3. Consolidate build_hail_2.sh and build_hail_bcftools.sh. 

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
17. https://hail.is/docs/0.2/api.html
18. https://discuss.hail.is/


