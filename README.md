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

</details>

## The Pipeline 

:arrow_down: **Click on the drop down arrows to view a module.** :arrow_down:

<details>
<summary>Module 1: Re-alignment to GRCh38</summary>

### Module 1: Re-alignment 

Re-alignment generally adheres to GATK best-practices guidelines (1,2,3,4). You pass your single original BAM or batch file of BAMs (as describe in **Inputs**) to submit-pre-process-pipe.sh. This will then either submit an sbatch job for running pre-process-pipe on your single BAM or generate a swarm (5) of pre-process-pipe.sh jobs. A swarm just deploys multiple similar jobs in parallel; in this case the same program with different input BAMs. 

Pre-process-pipe.sh starts by using RevertSam to extract reads from the GRCh37 aligned bam, and removes alignment features, producing an unaligned BAM (uBAM). We use an option in the RevertSam tool to output uBAMs for each read-group. This allows for parallelization of the next steps, for efficiency. The per-read-group uBAMs are then flowed through MarkIlluminaAdapters, SamToFastq and bwa-mem2, producing a GRCh38 aligned BAM (for a given read group). In the last step, I merge all the re-aligned per-read-group BAMs back into a single output; this will be named "[SAMPLE].fullmerged.bam"

<div align="center">
  <img width="3498" height="1414" src="https://github.com/user-attachments/assets/35527034-b8fd-4f88-ba1d-9dbbded07cdd">
</div>

###### Figure 2. Flowchart for re-alignment. Lighter yellow circles are inputs/intermediates. The dark yellow circle is the output. The yellow rectangles are the wrapper scripts I designed. Note that for a batch, the pre-process-pipe.sh will be running many times in unison, via a swarm job (5). 

```
sbatch --mem=5g --gres=lscratch:25 submit-pre-process-pipe.sh -b [original bam]

  -h help
```

Note that computational resource allocations are hard-coded into pre-process-pipe.sh. **Memory, CPUs and lscratch are not especially necessary for the submit script, but it does need to be run with sbatch**.

</details>

<details>
<summary>Module 2: Alignment cleanup</summary>

### Module 2: Alignment cleanup 

This module follows instructions from GATK's pre-variant calling recommendations (4). The steps are:
1) MarkDuplicatesSpark
2) BaseRecalibrator (1)
3) ApplyBQSR
4) BaseRecalibrator (2)
5) AnalyzeCovariates

The goal of alignment cleanup is to optimize an aligned BAM for variant calling by mitigating artifacts. Duplicate read pairs can arise from the same DNA fragments in sequencing, and unless marked variant callers may consider them as independent sequences, which would be inaccurate. __We fully remove duplicate read pairs.__

Base recalibration models the base quality scores using prior knowledge of variant sites, to identify patterns of sequencing bias. The model is applied to the bam to correct these biases. GATK documentation then recommends a second round of BQSR, which is compared to the first BQSR model using CovariateAnalysis, for quality assurance. 

MarkDuplicatesSpark is run on the whole BAM because it needs the whole file for accuracy, but the remaining steps are split by chromosome and run in parallel for efficiency. The per-chromosome base recalibrated BAMs are aggregated into a **final output with the naming convention {sample}_bqsr_1.bam**

To run this module, pass the path for an __original BAM__ to the submit script, or pass a __text file containing one path per line to run a batch__, as such: 

```
sbatch --time=days-hours:minutes:seconds submit_alignment_cleanup.sh [arguments]

  -b original BAM path
  -f batch text file, one original BAM path per line 
  -h help
```

</details>

<details>
<summary>Module 3: Variant calling with GATK</summary>

### Module 3: Variant calling with GATK 

GATK variant calling mostly subscribes to the recommendations in GATK's documentation and tutorials (5,6). The script, variant_calling_GATK.sh, generates swarms to run HaplotypeCaller on each chromosome of each sample in parallel. A second swarm combines gVCFs produced by HaplotypeCaller on a per-chromosome basis across samples with GATK's CombineGVCFs, resulting in 24 gVCFs. The chromosome-combined gVCFs are then genotyped. 

After genotyping, we then need to filter the called variants. GATK has a program for this called Variant Quality Score Recalibration, which is akin to BQSR. Again using prior knowledge from curated datasets like dnsnp and 1000Genomes, and annotations in our VCF, VariantRecalibrator tries to model variant scores that are likely to be true variants. All genotyped per-chromosome VCFs are combined before being subject to recalibration. Then we generate one VariantRecalibrator model for SNPs and another for indels. The models are applied back to the VCF to organize variants into tranches, effectively filtering them. The SNP model is applied first result in a fully variant quality score recalibrated VCF of SNPs and indels in our batch. 

<div align="center">
  <img width="3518" height="976" src="https://github.com/user-attachments/assets/23701839-efda-4d31-9b2a-fa5ca9493e3b">
</div>

<br />

```
sbatch --mem=[] --cpus-per-task=[] --gres=lscratch:[] variant_calling_GATK.sh -o [batchfile]

  -o This is a textfile of the original bam paths, assuming you have bams in locations other than the working directory and created this file for use in pre-process-pipe.sh. One path per line.
  -h help

```

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

Even though the MatrixTable can accomodate much larger quanities of data than we store in the Kastner lab, I chose to split the database into four groups of chromosomes: chr1-4, chr5-10, chr11-17 and chr18-22,X,Y. This was mostly to improve the speed for querying the database because it allows the query to be parallelized. 

Batch VCFs need to be modified before building these four hail databases. To perform these modifications, I have written build_hail_2.sh and build_hail_bcftools.sh. These are a little redundant, and will be consolidated at some point, per my laundry list. Bcftools and GATK variants are stored in separate hail databases because they have discongruant annotations. When performing queries with either submit-sample-query.sh or hail-gene-query.py, variants from both callers are consolidated in the outputs. The scripts:
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
    db1.*.[rundate].vcf.gz
    db2.*.[rundate].vcf.gz
    db3.*.[rundate].vcf.gz
    db4.*.[rundate].vcf.gz 

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
1. Make a new directory in /data/Kastner_PFS/WGS/cohort_db of the name version_[rundate]. You can see version_091725 as an example. In that example you will see the 4 VCFs as well as 4 folders with the same name but a '.mt' suffix instead of '.vcf.gz'. These are the matrix tables (sub-databases) constructed from the four VCFs.
2. Move the VCFs to your new directory from buildhail_*.
3. You need to combine the new *concat VCFs with the combined VCFs from the previous update. I use bcftools merge, but this will send you an error if there are duplicate samples. There may be duplicate runs on a given sample if say, a new family member was sequenced, and we wanted to variant call them together. Start by identifying duplicate samples with something like:

```
bcftools query -l [NEW DB1] > samples.txt
bcftools query -l [OLD COMBINED DB1] >> samples.txt
sort samples.txt | uniq -cd
(or some variation of this)
```
4. Then, remove the samples from the OLD VCF. 
   
6. Assuming the batches you passed to build_hail_2.sh and build_hail_bcftools.sh have not been included in previous versions of the database, merge the previous versions' VCFs with yours. Once merged, you will have your four cohort representive VCFs from which you can finally build your matrix tables. For example:

```
bcftools merge -m none -Oz -o db1.gatk.[DATE].vcf.gz [previous version db1 vcf] [new db1 vcf]
bcftools merge -m none -Oz -o db1.gatk.[RUNDATE].vcf.gz db1.gatk.091725.vcf.gz  db1.concat.[rundate].vcf.gz
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

<br />

<div align="center">
  <img width="3162" height="980" alt="Decision Tree Team Whiteboard in Burgundy Peach Simple Colorful Style (5)" src="https://github.com/user-attachments/assets/2f7608a8-9d85-4718-99a8-71f10a5f707c" />
</div>

<br />

</details>

<details>
<summary>Module 8: Transferring completed VCFs and BAMs to long-term storage</summary>
  
### Module 8: Transferring completed VCFs and BAMs to long-term storage 

Storing our re-aligned samples and variants is an essential final step in this workflow. All data from modules 1-5 should be in your working directory (i.e: /data/[USERNAME]/working_dir/). Some of the files you produced need to be retained. See the table below for a list of retained files and their storage destination. 

| File | Destination | Path |
|------|-------------|------|
|[SAMPLE]_bqsr_1.bam | ketu | /mnt/[USERNAME]/ketu_labs/Kastner/WGS/[YEAR ORIGINAL DATA PRODUCED]/[SAMPLE] |
|original BAM | ketu | /mnt/[USERNAME]/ketu_labs/Kastner/WGS/[YEAR ORIGINAL DATA PRODUCED]/[SAMPLE] |
| sample GATK VCF (unannotated) | biowulf | /data/Kastner_PFS/WGS/[YEAR ORIGINAL DATA PRODUCED]/[SAMPLE] |
| sample bcftools VCF (unannotated) | biowulf | /data/Kastner_PFS/WGS/[YEAR ORIGINAL DATA PRODUCED]/[SAMPLE] |
| original VCFs | biowulf | /data/Kastner_PFS/WGS/[YEAR ORIGINAL DATA PRODUCED]/[SAMPLE] |

To perform these transfers for a given batch use data_organization.sh. You must pass the batch-[DATE].txt file, and you can optionally pass the __unannotated__ batch GATK and bcftools VCFs. The tool will run if you pass just the batch text file, or the batch file and just one of the batch VCFs. It will generate appropriate outputs from whatever files are provded. I use bcftools view with the -s flag to extract a sample VCF from a given batch VCF. The program produces 1) a globus transfer file with the naming scheme globus-transfer-[CURRENT DATE]-[BATCH NAME].txt, 2) a swarm file with the name vcf_transfer_[BATCH NAME].swarm, and 3) a text file named orig_vcf_transfer_[DATE]-[BATCHNAME].txt. To avoid data loss, I do not initite these transfers from within the program. Instead, I __strongly encourage__ the user to review all three transfer files for errors before manually submitting them. Run the program and initiate the subsequent transfers as below: 

```

sbatch [parameters] globus_transfer.sh -o [batch] -g [GATK VCF (unannotated)] -b [bcftools VCF (unannotated)]

-o batch file of paths to original bams 
-g gatk batch VCF
-b bcftools batch VCF
-h help

```

Use the globus endpoint search to find the globus origin and destination IDs for biowulf and ketu. The origin endpoint ID (biowulf) should be e2620047-6d04-11e5-ba46-22000b92c6ec. Use the flag --dry-run to print globus submission output, as an additional check before actually submitting. Remove the flag to actually submitt the transfer. 

Submit the swarm to transfer new and original sample VCFs from your working director to Kastner_PFS. I use the -b flag to group a few bcftools view commands onto the same swarm job, because a single command is very quick, and grouping commands allows for more efficient use of resources. 

```

globus endpoint search NIH
globus endpoint search nhgrishell02
globus transfer --no-verify-checksum --batch globus-transfer-${rundate}-${batch_prefix}.txt {origin endpoint globus ID} {destination endpoint globus ID} --dry-run

swarm --module bcftools -b 3 vcf_transfer_${batch_prefix}.swarm

```
The final file, the text file, contains just a series of mv commands. These commands will take the original VCFs for a given sample and move them into their prescribed directory (i.e:/data/Kastner_PFS/WGS/[YEAR OF ORIGN]/[SAMPLE]). These VCFs are mostly stored in seq_data_storage, but some are stored elsewhere, so in an effort to consolidate everything in one directory with consistent naming, we are performing these moves. The mv command is very fast, so if you want to make this txt into a swarm file you can, but I recommend just running them through a loop; it will take just a few seconds. 

```
while read command; do ${command} ; done < [TEXT FILE]
```

I plan to add another module to this progran to move the __annotated__ batch VCFs to their destination in /data/Kastner_PFS/WGS/cohort_db/batch_VCFs/. 

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
4. Transfer annotated batch VCFs to Kastner_PFS with data_organization. 

## References 

1. https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
2. https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM
3. https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
4. https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
5. https://hpc.nih.gov/apps/swarm.html
6. https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels
7. https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR
8. https://samtools.github.io/bcftools/bcftools.html
9. https://www.htslib.org/workflow/filter.html
10. https://speciationgenomics.github.io/filtering_vcfs/
11. https://rstudio.github.io/renv/articles/renv.html
12. https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic
13. https://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html
14. https://bioinformatics.stackexchange.com/questions/22124/variants-from-multiple-tools-normalization-before-or-after-annotation-with-vep
15. https://www.ensembl.info/2020/05/26/normalising-variants-to-standardise-ensembl-vep-output/
16. https://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options
17. https://blog.hail.is/introtohail/
18. https://hail.is/docs/0.2/api.html
19. https://discuss.hail.is/


