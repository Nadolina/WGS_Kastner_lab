# Whole genome variant calling pipeline 

This README details the four modules of the WGS pipeline developed for rare variant calling, the choices made in development, and how to use the various scripts. The modules are as such:
1. Revert the original BAM and re-align to GRCh38.
2. Cleanup the alignment to prepare it for variant calling.
3. Variant calling and filtering with GATK.
4. Variant calling and filtering with bcftools.

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

## The pipeline 

### Module 1: Re-alignment 

The process of extract reads from an aligned BAM mostly adheres to the instructions in (2). This extracts the sequences, and removes alignment features, producing an unaligned BAM (uBAM). Then, we mark any adapters which may be artifacts from sequencing and can interfere with alignment. MarkIlluminaAdapters produces a new BAM with these demarcations, as well as a metrics file describing the adapters found, but for space, the demarcated BAM is deleted, and the metrics are retained. In the last step, I pipeline the adapter-marked bam through SamToFastq, bwa-mem2 alignment to GRCh38 with decoys, and MergeBamAlignment to concatenate the unmapped reads back to the final bam. 

Samtools stats is run on both the original and final bam. Fastqc is run on the reverted bam to look at the quality of the sequence reads. 

<div align="center">
  <img src="https://github.com/user-attachments/assets/5df34f9f-9167-4ac3-8ec9-7742378c3a49">
</div>

## Cumulative outputs

```
[sample]/
  [sample]_out/
    [sample].reverted.bam
    [sample].mergedaln.bam
    [sample].mkAdapter.metrics.txt
  qc_out/
    [sample].reverted_fastqc.html
    [sample].reverted_fastqc.zip
    [sample].originalaln.stats
    [sample].mergedaln.stats
    
```


## References 

1. https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
2. https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM


