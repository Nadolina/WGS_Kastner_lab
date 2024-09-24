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

## The pipeline 

### Module 1: Re-alignment 

<div align="center">
  <img src="https://github.com/user-attachments/assets/5df34f9f-9167-4ac3-8ec9-7742378c3a49">
</div>

