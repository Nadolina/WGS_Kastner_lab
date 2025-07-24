#!/usr/bin/env python3

import argparse
import hail as hl 
import pandas as pd
import time
from datetime import date
from pathlib import Path
# from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor

## Parsing gene symbol from command line input 
p=argparse.ArgumentParser()
p.add_argument("-g","--gene", help = "Pass gene symbol to query.")
args = p.parse_args()
if args.gene:
    print ("Gene to query: % s" % args.gene)
    gene=str(args.gene)

## Initializing hail and collecting pre-built matrix tables for querying
logname="".join(['hail-',gene,'-',date.today().strftime('%m%d%y'),'.log'])
hl.init(spark_conf={'spark.driver.memory': '25g'}, log=logname)

constraint=hl.read_table('/data/Kastner_PFS/references/gnomad_v4.1/gnomad.v4.1.constraint_metrics.ht')
constraint=constraint.key_by('gene','gene_id','transcript')

## FUNCTIONS --------------------------------

def gene_search(db):

    print ("\nInitializing hail and loading matrix table %s.\n" % db)
    mt=hl.read_matrix_table(db)

    query=mt.filter_rows(mt.info['SYMBOL'].contains(gene))
    if query.rows().count()==0:
        print ('%s not found in %s' % (gene,db))
        return 
    else:
        query=query.persist()

    dim=query.count()
    print ("\nNumber of rows matching %s in %s: " % (gene,db), dim[0], "\n")

    ## Adding some aggregate statistics from hail's built in variant_qc function 
    query=hl.variant_qc(query)

    ## Collecting heterozygous and alt-homozygous samples 
    print ("\nCollecting samples with variants observed in % s" % gene)
    query=query.annotate_rows(n_het=query.variant_qc['n_het'],
                        n_hom_alt=query.variant_qc['homozygote_count'][1],
                        het_samples=hl.agg.filter(query.GT.is_het(),hl.agg.collect(query.s)),
                        hom_alt_samples=hl.agg.filter(query.GT.is_hom_var(),hl.agg.collect(query.s)))
    
    tb=query.localize_entries()
    query=query.unpersist()
    return (tb)

# Function to safely extract string from single-item list
def extract_single_string(cell):
    if isinstance(cell, list) and len(cell) == 1:
        return str(cell[0])
    elif isinstance(cell, list) and len(cell) == 0:
        return None  # or '' if you prefer an empty string
    return cell  # fallback

## QUERY -------------------------------------

start=time.time()

## Querying the GATK variant calls in parallel
if __name__ == '__main__':
    databases=['/data/Kastner_PFS/WGS/cohort_db/version_070925/db1.concat.070925.mt',
            '/data/Kastner_PFS/WGS/cohort_db/version_070925/db2.concat.070925.mt',
            '/data/Kastner_PFS/WGS/cohort_db/version_070925/db3.concat.070925.mt',
            '/data/Kastner_PFS/WGS/cohort_db/version_070925/db4.concat.070925.mt']
    with ThreadPoolExecutor(max_workers=4) as executor:
        queries=list(executor.map(gene_search,databases))
    
    for query in queries:
        if query is not None:
            g_tb=query

## Querying the bcftools calls in parallel
if __name__ == '__main__':
    databases=['/data/Kastner_PFS/WGS/cohort_db/bcftools_071425/db1.concat.071425.mt',
        '/data/Kastner_PFS/WGS/cohort_db/bcftools_071425/db2.concat.071425.mt',
        '/data/Kastner_PFS/WGS/cohort_db/bcftools_071425/db3.concat.071425.mt',
        '/data/Kastner_PFS/WGS/cohort_db/bcftools_071425/db4.concat.071425.mt'
        ]
    with ThreadPoolExecutor(max_workers=4) as executor:
        queries=list(executor.map(gene_search,databases))
    
    for query in queries:
        if query is not None:
            b_tb=query


btb_semijoin=b_tb.semi_join(g_tb).flatten().key_by('locus','alleles') ## Shared loci+alleles between gatk and bcftools calls 
btb_antijoin=b_tb.anti_join(g_tb).flatten().key_by('locus','alleles') ## Loci/alleles seen in bcftools and not GATK 

# # btb_flat=btb_antijoin.flatten()

#Renaming bcftools columns to distinguish them when antijoin unioned to gtb 
btb_antijoin=btb_antijoin.transmute(b_n_hom_alt=btb_antijoin['n_hom_alt'],
                        b_n_het=btb_antijoin['n_het'],
                        b_het_samples=btb_antijoin['het_samples'],
                        b_hom_alt_samples=btb_antijoin['hom_alt_samples'])

## For shared loci/alleles, I am just annotating the GATK call rows with bcftools call info
g_tb=g_tb.annotate(b_het_samples=btb_semijoin[g_tb.key].het_samples,
        b_hom_alt_samples=btb_semijoin[g_tb.key].hom_alt_samples,
        b_n_hom_alt=btb_semijoin[g_tb.key]['n_hom_alt'],
        b_n_het=btb_semijoin[g_tb.key]['n_het'],
        b_DP=btb_semijoin[g_tb.key]['info.DP'],
        b_AF=btb_semijoin[g_tb.key]['variant_qc.AF'],
        b_AC=btb_semijoin[g_tb.key]['variant_qc.AC'],
        b_AN=btb_semijoin[g_tb.key]['variant_qc.AN'])
gtb_flat=g_tb.flatten().key_by('locus','alleles')

## Union combines tables that do not have shared keys - so it retains unique locus/allele pairs from GATK and bcftools 
union_tb=gtb_flat.union(btb_antijoin,unify=True)

#temporarily re-key the table to add-in constraint information
union_tb=union_tb.annotate(gene=union_tb['info.SYMBOL'][0], 
                gene_id=union_tb['info.Gene'][0], 
                transcript=union_tb['info.Feature'][0])
union_tb=union_tb.key_by('gene','gene_id','transcript')
union_tb=union_tb.annotate(pLI=constraint[union_tb.key].lof['pLI'], 
                LOEUF=constraint[union_tb.key].lof.oe_ci.upper, 
                mis_z_score=constraint[union_tb.key].mis.z_score)
union_tb=union_tb.key_by('locus','alleles')

## Transform to pandas dataframe and re-order columns for output 
df=union_tb.to_pandas()
cols=df.columns.tolist()
for col in cols:
    df[col] = df[col].apply(extract_single_string)

colorder=['locus','alleles','n_het','n_hom_alt','het_samples','hom_alt_samples','b_n_het','b_n_hom_alt',
            'b_het_samples','b_hom_alt_samples','b_DP','b_AF','b_AC','b_AN','rsid','qual','info.DP',
            'variant_qc.AF','variant_qc.AN','variant_qc.dp_stats.mean','variant_qc.dp_stats.min',
            'variant_qc.dp_stats.max','variant_qc.gq_stats.mean','info.gnomADg',
            'info.gnomADg_AC', 'info.gnomADg_AN','info.gnomADg_nhomalt_joint',
            'info.gnomADg_AF_joint','info.gnomADg_AF_grpmax','info.gnomADg_AC_XY','info.gnomADg_non_par',
            'info.gnomADg_sift_max','info.gnomADg_polyphen_max','info.ClinVar',
            'info.ClinVar_CLNSIG','info.GERP','info.Consequence','info.IMPACT','info.SYMBOL',
            'info.Gene','info.Feature_type','info.Feature','info.BIOTYPE','info.EXON','info.INTRON',
            'info.HGVSc','info.HGVSp','info.cDNA_position','info.CDS_position','info.Protein_position',
            'info.Amino_acids','info.Codons','info.Existing_variation','info.DISTANCE','info.STRAND',
            'info.FLAGS','info.FLAGS','info.PICK','info.SYMBOL_SOURCE','info.HGNC_ID','info.CANONICAL',
            'info.SOURCE','info.SOMATIC','info.PHENO','info.CADD_PHRED','info.SpliceAI_pred_DP_AG',
            'info.SpliceAI_pred_DP_AL','info.SpliceAI_pred_DP_DG','info.SpliceAI_pred_DP_DL',
            'info.SpliceAI_pred_DS_AG','info.SpliceAI_pred_DS_AL','info.SpliceAI_pred_DS_DG',
            'info.SpliceAI_pred_DS_DL','info.SpliceAI_pred_SYMBOL','info.MaxEntScan_alt','info.MaxEntScan_diff',
            'info.MaxEntScan_ref','info.REVEL','info.am_class','info.am_pathogenicity','info.am_protein_variant',
            'info.am_transcript_id','info.am_uniprot_id','pLI','LOEUF','mis_z_score']
reorder_df=df[colorder]      

strlist=["hail-",gene,"-",date.today().strftime('%m%d%y'),".tsv"]
outname="".join(strlist)
print ("\nOutput query table saved to : %s " % outname)
reorder_df.to_csv(outname, sep='\t')

end=time.time()

print (f"Time taken to perform query: {end-start} seconds")

