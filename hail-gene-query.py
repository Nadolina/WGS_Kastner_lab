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
p.add_argument(
    "-g","--genes", 
    help = "Comma-separated list of up to 100 gene symbols to query.")
p.add_argument(
    "-o", "--output",
    help="Name of output file for resulting tsv (ex. 'MEFV-query-032026.tsv')"
)

args = p.parse_args()
if args.genes:
    gene_list = [g.strip() for g in args.genes.split(",") if g.strip()]
    if len(gene_list) > 100:
        raise ValueError("Maximum 100 genes allowed.")
    print("Genes to query:", gene_list)
else:
    raise ValueError("Please provide at least one gene with -g or --genes.")
gene_str = "_".join(gene_list[:3]) + ("_etc" if len(gene_list) > 3 else "")
logname = f"hail-{gene_str}-{date.today().strftime('%m%d%y')}.log"

## REFERENCE FILES --------------------------------
GTF = '/data/Kastner_PFS/references/HG38/gencode.v49.annotation.gtf.gz'
CONSTRAINT_HT = '/data/Kastner_PFS/references/gnomad_v4.1/gnomad.v4.1.constraint_metrics.ht'

## QUERY GENE COORDINATES ---------------------------------

gtf = pd.read_csv(
    GTF, sep='\t', comment='#', header=None, usecols=[0,2,3,4,8],
    names=['chrom','feature','start','end','attrs']
) 

## Resolve gene coordinates from GENCODE v49 GTF before starting Spark
gtf = gtf[gtf['feature'] == 'gene'].copy() ## Only keep gene features to simplify parsing
gtf['gene_name'] = gtf['attrs'].str.extract(r'gene_name "([^"]+)"') ## creating a "gene_name" column 
gene_coords = gtf[gtf['gene_name'].isin(gene_list)][['gene_name','chrom','start','end']]

missing = set(gene_list) - set(gene_coords['gene_name'])
if missing:
    print(f"Warning: no GENCODE coordinates found for: {missing}")
if gene_coords.empty:
    raise ValueError("Could not resolve genomic coordinates for any gene from GENCODE GTF.")
gene_intervals_raw = list(gene_coords.itertuples(index=False, name=None))  # creating list of gene coordinates 

hl.init(spark_conf={'spark.driver.memory': '30g'}, log=logname)

## Build Hail locus intervals now that hl is initialized
gene_intervals = []
for gene, chrom, start, end in gene_intervals_raw:
    try:
        # GENCODE uses chrN naming and 1-based coords; Hail intervals are half-open
        gene_intervals.append(
            hl.parse_locus_interval(f'{chrom}:{start}-{end}', reference_genome='GRCh38')
        )
    except Exception as e:
        print(f"Warning: could not parse interval for '{gene}' ({chrom}:{start}-{end}): {e}")
if not gene_intervals:
    raise ValueError("No valid Hail intervals could be constructed from gene coordinates.")

## FUNCTIONS --------------------------------

def gene_search(mt, source_name):

    print ("\nRunning query in %s.\n" % source_name)

    mt=hl.filter_intervals(mt, gene_intervals)

    query=mt.filter_rows(hl.literal(gene_list).contains(mt.info['SYMBOL'][0]))
    if query.rows().count()==0:
        print ('%s not found in %s' % (gene_list,source_name))
        return 
    # else:
    #     query=query.persist()

    dim=query.count()
    print ("\nNumber of rows matching %s in %s: " % (gene_list,source_name), dim[0], "\n")

    ## Adding some aggregate statistics from hail's built in variant_qc function 
    query=hl.variant_qc(query)

    ## Collecting heterozygous and alt-homozygous samples 
    print ("\nCollecting samples with variants observed in % s" % gene_list)
    query=query.annotate_rows(
        n_het=query.variant_qc['n_het'],
        n_hom_alt=query.variant_qc['homozygote_count'][1],
        het_samples=hl.agg.filter(query.GT.is_het(),hl.agg.collect(query.s)),
        hom_alt_samples=hl.agg.filter(query.GT.is_hom_var(),hl.agg.collect(query.s))
    )
    
    tb=query.localize_entries()
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

## Union GATK databases first, then query once
if __name__ == '__main__':
    databases=['/data/Kastner_PFS/WGS/cohort_db/version_042426/db1.combined.041726.mt',
            '/data/Kastner_PFS/WGS/cohort_db/version_042426/db2.combined.041726.mt',
            '/data/Kastner_PFS/WGS/cohort_db/version_042426/db3.combined.041726.mt',
            '/data/Kastner_PFS/WGS/cohort_db/version_042426/db4.combined.041726.mt']
    gatk_mts = [hl.read_matrix_table(db) for db in databases]
    g_mt = hl.MatrixTable.union_rows(*gatk_mts)
    g_tb = gene_search(g_mt, 'GATK (union of 4 MTs)')

## Union bcftools databases first, then query once
    databases=['/data/Kastner_PFS/WGS/cohort_db/bcftools_042426/db1.combined.042226.mt',
        '/data/Kastner_PFS/WGS/cohort_db/bcftools_042426/db2.combined.042226.mt',
        '/data/Kastner_PFS/WGS/cohort_db/bcftools_042426/db3.combined.042226.mt',
        '/data/Kastner_PFS/WGS/cohort_db/bcftools_042426/db4.combined.042226.mt'
        ]
    bcftools_mts = [hl.read_matrix_table(db) for db in databases]
    b_mt = hl.MatrixTable.union_rows(*bcftools_mts)
    b_tb = gene_search(b_mt, 'bcftools (union of 4 MTs)')


btb_semijoin=b_tb.semi_join(g_tb).flatten().key_by('locus','alleles') ## Shared loci+alleles between gatk and bcftools calls 
btb_antijoin=b_tb.anti_join(g_tb).flatten().key_by('locus','alleles') ## Loci/alleles seen in bcftools and not GATK 

#Renaming bcftools columns to distinguish them when antijoin unioned to gtb 
btb_antijoin=btb_antijoin.transmute(
    b_n_hom_alt=btb_antijoin['n_hom_alt'],
    b_n_het=btb_antijoin['n_het'],
    b_het_samples=btb_antijoin['het_samples'],
    b_hom_alt_samples=btb_antijoin['hom_alt_samples']
    )

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

## Loading consraint table and annotating union table with constraint metrics
constraint = hl.read_table(CONSTRAINT_HT)
constraint = constraint.key_by('gene','gene_id','transcript')

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

print ("\nOutput query table saved to : %s " % args.output)
reorder_df.to_csv(args.output, sep='\t')

end=time.time()

print (f"Time taken to perform query: {end-start} seconds")

