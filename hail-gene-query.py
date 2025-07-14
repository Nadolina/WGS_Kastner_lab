#!/usr/bin/env python3

import argparse
import hail as hl 
import pandas as pd
import time
from datetime import date
from pathlib import Path

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

## Loading constraint matrix table from gnomad v4.1
constraint=hl.read_table('/data/Kastner_PFS/references/gnomad_v4.1/gnomad.v4.1.constraint_metrics.ht')
constraint=constraint.key_by('gene','gene_id','transcript')

## FUNCTIONS ---------------------

def gene_search(db):

    print ("\nInitializing hail and loading matrix tables.\n")
    mt=hl.read_matrix_table(db)

    print ("\nStarting query of % s in % s." % (args.gene,db))

    query=mt.filter_rows(mt.info['SYMBOL'].contains(gene))
    if query.rows().count()==0:
        print ('%s not found in %s' % (gene,db))
        return None

    print ("\nStructure of matrix table:\n")
    query.describe()

    dim=query.count()
    print ("\nNumber of rows matching %s: " % gene, dim[0], "\n")

    ## Adding some aggregate statistics from hail's built in variant_qc function 
    query=hl.variant_qc(query)

    ## Collecting heterozygous and alt-homozygous samples 
    print ("\nCollecting samples with variants observed in % s" % gene)
    query=query.annotate_rows(het_samples=hl.agg.filter(query.GT.is_het(),hl.agg.collect(query.s)),
                         hom_alt_samples=hl.agg.filter(query.GT.is_hom_var(),hl.agg.collect(query.s)))

    tb=query.localize_entries()

    return (tb)


def build_df(g_tb,b_tb):

    ## Storing table in chache for improved speed 
    g_tb=g_tb.persist()
    b_tb=b_tb.persist()

    ## adding pLI, LOEUF and mis_z_score from gnomad's pre-built constraint hail table to GATK db - will add to remaining bcftools variants later.
    ## https://gnomad.broadinstitute.org/help/constraint
    ## https://gnomad.broadinstitute.org/data#v4-constraint
    print ("\nAdding pLI, LOEUF and missense z-score annotations from gnomad constraint table to %s." % g_tb)
    g_tb=g_tb.annotate(gene=g_tb.info['SYMBOL'][0], gene_id=g_tb.info['Gene'][0], transcript=g_tb.info['Feature'][0])
    g_tb=g_tb.key_by('gene','gene_id','transcript') 
    ## need to key-by same columns as constraint table for annotation (no locus-allele in constraint)
    g_tb=g_tb.annotate(pLI=constraint[g_tb.key].lof['pLI'], LOEUF=constraint[g_tb.key].lof.oe_ci.upper, mis_z_score=constraint[g_tb.key].mis.z_score)

    ## re-keying our table by locus allele 
    g_tb=g_tb.key_by('locus','alleles')

    ## Annotating the GATK table to show variants also called by BCFTOOLS
    g_tb=g_tb.annotate(bcftools_het_samples=b_tb[g_tb.key].het_samples,
        bcftools_hom_samples=b_tb[g_tb.key].hom_alt_samples,
        bcftools_DP=b_tb[g_tb.key].info['DP'],
        bcftools_AF=b_tb[g_tb.key].variant_qc['AF'][1],
        bcftools_AC=b_tb[g_tb.key].variant_qc['AC'][1],
        bcftools_AN=b_tb[g_tb.key].variant_qc['AN'])

    ## Initializing dictionary to build output table from 
    print ("\nBuilding output table.")
    hail_dict={
        "locus": g_tb.locus.collect(),
        "alleles": g_tb.alleles.collect(),
        "rsid": g_tb.rsid.collect(),
        "qual": g_tb.qual.collect(),
        "AN": g_tb.variant_qc['AN'].collect(),
        "AC" : g_tb.variant_qc['AC'][1].collect(),
        "AF": g_tb.variant_qc['AF'][1].collect(),
        "DP": g_tb.info['DP'].collect(),
        "n_het": g_tb.variant_qc['n_het'].collect(),
        "n_hom_alt": g_tb.variant_qc['homozygote_count'][1].collect(),
        "het_samples": g_tb.het_samples.collect(),
        "hom_alt_samples": g_tb.hom_alt_samples.collect(),
        "bcftools_het_samples": g_tb.bcftools_het_samples.collect(),
        "bcftools_hom_alt_samples": g_tb.bcftools_hom_samples.collect(),
        "bcftools_AN": g_tb.bcftools_AN.collect(),
        "bcftools_AC": g_tb.bcftools_AC.collect(),
        "bcftools_AF": g_tb.bcftools_AF.collect(),
        "bcftools_DP": g_tb.bcftools_DP.collect(),
        "pLI": g_tb.pLI.collect(),
        'LOEUF': g_tb.LOEUF.collect(),
        'mis_z_score': g_tb.mis_z_score.collect()
        }


    ## Can't modify "Number" in info lines describing annotations in VCF in convenient way
    ## So, hail expects most of these fields to be arrays. This function collects the first (and only) value if the class of the field is Array. 
    def gettype(name):
        if "Array" in str(g_tb.info[name].__class__):
            hail_dict["{0}".format(name)]=g_tb.info[name][0].collect()
        else:
            hail_dict["{0}".format(name)]=g_tb.info[name].collect()

    ## Applying gettype function to subfields in VCF's INFO field 
    rowlist=list(g_tb.row.info)
    list(map(gettype,rowlist))

    ## Generating a pandas dataframe 
    return(pd.DataFrame(hail_dict))

    # Removing table from cache
    g_tb=g_tb.unpersist()
    b_tb=b_tb.unpersist()


## RUNNING QUERY --------------------

# Tracking query time 
start=time.time() 

g_tb_1=gene_search('/data/Kastner_PFS/WGS/cohort_db/version_042825/db1.concat.042825.mt')
b_tb_1=gene_search('/data/Kastner_PFS/WGS/cohort_db/bcftools_043025/db1.concat.050125.mt')

g_tb_2=gene_search('/data/Kastner_PFS/WGS/cohort_db/version_042825/db2.concat.042825.mt')
b_tb_2=gene_search('/data/Kastner_PFS/WGS/cohort_db/bcftools_043025/db2.concat.050125.mt')

g_tb_3=gene_search('/data/Kastner_PFS/WGS/cohort_db/version_042825/db3.concat.042825.mt')
b_tb_3=gene_search('/data/Kastner_PFS/WGS/cohort_db/bcftools_043025/db3.concat.050125.mt')

g_tb_4=gene_search('/data/Kastner_PFS/WGS/cohort_db/version_042825/db4.concat.042825.mt')
b_tb_4=gene_search('/data/Kastner_PFS/WGS/cohort_db/bcftools_043025/db4.concat.050125.mt')

if g_tb_1 is not None:
    db_df=build_df(g_tb_1,b_tb_1)
elif g_tb_2 is not None:
    db_df=build_df(g_tb_2,b_tb_2)
elif g_tb_3 is not None:
    db_df=build_df(g_tb_3,b_tb_3)
elif g_tb_4 is not None:
    db_df=build_df(g_tb_4,b_tb_4)

end=time.time()
print (f"Time taken to perform query: {end-start} seconds")

# db_all_df=pd.concat([db1_df,db2_df,db3_df,db4_df])

# Using current date and input gene to construct output file name 
strlist=["hail-",gene,"-",date.today().strftime('%m%d%y'),".tsv"]
outname="".join(strlist)
print ("\nOutput query table saved to : %s " % outname)
db_df.to_csv(outname,sep='\t')
