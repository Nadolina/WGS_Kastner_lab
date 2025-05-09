#!/usr/bin/env python3

import argparse
import hail as hl 
import pandas as pd
import time
from datetime import date


## INPUTS ------------------------------

p=argparse.ArgumentParser()
p.add_argument("-s","--samples", 
    nargs = "*",
    help = "Pass sample(s) to query.")
p.add_argument("-n","--num_non_ref",
    type = int,
    help = "Number of non-reference genotypes to look for per variant. Used to filter variants that have no calls or no alt genotypes. Should be less than or equal to the number of samples submitted. [Default 1]")
args = p.parse_args()

if args.samples:
    print ("\nSample(s) to query: % r" % args.samples)
    samples=args.samples
if args.num_non_ref:
    print ("Number of non-ref genotypes to filter for: %i\n" % args.num_non_ref)
    num_non_ref=args.num_non_ref

hl.init(spark_conf={'spark.driver.memory': '25g'}, log='/data/Kastner_PFS/WGS/cohort_db/version_031925/query_logs/testing_032525.log', append=True)
constraint=hl.read_table('/data/brajukanm/hail_testing/constraint/gnomad.v4.1.constraint_metrics.ht')
constraint=constraint.key_by('gene','gene_id','transcript')

## FUNCTIONS ------------------------------

def sample_search(db):

    print ("\nInitializing hail and loading matrix tables.\n")
    mt=hl.read_matrix_table(db)

    ## Removing all samples that are not in the list passed 

    set_to_remove=hl.literal(set(samples))
    query=mt.filter_cols(set_to_remove.contains(mt['s']))

    print ("\nConfirming correct samples were subsetted from %s:" % db)
    query.s.show()

    ## Removing any variants that have no calls or do not contain 2 non-ref genotypes
    query=hl.variant_qc(query)
    query=query.filter_rows(query.variant_qc['n_non_ref']>=num_non_ref) 

    query=query.filter_rows((query.info['IMPACT'][0] == 'HIGH') | (query.info['IMPACT'][0] == 'MODERATE'))

    query=query.annotate_rows(het_samples=hl.agg.filter(query.GT.is_het(),hl.agg.collect(query.s)),
                        hom_alt_samples=hl.agg.filter(query.GT.is_hom_var(),hl.agg.collect(query.s)))

    tb=query.localize_entries()

    return (tb)


def build_df(g_tb,b_tb):

    ## Storing table in chache for improved speed 
    g_tb=g_tb.persist()
    b_tb=b_tb.persist()

    ## Loading constraint matrix table
    constraint=hl.read_table('/data/brajukanm/hail_testing/constraint/gnomad.v4.1.constraint_metrics.ht')

    print ("\nAdding pLI, LOEUF and missense z-score annotations from gnomad constraint table to %s." % g_tb)
    g_tb=g_tb.annotate(gene=g_tb.info['SYMBOL'][0], gene_id=g_tb.info['Gene'][0], transcript=g_tb.info['Feature'][0])
    g_tb=g_tb.key_by('gene','gene_id','transcript') 
    ## need to key-by same columns for annotation (no locus-allele in constraint)
    constraint=constraint.key_by('gene','gene_id','transcript')
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

    ## Removing table from cache
    g_tb=g_tb.unpersist()
    b_tb=b_tb.unpersist()

    ## Generating a pandas dataframe 
    return(pd.DataFrame(hail_dict))


## RUNNING QUERY --------------------

# Tracking query time 
start=time.time() 

g_tb_1=sample_search('/data/Kastner_PFS/WGS/cohort_db/version_042825/db1.concat.042825.mt')
b_tb_1=sample_search('/data/Kastner_PFS/WGS/cohort_db/bcftools_043025/db1.concat.050125.mt')

g_tb_2=sample_search('/data/Kastner_PFS/WGS/cohort_db/version_042825/db2.concat.042825.mt')
b_tb_2=sample_search('/data/Kastner_PFS/WGS/cohort_db/bcftools_043025/db2.concat.050125.mt')

g_tb_3=sample_search('/data/Kastner_PFS/WGS/cohort_db/version_042825/db3.concat.042825.mt')
b_tb_3=sample_search('/data/Kastner_PFS/WGS/cohort_db/bcftools_043025/db3.concat.050125.mt')

g_tb_4=sample_search('/data/Kastner_PFS/WGS/cohort_db/version_042825/db4.concat.042825.mt')
b_tb_4=sample_search('/data/Kastner_PFS/WGS/cohort_db/bcftools_043025/db4.concat.050125.mt')

db1_df=build_df(g_tb_1,b_tb_1)
db2_df=build_df(g_tb_2,b_tb_2)
db3_df=build_df(g_tb_3,b_tb_3)
db4_df=build_df(g_tb_4,b_tb_4)

end=time.time()
print (f"Time taken to perform query: {end-start} seconds")

db_all_df=pd.concat([db1_df,db2_df,db3_df,db4_df])

strlist=["hail-","-".join(samples),"-",date.today().strftime('%m%d%y'),".tsv"]
outname="".join(strlist)
print ("\nOutput query table saved to : %s " % outname)
db_all_df.to_csv(outname,sep='\t')