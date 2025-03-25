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

def query(db):

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

    query=query.drop('AD','GT','DP','GQ','MIN_DP','PGT','PID','PL','PS','RGQ','SB')

    tb=query.localize_entries()

    tb=tb.persist()

    print ("\nAdding pLI, LOEUF and missense z-score annotations from gnomad constraint table.")
    tb=tb.annotate(gene=tb.info['SYMBOL'][0], gene_id=tb.info['Gene'][0], transcript=tb.info['Feature'][0])
    tb=tb.key_by('gene','gene_id','transcript') 
    
    tb=tb.annotate(pLI=constraint[tb.key].lof['pLI'], LOEUF=constraint[tb.key].lof.oe_ci.upper, mis_z_score=constraint[tb.key].mis.z_score)

    ## re-keying our table by locus allele 
    tb=tb.key_by('locus','alleles')

    print ("\nBuilding output table.")
    hail_dict={
        "locus": tb.locus.collect(),
        "alleles": tb.alleles.collect(),
        "rsid": tb.rsid.collect(),
        "qual": tb.qual.collect(),
        "AN": tb.variant_qc['AN'].collect(),
        "AC" : tb.variant_qc['AC'][1].collect(),
        "AF": tb.variant_qc['AF'][1].collect(),
        "DP": tb.info['DP'].collect(),
        "n_het": tb.variant_qc['n_het'].collect(),
        "n_hom_alt": tb.variant_qc['homozygote_count'][1].collect(),
        "het_samples": tb.het_samples.collect(),
        "hom_alt_samples": tb.hom_alt_samples.collect(),
        "pLI": tb.pLI.collect(),
        'LOEUF': tb.LOEUF.collect(),
        'mis_z_score': tb.mis_z_score.collect()
        }

    ## Can't modify "Number" in info lines describing annotations in VCF in convenient way
    ## So, hail expects most of these fields to be arrays. This function collects the first (and only) value if the class of the field is Array. 
    def gettype(name):
        if "Array" in str(tb.info[name].__class__):
            hail_dict["{0}".format(name)]=tb.info[name][0].collect()
        else:
            hail_dict["{0}".format(name)]=tb.info[name].collect()

    ## Applying gettype function to subfields in VCF's INFO field 
    rowlist=list(tb.row.info)
    list(map(gettype,rowlist))

    ## Removing table from cache
    tb=tb.unpersist()

    ## Generating a pandas dataframe 
    return(pd.DataFrame(hail_dict))


## RUNNING QUERY ------------------------------

start=time.time()

db1_df=query('/data/Kastner_PFS/WGS/cohort_db/version_031925/db1.concat.031725.mt')
db2_df=query('/data/Kastner_PFS/WGS/cohort_db/version_031925/db2.concat.031725.mt')
db3_df=query('/data/Kastner_PFS/WGS/cohort_db/version_031925/db3.concat.031725.mt')
db4_df=query('/data/Kastner_PFS/WGS/cohort_db/version_031925/db4.concat.031725.mt')

end=time.time()
print (f"Time taken to perform query: {end-start} seconds")

db_all_df=pd.concat([db1_df,db2_df,db3_df,db4_df])


strlist=["hail-","-".join(samples),"-",date.today().strftime('%m%d%y'),".tsv"]
outname="".join(strlist)
print ("\nOutput query table saved to : %s " % outname)
db_all_df.to_csv(outname,sep='\t')