#!/usr/bin/env python3

import argparse
import hail as hl 
import pandas as pd
import time
from datetime import date

## Parsing gene symbol from command line input 
p=argparse.ArgumentParser()
p.add_argument("-g","--gene", help = "Pass gene symbol to query.")
args = p.parse_args()
if args.gene:
    print ("Gene to query: % s" % args.gene)
    gene=str(args.gene)

# Tracking query time 
start=time.time() 

## Initializing hail and collecting pre-built matrix tables for querying
print ("\nInitializing hail and loading matrix tables.\n")
hl.init()
mt=hl.read_matrix_table('subset.gatk.022725.mt')
constraint=hl.read_table('/data/brajukanm/hail_testing/constraint/gnomad.v4.1.constraint_metrics.ht')

## Performing query 
print ("\nStarting query of % s." % args.gene)
query=mt.filter_rows(mt.info['SYMBOL'].contains(gene))

print ("\nStructure of matrix table:\n")
query.describe()

dim=query.count()
print ("\nNumber of rows matching %s: " % gene, dim[0], "\n")

## Adding some aggregate statistics from hail's built in variant_qc function 
query=hl.variant_qc(query)

## Collecting heterozygous and alt-homozygous samples 
print ("\nCollecting samples with variants observed in % s" % gene)
annot=query.annotate_rows(het_samples=hl.agg.filter(query.GT.is_het(),hl.agg.collect(query.s)),
                    hom_alt_samples=hl.agg.filter(query.GT.is_hom_var(),hl.agg.collect(query.s)))

annot=annot.drop('AD','GT','DP','GQ','MIN_DP','PGT','PID','PL','PS','RGQ','SB')

tb=annot.localize_entries()

## Stoing table in chache for improved speed 
tb=tb.persist()

## adding pLI, LOEUF and mis_z_score from gnomad's pre-built constraint hail table 
## https://gnomad.broadinstitute.org/help/constraint
## https://gnomad.broadinstitute.org/data#v4-constraint
print ("\nAdding pLI, LOEUF and missense z-score annotations from gnomad constraint table.")
tb=tb.annotate(gene=tb.info['SYMBOL'][0], gene_id=tb.info['Gene'][0], transcript=tb.info['Feature'][0])
tb=tb.key_by('gene','gene_id','transcript') 
## need to key-by same columns for annotation (no locus-allele in constraint)
constraint=constraint.key_by('gene','gene_id','transcript')
tb=tb.annotate(pLI=constraint[tb.key].lof['pLI'], LOEUF=constraint[tb.key].lof.oe_ci.upper, mis_z_score=constraint[tb.key].mis.z_score)

## re-keying our table by locus allele 
tb=tb.key_by('locus','alleles')

## Initializing dictionary to build output table from 
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

## Using current date and input gene to construct output file name 
strlist=["hail-",gene,"-",date.today().strftime('%m%d%y'),".tsv"]
outname="".join(strlist)
print ("\nOutput query table saved to : %s " % outname)

## Generating a pandas dataframe and then an output table
hail_df=pd.DataFrame(hail_dict)
hail_df.to_csv(outname,sep='\t')

end=time.time()
print (f"Time taken to perform query: {end-start} seconds")

## Removing table from cache
tb=tb.unpersist()

