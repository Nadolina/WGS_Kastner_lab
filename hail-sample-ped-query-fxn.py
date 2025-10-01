#!/usr/bin/env python3

import argparse
import hail as hl 
import os
from pathlib import Path
from datetime import date
import traceback
import time 
import pandas as pd


def query_filter(p, m, s):

    print ("\nRetaining rows where probands have non-ref genotypes...")
    mt_filtered = m.filter_rows(
        hl.agg.count_where(
            hl.literal(p).contains(m.s) & m.GT.is_non_ref()
        ) > 0
    )

    print ("\nFiltering to genotypes in the provided sample list and to HIGH or MODERATE impact variants...")
    query = mt_filtered.filter_cols(hl.literal(set(s)).contains(mt_filtered['s']))
    query = query.filter_rows((query.info['IMPACT'][0] == 'HIGH') | (query.info['IMPACT'][0] == 'MODERATE'))

    query = query.annotate_rows(
        het_samples=hl.agg.filter(query.GT.is_het(), hl.agg.collect(query.s)),
        hom_alt_samples=hl.agg.filter(query.GT.is_hom_var(), hl.agg.collect(query.s))
    )

    return query

def common_annots(mt,q):

    ## mt is the original matrix table, q is the table up to this point 
    ## using variant_qc on the rows shared with the whole cohort to get n_het and n_hom_alt in the cohort

    ht=q.rows()
    constraint=hl.read_table('/data/Kastner_PFS/references/gnomad_v4.1/gnomad.v4.1.constraint_metrics.ht')
    constraint=constraint.key_by('gene','gene_id','transcript')

    mt_semijoin=hl.variant_qc(mt.semi_join_rows(ht))

    tb_semijoin=mt_semijoin.localize_entries().key_by('locus','alleles')
    print ("\nAnnotating variants with cohort-level counts of heterozygous and homozygous alternate carriers...")
    ht = ht.annotate(
        cohort_n_het=tb_semijoin.index(ht.key).variant_qc['n_het'],
        cohort_n_hom_alt = tb_semijoin.index(ht.key).variant_qc['homozygote_count'][1]
    )

    ## Flattening the structs to make it easier to view in the output table
    ht=ht.flatten()

    ## Modifying the keys in the table so I can annotate it with constraint metrics

    print ("\nAnnotating variants with constraint metrics...")
    ht_annot=ht.annotate(gene=ht['info.SYMBOL'][0],gene_id=ht['info.Gene'][0],transcript=ht['info.Feature'][0]).key_by('gene','gene_id','transcript')
    ht_annot = ht_annot.annotate(
        pLI=constraint.index(ht_annot.key).lof['pLI'],
        LOEUF=constraint.index(ht_annot.key).lof.oe_ci.upper,
        mis_z_score=constraint.index(ht_annot.key).mis.z_score
    )

    return ht_annot


def format_query(q, db, samples, outdir):

    dbname=os.path.basename(db)

    print ("\nPreparing table for output...")
    for field in q.row:
        if field in set(q.key): continue
        if field in {"het_samples", "hom_alt_samples"}: continue
        dtype = q[field].dtype
        if dtype == hl.tarray(hl.tstr) or dtype == hl.tarray(hl.tint32):
                q = q.annotate(**{field: hl.delimit(q[field], ',')})

    dbprefix=dbname.replace(".mt","")
    ht_out="".join(["hail-query-","-".join(samples),"-out/",dbprefix,".query.ht"])
    print (f"Writing output to {ht_out}...")
    q.write(ht_out, overwrite=True)

def sample_search_nontrio(db, samples, locus=None):
    try:
        print ("Incomplete trio identified, performing query without inheritance analysis...") 

        mt = hl.read_matrix_table(db)

        # Filter to specific locus if provided
        if locus is not None:
            print(f"Filtering to locus: {locus}")
            interval = hl.parse_locus_interval(f"{locus}", reference_genome="GRCh38")
            mt = mt.filter_rows(interval.contains(mt.locus))

        pedtsv=pd.read_csv(args.pedigree, delim_whitespace=True, header=None)
        probands = [str(pro) for pro in pedtsv.loc[pedtsv[5] == 2, 1].tolist()]
        if len(probands) == 0:
            print ("ERROR: No affected individuals found in the pedigree.")
            return
        print ("Affected individuals in the pedigree:", probands)

        query = query_filter(probands, mt, samples)

        ht_annot = common_annots(mt, query)

        format_query(ht_annot, db, samples, args.outdir)

    except Exception as e:
        print(f"[{db}] ERROR: {e}")
        traceback.print_exc()


def sample_search_trio(db, samples, locus=None):
    try: 
        print ("Complete trio identified, performing query with inheritance analysis...")

        mt = hl.read_matrix_table(db)

        ped= hl.Pedigree.read(args.pedigree) 
        fam=hl.import_fam(args.pedigree)
        pedtsv = pd.read_csv(args.pedigree, delim_whitespace=True, header=None) ## Reading in as a tsv so I can use the affection status
        affected_samples = set(pedtsv.loc[pedtsv[5] == 2, 1].astype(str))

        # Filter to specific locus if provided
        if locus is not None:
            print(f"Filtering to locus: {locus}")
            interval = hl.parse_locus_interval(f"{locus}", reference_genome="GRCh38")
            mt = mt.filter_rows(interval.contains(mt.locus))
        
        probands = [trio.s for trio in ped.complete_trios() if trio.s is not None and str(trio.s) in affected_samples]

        if len(probands) == 0: 
            print ("ERROR: No affected probands found in the pedigree.")
        else: 
            
            query= query_filter(probands, mt, samples)

            query_trio = hl.trio_matrix(query,ped,complete_trios=True)

            print ("\nAnnotating variants with inheritance patterns in probands...")
            trio_annots = [
                hl.struct(
                    proband=t.s,
                    father=t.pat_id,
                    mother=t.mat_id,
                    denovo=(
                        (hl.agg.count_where((query.s == t.s) & query.GT.is_non_ref()) > 0)
                        & (hl.agg.count_where((query.s == t.pat_id) & query.GT.is_hom_ref()) == 0)
                        & (hl.agg.count_where((query.s == t.mat_id) & query.GT.is_hom_ref()) == 0)
                    ),
                    ar_inherited=(
                        (hl.agg.count_where((query.s == t.s) & query.GT.is_hom_var()) > 0)
                        & (hl.agg.count_where((query.s == t.pat_id) & query.GT.is_het()) > 0)
                        & (hl.agg.count_where((query.s == t.mat_id) & query.GT.is_het()) > 0)
                    ),
                    ad_inherited=hl.if_else(
                        (
                            (hl.agg.count_where((query.s == t.s) & query.GT.is_non_ref()) > 0)
                            & (
                                (hl.agg.count_where((query.s == t.pat_id) & query.GT.is_non_ref()) > 0)
                                | (hl.agg.count_where((query.s == t.mat_id) & query.GT.is_non_ref()) > 0)
                            )
                        )
                        & ~(
                            (hl.agg.count_where((query.s == t.s) & query.GT.is_hom_var()) > 0)
                            & (hl.agg.count_where((query.s == t.pat_id) & query.GT.is_het()) > 0)
                            & (hl.agg.count_where((query.s == t.mat_id) & query.GT.is_het()) > 0)
                        ),
                        True,
                        False
                    ),
                    x_recessive=(
                        (query.locus.contig == "X") & (
                            # Male proband: het or hom-alt, mother is carrier
                            (t.is_female==False) &
                            (hl.agg.count_where((query.s == t.s) & query.GT.is_non_ref()) > 0) &
                            (hl.agg.count_where((query.s == t.mat_id) & query.GT.is_non_ref()) > 0)
                            )
                            |
                            # Female proband: hom-alt, both parents are carriers
                            ((t.is_female==True) &
                            (hl.agg.count_where((query.s == t.s) & query.GT.is_hom_var()) > 0) &
                            (hl.agg.count_where((query.s == t.pat_id) & query.GT.is_non_ref()) > 0) &
                            (hl.agg.count_where((query.s == t.mat_id) & query.GT.is_non_ref()) > 0)
                            )
                        ),
                    y_recessive=(
                        (query.locus.contig == "Y") &
                        (t.is_female==False) &  # proband is male
                        (hl.agg.count_where((query.s == t.s) & query.GT.is_non_ref()) > 0) &
                        (hl.agg.count_where((query.s == t.pat_id) & query.GT.is_non_ref()) > 0)
                    )
                )
                for t in ped.complete_trios()
                if str(t.s) in affected_samples
            ]

            query = query.annotate_rows(trio_inheritance=hl.array(trio_annots))

        ht_annot = common_annots(mt, query)

        ht_exploded = ht_annot.explode('trio_inheritance')

        ht_exploded = ht_exploded.annotate(
            proband=ht_exploded.trio_inheritance.proband,
            father=ht_exploded.trio_inheritance.father,
            mother=ht_exploded.trio_inheritance.mother,
            ad_inherited=ht_exploded.trio_inheritance.ad_inherited,
            denovo=ht_exploded.trio_inheritance.denovo,
            ar_inherited=ht_exploded.trio_inheritance.ar_inherited,
            x_recessive=ht_exploded.trio_inheritance.x_recessive,
            y_recessive=ht_exploded.trio_inheritance.y_recessive
        ).drop('trio_inheritance')
        

        format_query(ht_exploded, db, samples, args.outdir)

    except Exception as e:
        print(f"[{db}] ERROR: {e}")
        traceback.print_exc()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--samples", nargs='+', required=True)
    p.add_argument("--pedigree", required=False)
    p.add_argument("--locus", required=False)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    hl.init(
            tmp_dir=f"/lscratch/{os.environ.get('SLURM_JOB_ID', 'tmp')}",
            log=f"hail_{os.getpid()}.log",
            spark_conf={
                'spark.driver.memory': '20g',
                'spark.sql.shuffle.partitions': '8',
                'spark.default.parallelism': '8'
            })

    if args.pedigree is not None:
        ped=hl.Pedigree.read(args.pedigree) 
        trios=ped.complete_trios()

        if len(trios)==0:
            sample_search_nontrio(args.db, args.samples, args.locus)
        else:
            sample_search_trio(args.db, args.samples, args.locus)

