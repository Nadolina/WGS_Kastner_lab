#!/usr/bin/env python3

import argparse
import hail as hl 
import os
from pathlib import Path
from datetime import date
import traceback
import time 

def query_pair_combine(g_ht_path,b_ht_path,outdir):
    try: 
        hl.init(
            tmp_dir=f"/lscratch/{os.environ.get('SLURM_JOB_ID', 'tmp')}",
            log=f"hail_{os.getpid()}.log",
            spark_conf={
                'spark.driver.memory': '20g',
                'spark.sql.shuffle.partitions': '8',
                'spark.default.parallelism': '8'
            })

        print (f"\nReading hail tables...\nGATK table: {g_ht_path}\nbcftools table: {b_ht_path}\n")
        g_ht=hl.read_table(g_ht_path)
        b_ht=hl.read_table(b_ht_path)
        g_ht_prefix=os.path.basename(g_ht_path).replace('.ht','')
        b_ht_prefix=os.path.basename(b_ht_path).replace('.ht','')

        ## key tables
        g_ht=g_ht.key_by('locus','alleles')
        b_ht=b_ht.key_by('locus','alleles')

        ## Identify shared variants between GATK and bcftools calls using locus-allele keys 
        b_semijoin = b_ht.semi_join(g_ht)
        b_antijoin = b_ht.anti_join(g_ht)
        b_semijoin = b_semijoin.annotate(bcftools_only="Both")

        ## Dropping fields that are specific to GATK or bcftools variant calling and cannot be unionized, running checkpoint 
        g_ht_drop = g_ht.drop('info.AF','info.AC','info.AN','info.AS_BaseQRankSum','info.AS_FS','info.AS_InbreedingCoeff','info.AS_MQ','info.AS_MQRankSum','info.AS_QD','info.END')
        g_ht_drop = g_ht_drop.annotate(
                bcftools_only = hl.if_else(
                hl.is_defined(b_semijoin[g_ht_drop.key]),
                "Both",
                "False"
            )
        )
        g_ht_drop = g_ht_drop.key_by('locus','alleles').checkpoint(f"{outdir}/{g_ht_prefix}_drop.ht", overwrite=True)
        
        b_antijoin_drop=b_antijoin.drop('info.INDEL','info.IMF','info.SGB','info.REF_AN','info.REF_AC','info.AC','info.AN','info.DP4','info.MQ')
        b_antijoin_drop = b_antijoin_drop.annotate(bcftools_only="True").checkpoint(f"{outdir}/{b_ht_prefix}_antijoin_drop.ht", overwrite=True)

        ht=g_ht_drop.union(b_antijoin_drop)
        dbprefix=(Path(g_ht_path).name).split(".")[0]
        ht_out="".join([outdir,"/",dbprefix,".query.ht"])
        print ("Results written to: ", ht_out)
        ht.write(ht_out, overwrite=True)

    except Exception as e:
        print(f"[{g_ht} {b_ht}] ERROR: {e}")
        traceback.print_exc()

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument("--g_ht_path", required=True)
    p.add_argument("--b_ht_path", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    query_pair_combine(args.g_ht_path, args.b_ht_path, args.outdir)
