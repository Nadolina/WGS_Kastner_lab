#!/usr/bin/env python3

import traceback
import hail as hl
import argparse
import os.path
import sys

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--db", required=False)
    p.add_argument("--samples",required=True)
    args = p.parse_args()

    try:

        SAMPLES=args.samples.split(",")
        JOINED_IDS="-".join(SAMPLES)

        if args.db:
            print("Database provided: db", args.db)
            DB="".join(["db", str(args.db),".query.ht"])

            if not os.path.exists("".join([args.outdir,"/",DB])):
                raise FileNotFoundError(f"The specified hail table does not exist: {DB}")
                sys.exit(1)
            else:
                print(f"The specified hail table exists: {DB}")

            hl.init(spark_conf={'spark.driver.memory': '10g'})

            db=hl.read_table("".join([args.outdir,"/",DB]))
            print ("Hail table read successfully.")

            OUTFILE=f"{args.outdir}/hail-query-{JOINED_IDS}-db{str(args.db)}.tsv"
            db.export(OUTFILE)
            print (f"Hail table exported to: {args.outdir}/hail-query-{JOINED_IDS}-db{str(args.db)}.tsv")

        elif args.db is None:
            print("No database provided, assuming all databases to be combined.")

            db1_path="".join([args.outdir,"/","db1.query.ht"])
            db2_path="".join([args.outdir,"/","db2.query.ht"])
            db3_path="".join([args.outdir,"/","db3.query.ht"])
            db4_path="".join([args.outdir,"/","db4.query.ht"])

            list_of_paths = [db1_path, db2_path, db3_path, db4_path]

            # Check if all hail tables exist before proceeding
            if not all(os.path.exists(p) for p in list_of_paths):
                # Find which tables are missing to include in the error message
                missing_paths = [p for p in list_of_paths if not os.path.exists(p)]
                raise FileNotFoundError(f"One or more required hail tables do not exist: {missing_paths}")
                exitsys.exit(1)
            else:
                print("All hail tables exist, generating final hail table...")
                
            ## Initating hail...
            hl.init(spark_conf={'spark.driver.memory': '20g'})

            ## Reading in hail tables
            db1=hl.read_table(db1_path)
            db2=hl.read_table(db2_path)
            db3=hl.read_table(db3_path)
            db4=hl.read_table(db4_path)
            print ("Hail tables read successfully.")
            
            ## Combining all four hail tables into one final hail table
            ht_final=db1.union(db2,db3,db4)
            print ("Number of variants in final table:", ht_final.count())

            ## Writing hail table to output TSV filels
            OUTFILE=f"{args.outdir}/hail-query-{JOINED_IDS}.tsv"
            ht_final.export(OUTFILE)
            print (f"Final combined hail table written to: {OUTFILE}")

    except Exception as e:
        print(f"ERROR: {e}")
        traceback.print_exc()