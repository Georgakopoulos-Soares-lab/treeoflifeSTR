from minditool import MindiTool
from pathlib import Path
import pandas as pd
import logging 
import json 
import gzip
import os
import pandas as pd

def load_bucket(bucket_id, schedule):
    """
    Load a bucket of accessions from a JSON file.
    """
    with open(schedule, 'r') as f:
        buckets = json.load(f)
    return buckets.get(str(bucket_id), [])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="MindiTool Extraction and Validation")
    parser.add_argument("--bucket_id", type=int, default=0, help="Bucket ID for parallel processing")
    parser.add_argument("--schedule", type=str, default="all", help="Schedule for processing (e.g., all, daily, weekly)")
    parser.add_argument("--pattern", type=str, default="STR")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory for extracted files")
    args = parser.parse_args()
    
    if args.outdir == ".":
        pattern__repr__ = "_".join(args.pattern.split(","))
        args.outdir = f"extractions_{pattern__repr__}"
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(exist_ok=True)
    accessions = load_bucket(args.bucket_id, args.schedule)
    # Create logging directory
    logdir = outdir.joinpath("log_debug_nonbdna")
    logdir.mkdir(exist_ok=True)
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s', 
                        filename=f'{logdir}/mindi_tool_{args.bucket_id}.log', 
                        filemode='w'
                        )
    if not accessions:
        logging.info(f"No accessions found for bucket {args.bucket_id} in schedule {args.schedule}.")
        exit(0)
    else:
        logging.info(f"Loaded {len(accessions)} accessions for bucket {args.bucket_id} in schedule {args.schedule}.")
    status = outdir.joinpath("status")
    status.mkdir(exist_ok=True)
    total_accessions = len(accessions)
    subdir = outdir.joinpath("merged_bulk")
    subdir.mkdir(exist_ok=True)
    fout = {}
    for m in args.pattern.split(","):
        fout[m] = gzip.open(subdir.joinpath(f"bucket_{args.bucket_id}_mode_{m}.tsv.gz"), mode="wt", encoding="UTF-8"), \
                  open(subdir.joinpath(f"bucket_{args.bucket_id}_mode_{m}_empty.tsv"), mode="w", encoding="UTF-8")
    # Start
    for i, accession in enumerate(accessions, 1):
        logging.info(f"Processing accession: {accession}. Progress: {(i-1) * 1e2 / total_accessions:.2f}% (BucketID: {args.bucket_id}).")
        hunter = MindiTool(tempdir=outdir)
        result = hunter.extract(accession, 
                                pattern=args.pattern)
        if result is None:
            logging.error(f"Extraction failed for accession: {accession}.")
            continue
        for m in args.pattern.split(","):
            try:
                hunter.sanitize(accession, mode=m)
            except AssertionError:
                logging.error(f"Extraction failed for accession: {accession} (mode {m}).")
                os.remove(hunter.fn[m])
                continue
            try:
                result_df = pd.read_table(hunter.fnp[m])
                if result_df.shape[0] == 0:
                    fout[m][1].write(accession + "\n")
                # result_df.loc[:, "#assembly_accession"] = MindiTool.extract_id(accession)
                result_df.to_csv(fout[m][0], header=i==1, index=False, sep="\t")
            except pd.errors.EmptyDataError:
                fout[m][1].write(accession + "\n")
            os.remove(hunter.fn[m])
    for m in args.pattern.split(","):
        for i in range(2):
            fout[m][i].close()
    logging.info(f"Completed processing for bucket {args.bucket_id} in schedule {args.schedule}.")
    with status.joinpath(f"bucket_{args.bucket_id}_status.completed").open("w") as f:
        f.write(f"Bucket {args.bucket_id} has been completed\n")
