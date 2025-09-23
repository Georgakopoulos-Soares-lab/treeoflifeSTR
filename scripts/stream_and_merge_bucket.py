import json
from pathlib import Path
from tqdm import tqdm
from termcolor import colored
import gzip
from typing import Optional
import pandas as pd
from collections import defaultdict
from pybedtools import BedTool
from minditool import MindiTool
import re
from attr import field
import attr

@attr.s
class StreamAndMerge:
    indir: str = field(converter=str, init=True, repr=True)
    # log_indir: str = field(converter=str, init=True, repr=True)
    schedule: str = field(converter=str, init=True, repr=False)
    merge_outdir: str = field(init=False, repr=True)
    max_spacer: int = field(converter=int, default=7, init=True, repr=True)
    min_arm_length: int = field(converter=int, default=10, init=True, repr=True)
    min_sequence_length: int = field(converter=int, default=10, init=True, repr=True)
    min_consensus_repeats: int = field(converter=int, default=3, init=True, repr=True)
    GA_threshold: float = field(converter=float, default=0.9, init=True, repr=True)
    GT_threshold: float = field(converter=float, default=0.9, init=True, repr=True)
    AT_threshold: float = field(converter=float, default=0.8, init=True, repr=True)

    def __attrs_post_init__(self) -> None:
        self.indir = Path(self.indir).resolve()
        self.schedule = Path(self.schedule).resolve()
        if not self.schedule.is_file():
            raise FileNotFoundError()
        # self.log_indir = Path(self.log_indir).resolve()
        self.log_indir = self.indir.joinpath("log_debug_nonbdna")
        self.merge_outdir = self.indir.joinpath("merged")
        self.merge_outdir.mkdir(exist_ok=True)
        print(f"Chosen parameters for frozen instance {type(self).__name__}:")
        print(colored(f"Minimum consensus repeats allowed: {self.min_consensus_repeats}.", "green"))
        print(colored(f"Max spacer length: {self.max_spacer}.", "green"))
        print(colored(f"Min arm length: {self.min_arm_length}.", "green"))
        print(colored(f"Min sequence length: {self.min_sequence_length}.", "green"))
        print(colored(f"GA threshold: {self.GA_threshold:.2f}.", "green"))
        print(colored(f"GT threshold: {self.GT_threshold:.2f}.", "green"))
        print(colored(f"AT threshold: {self.AT_threshold:.2f}.", "green"))
        print(colored(f"Retrieving accessions from schedule {self.schedule}...", "green"))
        print(colored(f"Outsourcing results to directory: {self.merge_outdir}...", "green"))

    def load_bucket(self, bucket_id: int) -> list[str]:
        with self.schedule.open("r", encoding="UTF-8") as f:
            return json.load(f)[str(bucket_id)]

    @staticmethod
    def is_square_free(seq: str) -> bool:
        if seq == ".":
            return None
        return re.search(r"([agct]+)\1", seq) is None

    @staticmethod
    def is_cubic_free(seq: str) -> bool:
        if seq == ".":
            return None
        return re.search(r"([agct]+)\1\1", seq) is None

    @staticmethod
    def detect_STR_coverage(seq: str) -> float:
        """
        Returns the STR that coverages the most of sequence.
        """
        if seq == ".":
            return None
        motifs = r"[agct]"
        matches = re.finditer(r"(%s+)\1{2,}" % motifs, seq)
        covered = [0 for _ in range(len(seq))]
        for match in matches:
            # print(match)
            start, end = match.span()
            for i in range(start, end):
                covered[i] = 1
        return sum(covered) # / len(seq) if seq else 0.0
    
    def load_log(self, bucket_id: int, pattern: str) -> list[str]:
        extract_id = lambda accession: "_".join(Path(accession).name.split("_")[:2])
        log_file = self.log_indir.joinpath(f"mindi_tool_{bucket_id}.log")
        extracted_ids = set()
        failed_ids = set()
        bucket_is_complete = False
        with log_file.open("r", encoding="UTF-8") as f:
            for line in f:
                line = line.strip()
                if f"Extraction failed for accession" in line: # and f"mode {pattern}" in line:
                    grouped_id = re.search(r"(GC[AF]_.+)\.fna", line)
                    assembly_id = extract_id(grouped_id.group(1))
                    failed_ids.add(assembly_id)
                elif f"has passed all checks (mode {pattern})." in line:
                    grouped_id = re.search(r"(GC[AF]_.+)\.fna", line)
                    assembly_id = extract_id(grouped_id.group(1))
                    extracted_ids.add(assembly_id)
                elif f"Completed processing for bucket {bucket_id} in schedule" in line:
                    bucket_is_complete = True
        failed_ids -= extracted_ids
        extracted_ids = list(extracted_ids)
        print(colored(f"Total extracted ids found: {len(extracted_ids)} for pattern {pattern} (bucket {bucket_id}).", "green"))
        print(colored(f"Total failed ids found: {len(failed_ids)} for pattern {pattern} (bucket {bucket_id}).", "green"))
        destined_files = self.load_bucket(bucket_id)
        bucket_size = len(destined_files)
        files = []
        # extractions = [file for file in self.indir.resolve().glob(f"*_{pattern}.processed.tsv") if file.is_file()]
        for file in tqdm(destined_files):
            if extract_id(file) in extracted_ids:
                extracted_file = self.indir.joinpath(MindiTool.extract_name(file) + f"_{pattern}.processed.tsv")
                if extracted_file.is_file():
                    files.append(extracted_file)
        if bucket_is_complete:
            print(colored(f"Bucket {bucket_id} is Complete!", "green"))
            assert bucket_size == len(files) + len(failed_ids), f"Bucket {bucket_id} is complete, but there are missing files!"
        else:
            print(colored(f"Bucket {bucket_id} is In-Complete!", "yellow"))
            assert bucket_size > len(failed_ids) + len(files), f"Invalid number of files detected."
        print(colored(f"Total extracted files detected: {len(files)} (bucket {bucket_id}).", "green"))
        return files

    def load_empty(self, bucket_id: str, pattern: str) -> list[str]:
        empty_log = self.indir.joinpath("merged_bulk", f"bucket_{bucket_id}_mode_{pattern}_empty.tsv")
        empty_files = []
        with empty_log.open("r", encoding="UTF-8") as f:
            for line in f:
                line = line.strip()
                empty_files.append(MindiTool.extract_id(line))
        print(colored(f"Total empty files detected: {len(empty_files)} for pattern {pattern} (bucket_id {bucket_id}).", "green"))
        return empty_files

    def merge_bucket(self, bucket_id: str, 
                    pattern: str, 
                    partition_col: Optional[str] = None,
                    min_partition: Optional[int] = None,
                    max_partition: Optional[int] = None,
                    assembly_summary: Optional[str] = None,
                    multiplier: float = 1e3) -> None:
        """Merge coordinates and calculate density."""
        # files = self.load_bucket(bucket_id)
        files = self.load_log(bucket_id, pattern)
        merge_outdir = self.merge_outdir.joinpath(pattern)
        merge_outdir.mkdir(exist_ok=True)
        outfile_normal = merge_outdir.joinpath(f"bucket_{bucket_id}_{pattern}_raw.tsv.gz")
        outfile = merge_outdir.joinpath(f"bucket_{bucket_id}_{pattern}_merged.tsv.gz")
        outfile_density = merge_outdir.joinpath(f"bucket_{bucket_id}_{pattern}_density_merged.tsv.gz")
        outfile_empty = merge_outdir.joinpath(f"bucket_{bucket_id}_{pattern}_empty.tsv.gz")
        print(f"Outsourcing files to --> `{merge_outdir}` with pattern {pattern} and partition column {partition_col}...")
        empty_files = set(self.load_empty(bucket_id, pattern))
        found_empty = set()
        density_df = defaultdict(list)
        # Process sub-divisions of dataframes
        if pattern == "MR":
            subsets = [pattern, "HDNA", "GT"]
        else:
            subsets = [pattern]
        # START >
        with gzip.open(outfile_normal, "wt") as f1, \
             gzip.open(outfile, "wt") as f2, \
             gzip.open(outfile_empty, "wt") as f3:
            f3.write("#assembly_accession\n")
            if not partition_col:
                f2.write("seqID\tstart\tend\toverlapping\tpattern\taccession_id\n")
            else:
                f2.write("seqID\tstart\tend\toverlapping\tpattern\tpartition_col\taccession_id\n")
            for i, file in tqdm(enumerate(files, 1), total=len(files)):
                # dest_file = MindiTool.extract_name(file) + ".processed.merged.tsv"
                accession_id = MindiTool.extract_id(file)
                unique_partition = ["."]
                usecols = ["seqID", "start", "end"]
                df = pd.read_table(file)
                if partition_col:
                    usecols += [partition_col]
                    # unique_partition += list(df[partition_col].unique())
                    unique_partition += list(range(min_partition, max_partition))
                shape_before = df.shape[0]
                df_collection = dict()
                df.loc[:, "sequence"] = df["sequence"].str.lower()
                if pattern == "STR":
                    df = df[(df["consensus_repeats"] >= self.min_consensus_repeats) \
                                & (df["sequence_length"] >= self.min_sequence_length) \
                                & (df["sru"] <= 9) \
                                & (df["sru"] >= 1)
                            ].reset_index(drop=True)
                    df_collection["STR"] = df
                elif pattern == "IR" or pattern == "DR":
                    if isinstance(self.min_arm_length ,int) and isinstance(self.max_spacer, int):
                        df = df[(df["spacer_length"] <= self.max_spacer) & (df["arm_length"] >= self.min_arm_length)].reset_index(drop=True)
                    elif isinstance(self.min_arm_length ,int):
                        df = df[(df["arm_length"] >= self.min_arm_length)].reset_index(drop=True)   
                    elif isinstance(self.max_spacer, int):
                        df = df[(df["spacer_length"] <= self.max_spacer)].reset_index(drop=True)
                    df_collection[pattern] = df
                elif pattern == "MR":
                    df = df[(df["spacer_length"] <= self.max_spacer) & (df["arm_length"] >= self.min_arm_length)].reset_index(drop=True)
                else:
                    raise ValueError(f"Invalid pattern detected: {pattern}.")
                if df.shape[0] < shape_before:
                    print(colored(f"Warning! Df shape was altered from {shape_before} to {df.shape[0]}.", "red"))
                if df.shape[0] == 0:
                    # assert accession_id in empty_files, f"Accession {accession_id} was not found empty but it is?"
                    if accession_id not in empty_files:
                        print(colored(f"Warning! Accession {accession_id} was not found empty but it is?", "red"))
                    found_empty.add(accession_id)
                    f3.write(f"{accession_id}\n")
                    for subset in subsets:
                        for partition in unique_partition:
                            density_df["accession_id"].append(accession_id)
                            density_df["total_bp"].append(0)
                            density_df["pattern"].append(pattern)
                            density_df["bucket_id"].append(bucket_id) # irrelevant metadata
                            density_df["partition_col"].append(partition_col)
                            density_df["partition"].append(partition)
                    continue
                # In case of MR, calculate H-DNA
                if pattern == "MR":
                    assert df[df["arm_length"] != df["sequence_of_arm"].apply(len)].shape[0] == 0
                    df.loc[:, "arm_length"] = df["sequence_of_arm"].apply(len)
                    df.loc[:, "ga_proportion"] = (df["sequence_of_arm"].str.count("g|a")).div(df["arm_length"])
                    df.loc[:, "gt_proportion"] = (df["sequence_of_arm"].str.count("g|t")).div(df["arm_length"])
                    df.loc[:, "at_proportion"] = (df["sequence_of_arm"].str.count("a|t")).div(df["arm_length"])
                    df.loc[:, "square_free"] = df["sequence"].apply(StreamAndMerge.is_square_free).astype(int)
                    df.loc[:, "cubic_free"] = df["sequence"].apply(StreamAndMerge.is_cubic_free).astype(int)
                    # df.loc[: "sequence_STR_coverage"] = df["sequence"].apply(StreamAndMerge.detect_STR_coverage)
                    # df.loc[: "arm_STR_coverage"] = df["sequence_of_arm"].apply(StreamAndMerge.detect_STR_coverage)
                    # df.loc[: "spacer_STR_coverage"] = df["sequence_of_spacer"].apply(StreamAndMerge.detect_STR_coverage)
                    # Calculate H-DNA and GT threshold
                    # # # #
                    # H-DNA Calculation
                    df.loc[:, "is_HDNAmr"] = (((df["ga_proportion"] >= self.GA_threshold) | (df["ga_proportion"] < 1 - self.GA_threshold)) & (df["at_proportion"] < self.AT_threshold)).astype(int)
                    df.loc[:, "HDNAmr_strand"] = (df["ga_proportion"] >= 0.5).astype(int).apply(lambda x: "+" if x == 1 else "-")
                    # # # # 
                    # GT Calculation
                    df.loc[:, "is_GTmr"] = (((df["gt_proportion"] >= self.GT_threshold) | (df["gt_proportion"] < 1 - self.GT_threshold)) & (df["at_proportion"] < self.AT_threshold)).astype(int)
                    df.loc[:, "GTmr_strand"] = (df["gt_proportion"] >= 0.5).astype(int).apply(lambda x: "+" if x == 1 else "-")
                    # # # #
                    # # # # # # # ##
                    df_collection["MR"] = df
                    df_collection["HDNA"] = df[df["is_HDNAmr"] == 1].reset_index(drop=True)
                    df_collection["GT"] = df[df["is_GTmr"] == 1].reset_index(drop=True)
                elif pattern == "DR" or pattern == "IR":
                    pass
                    # df.loc[: "sequence_STR_coverage"] = df["sequence"].apply(StreamAndMerge.detect_STR_coverage)
                    # df.loc[: "arm_STR_coverage"] = df["sequence_of_arm"].apply(StreamAndMerge.detect_STR_coverage)
                    # df.loc[: "spacer_STR_coverage"] = df["sequence_of_spacer"].apply(StreamAndMerge.detect_STR_coverage)
                df.loc[:, "#assembly_accession"] = accession_id
                df.to_csv(f1, sep="\t", index=False, header=i==1)
                for subset, df in df_collection.items():
                    df = df[usecols]
                    for partition in unique_partition:
                        total_bp = 0
                        if partition != ".":
                            df_temp = df[df[partition_col] == partition].reset_index(drop=True)
                        else:
                            df_temp = df
                        if df_temp.shape[0] == 0:
                            density_df["accession_id"].append(accession_id)
                            density_df["total_bp"].append(total_bp)
                            density_df["pattern"].append(subset)
                            density_df["bucket_id"].append(bucket_id) # irrelevant metadata
                            density_df["partition_col"].append(partition_col)
                            density_df["partition"].append(partition)
                            continue
                        # save to disk and calculate total base pairs
                        with open(BedTool.from_dataframe(df_temp).sort().merge(c="4", o="count").fn, 
                                mode="r",
                                encoding="UTF-8") as g:
                            for line in g:
                                start, end = line.split("\t")[1:3]
                                start, end = int(start), int(end)
                                total_bp += end - start
                                f2.write(line.replace("\n", f"\t{subset}\t{partition}\t{accession_id}\n"))
                        # store density data
                        density_df["accession_id"].append(accession_id)
                        density_df["total_bp"].append(total_bp)
                        density_df["pattern"].append(subset)
                        density_df["bucket_id"].append(bucket_id) # irrelevant metadata
                        density_df["partition_col"].append(partition_col)
                        density_df["partition"].append(partition)
            density_df = pd.DataFrame(density_df)
            if assembly_summary is not None:
                assembly_summary = Path(assembly_summary).resolve()
                if assembly_summary.is_file():
                    # headers = list(pd.read_table("headers.txt").columns)
                    assembly_df = pd.read_table(assembly_summary)
                    density_df = density_df\
                                .merge(
                                            assembly_df,
                                            left_on="accession_id",
                                            right_on="#assembly_accession",
                                            how="left"
                                )
                    density_df.loc[:, f"density_ungapped_{pattern}"] = (density_df["total_bp"] * multiplier / density_df["genome_size_ungapped"]).round(3)
                    density_df.loc[:, f"density_{pattern}"] = (density_df["total_bp"] * multiplier / density_df["genome_size"]).round(3)
                else:
                    print(colored(f"Failed to merge with assembly summary. Reason: provided invalid assembly summary path {assembly_summary}.", "yellow"))
            density_df.to_csv(outfile_density, compression="gzip", mode="w", sep="\t", index=False, header=True)
            assert len(empty_files) <= len(found_empty), f"Empty list was not exhausted for pattern {pattern}. Remaining files: {len(empty_files)} (bucket {bucket_id})."
            assert empty_files.issubset(found_empty), f"Empty files were not found in the log for pattern {pattern}. Remaining files: {empty_files - found_empty} (bucket {bucket_id})."
            print(colored(f"Bucket merge {bucket_id} for pattern {pattern} has completed succesfully!", "green"))
            return

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("bucket_id", type=int, default=0)
    parser.add_argument("--schedule", type=str, default="schedule_tandem_extractions.json")
    parser.add_argument("--indir", type=str, default="extractions_STR_MR")
    # parser.add_argument("--log_indir", type=str, default="log_debug_nonbdna")
    parser.add_argument("--pattern", type=str, default="STR")
    parser.add_argument("--partition_col", type=str, default="sru")
    parser.add_argument("--min_partition", type=int, default=1)
    parser.add_argument("--max_partition", type=int, default=10)
    parser.add_argument("--assembly_summary", type=str, default="assembly_summary.txt.gz")
    parser.add_argument("--min_sequence_length", type=int, default=10)
    parser.add_argument("--max_spacer", type=int, default=7)
    parser.add_argument("--min_arm_length", type=int, default=10)
    parser.add_argument("--min_consensus_repeats", type=int, default=3)
    parser.add_argument("--GA_threshold", type=float, default=0.9)
    parser.add_argument("--GT_threshold", type=float, default=0.9)
    parser.add_argument("--AT_threshold", type=float, default=0.8)
    parser.add_argument("--multiplier", type=float, default=1e3)
    args = parser.parse_args()
    merger = StreamAndMerge(schedule=args.schedule, 
                                indir=args.indir, 
                                # log_indir=args.log_indir,
                                min_consensus_repeats=args.min_consensus_repeats,
                                min_arm_length=args.min_arm_length,
                                max_spacer=args.max_spacer,
                                min_sequence_length=args.min_sequence_length,
                                GA_threshold=args.GA_threshold,
                                GT_threshold=args.GT_threshold,
                                AT_threshold=args.AT_threshold)
    # logged_files = merger.load_log(args.bucket_id, args.pattern)
    # empty_files = merger.load_empty(args.bucket_id, args.pattern)
    merger.merge_bucket(bucket_id=args.bucket_id, 
                        pattern=args.pattern, 
                        partition_col=args.partition_col, 
                        min_partition=args.min_partition, 
                        max_partition=args.max_partition,
                        assembly_summary=args.assembly_summary,
                        multiplier=args.multiplier)
