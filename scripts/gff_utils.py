# Yet another GFF parser
# Coverage Pipeline - IGS Lab
# Nikol <3

# CTTCGGGGAAAGGGAAAAGGGGGGAAAGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGGCTCGGGGAAAGGGAAAAGGGGGACATGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGTTTTGGGGAAAGGGAAAAGGGGGTTTTGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGTTTGGGGGAAAGGGAAAAGGGGGACATGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGACGCGGGGAAAGGGAAAAGGGGGTCAAGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGGGGCGGGGAAAGGGAAAAGGGGGGAAAGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGTATAGGGGAAAGGGAAAAGGGGGACGCGGGGAAAGGGAAAAGGGGGTAAAGGGGAAAGGGAAAAGGGGGTTTTGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGTATCGGGGAAAGGGAAAAGGGGGTATAGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGTAAAGGGGAAAGGGAAAAGGGGGTTTGGGGGAAAGGGAAAAGGGGGCAGAGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGTTGTGGGGAAAGGGAAAAGGGGGTATAGGGGAAAGGGAAAAGGGGGGGGCGGGGAAAGGGAAAAGGGGGGTTGGGGGAAAGGGAAAAGGGGGTTTGGGGGAAAGGGAAAAGGGGGCAGAGGGGAAAGGGAAAAGGGGGTATCGGGGAAAGGGAAAAGGGGGTTTGGGGGAAAGGGAAAAGGGGGCCCTGGGGAAAGGGAAAAGGGGGAGCTGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGTCAAGGGGAAAGGGAAAAGGGGGACGCGGGGAAAGGGAAAAGGGGGTTGTGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGGGGCGGGGAAAGGGAAAAGGGGGGAAAGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGGGTAGGGGAAAGGGAAAAGGGGGGTTGGGGGAAAGGGAAAAGGGGGTTTGGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGCCCAGGGGAAAGGGAAAAGGGGGGAAAGGGGAAAGGGAAAAGGGGGACGCGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGTAAAGGGGAAAGGGAAAAGGGGGTTTGGGGGAAAGGGAAAAGGGGGCAGAGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGTTGTGGGGAAAGGGAAAAGGGGGTATAGGGGAAAGGGAAAAGGGGGGGGCGGGGAAAGGGAAAAGGGGGGTTGGGGGAAAGGGAAAAGGGGGTTTGGGGGAAAGGGAAAAGGGGGCAGAGGGGAAAGGGAAAAGGGGGTATAGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGTATCGGGGAAAGGGAAAAGGGGGTATAGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGACATGGGGAAAGGGAAAAGGGGGACGCGGGGAAAGGGAAAAGGGGGTTTGGGGGAAAGGGAAAAGGGGGTATAGGGGAAAGGGAAAAGGGGGACATGGGGAAAGGGAAAAGGGGGTATCGGGGAAAGGGAAAAGGGGGACGCGGGGAAAGGGAAAAGGGGGTAAAGGGGAAAGGGAAAAGGGGGTATAGGGGAAAGGGAAAAGGGGGCGCTGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGCAGAGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGTGTAGGGGAAAGGGAAAAGGGGGACGCGGGGAAAGGGAAAAGGGGGGGGCGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGCAGAGGGGAAAGGGAAAAGGGGGCGCTGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGGTTGGGGGAAAGGGAAAAGGGGGTTTGGGGGAAAGGGAAAAGGGGGCAGAGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGGTTGGGGGAAAGGGAAAAGGGGGTTTTGGGGAAAGGGAAAAGGGGGTTGTGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGGTTGGGGGAAAGGGAAAAGGGGGCAGAGGGGAAAGGGAAAAGGGGGTGCAGGGGAAAGGGAAAAGGGGGTATATATAAAGGGGAAAGGGAAAAGGGGGCCCTGGGGAAAGGGAAAAGGGGGACGCGGGGAAAGGGAAAAGGGGGCAGAGGGGAAAGGGAAAAGGGGGTTTTGGGGAAAGGGAAAAGGGGGTATCGGGGAAAGGGAAAAGGGGGAGTAGGGGAAAGGGAAAAGGGGGGGCTGGGGAAAGGGAAAAGGGGGAGAT

from pathlib import Path
import time
import argparse
import os
from tqdm import tqdm
import re
import csv
from functools import partial
import gzip
from termcolor import colored
import json
import logging
import polars as pl
from scheduling import MiniBucketScheduler
import threading
import pysam
import attr
from attr import field
from Bio import SeqIO
from Bio.Seq import Seq
from typing import Optional, ClassVar, Iterator, Iterable
from pybedtools import BedTool

class PolarityError(Exception):

    def __init__(self, msg: str) -> None:
        super().__init__(msg)

@attr.s(frozen=True, slots=True)
class Strands:

    strands: tuple[str, str] = field(init=False, factory=lambda : ("+", "-"))

    def __iter__(self) -> Iterator[str]:
        for strand in self.strands:
            yield strand

def parse_fasta(fasta: str) -> Iterator[tuple[str, str]]:
    fasta = Path(fasta).resolve()
    if fasta.name.endswith(".gz"):
        f = gzip.open(fasta, "rt")
    else:
        f = open(fasta, mode="r", encoding="utf-8")
    for record in SeqIO.parse(f, "fasta"):
        yield str(record.id), str(record.seq)
    f.close()

def extract_id(accession: str) -> str:
    accession = Path(accession).name
    if accession.count('_') > 2:
        return '_'.join(accession.split('_')[:2])
    return accession.split('.gff')[0]
    # return '.'.join(accession.split('.')[:2])

def extract_name(accession: str, suffix: str = '.tsv') -> str:
    return Path(accession).name.split(f'{suffix}')[0]

def evaluate_motif_strand(mode: str, sequence: str) -> str:
    if mode == "GC":
        # for GQuadruplex
        strand = "+" if sequence.count("G") >= sequence.count("C") else "-"
    elif mode == "GA":
        # for H-DNA
        strand = "+" if sequence.count("G") + sequence.count("A") >= len(sequence) else "-"
    else:
        raise ValueError(f"Invalid mode `{mode}`.")
    return strand

def merge(df: pl.DataFrame, 
            col: str,
            new_columns: Optional[str] = None,
            faidx: Optional[str] = None,
            merge_cols: Optional[list[str]] = None, 
            merged_ops: Optional[list[str]] = None, 
            delim: str = "*") -> pl.DataFrame:
    partition_list = set(df[col])
    merged_df: list[pl.DataFrame] = []

    if merged_cols is None:
        merged_cols = ["3"]
    if merged_ops is None:
        merged_ops = ["count"]
    if new_columns is None:
        new_columns = ["seqID", "start", "end", "counts"]

    for partition_col in partition_list:
        merged_df_bed = BedTool.from_dataframe(
                                        df.filter(pl.col(col) == partition_col)
                                                .select(["seqID", "start", "end"])
                                                .to_pandas()
                                        )
        if faidx is not None:
            merged_df_bed = merged_df_bed.sort(faidx=faidx)
        else:
            merged_df_bed = merged_df_bed.sort()
        df_temp = pl.read_csv(
                            merged_df_bed.merge(c=merged_cols, 
                                                o=merged_ops, 
                                                delim=delim)
                                        .fn,
                            has_header=False,
                            new_columns=new_columns,
                            separator="\t",
                        )\
                        .with_columns(
                            pl.lit(compartment).alias(col)
                        )
        merged_df.append(df_temp)
    merged_df = pl.concat(merged_df)
    return merged_df


def merge_polarity(df: pl.DataFrame, 
                    col: str,
                   default_col: str = "phase", 
                   delim: str = "*",
                   faidx: Optional[str] = None) -> pl.DataFrame:
    partition_list = set(df[col])
    merged_df: list[pl.DataFrame] = []
    # keep only valid polarity values
    df = df.filter((pl.col("strand") == "+") | (pl.col("strand") == "-"))

    for partition_col in partition_list:
        df_bed = BedTool.from_dataframe(
                                        df.filter(pl.col(col) == partition_col)
                                          .select(["seqID", "start", "end", "compartment", default_col, "strand"])
                                          .to_pandas()
                                    )
        if faidx is not None:
            df_bed = df_bed.sort(faidx=faidx)
        else:
            df_bed = df_bed.sort()
        df_temp = pl.read_csv(
                         df_bed.merge(c=["3", "4", "6"], 
                                             s=True,
                                             o=["count", "distinct", "distinct"], 
                                             delim=delim
                                             ).fn,
                         has_header=False,
                         new_columns=["seqID", "start", "end", "counts", col, "strand"],
                         separator="\t",
                        )
        merged_df.append(df_temp)
    merged_df = pl.concat(merged_df)
    return merged_df


@attr.s(slots=True, kw_only=True)
class GFFExtractor:
    
    cutoff: int = field(converter=int, default=0)
    compartments: Optional[list[str]] = field(factory=lambda : ["gene", 
                                                                "CDS", 
                                                                "intron", 
                                                                "five_prime_UTR", 
                                                                "three_prime_UTR", 
                                                                "exon", 
                                                                "region", 
                                                                "pseudogene", 
                                                                "genes", 
                                                                "pseudogenes"])
    names_mapping: dict[str, str] = field(factory=lambda : {
                                                            "gene": "Gene",
                                                            "exon": "Exon",
                                                            "region": "Genome",
                                                            "intron": "Intron",
                                                            "pseudogene": "Pseudogene",
                                                            })
    sequence_report: str = field(default=None)
    sequence_report_df: pl.DataFrame = field(init=False)
    promoter_kb: int = field(default=1_000, validator=attr.validators.instance_of(int), converter=int)
    terminator_kb: int = field(default=1_000, validator=attr.validators.instance_of(int), converter=int)
    delim: str = field(default="*", converter=str)
    GFF_FIELDS: ClassVar[list[str]] = ["seqID", "source", "compartment", "start", "end", "score", "strand", "phase", "attributes"]
    COVERAGE_FIELDS: ClassVar[list[str]] = ["totalHits", "overlappingBp", "compartmentLength", "coverage"]
    strands: Strands = field(factory=lambda : Strands()) 

    def __attrs_post_init__(self) -> None:
        if self.sequence_report is not None:
            self.sequence_report = Path(self.sequence_report).resolve()
            self.sequence_report_df = pl.read_csv(self.sequence_report, separator="\t")
        else:
            self.sequence_report_df = None

    @staticmethod
    def parse_attributes(attributes: str) -> dict[str, str]:
        attributes = attributes.split(";")
        attributes_map = {}
        for attr in attributes:
            key, value = attr.split("=")
            attributes_map[key] = value
        return attributes_map

    @staticmethod 
    def parse_biotype(attributes: str) -> str:
        attributes = GFFExtractor.parse_attributes(attributes)
        return attributes.get("biotype", attributes.get("gene_biotype", "."))

    @staticmethod 
    def parse_ID(attributes: str) -> str:
        attributes = GFFExtractor.parse_attributes(attributes)
        return attributes.get("ID", ".")

    @staticmethod
    def marshal(attributes: dict) -> str:
        return ";".join(f"{key}={value}" for key, value in attributes.items())
    
    @staticmethod
    def parse_gene_biotype(gff_table: pl.DataFrame, overload_biotype: bool = True) -> pl.DataFrame:
        return gff_table.with_columns(
                                pl.when(pl.col("compartment") == "gene")
                                .then(
                                    pl.when(pl.col("attributes").str.extract(r"gene_biotype=([^;]+)", 1) == "protein_coding")
                                    .then(pl.lit("protein_coding"))
                                    .otherwise(pl.lit("non_coding"))
                                    )
                                .otherwise(pl.lit("."))
                                .alias("biotype")
                            ) if overload_biotype else gff_table.with_columns(
                                                                            pl.when(pl.col("compartment") == "gene")
                                                                            .then(
                                                                                pl.col("attributes").str.extract(r"gene_biotype=([^;]+)", 1)
                                                                            )
                                                                            .otherwise(pl.lit("."))
                                                                            .alias("biotype")
                                                                    )

    @staticmethod
    def parse_compartment_ID(gff_table: pl.DataFrame) -> pl.DataFrame:
            return gff_table.with_columns(
                            pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("compartmentID")
                        )

    @staticmethod
    def parse_parent_ID(gff_table: pl.DataFrame) -> pl.DataFrame:
            return gff_table.with_columns(
                            pl.col("attributes").str.extract(r"Parent=([^;]+)", 1).alias("parentID")
                        )

    def read_gff(self, gff_file: str,
                        change_names: bool = True,
                        parse_biotype: bool = True,
                        overload_biotype: bool = True,
                        parse_ID: bool = False,
                        merge_with_parent: bool = False,
                        join_region: bool = False,
                        end_one_base: bool = False,
                        parse_length: bool = False,
                        parse_parentID: bool = False,
                        use_refseq: bool = False,
                        use_names: bool = False,
                        replace_pseudogene_with_gene: bool = True,
                        parse_promoters: bool = False,
                        parse_terminators: bool = False,
                 ) -> Optional[pl.DataFrame]:
        try:
            gff_table = pl.read_csv(gff_file, 
                                separator="\t",
                                has_header=False,
                                comment_prefix="#",
                                new_columns=GFFExtractor.GFF_FIELDS
                                )\
                        .with_columns(
                                        # 1-base --> 0-base
                                        pl.col("start") - 1,
                                        pl.col("compartment").replace({
                                                                       "genes": "gene", 
                                                                       "pseudogenes": "pseudogene"
                                                                       })
                            )
            if end_one_base:
                gff_table = gff_table.with_columns(
                                            pl.col("end") - 1
                                        )
        except pl.exceptions.NoDataError as e:
            logging.error(f"GFF file `{gff_file}` is empty.")
            logging.error(e)
            return

        if self.sequence_report_df is not None:
            if use_refseq :
                naming_column = "RefSeq seq accession"
            else:
                naming_column = "GenBank seq accession"

            seqID_mapping = dict(zip(self.sequence_report_df[naming_column],
                                     self.sequence_report_df["UCSC style name"]
                                     ))
            gff_table = gff_table.with_columns(
                                    pl.col("seqID").replace_strict(seqID_mapping).alias("seqID")
                                )
        if parse_length:
            gff_table = gff_table.with_columns(
                                    (pl.col("end") - pl.col("start")).alias("compartment_length")
                                    )

        if join_region:
            # WARNING
            # Here we use the `region` to dedude the total Genome Size
            # However, there are multiple viral assemblies that have more than one region 
            # corresponding to the same sequence ID.
            # This generates duplicate entries that can distort the data with the merge below.
            # If you are unsure, please use the genome size from the assembly summary
            # or calculate it directly from the fasta file.
            # Here we offer this alternative, but also we deduplicate.
            regions_df = gff_table.filter(pl.col("compartment") == "region")\
                                  .select(["seqID", "start", "end"])\
                                  .rename(mapping={
                                      "start": "region_start",
                                      "end": "region_end"
                                      })\
                                  .unique(["seqID"], keep="first")
            if not regions_df["seqID"].value_counts().filter(pl.col("count") > 1).is_empty():
                print(colored("WARNING! Multiple regions detected when deriving the total genome size.\nPossibly the GFF file contains multiple regions for the same sequence ID.", "red"))
            gff_table = gff_table.join(
                                        regions_df,
                                        on="seqID",
                                        how="left",
                                    )
                                 # .unique(["seqID", "start", "end", "compartment", "strand"], 
                                 #       maintain_order=True,
                                 #        keep="first")

        if replace_pseudogene_with_gene:
            gff_table = gff_table.with_columns(
                                pl.col("compartment").replace("pseudogene", "gene").alias("compartment")
                                )
        if self.compartments:
            compartments = set(self.compartments)
            gff_table = gff_table.filter(pl.col("compartment").is_in(compartments))
        if parse_biotype:
            gff_table = GFFExtractor.parse_gene_biotype(gff_table, overload_biotype=overload_biotype)
        if parse_ID:
            gff_table = GFFExtractor.parse_compartment_ID(gff_table)
        if parse_parentID:
            gff_table = GFFExtractor.parse_parent_ID(gff_table)
            additional_columns = ["seqID", "start", "end", "compartmentID"]
            if parse_biotype:
                additional_columns.append("biotype")
            if parse_ID and merge_with_parent:
                gff_table = gff_table.join(
                                            gff_table.select(additional_columns),
                                            left_on=["seqID", "parentID"],
                                            right_on=["seqID", "compartmentID"],
                                            how="inner"
                                            )
        if change_names:
            gff_table = gff_table.with_columns(
                                compartment=pl.col("compartment").replace(self.names_mapping)
                                    )
        if parse_promoters:
            gff_table = pl.concat([
                            gff_table,
                            self.parse_promoters(
                                    gff_table,
                                    promoter_kb=self.promoter_kb,
                                    filter_on="gene",
                                    gene_cutoff=self.cutoff,
                                    use_names=use_names,
                                    alias="Promoter",
                                )
                            ])
        if parse_terminators:
            gff_table = pl.concat([
                            gff_table,
                            self.parse_terminators(
                                    gff_table,
                                    terminator_kb=self.terminator_kb,
                                    filter_on="gene",
                                    gene_cutoff=self.cutoff,
                                    use_names=use_names,
                                    alias="Terminator",
                                )
                            ])
        return gff_table

    def parse_terminators(self, gff_table: pl.DataFrame, 
                                terminator_kb: int = 1_000, 
                                filter_on: str = "gene", 
                                gene_cutoff: int = 0, 
                                alias: str = "Terminator",
                                use_names: bool = False) -> pl.DataFrame:
        # gene_cutoff = 3, equals including the stop codon in prokaryotes
        selection_cols = ["seqID", "terminator_start", "terminator_end", "gene_length", "biotype", "region_start", "region_end"]
        if use_names:
            selection_cols.insert(5, "compartmentID")
        selection_cols.insert(5, "strand")
        if filter_on == "gene":
            gff_table = gff_table.filter(
                    (pl.col("compartment") == "gene") | (pl.col("compartment") == "Gene") | (pl.col("compartment") == "pseudogene") | (pl.col("compartment") == "Pseudogene")
                    )
        elif filter_on == "exon":
            gff_table = gff_table.filter(
                    (pl.col("compartment") == "exon") | (pl.col("compartment") == "Exon")
                    )
        return gff_table.filter(pl.col("end") <= pl.col("region_end"), 
                                (pl.col("strand") == "+") | (pl.col("strand") == "-")
                                )\
                        .with_columns(
                                    [
                                    pl.when(pl.col("strand") == "+")
                                        .then(pl.col("end") - gene_cutoff)
                                        .otherwise(pl.max_horizontal(0, pl.col("start") - terminator_kb))
                                        .alias("terminator_start"),

                                    pl.when(pl.col("strand") == "+")
                                        .then(pl.min_horizontal(pl.col("end") + terminator_kb, pl.col("region_end")))
                                        .otherwise(pl.col("start") + gene_cutoff)
                                        .alias("terminator_end"),

                                    (pl.col("end") - pl.col("start")).alias("gene_length"),
                                    ]
                                )\
                        .filter(pl.col("terminator_start") < pl.col("terminator_end"))\
                        .select(selection_cols)\
                        .rename({
                                 "terminator_start": "start",
                                 "terminator_end": "end"
                                 })\
                        .with_columns(
                                pl.lit(alias).alias("compartment")
                                )

    def parse_promoters(self, gff_table: pl.DataFrame, 
                              promoter_kb: int = 1_000, 
                              filter_on: str = "gene", 
                              gene_cutoff: int = 0, 
                              use_names: bool = False,
                              alias: str = "Promoter"
                        ) -> pl.DataFrame:
        # gene_cutoff = 3, equals including the start codon in prokaryotes
        selection_cols = ["seqID", "promoter_start", "promoter_end", "biotype", "region_start", "region_end"]
        if use_names:
            selection_cols.insert(5, "compartmentID")
        selection_cols.insert(5, "strand")
        if filter_on == "gene":
            gff_table = gff_table.filter(
                    (pl.col("compartment") == "gene") | (pl.col("compartment") == "Gene") | (pl.col("compartment") == "pseudogene") | (pl.col("compartment") == "Pseudogene")
                    )
        elif filter_on == "exon":
            gff_table = gff_table.filter(
                    (pl.col("compartment") == "exon") | (pl.col("compartment") == "Exon")
                    )
        return gff_table.filter(pl.col("end") <= pl.col("region_end"), 
                                (pl.col("strand") == "+") | (pl.col("strand") == "-")
                                )\
                             .with_columns(
                                            [
                                            pl.when(pl.col("strand") == "+")
                                                .then(pl.max_horizontal(0, pl.col("start") - promoter_kb))
                                                .otherwise(pl.col("end") - gene_cutoff)
                                                .alias("promoter_start"),

                                            pl.when(pl.col("strand") == "+")
                                                .then(pl.col("start") + gene_cutoff)
                                                .otherwise(pl.min_horizontal(pl.col("end") + promoter_kb, pl.col("region_end")))
                                                .alias("promoter_end"),

                                            (pl.col("end") - pl.col("start")).alias("gene_length"),
                                            ]).select(selection_cols)\
                        .rename({
                                 "promoter_start": "start",
                                 "promoter_end": "end"
                                 })\
                        .with_columns(
                                pl.lit(alias).alias("compartment")
                            )

    @staticmethod
    def find_stacking(sequence: str, stacker: str = "A", offset: int = 0, minimum_required: int = 4, reverse: bool = False) -> dict:
        """Detects polyN signals in a sequence"""
        if len(stacker) > 1:
            raise ValueError(f"Invalid stacker `{stacker}`. Accepts only a single nucleotide.")
        if reverse:
            sequence = sequence[::-1]
        sequence_length = len(sequence)
        total_observed = 0
        for i in range(sequence_length):
            if sequence[i] == stacker:
                total_observed += 1
            else:
                if total_observed >= minimum_required:
                    # AAATTTTTTAAAAA
                    # 0123456789
                    # AAAAATTTTTTAAA
                    # 01234567890123
                    # length = 14
                    # start = 3, end = 9
                    # 
                    if not reverse:
                        start_index = offset + i - total_observed
                        end_index = offset + i
                    else:
                        start_index = sequence_length - i + offset
                        end_index = sequence_length - i + offset + total_observed
                    return {
                            f"start_{stacker}": start_index,
                            f"end_{stacker}": end_index,
                            f"repetitions_{stacker}": total_observed,
                        }
                total_observed = 0
        return {} 

    @staticmethod
    def get_sequence(df: pl.DataFrame, 
                    genome: str, 
                    cutoff: int = 0, 
                    ascending: bool = True,
                    stackers: Optional[str] = None,
                    minimum_required: int = 4,
                    reverse: bool = False,
                    ) -> pl.DataFrame:
        if not Path(genome).is_file():
            raise FileNotFoundError(f"Fasta file `{genome}` not found.")
        def _create_index(genome: str) -> str | Path:
            pysam.faidx(genome)
            suffix = Path(genome).suffix 
            faidx = Path(genome).parent.joinpath(Path(genome).name.replace(suffix, f"{suffix}.fai"))
            return faidx
        faidx = _create_index(genome=genome)
        df = df.with_columns(
                    (
                        pl.col("gene_length").cast(pl.Utf8) + ":" + pl.col("biotype") + ":" + pl.col("region_end").cast(pl.Utf8)
                    ).alias("gene_length")
            )
        df_bed = BedTool.from_dataframe(df.to_pandas()).sort(faidx=faidx)
        sequences = df_bed.sequence(fi=genome, s=True, tab=True, name=True)
        sequence_df = []
        nucleotides = {"A", "G", "C", "T"}
        with open(sequences.seqfn, mode="r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t", fieldnames=["header", "sequence"])
            for row in reader:
                header = row["header"]
                sequence = row["sequence"]
                if any(n not in nucleotides for n in sequence):
                    continue
                gc_content = sequence.count("G") + sequence.count("C")
                sequence_length = len(sequence)
                strand = re.search(r"\(([^\)]+)\)", header).group(1)
                info = header.split("::")
                gene_length = info[0].split(":")[0]
                biotype = info[0].split(":")[1]
                region_end = info[0].split(":")[2]
                if "(" in info[1]:
                    info[1] = info[1].split("(")[0]
                seqID = info[1].split(":")[0]
                start, end = info[1].split(":")[-1].split("-")
                start = int(start)
                end = int(end)
                # if strand == "-":
                # sequence = str(Seq(sequence).reverse_complement())
                sequence_data = {
                                    "seqID": seqID,
                                    "start": start,
                                    "end": end,
                                    "gene_length": gene_length,
                                    "biotype": biotype,
                                    "strand": strand,
                                    "region_end": region_end,
                                    "sequence": sequence,
                                    "gc_content": gc_content,
                                    "sequence_length": sequence_length,
                                    "gc_percent": 1e2 * gc_content / sequence_length,
                                }
                if stackers:
                    for stacker in stackers:
                        tail_repeats = GFFExtractor.find_stacking(sequence=sequence,
                                                                  stacker=stacker,
                                                                  offset=0,
                                                                  minimum_required=minimum_required,
                                                                  reverse=reverse
                                                                )
                        sequence_data = sequence_data | tail_repeats
                if cutoff > 0:
                    if ascending:
                        # for instance, in case of promoters, fetch start codon
                        border_sequence = sequence[-cutoff:]
                    else:
                        # case for terminators, fetch stop codon
                        border_sequence = sequence[:cutoff]
                    sequence_data.update({"border_sequence": border_sequence})
                sequence_df.append(sequence_data)
        sequence_df = pl.DataFrame(sequence_df)
        return sequence_df


    def merge(self, gff_table: pl.DataFrame, faidx: Optional[str] = None) -> pl.DataFrame:
        unique_compartments = set(gff_table["compartment"])
        merged_gff: list[pl.DataFrame] = []
        for compartment in unique_compartments:
            gff_table_bed = BedTool.from_dataframe(
                                        gff_table.filter(pl.col("compartment") == compartment)
                                        .select(["seqID", "start", "end"])
                                        .to_pandas()
                                        )
            if faidx is not None:
                gff_table_bed = gff_table_bed.sort(faidx=faidx)
            else:
                gff_table_bed = gff_table_bed.sort()
            gff_temp = pl.read_csv(
                             gff_table_bed.merge(c="3", 
                                                 o="count", 
                                                 delim=self.delim).fn,
                             has_header=False,
                             new_columns=["seqID", "start", "end", "counts"],
                             separator="\t",
                            )\
                            .with_columns(
                                pl.lit(compartment).alias("compartment")
                            )
            merged_gff.append(gff_temp)
        merged_gff = pl.concat(merged_gff)
        return merged_gff

    def merge_polarity(self, gff_table: pl.DataFrame, 
                       default_col: str = "phase", 
                       faidx: Optional[str] = None) -> pl.DataFrame:
        unique_compartments = set(gff_table["compartment"])
        merged_gff: list[pl.DataFrame] = []
        # keep only valid polarity values
        gff_table = gff_table.filter((pl.col("strand") == "+") | (pl.col("strand") == "-"))
        for compartment in unique_compartments:
            gff_table_bed = BedTool.from_dataframe(
                                        gff_table.filter(pl.col("compartment") == compartment)
                                        .select(["seqID", "start", "end", "compartment", default_col, "strand"])
                                        .to_pandas()
                                        )
            if faidx is not None:
                gff_table_bed = gff_table_bed.sort(faidx=faidx)
            else:
                gff_table_bed = gff_table_bed.sort()
            gff_temp = pl.read_csv(
                             gff_table_bed.merge(c=["3", "4", "6"], 
                                                 s=True,
                                                 o=["count", "distinct", "distinct"], 
                                                 delim=self.delim
                                                 ).fn,
                             has_header=False,
                             new_columns=["seqID", "start", "end", "counts", "compartment", "strand"],
                             separator="\t",
                            )
            merged_gff.append(gff_temp)
        merged_gff = pl.concat(merged_gff)
        return merged_gff

    def parse_coverage(self, gff_table: pl.DataFrame, 
                            extraction_table: pl.DataFrame,
                            group: bool = True,
                            polarity: bool = False,
                            faidx: Optional[str] = None) -> pl.DataFrame:
        if isinstance(gff_table, pl.DataFrame):
            gff_table = gff_table.to_pandas()

        if polarity and ("strand" not in extraction_table and "motif_strand" not in extraction_table):
            raise PolarityError(f"Failure to detect polarity. Column `strand` is missing from the extraction table.")
        if polarity and "strand" not in gff_table:
            raise PolarityError(f"Failure to detect polarity. Column `strand` is missing from the GFF table.")
        if polarity and "strand" in extraction_table.columns:
            extraction_table = extraction_table.rename({"strand": "motif_strand"})

        gff_table_bed = BedTool.from_dataframe(gff_table)
        if faidx is not None:
            gff_table_bed = gff_table_bed.sort(faidx=faidx)
        else:
            gff_table_bed = gff_table_bed.sort()

        if not polarity:
            if isinstance(extraction_table, pl.DataFrame):
                extraction_table = extraction_table.to_pandas()
            extraction_table_bed = BedTool.from_dataframe(extraction_table)
            if faidx is not None:
                extraction_table_bed = extraction_table_bed.sort(faidx=faidx)
            else:
                extraction_table_bed = extraction_table_bed.sort()

            coverage_df = pl.read_csv(
                                gff_table_bed.coverage(extraction_table_bed).fn,
                                has_header=False,
                                separator="\t",
                                new_columns=gff_table.columns.tolist() + GFFExtractor.COVERAGE_FIELDS
                                )\
                        .with_columns(
                                coverage=(1e6 * pl.col("coverage")),
                                atLeastOne=(pl.col("totalHits") > 0).cast(pl.Int32)
                        )\
                        .with_columns(
                                (pl.col("atLeastOne") * pl.col("counts")).alias("atLeastOneUnmerged")
                        )
        else:
            extraction_pos_charge_bed = BedTool.from_dataframe(extraction_table.filter(pl.col("motif_strand") == "+").to_pandas())
            extraction_neg_charge_bed = BedTool.from_dataframe(extraction_table.filter(pl.col("motif_strand") == "-").to_pandas())
            if faidx is not None:
                extraction_pos_charge_bed = extraction_pos_charge_bed.sort(faidx=faidx)
                extraction_neg_charge_bed = extraction_neg_charge_bed.sort(faidx=faidx)
            else:
                extraction_pos_charge_bed = extraction_pos_charge_bed.sort()
                extraction_neg_charge_bed = extraction_neg_charge_bed.sort()
            coverage_pos_charge_df = pl.read_csv(
                                                gff_table_bed.coverage(extraction_pos_charge_bed).fn,
                                                has_header=False,
                                                separator="\t",
                                                new_columns=gff_table.columns.tolist() + GFFExtractor.COVERAGE_FIELDS
                                )\
                                .with_columns(
                                        pl.lit("+").alias("motif_strand")
                                )
            coverage_neg_charge_df = pl.read_csv(
                                                gff_table_bed.coverage(extraction_neg_charge_bed).fn,
                                                has_header=False,
                                                separator="\t",
                                                new_columns=gff_table.columns.tolist() + GFFExtractor.COVERAGE_FIELDS
                                )\
                                .with_columns(
                                        pl.lit("-").alias("motif_strand")
                                )
            coverage_df = pl.concat([
                                    coverage_pos_charge_df,
                                    coverage_neg_charge_df
                                    ])\
                            .with_columns(
                                        coverage=(1e6 * pl.col("coverage")),
                                        atLeastOne=(pl.col("totalHits") > 0).cast(pl.Int32)
                                )\
                            .with_columns(
                                        (pl.col("atLeastOne") * pl.col("counts")).alias("atLeastOneUnmerged"),
                                        
                                        pl.when(pl.col("strand") == pl.col("motif_strand"))
                                            .then(pl.lit("Non-Template"))
                                            .otherwise(pl.lit("Template"))
                                            .alias("strand_polarity")
                            )
        if not polarity and group:
            coverage_df = coverage_df.group_by("compartment",
                                               maintain_order=True)\
                                    .agg(
                                            pl.col("totalHits").sum(),
                                            pl.col("compartmentLength").count().alias("totalCompartments"),
                                            pl.col("counts").sum().alias("totalCompartmentsUnmerged"),
                                            pl.col("overlappingBp").sum(),
                                            pl.col("compartmentLength").sum(),
                                            pl.col("coverage").mean().alias("avg_coverage"),
                                            pl.col("coverage").median().alias("median_coverage"),
                                            pl.col("coverage").std().alias("std_coverage"),
                                            pl.col("coverage").quantile(0.975).alias("q975_coverage"),
                                            pl.col("coverage").quantile(0.025).alias("q025_coverage"),
                                            pl.col("coverage").min().alias("min_coverage"),
                                            pl.col("coverage").max().alias("max_coverage"),
                                            pl.col("atLeastOne").sum().alias("atLeastOne"),
                                            pl.col("atLeastOneUnmerged").sum().alias("atLeastOneUnmerged"),
                                            (1e2 * pl.col("atLeastOne").mean()).alias("perc_atLeastOne"),
                                        )\
                                    .with_columns(
                                            (1e6 * pl.col("overlappingBp") / pl.col("compartmentLength")).alias("coverage"),
                                        )
        elif polarity and group:
            coverage_df = coverage_df.group_by(["compartment", "strand_polarity"],
                                               maintain_order=True)\
                                    .agg(
                                            pl.col("totalHits").sum(),
                                            pl.col("compartmentLength").count().alias("totalCompartments"),
                                            pl.col("counts").sum().alias("totalCompartmentsUnmerged"),
                                            pl.col("overlappingBp").sum(),
                                            pl.col("compartmentLength").sum(),
                                            pl.col("coverage").mean().alias("avg_coverage"),
                                            pl.col("coverage").median().alias("median_coverage"),
                                            pl.col("coverage").std().alias("std_coverage"),
                                            pl.col("coverage").quantile(0.975).alias("q975_coverage"),
                                            pl.col("coverage").quantile(0.025).alias("q025_coverage"),
                                            pl.col("coverage").min().alias("min_coverage"),
                                            pl.col("coverage").max().alias("max_coverage"),
                                            pl.col("atLeastOne").sum().alias("atLeastOne"),
                                            pl.col("atLeastOneUnmerged").sum().alias("atLeastOneUnmerged"),
                                            (1e2 * pl.col("atLeastOne").mean()).alias("perc_atLeastOne"),
                                            (1e2 * pl.col("atLeastOneUnmerged").mean()).alias("perc_atLeastOneUnmerged"),
                                        )\
                                    .with_columns(
                                            (1e6 * pl.col("overlappingBp") / pl.col("compartmentLength")).alias("coverage"),
                                        )
        return coverage_df
    
    def parse_UTR(self, gff_table: pl.DataFrame) -> pl.DataFrame:
        raise NotImplementedError("UTR parsing utility has yet to be implemented but it's on my TODO list.")

    def parse_introns(self, gff_table: pl.DataFrame, maintain_order: bool = True) -> pl.DataFrame:
        if "compartmentID" not in gff_table.columns:
            gff_table = GFFExtractor.parse_compartment_ID(gff_table)
        if "parentID" not in gff_table.columns:
            gff_table = GFFExtractor.parse_parent_ID(gff_table)
            # I need parent to estimate the compartment size
        exons_table = gff_table.filter((pl.col("compartment") == "exon") | (pl.col("compartment") == "Exon"))
        exons_table = (exons_table
                             .join(
                                    gff_table.select(["seqID", "start", "end", "compartmentID"]),
                                    left_on=["seqID", "parentID"],
                                    right_on=["seqID", "compartmentID"],
                                    how="inner",
                                    suffix="_transcript",
                              )
                             .sort(
                                        by=["seqID", "parentID", "start"],
                                        descending=True
                                )
                )
        introns_table = []
        print("STARTED")
        group_id = 0
        for _, exon_group in tqdm(exons_table.group_by(["seqID", "parentID"], 
                                                        maintain_order=maintain_order), total=180_000):
            group_id += 1
            # it must be sorted here
            exons_df = sorted([row for row in exon_group.iter_rows(named=True)], 
                                key=lambda row: int(row['start']),
                                reverse=False)
            # invariant information (seqID, strand, parentID)
            seqID = exons_df[0]["seqID"]
            strand = exons_df[0]["strand"]
            parentID = exons_df[0]["parentID"]
            total_exons = len(exons_df)
            total_introns = 0
            start_transcript = exons_df[0]["start_transcript"]
            end_transcript = exons_df[0]["end_transcript"]
            transcript_length = end_transcript - start_transcript
            if strand == "-":
                intron_pos = "five_prime_intron_length"
                reverse_intron_pos = "three_prime_intron_length"
            else:
                intron_pos = "three_prime_intron_length"
                reverse_intron_pos = "five_prime_intron_length"

            # ordering within each transcript
            order = 1
            for exon_idx, exon in enumerate(exons_df, 0):
                exon_length = exon["end"] - exon["start"]
                exon = exon | {
                               "length": exon_length, 
                               "start_transcript": start_transcript,
                               "end_transcript": end_transcript,
                               "transcript_length": transcript_length,
                               "total_exons": total_exons,
                               }
                intron_id = exon["compartmentID"].replace("exon", "intron")
                if "attributes" in exon:
                    intron_attributes = exon["attributes"].replace("exon", "intron")
                else:
                    intron_attributes = ""

                # INTRON before first exon (?) RARE <
                if exon_idx == 0:
                    intron_start = exon["start_transcript"]
                    intron_end = exon["start"]
                    intron_length = intron_end - intron_start
                    if intron_start < intron_end:
                        intron = {
                                  "seqID": seqID,
                                  "source": exon["source"],
                                  "compartment": "Intron",
                                  "start": intron_start,
                                  "end": intron_end,
                                  "score": exon["score"],
                                  "strand": strand,
                                  "phase": exon["phase"],
                                  "attributes": intron_attributes,
                                  "length": intron_length,
                                  "order": order,
                                  "group_id": group_id,
                                  "compartmentID": f"Mindi:{intron_id}",
                                  "start_transcript": start_transcript,
                                  "end_transcript": end_transcript,
                                  "transcript_length": transcript_length,
                                  "parentID": parentID,
                                  "previous_exon_length": -1,
                                  "next_exon_length": exon["end"] - exon["start"],
                                  "total_exons": total_exons,
                                  }
                        introns_table.append(intron)
                        order += 1
                    introns_table.append(exon | {
                                                 "order": order, 
                                                 "group_id": group_id, 
                                                 })
                    if len(introns_table) > 1 and introns_table[-2]["compartment"] == "Intron":
                        introns_table[-1].update({intron_pos: intron_length})

                    order += 1
                # INTRON between previous exon and current exon
                else:
                    intron_start = exons_df[exon_idx-1]["end"]
                    intron_end = exon["start"]
                    intron_length = intron_end - intron_start
                    if intron_start < intron_end:
                        intron = {
                                    "seqID": seqID,
                                    "source": exon["source"],
                                    "compartment": "Intron",
                                    "start": intron_start,
                                    "end": intron_end,
                                    "score": exon["score"],
                                    "strand": strand,
                                    "phase": exon["phase"],
                                    "attributes": intron_attributes,
                                    "length": intron_length,
                                    "order": order,
                                    "group_id": group_id,
                                    "compartmentID": f"Mindi:{intron_id}",
                                    "start_transcript": start_transcript,
                                    "end_transcript": end_transcript,
                                    "transcript_length": transcript_length,
                                    "parentID": parentID,
                                    "previous_exon_length": exons_df[exon_idx-1]["end"] - exons_df[exon_idx-1]["start"],
                                    "next_exon_length": exon["end"] - exon["start"],
                                    "total_exons": total_exons,
                                }
                        introns_table.append(intron)
                        order += 1
                    # ADD INTRON first because its between the LAST and the CURRENT exon

                    if exon_idx == total_exons - 1:
                        order = -1
                    introns_table.append(exon | {
                                                 "order": order, 
                                                 "group_id": group_id,
                                                 })
                    if introns_table[-2]["compartment"] == "Intron":
                        introns_table[-1].update({intron_pos: intron_length})
                        introns_table[-3].update({reverse_intron_pos: intron_length})
                    order += 1

            # INTRON after last exon (?) RARE <
            intron_start = exon["end"]
            intron_end = exon["end_transcript"]
            # INSERT exon first
            intron_length = intron_end - intron_start
            if intron_start < intron_end:
                intron = {
                           "seqID": seqID,
                           "source": exon["source"],
                           "compartment": "Intron",
                           "start": intron_start,
                           "end": intron_end,
                           "score": exon["score"],
                           "strand": strand,
                           "phase": exon["phase"],
                           "attributes": intron_attributes,
                           "length": intron_length,
                           "order": order,
                           "group_id": group_id,
                           "compartment": "Intron",
                           "compartmentID": f"Mindi:{intron_id}",
                           "start_transcript": start_transcript,
                           "end_transcript": end_transcript,
                           "transcript_length": transcript_length,
                           "parentID": parentID,
                           "previous_exon_length": exon["end"] - exon["start"],
                           "next_exon_length": -1,
                           "total_exons": total_exons,
                        }
                introns_table.append(intron)
        print("FINISHED!")
        introns_table = pl.DataFrame(introns_table)
        if introns_table.shape[0] == 0:
            print(colored(f"Introns table is empty!", "red"))
        return introns_table

    def drop_first_last(self, gff_table: pl.DataFrame, factorize: bool = False, compartment: str = "exon") -> pl.DataFrame:
        if "parentID" not in gff_table.columns:
            gff_table = GFFExtractor.parse_parent_ID(gff_table)
        gff_table = ( 
                     gff_table.filter((pl.col("compartment") == compartment) | (pl.col("compartment") == compartment.capitalize()))
                              .group_by("parentID", maintain_order=True)
                              .map_groups(lambda g: g.sort("start", descending=False)[1:-1])
                    )
        if factorize:
            gff_table = gff_table.with_columns(
                                    pl.col("parentID").rank("dense").alias("group_id") - 1
                                        )
        return gff_table

    def filter_TSS_TES_from_exons(self, gff_table: pl.DataFrame) -> pl.DataFrame:
        if "compartmentID" not in gff_table.columns:
            gff_table = GFFExtractor.parse_compartment_ID(gff_table)
        if "parentID" not in gff_table.columns:
            gff_table = GFFExtractor.parse_parent_ID(gff_table)
        TSS_TES = gff_table.filter((pl.col("compartment") == "gene") | (pl.col("compartment") == "Gene"))\
                           .select(["seqID", "start", "end", "compartmentID"])
        gff_table = (
                        gff_table.filter((pl.col("compartment") == "exon") | (pl.col("compartment") == "Exon"))
                                 .join(
                                        gff_table,
                                        left_on=["seqID", "parentID"],
                                        right_on=["seqID", "compartmentID"],
                                        how="left",
                                        suffix="_transcript",
                                    )
                                 .join(
                                       TSS_TES,
                                       left_on=["seqID", "parentID_transcript"],
                                       right_on=["seqID", "compartmentID"],
                                       how="left",
                                       suffix="_gene",
                                    )
                                 .sort(
                                        by=["seqID", "start"], descending=True
                                  )
                    )
        gff_table = gff_table.filter(
                        (pl.col("start") != pl.col("start_gene")) | (pl.col("end") != pl.col("end_gene"))
                )
        return gff_table

class CoverageExtractor(GFFExtractor):

    def __init__(self, out: Optional[str] = None,
                       schedule: Optional[str] = None, 
                       design: Optional[str] = None, 
                       float_precision: int = 2,
                       sleeping_time: float = 200,
                       biotypes: Optional[list[str]] = None,
                       compartments: Optional[list[str]] = None,
                       promoter_kb: int = 100,
                       terminator_kb: int = 100,
                       faidx: Optional[str] = None) -> None:
        super().__init__(compartments=compartments,
                         promoter_kb=promoter_kb,
                         terminator_kb=terminator_kb)
        if out:
            self.out = Path(out).resolve()
        else:
            self.out = out

        if schedule:
            self.schedule = Path(schedule).resolve()
            if not self.schedule.is_file():
                raise FileNotFoundError(f"Could not detect schedule file at `{self.schedule}`.")
        else:
            self.schedule = schedule
        self.faidx = faidx
        self.float_precision = float_precision

        if compartments is None:
            self.compartments = ["gene", "CDS", "exon", "region"]
        else:
            self.compartments = compartments
        if biotypes is None:
            # protein coding biotype --> parses only protein coding genes coverage
            # non coding biotype --> parses only non coding genges coverage
            # `.` biotype --> refers to generic information: parses genes, exons, CDS, etc. for total coverage
            self.biotypes = ["protein_coding", "non_coding", "."]
        else:
            self.biotypes = biotypes

        self.extractions = dict()
        self.sleeping_time = sleeping_time
        def _sniff_delimiter(design: str) -> str:
            """Determines if CSV file is tab or comma delimited."""
            if not Path(design).is_file():
                raise FileNotFoundError(f"File `{design}` was not found.")
            with open(design, mode="r", encoding="utf-8") as f:
                for line in f:
                    return "\t" if line.count("\t") > line.count(",") else ","
            print(colored(f"File `{design}` is empty! Returning empty delimiter.", "red"))
            return ""

        if design is None:
            return

        with open(design, mode="r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=_sniff_delimiter(design))
            for row in reader:
                accession_id = row['accession_id']
                extraction_file = row['extraction']
                self.extractions[accession_id] = extraction_file
        total_extractions = len(self.extractions)
        color = "red" if total_extractions == 0 else "green"
        print(colored(f"Total extractions {total_extractions} have been loaded.", color))
        if total_extractions == 0:
            raise ValueError(f"Could not detected any extractions from design file `{design}`.")

    def load_bucket(self, bucket_id: int) -> dict[str, list[str]]:
        with open(self.schedule, mode="r", encoding="utf-8") as f:
            bucket = json.load(f)[str(bucket_id)]
            print(colored(f"Total {len(bucket)} records have been loaded (bucket {bucket_id}).", "green"))
            return bucket

    @staticmethod
    def _sniff_delimiter(file: str) -> str:
        if not Path(file).is_file():
            raise FileNotFoundError(f"File `{file}` was not found.")
        if Path(file).name.endswith(".gz"):
            f = gzip.open(file, "rt")
        else:
            f = open(file, mode="r", encoding="utf-8")
        header = None
        for line in f:
            header = line
        f.close()
        if header is None:
            print(colored(f"File `{file}` is empty! Returning empty delimiter.", "red"))
            return ""
        return "\t" if header.count("\t") > header.count(",") else ","

    class _TrackProgress:
        def __init__(self, bucket_id: int, 
                           total_records: int,
                           sleeping_time: float) -> None:
            self.track = 0
            self.total_records = total_records
            self.bucket_id = bucket_id
            self.sleeping_time = sleeping_time

        def start(self) -> None:
            while True:
                progress = self.track * 1e2 / self.total_records
                logging.info(f"Current progress for bucket `{self.bucket_id}`: {progress:.2f}%.")
                time.sleep(self.sleeping_time)

    def create_schedule(self, input_file: str, total_buckets: int = 100) -> dict[str, list[str]]:
        scheduler = MiniBucketScheduler()
        schedule = scheduler.schedule_from_file(input_file, total_buckets=total_buckets)
        scheduler.saveas(schedule, dest=self.schedule)
        return schedule

    def process(self, gff_table: pl.DataFrame,
                      extraction_table: pl.DataFrame,
                      partition_col: Optional[str] = None,
                      unique_partitions: Optional[Iterable[str]] = None,
                      group: bool = True,
                      default_col: str = "phase",
                      polarity: bool = False,
                      return_df: bool = True,
                      accession_id: Optional[str] = None,
                    ) -> Optional[list[pl.DataFrame] | pl.DataFrame]:
        """"""
        coverage_df = []
        # split GFF file to process various biotypes
        for biotype in self.biotypes:
            if biotype != ".":
                gff_table_temp = gff_table.filter(pl.col("biotype") == biotype)
            else:
                gff_table_temp = gff_table
            if gff_table_temp.shape[0] == 0:
                continue
            
            # handle strand
            if not polarity:
                gff_table_merged = self.merge(gff_table=gff_table_temp, 
                                                faidx=self.faidx)
            else:
                gff_table_merged = self.merge_polarity(gff_table=gff_table_temp, 
                                                        default_col=default_col, 
                                                        faidx=self.faidx)

            # if partition column has been specified
            # then process the GFF coverage individually against the GFF features
            if partition_col:
                assert isinstance(unique_partitions, set)
                for partition in unique_partitions:
                    extraction_table_temp = extraction_table.filter(pl.col(partition_col) == partition)
                    coverage_table = self.parse_coverage(
                                            gff_table_merged,
                                            extraction_table_temp,
                                            group=group,
                                            polarity=polarity,
                                            faidx=self.faidx)\
                                        .with_columns(
                                                    pl.lit(biotype).alias("biotype"),
                                                    pl.lit(str(partition)).alias(partition_col)
                                        )
                    if accession_id:
                        coverage_table = coverage_table.with_columns(
                                                    pl.lit(accession_id).alias("#assembly_accession"),
                                                    )
                    coverage_df.append(coverage_table)

            # calculate the generic coverage against GFF features
            coverage_table = self.parse_coverage(
                                                gff_table_merged, 
                                                extraction_table, 
                                                group=group, 
                                                polarity=polarity,
                                                faidx=self.faidx)\
                                .with_columns(
                                        pl.lit(biotype).alias("biotype")
                                    )
            if accession_id:
                coverage_table = coverage_table.with_columns(
                                        pl.lit(accession_id).alias("#assembly_accession")
                                        )
            if partition_col:
                coverage_table = coverage_table.with_columns(
                                                    pl.lit("generic").alias(partition_col)
                                                )
            coverage_df.append(coverage_table)

        if return_df:
            coverage_df = pl.concat(coverage_df)
            if accession_id:
                coverage_df = coverage_df.sort(by=["#assembly_accession", "compartment", "biotype", "coverage"], 
                                                  descending=True)
            else:
                coverage_df = coverage_df.sort(by=["compartment", "biotype", "coverage"], descending=True)
        return coverage_df

    def process_bucket(self, bucket_id: int, 
                            partition_col: Optional[str] = None, 
                            group: bool = True,
                            polarity: bool = False,
                            default_col: str = "phase",
                            mode: str = "GC",
                            overload_biotype: bool = True,
                            replace_pseudogene_with_gene: bool = True,
                            use_terminators: bool = False,
                            use_promoters: bool = False,
                            sleeping_time: float = 200) -> None:
        bucket = self.load_bucket(bucket_id=bucket_id)
        logging.info(f"Processing bucket `{bucket_id}`...")
        print(f"Processing bucket `{bucket_id}`...")
        coverage_df = []
        unique_partitions = None
        tracker = CoverageExtractor._TrackProgress(bucket_id=bucket_id,
                                                    total_records=len(bucket),
                                                    sleeping_time=sleeping_time)
        daemon = threading.Thread(target=tracker.start, daemon=True, name="LoggingCoverageDaemon")
        daemon.start()
        empty_accessions = []
        for gff_file in bucket:
            print(gff_file)
            selection_items = ["seqID", "start", "end"]
            tracker.track += 1
            accession_id = extract_id(gff_file)
            extraction_filename = self.extractions.get(accession_id)
            if extraction_filename is None:
                logging.info(f"Failed to find extraction file for accession id `{accession_id}`.")
                continue
            delimiter = CoverageExtractor._sniff_delimiter(extraction_filename)
            try:
                extraction_table = pl.read_csv(extraction_filename, separator=delimiter)
                if extraction_table.shape[0] == 0:
                    is_empty = True
                else:
                    is_empty = False
            except pl.exceptions.NoDataError as e:
                logging.info(f"Accession ID `{accession_id}` was found empty.\n{e}")
                logging.warning(f"Failed to process gff file `{gff_file}`. Reason: empty dataset.")
                empty_accessions.append(accession_id)
                schema = {
                            "seqID": pl.Utf8,
                            "start": pl.Int32,
                            "end": pl.Int32
                        }
                if partition_col is not None:
                    schema.update({partition_col: pl.Utf8})
                if polarity:
                    schema.update({"motif_strand": pl.Utf8})
                extraction_table = pl.DataFrame([], schema=schema)
                is_empty = True

            if is_empty:
                empty_accessions.append(accession_id)
            extraction_table = extraction_table.rename({col: col[:1].lower() + col[1:] for col in extraction_table.columns})
            
            # check if chromosome column is set correctly
            if "seqID" not in extraction_table.columns and "chromosome" not in extraction_table.columns and "sequence_name" not in extraction_table:
                raise KeyError(f"Invalid column for chromosome ID.")
            if "seqID" not in extraction_table.columns and "chromosome" in extraction_table:
                extraction_table = extraction_table.rename({"chromosome": "seqID"})
            if "seqID" not in extraction_table.columns and "sequence_name" in extraction_table:
                extraction_table = extraction_table.rename({"sequence_name": "seqID"})

            # attempt to evaluate motif strand if the column does not exist
            if polarity and "motif_strand" not in extraction_table.columns \
                    and "strand" not in extraction_table.columns \
                    and "sequence" not in extraction_table.columns:
                raise KeyError(f"Strand column does not exist in extraction table, neither there is a logical way to derive it.")
            elif polarity and "strand" in extraction_table.columns:
                extraction_table = extraction_table.rename({"strand": "motif_strand"})
            elif polarity and "motif_strand" not in extraction_table.columns:
                extraction_table = extraction_table.with_columns(
                                                pl.col("sequence").str.to_uppercase()
                                                .map_elements(
                                                            partial(evaluate_motif_strand, mode), return_dtype=str
                                                            )
                                                .alias("motif_strand")
                                            )
                selection_items.append("motif_strand")

            if partition_col is not None:
                selection_items.append(partition_col)
                if partition_col not in extraction_table.columns:
                    raise KeyError(f"Invalid partition column. `{partition_col}` was not found in the dataframe.")
                unique_partitions = set(extraction_table[partition_col])
            extraction_table = extraction_table.select(selection_items)
            gff_table = self.read_gff(gff_file, 
                                      parse_biotype=True, 
                                      change_names=True,
                                      join_region=False,
                                      overload_biotype=overload_biotype,
                                      replace_pseudogene_with_gene=replace_pseudogene_with_gene,
                                      parse_promoters=use_promoters, 
                                      parse_terminators=use_terminators)
            if gff_table is None or gff_table.shape[0] == 0:
                logging.warning(f"Failed to process gff file `{gff_file}`. Reason: no records were found in the annotation file.")
                if is_empty:
                    empty_accessions = empty_accessions[:-1]
                continue

            coverage_table = self.process(
                                        gff_table, 
                                        extraction_table=extraction_table, 
                                        partition_col=partition_col,
                                        default_col=default_col,
                                        group=group,
                                        unique_partitions=unique_partitions,
                                        polarity=polarity,
                                        accession_id=accession_id
                                        )
            if isinstance(coverage_table, pl.DataFrame):
                coverage_df.append(coverage_table)
            elif isinstance(coverage_table, list):
                coverage_df.extend(coverage_table)
        coverage_df = pl.concat(coverage_df)
        if polarity:
            dest = Path(f"{out}/coverage_bucket_{bucket_id}_{polarity}_{mode}.txt").resolve()
            dest_empty = Path(f"{out}/empty_accessions_coverage_bucket_{bucket_id}_{polarity}_{mode}.txt").resolve()
        else:
            dest = Path(f"{out}/coverage_bucket_{bucket_id}_{polarity}.txt").resolve()
            dest_empty = Path(f"{out}/empty_accessions_coverage_bucket_{bucket_id}_{polarity}.txt").resolve()

        logging.info(f"Saving coverage output at `{dest}`... (bucket {bucket_id}).")
        print(colored(f"Saving coverage output at `{dest}`... (bucket {bucket_id}).", "green"))
        coverage_df.write_csv(dest,
                              separator="\t",
                              include_header=True,
                              float_precision=self.float_precision
                              )
        with open(dest_empty, mode="w", encoding="utf-8") as f:
            for accession_id in empty_accessions:
                f.write(accession_id + "\n")
        print(colored(f"Output has been saved succesfully (bucket {bucket_id}).", "green"))
        print(colored(f"Bucket `{bucket_id}` has been processed succesfully.", "green"))
        logging.info(f"Output has been saved succesfully (bucket {bucket_id}).")
        logging.info(f"Bucket `{bucket_id}` has been processed succesfully.")
        return

@attr.s(slots=True)
class Expander:

    window_size: int = field(converter=int, 
                             validator=attr.validators.instance_of(int), 
                             default=1_000)

    def expand_windows(self, df: pl.DataFrame, 
                            loci: str = "start", 
                            expand_below: bool = True, 
                            expand_above: bool = True,
                            maintain_start_end: bool = False
                        ) -> pl.DataFrame:
        if loci != "start" and loci != "end" and loci != "mid":
            raise ValueError(f"Invalid expanding location `{loci}`. Please use: `start`, `end` or `mid`.")
        df = df.filter((pl.col("strand") == "+") | (pl.col("strand") == "-"))
        if loci == "mid":
            df = df.with_columns(
                        ((pl.col("start") + pl.col("end")) // 2).alias("mid")
                    )
            if expand_below:
                df = df.drop("start").with_columns(
                        (pl.max_horizontal(pl.col("mid") - self.window_size, 0)).alias("start"),
                    )
            if expand_above:
                df = df.drop("end").with_columns(
                        (pl.col("mid") + self.window_size + 1).alias("end")
                )
        else:
            opposite_loci = "end" if loci == "start" else "start"
            df = df.with_columns(
                        pl.when(pl.col("strand") == "+")
                          .then(pl.max_horizontal(pl.col(loci) - self.window_size, 0))
                          .otherwise(pl.max_horizontal(pl.col(opposite_loci) - self.window_size, 0))
                          .alias("expanded_start"),

                        pl.when(pl.col("strand") == "+")
                          .then(pl.col(loci) + self.window_size + 1)
                          .otherwise(pl.col(opposite_loci) + self.window_size + 1)
                          .alias("expanded_end")
                          )
            if not maintain_start_end:
                df = df.drop(["start", "end"])\
                            .rename({
                                    "expanded_start": "start",
                                    "expanded_end": "end"
                                })
        df = df.select(["seqID", "start", "end"] + [col for col in df.columns if col != "seqID" and col != "start" and col != "end"])
        return df

    def expand_half_windows(self, df: pl.DataFrame, loci: str = "start", expand_upstream: bool = True) -> pl.DataFrame:
        if loci != "start" and loci != "end" and loci != "mid":
            raise ValueError(f"Invalid expanding location `{loci}`. Please use: `start`, `end` or `mid`.")
        df = df.filter((pl.col("strand") == "+") | (pl.col("strand") == "-"))
        if loci == "mid":
            df = df.with_columns(
                        ((pl.col("start") + pl.col("end")) // 2).alias("mid")
                    )
            if expand_upstream:
                df = df.drop("start").with_columns(
                                    (pl.max_horizontal(pl.col("mid") - self.window_size, 0)).alias("start"),
                            )\
                                .rename({"mid": "end"})
            else:
                df = df.drop("end").with_columns(
                                (pl.col("mid") + self.window_size).alias("end")
                        )\
                            .rename({"mid": "start"})
        else:
            opposite_loci = "end" if loci == "start" else "start"
            if expand_upstream:
                df = df.with_columns(
                            pl.when(pl.col("strand") == "+")
                              .then(pl.max_horizontal(pl.col(loci) - self.window_size, 0))
                              .otherwise(pl.col(opposite_loci))
                              .alias("expanded_start"),

                            pl.when(pl.col("strand") == "+")
                                .then(pl.col(loci))
                                .otherwise(pl.col(opposite_loci) + self.window_size)
                                .alias("expanded_end")

                          )\
                    .rename({
                            "start": "start_contracted",
                            "end": "end_contracted"
                            })\
                    .rename({
                            "expanded_start": "start",
                            "expanded_end": "end"
                        })
            else:
                df = df.with_columns(
                            pl.when(pl.col("strand") == "+")
                              .then(pl.col(loci))
                              .otherwise(pl.max_horizontal(pl.col(opposite_loci) - self.window_size, 0))
                              .alias("expanded_start"),

                            pl.when(pl.col("strand") == "+")
                                .then(pl.col(loci) + self.window_size)
                                .otherwise(pl.col(opposite_loci))
                                .alias("expanded_end")

                          )\
                    .rename({
                            "start": "start_contracted",
                            "end": "end_contracted"
                            })\
                    .rename({
                            "expanded_start": "start",
                            "expanded_end": "end"
                        })
            first_cols = ["seqID", "start", "end", "start_contracted", "end_contracted"]
            df = df.select(first_cols + [col for col in df.columns if col not in first_cols])
        return df


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Calculates the coverage of motifs across genomic subcompartments of interest.""")
    parser.add_argument("schedule", type=str, default="schedule.json")
    parser.add_argument("--out", type=str, default="coverage_out")
    parser.add_argument("--bucket_id", type=int, default=0)
    parser.add_argument("--design", type=str, default="design.csv")
    parser.add_argument("--genome", type=str, default=None)
    parser.add_argument("--faidx", type=str, default=None)
    parser.add_argument("--sleeping_time", type=float, default=60)
    parser.add_argument("--group", type=int, default=1, choices=[0, 1])
    parser.add_argument("--replace_pseudogene_with_gene", type=int, default=1, choices=[0, 1])
    parser.add_argument("--partition_col", type=str, default=None)
    parser.add_argument("--partition_groups", nargs="+", default=list(range(0, 9)), type=list)
    parser.add_argument("--overload_biotype", type=int, choices=[0, 1], default=1)
    parser.add_argument("--polarity", type=int, default=0, choices=[0, 1])
    parser.add_argument("--default_col", type=str, default="phase", choices=["score", "phase"])
    parser.add_argument("--strand_mode", type=str, default="GC", choices=["GC", "GA"])
    parser.add_argument("--use_promoters", type=int, choices=[0, 1], default=0)
    parser.add_argument("--use_terminators", type=int, choices=[0, 1], default=0)
    parser.add_argument("--promoter_kb", type=int, default=100)
    parser.add_argument("--terminator_kb", type=int, default=100)

    args = parser.parse_args()
    out = Path(args.out).resolve()
    if not out.is_dir():
        raise ValueError(f'Failure to create destination directory at `{out}`. Is it a nested structure?')

    schedule = args.schedule
    sleeping_time = args.sleeping_time
    group = bool(args.group)
    bucket_id = args.bucket_id
    design = args.design
    partition_col = args.partition_col
    partition_groups = args.partition_groups
    genome = args.genome
    overload_biotype = bool(args.overload_biotype)
    replace_pseudogene_with_gene = bool(args.replace_pseudogene_with_gene)
    faidx = args.faidx
    polarity = bool(args.polarity)
    default_col = args.default_col
    strand_mode = args.strand_mode
    use_promoters = bool(args.use_promoters)
    use_terminators = bool(args.use_terminators)
    promoter_kb = args.promoter_kb
    terminator_kb = args.terminator_kb

    if args.partition_col:
        out = out.joinpath(f"partition_{partition_col}")
        biolog_file = out.joinpath(f"biologs/coverage_bucket_{bucket_id}_polarity_{polarity}_partition_{partition_col}.log")
    else:
        biolog_file = out.joinpath(f"biologs/coverage_bucket_{bucket_id}_polarity_{polarity}.log")

    out.mkdir(exist_ok=True, parents=True)
    if genome and faidx is None:
        if not os.path.exists(genome):
            raise FileNotFoundError(f"Genome fasta file `{genome}` not found.")
        if not genome.endswith(".fa") or not genome.endswith(".fna"):
            raise ValueError(f"Genome fasta file `{genome}` has invalid suffix. It should either be `fna` or `fa`.")
        pysam.faidx(genome)
        # replace faidx with the index
        faidx = genome.replace(".fna", ".faidx.fna").replace(".fa", ".faidx.fa")
    
    # log destination
    biolog_file.parent.mkdir(exist_ok=True)
    logging.basicConfig(
                        level=logging.INFO,
                        filemode="a+",
                        format="%(asctime)s:%(levelname)s:%(message)s",
                        filename=biolog_file
                        )

    print(colored(f"Initializing extraction process.", "blue"))
    extractor = CoverageExtractor(out=out,
                                  schedule=schedule, 
                                  design=design,
                                  sleeping_time=sleeping_time,
                                  faidx=faidx)
    extractor.process_bucket(bucket_id=bucket_id, 
                             sleeping_time=sleeping_time, 
                             group=group, 
                             polarity=polarity,
                             overload_biotype=overload_biotype,
                             mode=strand_mode,
                             replace_pseudogene_with_gene=replace_pseudogene_with_gene,
                             default_col=default_col,
                             partition_col=partition_col)
    print(colored(f"Process has been completed succesfully.", "green"))
