# Mindi: A Non B-DNA extractor tool for data analysis

__author__ = "Nikol Chantzi"
__version__ = "1.0.1"
__email__ = "nmc6088@psu.edu"

import sys
import os
import shutil
import gzip
import logging
import uuid
import csv
from termcolor import colored
from pathlib import Path
import json
import tempfile
import subprocess
from dotenv import load_dotenv
from typing import ClassVar, Optional
import pandas as pd
from utils import parse_fasta
import attr 
from attr import field

csv.field_size_limit(sys.maxsize)

@attr.s
class MindiTool:

    TAIL_FIELDS: ClassVar[list[str]] = ["seqID", 
                                        "start", 
                                        "end", 
                                        "sequence", 
                                        "length", 
                                        "consensus", 
                                        "sru", 
                                        "consensus_repeats"]

    max_spacer: int = field(converter=int, default=100)
    fn: dict = field(factory=dict, init=False)
    fnp: dict = field(factory=dict, init=False)
    nonBDNA: str = field(init=False, default=None)
    nucleotides: ClassVar[set[str]] = {"a", "t", "g", "c"}
    tempdir: Path = field(default=Path().cwd(), init=True)
    HDNA_max_at_content: ClassVar[float] = 0.8
    HDNA_min_pyrine: ClassVar[float] = 0.9
    HDNA_max_pyrine: ClassVar[float] = 0.9

    def __attrs_post_init__(self) -> None:
        load_dotenv()
        if self.nonBDNA is None:
            self.nonBDNA = os.getenv('nonBDNA')
        if self.tempdir is None:
            self.tempdir = Path().cwd()
        else:
            self.tempdir = Path(self.tempdir).resolve()
            self.tempdir.mkdir(exist_ok=True)
        self.nonBDNA = self.nonBDNA
        
        with open(Path(__file__).resolve().parent.joinpath("defaults.json"),
                  mode="r",
                  encoding="utf-8") as f:
            data = json.load(f)
        self.defaults = data["params"]
        self.MINDI_FIELDS = data["headers"]
        if not Path(self.nonBDNA).is_file():
            raise FileNotFoundError(f"Invalid executable path {self.nonBDNA}.")

    @staticmethod
    def extract_name(accession: str) -> str:
        accession = Path(accession).resolve()
        accession_name = accession.name
        if accession_name.endswith('.fna.gz'):
            return accession_name.split('.fna.gz')[0]
        if accession_name.endswith('.fa.gz'):
            return accession_name.split('.fa.gz')[0]
        if accession_name.endswith('.fna'):
            return accession_name.split('.fna')[0]
        if accession_name.endswith('.fa'):
            return accession_name.split('.fa')[0]
        return accession_name
        
    @staticmethod
    def extract_id(accession: str) -> str:
        return '_'.join(Path(accession).name.split("_")[:2])

    def _generate_repeats(self, accession: str, pattern: list[str], **kwargs) -> "MindiTool":
        self.reset()
        accession = Path(accession).resolve()
        initial_accession = accession
        if not accession.is_file():
            raise FileNotFoundError(f'Unable to detect accession {accession}.')
        accession_name = MindiTool.extract_name(accession)
        # accession_tmp_dir = self.tempdir.joinpath(accession_name)
        accession_tmp_dir = tempfile.TemporaryDirectory(dir=self.tempdir,
                                                        prefix=f"{accession_name}_{'_'.join(pattern)}.",
                                                        )
        accession_tmp_dir_path = self.tempdir.joinpath(accession_tmp_dir.name)
        # if accession_tmp_dir.is_dir():
        #    shutil.rmtree(accession_tmp_dir)
        cur_dir = os.getcwd()
        # handle zipped files
        if accession.name.endswith(".gz"):
            # gzip compression strategy
            with tempfile.NamedTemporaryFile(dir=accession_tmp_dir_path,
                                             delete=False,
                                             suffix='.fna') as unzipped_tmp:
                with gzip.open(accession, 'rb') as handler:
                    unzipped_tmp.write(handler.read())
                accession = Path(unzipped_tmp.name).name
        os.chdir(accession_tmp_dir_path)
        rand_accession_name = str(uuid.uuid4())
        # unfortunately non b-gfa breaks with large names 
        if len(Path(accession).name) > 60:
            shutil.copy(accession, accession.name)
            accession = accession.name
        command = f"{self.nonBDNA} -seq {accession} -out {rand_accession_name} -skipAPR -skipSlipped -skipCruciform -skipTriplex -skipWGET"
        seen_modes = set()
        skipped = set()
        for m in pattern:
            if m == "IR":
                command += f" -minIRrep 10 -maxIRspacer 8"
            elif m == "MR":
                command += f" -minMRrep {kwargs['MR']['minrep']} -maxMRspacer {kwargs['MR']['maxspacer']}"
            elif m == "DR":
                command += f" -minDRrep {kwargs['DR']['minrep']} -maxDRrep {kwargs['DR']['maxrep']} -maxDRspacer {kwargs['DR']['maxspacer']}"
            seen_modes.add(m)
        for m in ["IR", "MR", "DR", "Z", "GQ", "STR"]:
            if m not in seen_modes:
                if m == "IR":
                    command += " -skipIR"
                elif m == "MR":
                    command += " -skipMR"
                elif m == "DR":
                    command += " -skipDR"
                elif m == "STR":
                    command += " -skipSTR"
                elif m == "Z":
                    command += " -skipZ"
                elif m == "GQ":
                    command += " -skipGQ"
                skipped.add(m)
        if not seen_modes:
            raise ValueError(f'Unknown pattern match `{pattern}`.')
        print(f"Input pattern detected: `{pattern}`.\nRunning command `{command}`...")
        try:
            _ = subprocess.run(command, shell=True,
                                    check=True,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL,
                            )
        except subprocess.CalledProcessError as e:
            print(f"Encountered a process error while processing accession `{initial_accession}` for pattern `{pattern}` using the command `{command}`.")
            print(e)
            os.chdir(cur_dir)
            accession_tmp_dir.cleanup()
            return None
            # raise e

        # check if operation was succesful
        for mode in pattern:
            if not Path(rand_accession_name + f"_{mode}.tsv").is_file():
                raise FileNotFoundError(f"Failed to extract {mode} for {accession}.")
        # continue with processing
        for mode in pattern:
            shutil.move(rand_accession_name + f"_{mode}.tsv", accession_name + f"_{mode}.tsv")
            shutil.move(rand_accession_name + f"_{mode}.gff", accession_name + f"_{mode}.gff")
            destination = self.tempdir.joinpath(accession_name + f'_{mode}.tsv')
            if destination.is_file():
                os.remove(destination)
            out_tsv = accession_tmp_dir_path.joinpath(accession_name + f'_{mode}.tsv')
            out_gff = accession_tmp_dir_path.joinpath(accession_name + f'_{mode}.gff')
            shutil.move(out_tsv, destination)
            # remove redundant files
            os.unlink(out_gff)
            # modify tmp file pointer
            self.fn[mode] = self.tempdir.joinpath(accession_name + f'_{mode}.tsv')
        accession_tmp_dir.cleanup()
        os.chdir(cur_dir)
        return self

    def get_right_arm(self, left_arm: str, mode: str) -> str:
        """
        Returns the right arm of a repeat given the left arm.
        :param left_arm: The left arm of the repeat.
        :param mode: The type of repeat (IR, DR, MR).
        :return: The right arm of the repeat.
        """
        if mode == "IR":
            return MindiTool.reverse(left_arm)
        elif mode == "DR":
            return left_arm
        elif mode == "MR":
            return left_arm[::-1]
        else:
            raise ValueError(f"Unknown mode {mode} for getting right arm.")

    @property
    def extractions(self) -> list[str]:
        return list(self.fn.keys())

    def reset(self) -> None:
        self.fn = {}
        self.fnp = {}

    def moveto(self, dest: str) -> None:
        dest = Path(dest).resolve()
        for m in self.extractions:
            shutil.move(self.fnp[m], dest)
            self.fnp[m] = dest.joinpath(self.fnp[m].name)
            if m == "IR" or m == "DR" or m == "MR":
                shutil.move(self.fnp[m], dest)
            else:
                shutil.move(self.fn[m], dest)

    def to_dataframe(self, mode: str, usecols: bool = True) -> pd.DataFrame:
        # STILL EXPERIMENTAL
        if mode == 'RE' or mode == "TAIL":
            tmp_file = self.fn[mode]
        else:
            tmp_file = self.fnp[mode]
        if usecols:
            mindi_frame = pd.read_table(tmp_file)
            if mode == "IR" or mode == "DR" or mode == "MR":
                mindi_frame.loc[:, "sequence_of_spacer"] = mindi_frame["sequence_of_spacer"].fillna(".")
        else:
            mindi_frame = pd.read_table(tmp_file, header=None, skiprows=1)
        return mindi_frame

    def set_tempdir(self, tempdir: os.PathLike[str]) -> None:
        print(f"Resetting tempdir to --> `{tempdir}`.")
        self.tempdir = Path(tempdir).resolve()
        self.tempdir.mkdir(exist_ok=True)

    def extract(self, accession: str, pattern: list[str] | str, **kwargs):
        # filter duplicates
        if isinstance(pattern, str):
            pattern = list(set(pattern.split(",")))
        if not isinstance(pattern, list):
            raise TypeError(f"Unknown type {type(pattern)}.")

        # process arguments/parameters for each mode
        for m in pattern:
            if m == 'IR' or m == 'DR' or m == 'MR':
                if m in kwargs:
                    kwargs[m]['minrep'] = kwargs[m].get('minrep', self.defaults[m]['minrep'])
                    kwargs[m]['maxspacer'] = kwargs[m].get('maxspacer', self.defaults[m]['maxspacer'])
                else:
                    kwargs[m] = {}
                    for attr, val in self.defaults[m].items():
                        kwargs[m]['maxspacer'] = self.defaults[m]['maxspacer']
                        kwargs[m]['minrep'] = self.defaults[m]['minrep']
            if m == 'DR':
                if m in kwargs:
                    kwargs[m]['maxrep'] = kwargs[m].get('maxrep', self.defaults[m]['maxrep'])
                else:
                    kwargs[m]['maxrep'] = self.defaults[m]['maxrep']
            if m == 'Z' or m == 'STR' or m == 'GQ':
                kwargs[m] = {}

        labels = {
                  "MR": "Mirror Repeats",
                  "IR": "Inverted Repeats",
                  "DR": "Direct Repeats",
                  "STR": "Short Tandem Repeats",
                  "GQ": "GQuadruplex",
                  "APR": "Poly-A Tract",
                  "Z": "Z-DNA",
                  }
        stars = "=" * 10
        for m in pattern:
           print(f"{stars}<{labels[m]}> parameters{stars}")
           for attr, val in kwargs[m].items():
               print(f"{attr}={val}")
           print(stars)

        result = self._generate_repeats(accession=accession, pattern=pattern, **kwargs)
        if result is None:
            return None
        return self
    
    @staticmethod
    def filter_HDNA(mindi_table: pd.DataFrame, threshold: float = 0.9) -> pd.DataFrame:
        mindi_table = mindi_table.copy()
        mindi_table["sequence"] = mindi_table["sequence"].str.lower()
        mindi_table.loc[:, "pyrine_proportion"] = (mindi_table["sequence"].str.count("a|g")).div(mindi_table["sequenceLength"])
        mindi_table.loc[:, "pyrimidine_proportion"] = 1.0 - mindi_table["pyrine_proportion"]
        # mindi_table.loc[:, "pyrimidine"] = (mindi_table["sequence"].str.count("c|t")).div(mindi_table["sequenceLength"])
        mindi_table.loc[:, "at_proportion"] = (mindi_table["sequence"].str.count("a|t")).div(mindi_table["sequenceLength"])
        mindi_table = mindi_table[(mindi_table["pyrimidine_proportion"] >= threshold) | (mindi_table["pyrine_proportion"] >= threshold)]\
                            .reset_index(drop=True)
        mindi_table.loc[:, "motif_strand"] = (mindi_table["pyrine_proportion"] >= 0.5).astype(int).apply(lambda x: "+" if x == 1 else "-")
        # H-DNA Template
        # non template --> ag rich strand
        # template --> ct rich strand
        return mindi_table

    @staticmethod
    def filter_GT(mindi_table: pd.DataFrame, threshold: float = 0.9) -> pd.DataFrame:
        mindi_table = mindi_table.copy()
        mindi_table["sequence"] = mindi_table["sequence"].str.lower()
        mindi_table.loc[:, "gt_proportion"] = (mindi_table["sequence"].str.count("t|g")).div(mindi_table["sequenceLength"])
        mindi_table.loc[:, "ca_proportion"] = 1.0 - mindi_table["gt_proportion"]
        # mindi_table.loc[:, "at_proportion"] = (mindi_table["sequence"].str.count("a|t")).div(mindi_table["sequenceLength"])
        mindi_table.loc[:, "at_proportion"] = (mindi_table["sequence"].str.count("a|t")).div(mindi_table["sequenceLength"])
        mindi_table = mindi_table[(mindi_table["gt_proportion"] >= threshold) | (mindi_table["ca_proportion"] >= threshold)]\
                            .reset_index(drop=True)
        mindi_table.loc[:, "motif_strand"] = (mindi_table["gt_proportion"] >= 0.5).astype(int).apply(lambda x: "+" if x == 1 else "-")
        return mindi_table

    def extract_HDNA(self) -> pd.DataFrame:
        # assert current extraction is of mode 'Mirror Repeat'
        if 'MR' not in self.fnp:
            raise ValueError("Cannot proceed with H-DNA filtering. Mirror Repeats have not been extracted.")
        mindi_table = pd.read_table(self.fnp["MR"])
        return MindiTool.get_HDNA(mindi_table)

    def cleanup(self) -> None:
        if self.fnp and Path(self.fnp).is_file():
            os.remove(self.fnp)
        else:
            print(colored(f"WARNING! Processed file {self.fnp} does not exist.", "red"))

    @staticmethod
    def complement(nucleotide: str) -> str:
        if nucleotide == "a":
            return "t"
        if nucleotide == "t":
            return "a"
        if nucleotide == "g":
            return "c"
        if nucleotide == "c":
            return "g"
        if nucleotide == "n":
            return "n"
        raise ValueError(f"Unknown nucleotide {nucleotide}.")

    @staticmethod
    def reverse(kmer: str) -> str:
        return ''.join(MindiTool.complement(c) for c in kmer)[::-1]

    def retrieve_spacer(self, sequence: str, 
                        sequence_of_arm: str,
                        sequence_length: int, 
                        repeat: int, mode: str) -> Optional[str]:
        """"""
        if mode == "IR" or mode == "DR" or mode == "MR":
            true_spacer_length = sequence_length - 2 * repeat
            right_arm = sequence[repeat+true_spacer_length:]
            if any(n not in MindiTool.nucleotides for n in right_arm) or any(n not in MindiTool.nucleotides for n in sequence_of_arm):
                logging.warning("Invalid record detected with non-AGCT nucleotide present in sequence of arm.")
                return 

            true_spacer = sequence[repeat:repeat+true_spacer_length]
            if len(true_spacer) == 0:
                true_spacer = "."
            # skip maximum spacer length
            if (isinstance(self.max_spacer, int) and true_spacer_length > self.max_spacer):
                # invalid record?
                raise ValueError(f"Invalid record detected. Max spacer was found {true_spacer_length} but cannot exceed {self.max_spacer}.")
            right_arm = self.get_right_arm(sequence_of_arm, mode=mode)
            if right_arm != sequence[repeat+true_spacer_length:]:
                logging.warning(f"Invalid record detected with arm {sequence_of_arm} not equal to the right arm {right_arm}.")
                return 
        else:
            # set null spacer for all other NON BDNA modes
            true_spacer = "."
        return true_spacer

    def sanitize(self, accession: os.PathLike[str], mode: str) -> None:
        print(colored(f"Validating accession `{accession}` for {mode=}...", "blue"))
        filtered_file = (
                    Path(self.fn[mode]).resolve()
                    .parent
                    .joinpath(Path(self.fn[mode])
                    .name.split(".tsv")[0] + ".processed.tsv")
                    )
        fout = filtered_file.open("w", encoding="UTF-8")
        fieldnames = ["seqID", "start", "end", "sequence_of_arm", "arm_length", 
                        "sequence_of_spacer", "spacer_length", 
                        "sequence", "sequence_length", 
                    ]
        if mode == "STR":
            fieldnames += ["sru", "consensus_repeats"]
        fieldnames += ["#assembly_accession"]
        writer = csv.DictWriter(fout, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        try:
            df = pd.read_table(self.fn[mode])
            if df.shape[0] == 0:
                raise pd.errors.EmptyDataError()
        except pd.errors.EmptyDataError:
            fout.close()
            self.fnp[mode] = filtered_file
            logging.info(f'Accession {accession} has passed all checks (mode {mode}).')
            logging.info(f"Accession {accession} is EMPTY!")
            return self.fnp[mode]
        total = 0
        validated = 0
        new_row = {}
        for seqID, seq in parse_fasta(accession):
            temp = df[df['Sequence_name'] == seqID]
            total += temp.shape[0]
            fasta_length = len(seq)
            for _, row in temp.iterrows():
                # make this zero base
                start = row['Start'] - 1
                end = row['Stop']
                sequence = row['Sequence']
                arm_length = int(row["Repeat"])
                sequence_of_arm = sequence[:arm_length]
                length = len(sequence)
                # 
                # # # # # # 
                if end > fasta_length:
                    diff = end - fasta_length
                    sequence = sequence[:len(sequence)-diff]
                    end = min(end, fasta_length)
                    length = len(sequence)
                    logging.warning(f"There was a change of record for {accession=} due to end of sequence {end} > {fasta_length}.")
                    # print(f"There was a change of record for {accession=} due to end of sequence {end} > {fasta_length} (mode {mode}).")
                # Process STR
                if mode == "STR":
                    # consensus_repeats, remainder = row['Repeated'][1:].split('+')
                    # consensus_repeats = int(consensus_repeats)
                    # remainder = int(remainder)
                    sru = len(sequence_of_arm)
                    sequence = sequence[:length - length%sru]
                    length -= length % sru
                    end = start + length
                    consensus_repeats = sequence.count(sequence_of_arm)
                    if sru * consensus_repeats < length:
                        sequence = sequence[:consensus_repeats * sru]
                        length = consensus_repeats * sru
                        end = start + length
                        # print(f"There was a change of record for {accession=} due to STR processing (mode {mode}).")
                        logging.warning(f"There was a change of record for {accession=} due to STR processing (mode {mode}).")
                    new_row["consensus_repeats"] = consensus_repeats
                    new_row["sru"] = sru
                # Process STR
                fasta_seq = seq[start: end].lower()
                # general assertions
                assert end - start == length, f"Length mismatch for {seqID} in {accession}. Expected {length} but got {end - start}."
                assert sequence == fasta_seq and len(fasta_seq) == length, f"Sequence mismatch for {seqID} in {accession}. Expected {sequence} but got {fasta_seq}."
                # update new row
                new_row["seqID"] = seqID
                new_row["start"] = start
                new_row["end"] = end
                new_row["sequence_of_arm"] = sequence_of_arm
                new_row["arm_length"] = arm_length
                # Spacer retrieval 
                sequence_of_spacer = self.retrieve_spacer(sequence=sequence, 
                                                            sequence_of_arm=sequence_of_arm, 
                                                            sequence_length=length, 
                                                            repeat=arm_length, 
                                                            mode=mode)
                assert sequence_of_spacer is not None, f"Invalid sequence of spacer for {seqID} in {accession}. Critical {row}."
                if sequence_of_spacer == ".":
                    spacer_length = 0
                else:
                    spacer_length = len(sequence_of_spacer)
                new_row["sequence_of_spacer"] = sequence_of_spacer
                new_row["spacer_length"] = spacer_length
                new_row["sequence"] = sequence
                new_row["sequence_length"] = length
                new_row["#assembly_accession"] = MindiTool.extract_id(accession)
                spacer_seq = new_row['sequence_of_spacer']
                if spacer_seq == ".":
                    spacer_seq = ""
                # mode specific checks
                if mode == "MR":
                    assert len(sequence_of_arm) * 2 + len(spacer_seq) == len(sequence) == end - start
                    assert sequence_of_arm + spacer_seq + sequence_of_arm[::-1] == fasta_seq
                elif mode == "DR":
                    assert len(sequence_of_arm) * 2 + len(spacer_seq) == len(sequence) == end - start
                    assert sequence_of_arm + spacer_seq + sequence_of_arm == fasta_seq 
                elif mode == "IR":
                    assert len(sequence_of_arm) * 2 + len(spacer_seq) == len(sequence) == end - start
                    assert sequence_of_arm + spacer_seq + MindiTool.reverse(sequence_of_arm) == fasta_seq
                elif mode == "STR":
                    consensus_repeats = new_row['consensus_repeats']
                    sru = new_row['sru']
                    assert (end - start) % sru == 0
                    assert sru == len(sequence_of_arm) == new_row['arm_length']
                    assert sequence_of_arm * consensus_repeats == fasta_seq == sequence
                    assert sru * consensus_repeats == len(fasta_seq) == length == end - start

                    # non extendable STRs
                    if start > sru:
                        assert sequence_of_arm * (consensus_repeats + 1) != fasta_seq[start-sru:end]
                    if end + sru < fasta_length:
                        assert sequence_of_arm * (consensus_repeats + 1) != fasta_seq[start:end+sru]
                    # if sequence_of_arm * consensus_repeats != fasta_seq:
                    #   logging.error(f"Skipping bad record {row}.")
                    #    continue
                    # elif sru * consensus_repeats != len(fasta_seq):
                    #   logging.error(f"Skipping bad record {row}.")
                    #    continue
                writer.writerow(dict(new_row))
                validated += 1
            # assert validated == total
        if validated == total:
            logging.info(f'Accession {accession} has passed all checks (mode {mode}).')
        else:
            logging.info(f'Accession {accession} has FAILED checks {total},VALIDATED={validated},FAILED={total-validated} (mode {mode}).')
        fout.close()
        self.fnp[mode] = filtered_file
        return self.fnp[mode]

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Mindi: A Non B-DNA extractor tool for data analysis.")
    parser.add_argument("fasta", type=str, help="Path to the fasta file to extract from.")
    parser.add_argument("-p", "--pattern", type=str, default="STR",
                        help="Comma separated list of patterns to extract. "
                             "Available patterns: IR, DR, MR, STR, GQ, Z.")
    parser.add_argument("-d", "--destination", type=str, default="extractions",
                        help="Destination directory to save the extracted files.")
    args = parser.parse_args()
    fasta = Path(args.fasta).resolve()
    if not fasta.is_file():
        raise FileNotFoundError(f"Invalid fasta file {fasta}.")

    mindi = MindiTool(tempdir=args.destination)
    mindi.extract(accession=fasta, 
                  pattern=args.pattern.split(","))
    for m in args.pattern.split(","):
        mindi.sanitize(accession=fasta, mode=m)
