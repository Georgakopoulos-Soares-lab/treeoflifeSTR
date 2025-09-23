from ushuffle import shuffle, Shuffler
from Bio.SeqIO.FastaIO import SimpleFastaParser
import json
from pathlib import Path
import gzip

extract_name = lambda accession: Path(accession).name.split('.fna')[0]
extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])

def parse_fasta(accession):
    accession = Path(accession).resolve()
    if accession.name.endswith(".gz"):
        file = gzip.open(accession, 'rt')
    else:
        file = open(accession, mode='r', encoding='UTF-8')
    for record in SimpleFastaParser(file):
        yield record[0].split(" ")[0], record[1].lower()
    file.close()

def shuffle_genome(fasta, outdir, level = 2):
    name = extract_name(fasta)
    MAX = 10_000_000
    shuffled_genome = open(f"{outdir}/{name}_level_{level}.shuffled.fna", mode="wb")
    for seqID, seq in parse_fasta(fasta):
        shuffled_seq = b""
        total_len = len(seq)
        for chunk_idx in range(0, total_len, MAX):
            chunk = seq[chunk_idx: chunk_idx + MAX]
            chunk_encoded = chunk.encode("utf-8")
            shuffler = Shuffler(chunk_encoded, level)
            shuffled = shuffler.shuffle()
            shuffled_seq += shuffled
        assert len(shuffled_seq) == total_len
         # Wrap lines at 60 chars
        wrapped = b"\n".join([shuffled_seq[i:i+60] for i in range(0, len(shuffled_seq), 60)])
        shuffled_genome.write(b">%s\n%s\n" % (seqID.encode("utf-8"), wrapped))
        # shuffled_genome.write(b">%s\n%s\n" % (seqID.encode("utf-8"), shuffled_seq))
    shuffled_genome.close()

def load_bucket(bucket_id, schedule_path):
    with open(schedule_path, mode="r", encoding="UTF-8") as f:
        return json.load(f)[str(bucket_id)]

def process_bucket(outdir, bucket_id, schedule, level = 2):
    bucket = load_bucket(bucket_id, schedule)
    outdir = Path(outdir).resolve()
    design_dest = open(f"design_bucket_{bucket_id}.csv", mode="w")
    for fasta in bucket:
        shuffle_genome(fasta, outdir=outdir, level=level)
        name = extract_name(fasta)
        accession_id = extract_id(fasta)
        target = Path(f"{outdir}/{name}_level_{level}.shuffled.fna")
        print(target)
        assert target.is_file()
        design_dest.write(f"{accession_id},{fasta},{target},{bucket_id}\n")
    design_dest.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""Genome shuffler for creation of negative controls.""")
    parser.add_argument("--schedule", type=str, default="schedule_shuffling.json")
    parser.add_argument("--bucket_id", type=int, default=0)
    parser.add_argument("--outdir", type=str, default="scratch/nmc6088/zimin_shuffling")
    parser.add_argument("--level", type=int, default=2)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True, parents=True)
    bucket_id = args.bucket_id 
    schedule = args.schedule
    level = args.level
    process_bucket(outdir=outdir, 
                   bucket_id=bucket_id, 
                   schedule=schedule, 
                   level=level)
