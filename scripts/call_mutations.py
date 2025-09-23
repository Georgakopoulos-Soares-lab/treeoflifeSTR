#!/usr/bin/env python3
import os
from tqdm import tqdm
import argparse
from pathlib import Path
from subprocess import run
from termcolor import colored

def run_cm(command):
    print(colored(f"Running command: {command}", "yellow"))
    result = run(command, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with return code {result.returncode}")
minimap="/home/dollzeta/biolab/pangenomics/minimap2-2.29_x64-linux/minimap2"
if __name__ == "__main__":
    # --- ARGPARSE ---
    parser = argparse.ArgumentParser(description="Align a genome to CHM13v2.0 and call variants with bcftools.")
    parser.add_argument("--indir", type=str, default="genomes", help="Input directory containing genome subdirectories.")
    parser.add_argument("--species_taxid", type=int, default=28450)
    parser.add_argument("--k", type=int, default=19)
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use.")
    parser.add_argument("--cluster", type=str, default="cluster_species.txt")
    parser.add_argument("--secondary", default="no", type=str, help="Include secondary alignments (yes/no).")
    parser.add_argument("--ref", type=str)
    args = parser.parse_args()
    ref = args.ref
    indir = Path(args.indir).resolve()
    accessions = []
    with open(args.cluster) as f:
        for line in f:
            line = line.strip().split("\t")
            if args.species_taxid == int(line[0]):
                organism = line[1]
                accessions = line[2].split(";")
                print(colored(f"Using reference {ref} for taxid {args.species_taxid}.", "green"))
                break

    if not Path("genomes").is_dir():
        with open("input.txt", "w") as f:
            for accession in accessions:
                f.write(f"{accession}\n")
        os.system("datasets download genome accession --inputfile input.txt --filename genomes.zip --include genome")
        if Path("genomes.zip").exists():
            os.system("unzip -o genomes.zip -d genomes")
            os.remove("genomes.zip")
        if ref is None:
            print(colored("No reference provided, using default.", "yellow"))
        else:
            ref = None
    
    extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
    files = {extract_id(file): file for file in Path("genomes").glob("**/*.fna")}
    print(len(files), "genomes found.")

    # --- PATHS ---
    VCF_locations = []
    ref_loc = files.get(ref, ref)
    for accession in tqdm(accessions):
        genome = files.get(accession)
        bam = f"{accession}_aligned_with_{ref}.bam"
        sorted_bam = f"{accession}_aligned_with_{ref}.sorted.bam"
        raw_vcf = f"{accession}_variants.raw.vcf.gz"
        normalized_vcf = f"{accession}_variants.norm.vcf.gz"
        # genome = indir.joinpath(accession, "genome.fasta")
        run_cm(f"{minimap} --version")
        run_cm(f"{minimap} -t {args.threads} -ax asm20 --secondary={args.secondary} --eqx {ref_loc} {genome}"
                f"| samtools sort -@{args.threads} -m4G -o {sorted_bam}")
        run_cm(f"samtools index {sorted_bam}")
        run_cm(f"bcftools mpileup -Ou -f {ref_loc} -a DP,AD,ADF,ADR,SP,DP4,INFO/AD,INFO/ADF,INFO/ADR {sorted_bam} "
                f"| bcftools call --ploidy 1 -mv -Oz -o {raw_vcf}")
        run_cm(f"bcftools norm -m-any -f {ref_loc} --check-ref w -Oz -o {normalized_vcf} {raw_vcf}")
        run_cm(f"tabix -p vcf {normalized_vcf}")
        print(colored(f"Done! Normalized, indexed VCF saved to: {normalized_vcf}", "cyan"))
        VCF_locations.append(normalized_vcf)
        break
