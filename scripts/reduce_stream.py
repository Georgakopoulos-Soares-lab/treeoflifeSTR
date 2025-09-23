if __name__ == "__main__":
    from pathlib import Path
    import argparse
    import gzip
    import pandas as pd
    from termcolor import colored
    from tqdm import tqdm
    # import dask.dataframe as dd

    parser = argparse.ArgumentParser(description="Merge final buckets of extractions.")
    parser.add_argument("--indir", type=str, default="extractions_STR_MR")
    parser.add_argument("--pattern", type=str, default="STR", choices=["STR", "MR", "DR", "IR"])
    parser.add_argument("--tree", type=str, default="assembly_summary_tree.txt.gz", help="Tree file for the assembly summary.")
    args = parser.parse_args()
    indir = Path(args.indir).resolve()
    tree_file = Path(args.tree).resolve()
    pattern = args.pattern
    indir = indir.joinpath("merged", pattern)
    assembly_df = pd.read_table(tree_file)
    if "subkingdom" in assembly_df:
        assembly_df.drop(columns=["subkingdom"], inplace=True)
    print(f"Sourcing files from `{indir}` <---")
    # Load files
    empty_files = [file for file in indir.glob(f"*_{pattern}_empty.tsv.gz") if file.is_file()]
    density_files = [file for file in indir.glob(f"*_{pattern}_density_merged.tsv.gz") if file.is_file()]
    merged_files = [file for file in indir.glob(f"*_{pattern}_merged.tsv.gz") if file.is_file()]
    raw_files = [file for file in indir.glob(f"*_{pattern}_raw.tsv.gz") if file.is_file()]

    # Print summary
    print(colored(f"Total empty files < {len(empty_files)}.", "green"))
    print(colored(f"Total merged files < {len(merged_files)}.", "green"))
    print(colored(f"Total raw files < {len(raw_files)}.", "green"))
    print(colored(f"Total density files < {len(density_files)}.", "green"))

    # Create output files
    outfile_density = indir.joinpath(f"extractions_{pattern}_density_merged.tsv.gz")
    outfile_merged = indir.joinpath(f"extractions_{pattern}_merged.tsv.gz")
    outfile_raw = indir.joinpath(f"extractions_{pattern}_raw.tsv.gz")
    outfile_empty = indir.joinpath(f"extractions_{pattern}_empty.tsv.gz")
    def write_output(outfile, files):
        with gzip.open(outfile, "wt", encoding="UTF-8") as f:
            for i, file in tqdm(enumerate(files, 1), total=len(files)):
                try:
                    df = pd.read_table(file).rename(columns={"accession_id": "#assembly_accession"})
                    # if "domain" not in df.columns:
                    #    df = df.merge(assembly_df,
                    #               how="left",
                    #               on="#assembly_accession")
                except pd.errors.EmptyDataError:
                    print(colored(f"File {file} is empty, skipping.", "yellow"))
                    continue
                df.to_csv(f, sep="\t", index=False, header=i==1)
    write_output(outfile_empty, empty_files)
    write_output(outfile_density, density_files)
    write_output(outfile_merged, merged_files)
    write_output(outfile_raw, raw_files)
