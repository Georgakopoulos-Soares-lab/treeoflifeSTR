import subprocess
import json
from tqdm import tqdm
from utils import load_bucket
import tempfile
import shutil
import os
from pathlib import Path

def agatify(gff_files: list[str], destination: os.PathLike[str]) -> None:
    extract_name = lambda accession: Path(accession).name.split('.gff')[0]
    total_gff_files = len(gff_files)
    destination = Path(destination).resolve()
    destination.mkdir(exist_ok=True)

    # remove potential duplicates
    gff_files = list(set(gff_files))

    # start processing
    for file in tqdm(gff_files, total=total_gff_files, leave=True, position=0):
        file = Path(file).resolve()
        dest_name = extract_name(file) + ".agat.gff"

        destination_file = destination.joinpath(dest_name)
        if destination_file.is_file():
            print(f"Destination file '{destination_file}' already exists!")
            os.remove(destination_file)
            # raise ValueError(f"File {destination_file} already exists!")

        tempdir = tempfile.TemporaryDirectory(prefix=extract_name(file) + ".")

        cur_dir = os.getcwd()
        os.chdir(tempdir.name)
        command = f"agat_convert_sp_gxf2gxf.pl -g {file} -o {dest_name}"
        subprocess.run(
                command,
                shell=True,
                check=True,
                # stdout=subprocess.DEVNULL,
                # stderr=subprocess.DEVNULL
                )

        shutil.move(dest_name, destination)
        os.chdir(cur_dir)
        tempdir.cleanup()


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("schedule_path", type=str, default="schedule_100.json")
    parser.add_argument("--bucket_id", type=int, default=0)
    parser.add_argument("--destination", type=str, default="gff_agatified")

    args = parser.parse_args()
    destination = Path(args.destination).resolve()
    destination.mkdir(exist_ok=True)
    bucket_id = args.bucket_id
    schedule_path = args.schedule_path

    print(f"Retrieving schedule from path {schedule_path}...")
    print(f"Processing bucket {bucket_id}...")
    with open(schedule_path) as f:
        data = json.load(f)
    gff_files = data[str(bucket_id)]
    print(len(gff_files))
    agatify(gff_files, destination=destination)

    print(f"Bucket {bucket_id} has been processed succesfully.")

