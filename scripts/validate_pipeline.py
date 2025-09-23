from pathlib import Path
from collections import defaultdict
from tqdm import tqdm
from termcolor import colored
import json
from minditool import MindiTool

def load_schedule(schedule) -> dict[str, list[str]]:
    with open(schedule, mode='r') as f:
        return json.load(f)

# Validate pipeline
def validate_pipeline(schedule_name, indir, pattern) -> None:
  schedule = load_schedule(schedule_name)
  # total_buckets = len(schedule)
  total_files_missing = 0
  missing_per_bucket = defaultdict(int)
  missing_files = []
  for bucket_id, files in tqdm(schedule.items(), total=len(schedule)):
      for file in files:
          file_name = MindiTool.extract_name(file)
          extracted_file = indir / Path(file_name + f"_{pattern}.processed.tsv")
          if not extracted_file.is_file():
              missing_files.append(file)
              total_files_missing += 1
              missing_per_bucket[bucket_id] += 1
  if total_files_missing > 0:
      print(colored(f"Missing {total_files_missing} files in the schedule: `{schedule_name}`.", "red"))
      with open(indir / "log_debug_nonbdna/missing_files.log", mode="w", encoding="UTF-8") as f:
          for file in missing_files:
              f.write(f"{file}\n")
      with open(indir / "log_debug_nonbdna/missing_files_per_bucket.log", mode="w", encoding="UTF-8") as f:
          for bucket_id, count in missing_per_bucket.items():
            f.write(f"{bucket_id}\t{count}\t{count * 1e2 / len(schedule[bucket_id]): .2f}\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Validate pipeline by checking if all files in the schedule are processed.")
    parser.add_argument("--schedule", type=str, required=True, help="Path to the schedule file.")
    parser.add_argument("--indir", type=str, required=True, help="Input directory containing the processed files.")
    parser.add_argument("--pattern", type=str, default="bdna", help="Pattern to match the processed files.")
    args = parser.parse_args()
    validate_pipeline(args.schedule, Path(args.indir).resolve(), args.pattern)
