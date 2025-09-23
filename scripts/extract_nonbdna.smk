from pathlib import Path
import json
from collections import defaultdict
from minditool import MindiTool
from validate_pipeline import validate_pipeline, load_schedule

configfile: "config_IR.yaml"
schedule = Path(config["schedule"]).resolve()
if not schedule.is_file():
    raise FileNotFoundError(f"Schedule file {schedule} does not exist.")
indir = Path(config["indir"]).resolve()
if not indir.is_dir():
    raise FileNotFoundError(f"Input directory {indir} does not exist or is not a directory.")

densities_dir = indir.joinpath("densities")
pattern = config["pattern"]
partition_col = config["partition_col"]
min_partition = int(config["min_partition"])
max_partition = int(config["max_partition"])
max_spacer = int(config["max_spacer"])
min_arm_length = int(config["min_arm_length"])
min_sequence_length = int(config["min_sequence_length"])
min_consensus_repeats = int(config["min_consensus_repeats"])
assembly_summary = Path(config["assembly_summary"]).resolve()

# Tree Of Life
ranks = ["domain",
         "kingdom",
         "phylum",
         "order",
         "class",
         "family"]
total_buckets = len(load_schedule(schedule))
print(f"Total buckets detected: {total_buckets}")
validate_pipeline(schedule, indir, pattern)

rule all:
    input:
        "%s/merged/%s/extractions_%s_merged.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/extractions_%s_raw.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/extractions_%s_empty.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/extractions_%s_density_merged.tsv.gz" % (indir, pattern, pattern),
        "%s/density/density_data_merged_1_%s.txt" % (indir, pattern),
        expand("%s/density/density_data_merged_1_%s_{rank}.txt" % (indir, pattern),
               rank=ranks)

rule merge_extractions:
    input:
        expand([
            "%s/status/bucket_{bucket_id}_status.completed" % indir,
            "%s/log_debug_nonbdna/mindi_tool_{bucket_id}.log" % indir], 
               bucket_id=range(total_buckets))
    output:
        "%s/merged/%s/bucket_{bucket_id}_%s_merged.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/bucket_{bucket_id}_%s_raw.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/bucket_{bucket_id}_%s_empty.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/bucket_{bucket_id}_%s_density_merged.tsv.gz" % (indir, pattern, pattern)
    shell:
        """
        python stream_and_merge_bucket.py \
            {wildcards.bucket_id} \
            --schedule {schedule} \
            --indir {indir} \
            --indir {indir} \
            --pattern {pattern} \
            --assembly_summary {assembly_summary} \
            --partition_col {partition_col} \
            --min_sequence_length {min_sequence_length} \
            --min_partition {min_partition} \
            --max_partition {max_partition} \
            --max_spacer {max_spacer} \
            --min_arm_length {min_arm_length} \
            --min_consensus_repeats {min_consensus_repeats}
        """

rule reduce_stream:
    input:
        expand("%s/merged/%s/bucket_{bucket_id}_%s_merged.tsv.gz" % (indir, pattern, pattern), bucket_id=range(total_buckets)),
        expand("%s/merged/%s/bucket_{bucket_id}_%s_raw.tsv.gz" % (indir, pattern, pattern), bucket_id=range(total_buckets)),
        expand("%s/merged/%s/bucket_{bucket_id}_%s_empty.tsv.gz" % (indir, pattern, pattern), bucket_id=range(total_buckets)),
        expand("%s/merged/%s/bucket_{bucket_id}_%s_density_merged.tsv.gz" % (indir, pattern, pattern), bucket_id=range(total_buckets))
    output:
        "%s/merged/%s/extractions_%s_merged.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/extractions_%s_raw.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/extractions_%s_empty.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/extractions_%s_density_merged.tsv.gz" % (indir, pattern, pattern),
    shell:
        """
        python reduce_stream.py \
                --indir {indir} \
                --pattern {pattern} \
                --tree {assembly_summary}
        """

rule calculate_density:
    input:
        "%s/merged/%s/extractions_%s_merged.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/extractions_%s_raw.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/extractions_%s_empty.tsv.gz" % (indir, pattern, pattern),
        "%s/merged/%s/extractions_%s_density_merged.tsv.gz" % (indir, pattern, pattern),
    output:
        expand("%s/density/density_data_merged_1_%s_{rank}.txt" % (indir, pattern),
               rank=ranks),
        "%s/density/density_data_merged_1_%s.txt" % (indir, pattern),
    params:
        step_threshold=float(config["step_threshold"]),
        seq_col=config["seq_col"],
    shell:
        """
        python calc_densities.py \
                  --dest {indir} \
                  --source {input[0]} \
                  --empty {input[2]} \
                  --assembly_summary {assembly_summary} \
                  --seq_col {params.seq_col} \
                  --step_threshold {params.step_threshold} \
                  --mode {pattern} \
                  --merged 1 \
                  --calculate_gc 0
            """
