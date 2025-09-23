# Code for "Repertoire of Short Tandem Repeats across the Tree of Life"


## Extractions


Directory structure:

```
.
├── notebooks
└── scripts
└── data
``` 

The notebooks directory contains Jupyter notebooks used for data analysis. The scripts directory contains Python scripts for data processing and extractions of STRs from genomic data.

We will showcase an example of how to extract STRs from a given collection FASTA file using the provided scripts.

The extractions were orchestrated and executed on HPC cluster nodes using SLURM job scheduling.


## Dataset

The data/filtered_assemblies_.txt file, contains the complete list of ~118 thousand organism genomes of RefSeq and GenBank assemblies that were used for the analysis.

For instance,

```
GCF_000002515.2 fungi   1
GCF_000002725.2 protozoa        1
GCF_000002765.6 protozoa        1
GCF_000002985.6 invertebrate    1
GCF_000005825.2 bacteria        1
GCF_000005845.2 bacteria        1
```

Also, the extracted STR densities are provided in the data/extractions_STR/ for the original genomes, and in the data/extractions_STR_shuffled/ for the shuffled genomes.
These tables were used to calculate the enrichment between the shuffled and real genomic sequences.

## Requirements

The following python libraries are required to succesfully run the notebooks:

```
statsmodels
numpy
pandas
biopython
polars
ushuffle
pybedtools
matplotlib
seaborn
bokeh
snakemake
```

## Extraction

We will perform a step-by-step extraction of STRs for a given collection of fasta files.

1. **Prepare the environment**: Ensure you have Python 3.x installed along with the required libraries. You can create a virtual environment and install dependencies using pip. 
Python can be installed using micromamba. We will need to install NCBI datasets tool as well.

```bash
micromamba create -n str_env python=3.8
micromamba activate str_env
pip install biopython pandas numpy
micromamba install -c conda-forge ncbi-datasets-cli
```

2. Run the following script:

```bash
bash download_sample.sh
```

This command will download the fasta files and schedule them into buckets that will be processed in parallel.

By default download_sample.sh will download 20 samples from the list assemblies_sample.txt which will be split into 5 buckets. You can modify the number of samples and buckets by changing the parameters in the script.
For instance, we can run 100 samples in 10 buckets by changing the parameters as follows:

Additionally, we need to clone the repository for non-Bgfa and compile the C files.

We can run:

```
git clone git@github.com:abcsFrederick/non-B_gfa.git
```

to clone the repository and then compile the C files:
```
cd non_B_gfa/
make
```
Note that you need to have gcc and make installed in your system.

After compiling, we have to modify the provided .env file as follows:

```
nonBDNA=<gfa executable path>
```

Now, we are ready to run the extraction script.

```bash
bash download_sample.sh 10
```
to download and split the files into 10 buckets (note that this will redownload the files).

Now, we can run the extraction script:

```bash
SCHEDULE=$1
TOTAL_BUCKETS=5
bash run_main.sh $SCHEDULE STR $TOTAL_BUCKETS
```

After the script finishes, you should be able to see the extractions_STR/ directory which contains the raw and the processed files of perfect STRs.

In the SLURM cluster we need to use the submission script to parallelize the jobs.

## Negative Controls

Negative controls preserving approximately the dinucleotide composition were performed using Ushuffle (Jiang et al., 2008) Python wrapper.

You need to install the Python ushuffle package:

```bash
pip install ushuffle
```

For each fasta file, 5Gb of shuffled dinucleotide preserving sequences were generated and concatenated into a single FASTA file.

For instance, if we want to shuffle the sequences in bucket 0, we can run:

```bash
BUCKET_ID=0
bash shuffle.sh <schedule_name>.json shuffled_controls/ 2 $BUCKET_ID
```

Again, because of large number of assemblies, we need to use SLURM submission script to parallelize the jobs. Then, we will repeat the extraction process on the shuffle genomes, as was previously described.


## GFF Sanitization

The GFF files were sanitized using AGAT. AGAT can be installed using conda or micromamba.

```bash
micromamba install -c conda-forge agat
```

The python wrappers are provided in the scripts/ directory below:

```
python agatify_gff.py
```

or 
```
bash agatify_gff.sh
```

Note that the scripts have been developed with the SLURM architecture in mind and require a schedule to be run. You can modify the scripts to run them in your local machine.

## Mutation Calling

The mutation calling across bacterial species were performed using minimap2, samtools, and bcftools.

For each bacterial reference genome we downloaded the strains corresponding to the same species_taxid, as provided from the assembly_summary.txt of RefSeq and GenBank databases.

We used the data/representative.txt (provided in the repository) to download the representative genomes of each species. The representative for each set of strains was chosen based on the NCBI classification. 
NCBI classifies the genomes into representative, reference, and other. We chose the representative genomes for each species. If a representative genome was not available, we chose the reference genome. If neither a representative nor a reference genome was available, we chose the first genome in the list of strains.

Subsequently, we used the script below to align each strain for each reference genome:

```
bash scripts/call_mutations.py --indir $INDIR --species_taxid $SPECIES_TAXID
```

## Notebooks

The notebooks directory contains Jupyter notebooks used for data analysis and visualization. You can open and run these notebooks using Jupyter Notebook or VSCode or another IDE that supports Jupyter notebooks.

Structure of the notebooks directory:

- The Simulation_STR.ipynb contains the statistical analysis of STR distributions within the negative controls (shuffled genomes) and the real genomes.
- The STR_dens.ipynb contains the analysis of STR densities across different taxa and the correlation with genome size, gc content, viral host distribution.
- The primates_STR.ipynb contains the analysis of STRs in primates and the correlation with generation time.
- The TSS_TES_enrichment.ipynb notebook contains the analysis of STR enrichment around TSS and TES regions across Eukaryotes, Bacteria, Archaea, and Viruses.
- The Human-Satellites-STR.ipynb notebook contains the analysis of STRs in human satellites and their comparison with known STR disease loci.
