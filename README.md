# Data‑Handling Pipelines for Training Tiberius

> **Status**: 🚧 *Work in progress* — contributions and bug reports are welcome. If something does not work or documentation is missing, please open an issue.


## Installation

The Pipelines require Nextflow, Singularity and gffcompare

### 1 · Install Nextflow

See Nextflow [Documentation](https://www.nextflow.io/docs/latest/install.html) for more information

```bash
# Download the launcher
curl -s https://get.nextflow.io | bash
chmod +x nextflow
# Put it somewhere on your $PATH
export PATH="$(pwd):$PATH"
# Test
nextflow -v
```

```bash
# Append this line to ~/.bashrc
echo 'export PATH="$PWD:$PATH"' >> ~/.bashrc
```

**Adapt the `config/nextflow.config` to the specifics of your cluster.** Currently, the configuration is set for the BRAIN cluster at University Greifswald

### 2 · Install gffcompare

Get the latest release from [GitHub](https://github.com/gpertea/gffcompare/releases/tag/v0.12.9) 

```bash
wget https://github.com/gpertea/gffcompare/releases/download/v0.12.9/gffcompare-0.12.9.Linux_x86_64.tar.gz
tar xzf gffcompare-0.12.9.Linux_x86_64.tar.gz
cd gffcompare-0.12.6.Linux_x86_64
export PATH="$(pwd):$PATH"
```

```bash
# Append this line to ~/.bashrc
echo 'export PATH="$PWD:$PATH"' >> ~/.bashrc
```


## Generating TfRecords from local data

Requirements: 
A set of matching genomic sequences (FASTA) files and annotations (gff3). 
Each match has to have an ID (e.g. latin species name) and the genome files have to be named <speciesID>.genome.fa and the annotations have to be named <speciesID>.gff3. All genome files have to be in one directory and all annnotation files have to be in one directory.

Preparation:
Configure the configuration file for the data at `config/config_dataprep.yaml`, it is filled in for the eudicot training: 
    - List all species IDs and split them into Training/Validation/Testing.
    - Change work_dir, annot_dir, genome_dir to appropriate location of your local files, output directory


Run:
```
nextflow run  --resume -c config/nextflow.config tiberius_dataprep/genome2tfrecords.nf --configFile config/config_dataprep.yaml
```

If you are using Slurm, you can also run above command as a simple job, it will spawn jobs on the required partitions.

Result:
At the specified working directory, two directories are the result:
- ``annot_gtf`` with the main result being ``<speciesID>_longest.gtf``, which will be used for training and evaluation
- ``tfrecords`` with the main result being ``train/<speciesID>_<number>.tfrecords`` and list of training species ``train/species.txt`` 


## Generating TFRecords from local data

### Input layout

Put matching pairs of genome FASTA and annotation GFF3 files in **two separate directories**:

| File            | Naming pattern                              |
| --------------- | ------------------------------------------- |
| Genome sequence | `<species_id>.genome.fa`                    |
| Gene annotation | `<species_id>.gff3`                         |

`<species_id>` is an arbitrary identifier (typically the Latin binomial—e.g. `arabidopsis_thaliana`) that must be identical for each matching pair. Every genome file must live in *one* directory (`genome_dir`) and every annotation file in *another* (`annot_dir`). 

### Configure the run

Edit `config/config_dataprep.yaml` to point the pipeline to your local data and to define the train/validation/test splits:

```yaml
species_split:
  train: [arabidopsis_thaliana, solanum_lycopersicum, capsicum_annuum]
  val: [brassica_napus]
  test:  [oryza_sativa]

work_dir:   /scratch/tiberius/work      # outputs go here
genome_dir: /data/genomes               # directory with *.genome.fa
annot_dir:  /data/annotations           # directory with *.gff3
```

### Launch

```bash
nextflow run tiberius_dataprep/genome2tfrecords.nf \
  -c config/nextflow.config \
  --configFile config/config_dataprep.yaml \
  --resume
```

> **On SLURM?** Submit above command as a batch job; the pipeline will distribute work across partitions automatically.

### Outputs

At `output_dir` the workflow produces two sub‑directories:

| Path         | Contents                                               | Use‑case                                                        |
| ------------ | ------------------------------------------------------ | --------------------------------------------------------------- |
| `annot_gtf/` | `<species_id>_longest.gtf`                             | Cleaned transcript models used during training/benchmarking |
| `tfrecords/` | `train/<species_id>_*.tfrecords` & `train/species.txt` | Sharded TFRecord files for training Tiberius

