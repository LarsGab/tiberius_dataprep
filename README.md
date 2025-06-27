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
cd gffcompare-0.12.9.Linux_x86_64
export PATH="$(pwd):$PATH"
```

```bash
# Append this line to ~/.bashrc
echo 'export PATH="$PWD:$PATH"' >> ~/.bashrc
```

### 3 · Install pyyaml
```bash
conda install anaconda::pyyaml
```

## Generating TFRecords from local data

### Required Input

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
min_seq_len: 500000 # minimum sequence length of genome FASTA files used for training
```

### Run the pipeline

```bash
nextflow run tiberius_dataprep/genome2tfrecords.nf \
  -c config/nextflow.config \
  --configYAML config/config_dataprep.yaml \
  --resume
```

> **On SLURM?** Submit above command as a batch job; the pipeline will distribute work across partitions automatically.

### Outputs

At `output_dir` the workflow produces two sub‑directories:

| Path         | Contents                                               | Use‑case                                                        |
| ------------ | ------------------------------------------------------ | --------------------------------------------------------------- |
| `annot_gtf/` | `<species_id>_longest.gtf`                             | Cleaned transcript models used during training/benchmarking |
| `tfrecords/` | `train/<species_id>_*.tfrecords` & `train/species.txt` | Sharded TFRecord files for training Tiberius


## Exon‑ and Transcript‑level Evaluation of Training Runs

The pipeline provides two evaluation phases:

1. **Per‑epoch validation** – run after every epoch on the **validation** species. These can be used for hyperparameter tuning and for choosing the final trainings weights.
2. **Final testing** – once you have a final training weights, evaluate its weights on the **test** species. 

### Required input

Supply a YAML file (example for eudicots at `config/config_train_eval.yaml`) with three top‑level sections:

| Key            | Purpose                                                             |
| -------------- | ------------------------------------------------------------------- |
| **Val\_Data**  | Map of *validation* species IDs → genome FASTA + reference GTF      |
| **Test\_Data** | Map of *test* species IDs → genome FASTA + reference GTF (optional) |
| **Training**   | One or more training runs to evaluate                               |

Example:

```yaml
Val_Data:
  arabidopsis_thaliana:
    Genome: arabidopsis_thaliana.genome.fa
    RefAnnot: arabidopsis_thaliana_longest.gtf

Test_Data:
  solanum_lycopersicum:
    Genome: solanum_lycopersicum.genome.fa
    RefAnnot: solanum_lycopersicum_longest.gtf

Training:
  train1:
    WeightsDir:
      0: train1/0/
      1: train1/1/ # restarted run 0, because of the time limit on our cluster
    HMM: False
    TestWeights:
      2_1: eudicots/train/train1/2/epoch_01/ # only used for testing
    EvalDir:
      eudicots/eval/train1 # output directorz
  train2:
    WeightsDir:
      0: eudicots/train/train2/0/
    HMM:
      True
    EvalDir:
      eudicots/eval/train2
    TestWeights:
      Null
```

### Run validation

```bash
nextflow run tiberius_dataprep/eval_training.nf \
  -c config/nextflow.config \
  --configFile config/config_train_eval.yaml \
  --resume
```

### Run testing (on the chosen epoch)

```bash
nextflow run tiberius_dataprep/eval_training.nf \
  -c config/nextflow.config \
  --configFile config/config_train_eval.yaml \
  --resume --use_test
```

> **SLURM users:** submit either command with `sbatch`; the workflow will fan out jobs across partitions automatically.


### Outputs

For each training run, the run produces the results for each evaluation (prediction as `.gtf` and accuracy as `.stats`) in the specified `EvalDir`:

### Plotting Evaluation Results

See `tiberius_dataprep/plot_acc.ipynb` for examples on how to plot the evaluation results.
