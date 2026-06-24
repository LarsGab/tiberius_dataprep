# Data‑Handling Pipelines for Training Tiberius

> **Status**: 🚧 *Work in progress* — contributions and bug reports are welcome. If something does not work or documentation is missing, please open an issue.


## Installation

The pipelines need three things on the host: Nextflow, a container runtime
(Singularity or Docker), and gffcompare (only used by `eval_training.nf`).
Everything else — Python, TensorFlow, Tiberius itself, `gffread` — lives inside
the pinned container.

> You do **not** need to `pip install .` to run the pipelines or open
> `plot_acc.ipynb`. A Python install is only required to run the unit tests
> (see *Running the tests* at the bottom of this file).

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

See *Choosing an execution profile* below for how to point the pipeline at your
cluster. If your site is not one of the bundled profiles, you can add your own
in a few lines.

### 2 · Install gffcompare

Get version 0.12.6 from [GitHub](https://github.com/gpertea/gffcompare/releases/tag/v0.12.6).
> **Important:**  Only use version *0.12.6*. Newer releases (e.g. v0.12.9) may report inaccurate results.


```bash
wget https://github.com/gpertea/gffcompare/releases/download/v0.12.6/gffcompare-0.12.6.Linux_x86_64.tar.gz
tar xzf gffcompare-0.12.6.Linux_x86_64.tar.gz
cd gffcompare-0.12.6.Linux_x86_64
export PATH="$(pwd):$PATH"
```

```bash
# Append this line to ~/.bashrc
echo 'export PATH="$PWD:$PATH"' >> ~/.bashrc
```

## Generating TFRecords from local data

### Required Input

Put matching pairs of genome FASTA and annotation files in **two separate directories**:

| File            | Naming pattern                              |
| --------------- | ------------------------------------------- |
| Genome sequence | `<species_id>.genome.fa`                    |
| Gene annotation | `<species_id>.gtf`, `.gff3`, or `.gff`      |

The pipeline picks the annotation file by extension in that order. Non-GTF files
are converted with `gffread -T` before further processing.

`<species_id>` is an arbitrary identifier (typically the Latin binomial—e.g. `arabidopsis_thaliana`) that must be identical for each matching pair. Every genome file must live in *one* directory (`genome_dir`) and every annotation file in *another* (`annot_dir`). 

### Configure the run

Both pipelines read the same `config/config.yaml`. The dataprep pipeline only
uses the top-level keys shown below; the `training:` block is added later for
the evaluation pipeline (see *Exon- and Transcript-level Evaluation* further
down).

```yaml
species_split:
  train: [arabidopsis_thaliana, solanum_lycopersicum, capsicum_annuum]
  val:   [brassica_napus]
  test:  [oryza_sativa]

work_dir:   /scratch/tiberius/work      # pipeline output root
genome_dir: /data/genomes               # directory with *.genome.fa
annot_dir:  /data/annotations           # directory with *.{gtf,gff3,gff}

min_seq_len: 500000 # minimum sequence length (bp) used for training; shorter sequences are filtered out
chunk_size:  9999   # size (bp) of each TFRecord chunk; forwarded to write_tfrecord_species.py as --wsize
```

### Choose a Tiberius version

The pipeline runs inside a pinned Singularity image. Use the version that matches
the weights you plan to train or evaluate. `gffread` is bundled in both images;
`gffcompare` (used only by `eval_training.nf`) stays a host dependency — see
*Installation* above.

| `--tiberius_version` | Image                                       |
| -------------------- | ------------------------------------------- |
| `2.x` *(default)*    | `docker://larsgabriel23/tiberius:2.0.6`     |
| `1.x`                | `docker://larsgabriel23/tiberius:1.1.8`     |

To use a custom image (for example, your own build), override the image directly:
`--container docker://your/image:tag`.

### Choosing an execution profile

The pipeline ships with a few ready-made profiles. Pick one with `-profile`:

| Profile              | Where it runs                                                                |
| -------------------- | ---------------------------------------------------------------------------- |
| `test`               | Local, no container. For `-stub-run`, CI, and tiny end-to-end tests.         |
| `local_docker`       | Your laptop, Docker.                                                         |
| `local_singularity`  | Your workstation, Singularity.                                               |
| `slurm_brain`        | BRAIN cluster (Univ. Greifswald), Singularity, GPU on `vision`.              |
| `slurm_generic`      | Any Slurm cluster, Singularity. No queue or `clusterOptions` assumed.        |

If your site is not one of these, copy `config/slurm_generic.config` to
`config/slurm_<yoursite>.config`, edit the executor/queue/cpus/memory/time to your
site's defaults, and add a matching entry to the `profiles { ... }` block in
`nextflow.config`. Then run with `-profile slurm_<yoursite>`. No `-c` flag is
needed — Nextflow auto-loads `nextflow.config` from the repo root.

### Run the pipeline

```bash
nextflow run tiberius_dataprep/genome2tfrecords.nf \
  -profile local_singularity \
  --config_yaml config/config.yaml \
  --tiberius_version 2.x \
  --resume
```

> **On SLURM?** Submit the same command as a batch job, replacing
> `local_singularity` with `slurm_brain` or `slurm_generic`. The pipeline will
> fan out jobs across the cluster automatically.

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

Add a `training:` block to the same `config/config.yaml` you used for dataprep.
The pipeline derives per-species genome and reference-annotation paths from the
existing `species_split`, `genome_dir`, and `work_dir` fields, so all that's
new is the training-run metadata:

```yaml
# (species_split, work_dir, genome_dir, ... already defined above)

training:
  train1:
    weights_dirs:          # validation: scan each for epoch_* subdirs
      - eudicots/train/train1/0/
      - eudicots/train/train1/1/
    test_weights:          # testing: run only these checkpoints
      - eudicots/train/train1/2/epoch_01/
    eval_dir: eudicots/eval/train1
    hmm: false             # true if the model includes the HMM head
  train2:
    weights_dirs:
      - eudicots/train/train2/0/
    test_weights: []
    eval_dir: eudicots/eval/train2
    hmm: true
```

Per species in `species_split.{val,test}`, the pipeline expects:

| Resource                | Auto-derived path                                       |
| ----------------------- | ------------------------------------------------------- |
| Genome FASTA            | `<genome_dir>/<species_id>.genome.fa`                   |
| Reference annotation    | `<work_dir>/annot_gtf/<species_id>_longest.gtf`         |

The reference annotation is the dataprep pipeline's `*_longest.gtf` output, so
run dataprep first.

Set `--tiberius_version` to the version that produced the weights (`1.x` for
weights trained with Tiberius 1.1.8, `2.x` for 2.0.6 — the default).

### Run validation

```bash
nextflow run tiberius_dataprep/eval_training.nf \
  -profile slurm_brain \
  --config_yaml config/config.yaml \
  --tiberius_version 2.x \
  --resume
```

### Run testing (on the chosen epoch)

```bash
nextflow run tiberius_dataprep/eval_training.nf \
  -profile slurm_brain \
  --config_yaml config/config.yaml \
  --tiberius_version 2.x \
  --resume --use_test
```

> **SLURM users:** submit either command with `sbatch`. Pick the profile that
> matches your site (`slurm_brain`, `slurm_generic`, or your own — see
> *Choosing an execution profile* above).


### Outputs

For each training run, the run produces the results for each evaluation (prediction as `.gtf` and accuracy as `.stats`) in the specified `EvalDir`:

### Plotting Evaluation Results

See `tiberius_dataprep/plot_acc.ipynb` for examples on how to plot the evaluation results.

## Running the tests

Unit tests cover the helpers in `tiberius_dataprep/util.py` (stats parsing,
best-epoch selection, plotting). Install the package with its test extras and
run `pytest`:

```bash
pip install -e ".[test]"
pytest
```

GitHub Actions runs the same suite on every push and pull request — see
`.github/workflows/ci.yml`.
