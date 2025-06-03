# Training protocol for the Training (on Eudicots)

This document describes the workflow used to create the Mammalian and Eudicot weights for Tiberius and serves as a default strategy for training a new clade.


Training consits of 3 steps:
1. Data preparation – convert genomes + annotations into TFRecords.

2. Model training – train the CNN + LSTM, then (optionally) fine-tune with an HMM head.

3. Evaluation – measure accuracy on hold-out validation and test sets.

## Required software


### Tiberius (**refactor_hmm_improve_testing branch**)

**Singularity**:
```
singularity build tiberius.sif docker://larsgabriel23/tiberius:dev
```

**GitHub**:
```
git pull https://github.com/Gaius-Augustus/Tiberius
cd Tiberius
git checkout refactor_hmm_improve_testing
pip install .
```

### Nextflow (optional)

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

### gffcompare

Get the latest release from [GitHub](https://github.com/gpertea/gffcompare/releases/tag/v0.12.9) 

```bash
wget https://github.com/gpertea/gffcompare/releases/download/v0.12.9/gffcompare-0.12.9.Linux_x86_64.tar.gz
tar xzf gffcompare-0.12.9.Linux_x86_64.tar.gz
cd gffcompare-0.12.6.Linux_x86_64
export PATH="$(pwd):$PATH"
```

## 1. Data Preparation 

This protocol assumes a set of species was already selected and genome and annotation files for each species was already prepared. Note that currently low confidence genes can not be filtered out from training and simply removeing them from an annotation can lead to incorrect training. Therefore, an annotation has to be used as whole or completely discarded before training. Additionally, it is assumed that the species were split into training, validation and testing.

It is assumed following data is present:

Matching pairs of genome FASTA and annotation GFF3 files in **two separate directories**:

| File            | Naming pattern                              |
| --------------- | ------------------------------------------- |
| Genome sequence | `<species_id>.genome.fa`                    |
| Gene annotation | `<species_id>.gff3`                         |

`<species_id>` is an arbitrary identifier (typically the Latin binomial—e.g. `arabidopsis_thaliana`) that must be identical for each matching pair. Every genome file must live in *one* directory (`genome_dir`) and every annotation file in *another* (`annot_dir`). 


### Preparing the tfRecords

**The recommended route is the `tiberius_dataprep` Nextflow pipeline**. For a manual run (shown for A. thaliana):


```bash
# convert gff3 to gtf
gffread -T athal.gff3 -o athal.gtf

# ensure correct format
tiberius/reformat_gtf.py --input athal.gtf \
        --out   athal_reformat.gtf \
        --prefix athal_

# get longest isoform per gene
get_longest_isoform.py  athal_reformat.gtf > athal_longest.gtf

# write 100 tfRecords per species
write_tfrecord_species.py \
        --wsize 9999 \
        --gtf   athal_longest.gtf \
        --fasta athal.genome.fa \
        --out   athal
```

Finally, list all training-set species—one per line—in `train_species.txt`.

## 2. Training

#### Inputs  
- `tfrecords/` – directory with TFRecord shards for every training species
- `train_species.txt` – species list
- `config.json`model configuration (default example [here](https://github.com/Gaius-Augustus/Tiberius/blob/main/docs/config.json))


Start a training (CNN+LSTM only):
```bash
tiberius/train.py  \
    --cfg config.json --out train_dir --data tfrecords/ \
    --train_species_file train_species.txt
```

*An epoch checkpoint is saved every 5 000 steps. On a single NVIDIA A100, a batch size of 500 is typical.*

To resume from a saved epoch:
```bash
tiberius/train.py  \
    --cfg config.json --out train_dir_2 --data tfrecords/ \
    --train_species_file train_species.txt --load previous_epoch
```

#### Finetuning with HMM
After training and *evaluating* the  CNN+LSTM Model of Tiberius, choose the epoch that has the highest accuracy on the validation data and restart the training with the HMM layer.

```bash
tiberius/train.py  \
    --cfg config.json --out train_dir_hmm --data tfrecords/ \
    --train_species_file train_species.txt --load_lstm previous_epoch --hmm
```
*Note that you have to reduce the batch size compared to the CNN-LStM training. For one A100 GPU, a batch size of 200 works.*


## 3. Evaluation 
Evaluation can be automated (see `README.md` in the Nextflow repo) or run manually as described below.

### Validation
Generate predictions for each saved epoch 

Running Tiberius with trainings epochs for CNN+LSTM model:

```bash
tiberius.py --genome athal.genome.fa \
        --model_lstm_old epoch_20 \
        --out train1.epoch20.athal.gtf
``` 

Running Tiberius with trainings epochs for CNN+LSTM+HMM model:

```bash
tiberius.py --genome athal.genome.fa \
        --model_old epoch_20 \
        --out train2.epoch20.athal.gtf
``` 


Computing Accuracy:
```bash
grep CDS athal.gtf > athal_cds.gtf
gffcompare \
    -r athal_cds.gtf \
    --strict-match -e 0 -T --no-merge \
    train2.epoch20.athal.gtf \
    -o train2.epoch20.athal.stats
```

### Plotting Evaluation Results

See `tiberius_dataprep/plot_acc.ipynb` for examples on how to plot the evaluation results and choosing a final epoch.