# BLISAcounts

A Snakemake pipeline of count enumeration for BLISA sequecing data.

## Description

With `umi_tools`, do:
- obtaining barcode whitelist
- barcode correction and extraction

With `extract_staggered_seq.py`, do:
- barcode extraction of multiplexed sequences with staggered primers

Subsequently, do:
- barcode deduplication and count enumeration

The final output is a count table with all combinations of well barcode - antibody barcode detected.
You may then use any data wrangling method of your choice to visualize the results.

## Usage

Install environment dependencies:
```
$ mamba env create -n BLISAcounts -f workflow/env/curated.yml
```

```
$ conda activate BLISAcounts
```

Setup project structure, in particular the input data, as indicated below.

In `config/config.yaml`, specify `sample` slot with sample metadata (e.g. `config/example.tsv`), 
and `plate` slot with antibody (e.g. `workflow/data/Zhong/HBV_example/AR005/metadata/ab_bc_ground_truth.csv`) and well (e.g. `workflow/data/Zhong/HBV_example/AR005/metadata/well_bc_ground_truth.csv`) ground truth metadta.


To give a concrete example, two toy data are provided to demonstrate the usage of BLISAcounts.

Please download the example data from the following link and put them under the BLISAcounts directory:
https://www.dropbox.com/scl/fo/pjsn5v20vqhsnw4m3pig1/h?rlkey=z51en6cp5w3ibwfzihpjkjtwk&dl=1

### Example 1: Standard sequences

In `config/config.yaml`, set `staggered: False`.

Run:
```
(BLISAcounts) $ snakemake --cores 1 --configfile=config/config.yaml
```

### Example 2: Multiplexed staggered sequences

In config/config.yaml, set `staggered: True`.

Run:
```
(BLISAcounts) $ snakemake --cores 1 --configfile=config/config_staggered.yaml
```

## Project structure

```
├── config
│   ├── config.yaml
│   └── example.tsv
├── LICENSE
├── README.md
└── workflow
    ├── data
    │   └── Zhong (user)
    │       ├── HBV_example (example 1: standard sequences)
    │       │   └── AR005 (library)
    │       │       ├── fastqs
    │       │       │   └── AR005_example.fastq.gz
    │       │       ├── logs
    │       │       ├── metadata
    │       │       └── outs
    │       └── SARS_CoV2_example (example 2: multiplex, staggered sequences)
    │           └── S1 (library)
    │       │       ├── fastqs
    │       │       │   └── S1_example.fastq.gz
    │       │       ├── metadata
    │       │       └── outs
    ├── envs
    │   └── BLISAcounts.yml
    ├── rules
    │   ├── common.smk
    │   ├── umi_tools.smk
    │   └── pipeline.smk
    ├── schemas
    │   ├── config.schema.yaml
    │   └── samples.schema.yaml
    ├── scripts
    │   ├── extract_staggered_seq.py
    │   └── wash_whitelist.py
    └── Snakefile
```
