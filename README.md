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
$ conda env create -f workflow/env/BLISAcounts.yml
```

```
$ conda activate BLISAcounts
```

Setup project structure, in particular the input data, as indicated below.

In `config/config.yaml`, set sample metadata under `samples`, and antibody(ab)/well 
ground truth metadta under `plate: well_settings`, `plate: ab_settings` in the config file.


Two toy examples are provided to demonstrate the usage of BLISAcounts.

Download example data from the following link:
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
