# Profiling metagenomes for containment of reference genomes

This repository contains code, including scripts, notebooks, and a workflow for profiling metagenomes for containment of refernece genomes. Specifically this workflow was used to profile fermented food metagenomes for containment of commonly used commercial or industrial fermenting strains. The workflow downloads sets of reference genomes and metagenomic samples, and uses the Sylph metagenome profiling tool to generate containment profiles.

## Workflow Usage

You will need to provide either lists of reference genome and metagenomic sample accessions that will be downloaded, or directories of pre-downloaded reference genomes and metagenomic samples. You can also mix and match these input options, for example providing a list of reference genomes for downloading and a directory of metagenomic samples.

If you are providing lists of reference genomes and/or metagenomic samples, you will need to provide two input files:

1. A list of reference genomes in a TSV file with `genome_name` and `accession` columns.
2. A list of metagenomic sample accessions in a TSV file with `sample_name` and `accession` columns.

Example of input files can be found in the `inputs` directory. 

If you are providing directories of pre-downloaded reference genomes and/or metagenomic samples, the genomes directory expects a set of FASTA files ending in .fna, .fa, or .fasta. For a directory of metagenomic samples, the workflow expects paired-end FASTQ files with the format `{sample_name}_{1,2}.fastq` or `{sample_name}_{1,2}.fastq.gz`.

To run the workflow, you will need to have Nextflow installed. This worklfow can be run using either conda or docker. Launch the workflow with: 
```
nextflow run main.nf \\
--ref_genomes_list inputs/ref_genomes.tsv \\
--samples_list inputs/samples.tsv \\
--ani_threshold 98 \\
--outdir <RESULTS_OUTPUT_DIRECTORY> \\
-profile <conda|docker>
```

Or for providing directories of pre-downloaded reference genomes and metagenomic samples:
```
nextflow run main.nf \\
--ref_genomes_dir <REFERENCE_GENOMES_DIRECTORY> \\
--samples_dir <METAGENOMIC_SAMPLES_DIRECTORY> \\
--ani_threshold 98 \\
--outdir <RESULTS_OUTPUT_DIRECTORY> \\
-profile <conda|docker>
```

The main output of the workflow is a TSV containing the results of `sylph profile` for each sammple and the collection of reference genomes. By default the `--ani_threshold` is set to 95 (species-level ANI threshold), but you can change this to a different threshold above 95.

## Analyzing Containment of Fermented Food Reference DB MAGs in Select Metagenomic Samples

We used this workflow to profile a set of ~1300 species-resolved genomes from fermented foods against a select set of fermented food metagenomic samples. Documentation for curating the set of reference genomes is available [here](https://github.com/MicrocosmFoods/fermentedfood_mags_curation), and documentation and scripts for curating the associated metadata is available [here](https://github.com/MicrocosmFoods/fermentedfood_mags_curation). The set of MAGs and associated metadata can also be found on Zenodo [TODO].

We selected a set of representative samples by selecting no more than 10 samples per food type. This selection is documented in `scripts/sample-selection.R`. We then profiled all ~1300 species-resolved genomes against these ~350 samples using the workflow described above with a 95% ANI threshold. The results from this run are available in `results/2025-05-28-mags-profiling/`.

The script `scripts/mags-sylph-analysis.R` was used to generate statistics and figures from this profiling run. Figures can be found in `figures`. 