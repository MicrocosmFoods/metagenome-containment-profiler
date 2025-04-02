# Profiling metagenomes for containment of reference genomes

This repository contains code, including scripts, notebooks, and a workflow for profiling metagenomes for containment of refernece genomes. Specifically this workflow was used to profile fermented food metagenomes for containment of commonly used commercial or industrial fermenting strains. The workflow downloads sets of reference genomes and metagenomic samples, and uses the Sylph metagenome profiling tool to generate containment profiles.

## Workflow Usage

You will need to provide two input files to run the workflow:

1. A list of reference genomes in a TSV file with `genome_name` and `accession` columns.
2. A list of metagenomic sample accessions in a TSV file with `sample_name` and `accession` columns.

Example of input files can be found in the `inputs` directory. 

To run the workflow, you will need to have Nextflow installed. This worklfow can be run using either conda or docker. Launch the workflow with: 
```
nextflow run main.nf \\
--ref_genomes_list inputs/ref_genomes.tsv \\
--samples_list inputs/samples.tsv \\
--outdir <RESULTS_OUTPUT_DIRECTORY> \\
-profile <conda|docker>
```

The main output of the workflow is a TSV containing the results of `sylph profile` for each sammple and the collection of reference genomes, which is filtered to only include hits with >= 98% ANI compared to the reference genome.

## Analyzing Containment of Industrial Fermenting Strains in Fermented Food Metagenomic Samples