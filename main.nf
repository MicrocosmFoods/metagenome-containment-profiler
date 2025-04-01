#! /usr/bin/env nextflow

// Description
// Profile metagenomes for containment of reference genomes
// Specifically, we are profiling fermented food metagenomic samples to identify the percent abundance that the sample might be made up by commonly used industrial strains

nextflow.enable.dsl=2

def date = new java.util.Date().format('yyyy-MM-dd')
params.outdir = "${date}-results"
params.threads=16

log.info """\

PROFILE METAGENOMES FOR CONTAINMENT OF REFERENCE GENOMES.
=================================================================
ref_genomes_list                : $params.ref_genomes_list
samples_list                    : $params.samples_list
outdir                          : $params.outdir
threads                         : $params.threads
"""

// Define input channels from TSV files
// Each TSV has two columns: name and accession
input_genomes = Channel
    .fromPath(params.ref_genomes_list)
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.genome_name, row.accession)}  // Creates tuple of [genome_name, accession]

fastq_samples = Channel
    .fromPath(params.samples_list)
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.sample_name, row.accession)}  // Creates tuple of [sample_name, accession]

// workflow steps
workflow {
    // Download reference genomes and samples
    ref_genomes_ch = download_ref_genomes(input_genomes)
    fastq_samples_ch = download_fastq_samples(fastq_samples)
    
    // QC the fastq samples
    qc_samples_ch = qc_fastq_samples(fastq_samples_ch)
    
    // Create sylph sketches
    ref_sketches_ch = sketch_references(ref_genomes_ch)
    
    // Run sylph profile
    profile_results = sylph_profile(
        ref_sketches_ch.collect(),
        qc_samples_ch)
    
    // combine profile results
    combined_profile_results = combine_profile_results(profile_results.collect())

}

process download_ref_genomes {
    tag "${accession}"
    conda "envs/ncbi_env.yml"
    container "quay.io/biocontainers/ncbi-genome-download:0.3.3--pyh7cba7a3_0"
    publishDir "${params.outdir}/reference_genomes", mode: 'copy'
    
    input:
    tuple val(genome_name), val(accession)
    
    output:
    tuple val(accession), path("*.fna")
    
    script:
    // Determine if the accession is from GenBank (GCA) or RefSeq (GCF) for section argument
    def section = accession.startsWith("GCA") ? "genbank" : "refseq"
    
    """
    # Download the genome
    ncbi-genome-download --section ${section} \
                         --assembly-accessions ${accession} \
                         --format fasta \
                         --output-folder temp \
                         --parallel ${params.threads} \
                         bacteria
    
    # Find the downloaded file (handles variable naming in the downloaded file)
    GENOME_FILE=\$(find temp/${section}/bacteria/${accession}/ -name "*.fna.gz")
    
    # Uncompress and rename to a simpler name
    gunzip -c \$GENOME_FILE > ${accession}.fna
    
    # Clean up temporary files
    rm -rf temp
    """
}

process download_fastq_samples {
    tag "${sample_name}"
    // conda "envs/sra_env.yml"
    container "quay.io/biocontainers/sra-tools:3.2.0--h4304569_0"
    input:
    tuple val(sample_name), val(accession)
    
    output:
    tuple val(accession), path("${accession}_*.fastq")
    
    script:
    """
    # Download the fastq files
    fasterq-dump ${accession} --split-files --outfile ${accession}
    """
}

process qc_fastq_samples {
    tag "${accession}"
    conda "envs/fastp.yml"
    container "quay.io/biocontainers/fastp:0.24.0--heae3180_1"
    
    input:
    tuple val(accession), path(reads)
    
    output: 
    tuple val(accession), path("${accession}_trimmed_1.fastq.gz"), path("${accession}_trimmed_2.fastq.gz")
    
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
          -o ${accession}_trimmed_1.fastq.gz \
          -O ${accession}_trimmed_2.fastq.gz \
          --json ${accession}_fastp.json \
          --html ${accession}_fastp.html \
          --thread ${params.threads}
    """
}

process sketch_references {
    tag "${accession}"
    conda "envs/sylph.yml"
    container "quay.io/biocontainers/sylph:0.8.1--ha6fb395_0"
    publishDir "${params.outdir}/reference_sketches", mode: 'copy'
    
    input:
    tuple val(accession_name), path(accession)
    
    output:
    path("${accession_name}.syldb")
    
    script:
    """
    sylph sketch ${accession} -o ${accession_name}
    """
}

process sylph_profile {
    tag "${accession}"
    conda "envs/sylph.yml"
    container "quay.io/biocontainers/sylph:0.8.1--ha6fb395_0"
    publishDir "${params.outdir}/profiles", mode: 'copy'
    
    input:
    path(reference_sketches)
    tuple val(accession), path(reads_1), path(reads_2)

    output:
    path("${accession}_profile.tsv")
    
    script:
    """
    sylph profile ${reference_sketches} -1 ${reads_1} -2 ${reads_2} -o ${accession}_profile.tsv
    """
}

process combine_profile_results {
    tag "combine_profile_results"
    conda "envs/pandas.yml"
    container "quay.io/biocontainers/pandas:1.5.3--pyh8642d48_0"
    publishDir "${params.outdir}/combined_profiles", mode: 'copy'
    
    input:
    path(tsv_files)

    output:
    path("combined_sylph_profiles.tsv")

    script:
    """
    python3 ${baseDir}/bin/process_profile_tsvs.py \
           . \
           combined_sylph_profiles.tsv
    """
}
