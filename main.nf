#! /usr/bin/env nextflow

// Description
// Profile metagenomes for containment of reference genomes
// Specifically, we are profiling fermented food metagenomic samples to identify the percent abundance that the sample might be made up by commonly used industrial strains

nextflow.enable.dsl=2

def date = new java.util.Date().format('yyyy-MM-dd')
params.outdir = "${date}-results"
params.threads=16

// Add the ANI threshold parameter with a default value of 95
params.ani_threshold = 95

log.info """\

PROFILE METAGENOMES FOR CONTAINMENT OF REFERENCE GENOMES.
=================================================================
ref_genomes_list                : $params.ref_genomes_list
ref_genomes_dir                 : $params.ref_genomes_dir
accessions_list                 : $params.accessions_list
fastq_dir                       : $params.fastq_dir
outdir                          : $params.outdir
threads                         : $params.threads
ani_threshold                   : $params.ani_threshold
"""

// input parameter validation
if (!params.accessions_list && !params.fastq_dir) {
    error "Either --accessions_list or --fastq_dir must be specified!"
}
if (params.accessions_list && params.fastq_dir) {
    error "Please specify either --accessions_list or --fastq_dir, not both!"
}

if (!params.ref_genomes_list && !params.ref_genomes_dir) {
    error "Either --ref_genomes_list or --ref_genomes_dir must be specified!"
}
if (params.ref_genomes_list && params.ref_genomes_dir) {
    error "Please specify either --ref_genomes_list or --ref_genomes_dir, not both!"
}

// input channels
// input genomes list
if (params.ref_genomes_list) {
    input_genomes = Channel
        .fromPath(params.ref_genomes_list)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row.containsKey('kingdom')) {
                error "Input file must contain 'kingdom' column"
            }
            if (!['bacteria', 'fungi'].contains(row.kingdom.toLowerCase())) {
                error "kingdom must be either 'bacteria' or 'fungi'"
            }
            return row
        }
} else {
    input_genomes = Channel
        .fromPath("${params.ref_genomes_dir}/*.{fna,fa,fasta}")  // accepts multiple extensions
        .map { file -> 
            def accession = file.baseName  // gets filename without extension
            tuple(accession, file)
        }
}


// fastq samples either from accession list to download or from fastq directory of pre-downloaded samples
if (params.accessions_list) {
    fastq_samples = Channel
        .fromPath(params.accessions_list)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_name, row.accession)}
} else {
    fastq_samples = Channel
        .fromFilePairs("${params.fastq_dir}/*_{1,2}.fastq{,.gz}", checkIfExists: true)
        .map { sample_name, files -> tuple(sample_name, files) }
}

// workflow steps
workflow {
    // parameter check for fastq samples
    if (params.accessions_list){
        // download the samples and QC
        downloaded_samples_ch = download_fastq_samples(fastq_samples)
        qc_samples_ch = qc_fastq_samples(downloaded_samples_ch)
    } else {
        // QC the fastq samples
        qc_samples_ch = qc_fastq_samples(fastq_samples)
    }
    
    if (params.ref_genomes_list) {
        // Download reference genomes
        ref_genomes_ch = download_ref_genomes(input_genomes)
    } else {
        ref_genomes_ch = input_genomes
    }
        
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
    memory "10G"
    errorStrategy 'ignore'
    
    input:
    tuple val(genome_name), val(accession), val(kingdom)
    
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
                         ${kingdom}
    
    # Find the downloaded file (handles variable naming in the downloaded file)
    GENOME_FILE=\$(find temp/${section}/${kingdom}/${accession}/ -name "*.fna.gz")
    
    # Uncompress and rename to a simpler name
    gunzip -c \$GENOME_FILE > ${accession}.fna
    
    # Clean up temporary files
    rm -rf temp
    """
}

process download_fastq_samples {
    tag "${sample_name}_${accession}"
    conda "envs/sra_env.yml"
    container "public.ecr.aws/biocontainers/sra-tools:3.1.1--h4304569_2"
    memory "20G"

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
    tag "${sample_name}"
    conda "envs/fastp.yml"
    container "quay.io/biocontainers/fastp:0.24.0--heae3180_1"
    memory "20G"
    errorStrategy { task.exitStatus in [1,143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries 1
    
    input:
    tuple val(sample_name), path(reads)
    
    output: 
    tuple val(sample_name), path("${sample_name}_trimmed_1.fastq.gz"), path("${sample_name}_trimmed_2.fastq.gz"), optional: true
    
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
          -o ${sample_name}_trimmed_1.fastq.gz \
          -O ${sample_name}_trimmed_2.fastq.gz \
          --json ${sample_name}_fastp.json \
          --html ${sample_name}_fastp.html \
          --thread ${params.threads}
    """
}

process sketch_references {
    tag "${accession}"
    conda "envs/sylph.yml"
    container "quay.io/biocontainers/sylph:0.8.1--ha6fb395_0"
    publishDir "${params.outdir}/reference_sketches", mode: 'copy'
    memory "10G"
    
    input:
    tuple val(accession), path(accession)
    
    output:
    path("${accession}.syldb")
    
    script:
    """
    sylph sketch ${accession} -o ${accession}
    """
}

process sylph_profile {
    tag "${accession}"
    conda "envs/sylph.yml"
    container "quay.io/biocontainers/sylph:0.8.1--ha6fb395_0"
    publishDir "${params.outdir}/profiles", mode: 'copy'
    memory "10G"
    
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
    container "public.ecr.aws/biocontainers/pandas:1.5.1_cv1"
    publishDir "${params.outdir}/combined_profiles", mode: 'copy'
    memory "10G"
    
    input:
    path(tsv_files)

    output:
    path("combined_sylph_profiles.tsv")

    script:
    """
    python3 ${baseDir}/bin/process_profile_tsvs.py \
           . \
           combined_sylph_profiles.tsv \
           --ani_threshold ${params.ani_threshold}
    """
}
