
from dataclasses import dataclass
import typing
import typing_extensions

from flytekit.core.annotation import FlyteAnnotation

from latch.types.metadata import NextflowParameter
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir

# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters

generated_parameters = {
    'ref_genomes_list': NextflowParameter(
        type=LatchFile,
        default=None,
        section_title='Input/output options',
        description='Path to TSV list of reference genomes to download, first column is the name of the genome and second column is the accession.',
    ),
    'ref_genomes_dir': NextflowParameter(
        type=LatchDir,
        default=None,
        section_title=None,
        description='Path to directory of pre-downloaded reference genomes in FASTA format and ending in .fna',
    ),
    'accessions_list': NextflowParameter(
        type=LatchFile,
        default=None,
        section_title=None,
        description='Path to TSV list of SRA fastq accessions to download, first column is the sample name and second column is the accession. Provide this if you did not pre-download FASTQ samples.',
    ),
    'fastq_dir': NextflowParameter(
        type=LatchDir,
        default=None,
        section_title=None,
        description='Path to directory of pre-downloaded FASTQ files, provide this if you already downloaded all sample accessions.',
    ),
    'ani_threshold': NextflowParameter(
        type=int,
        default=95,
        section_title=None,
        description='Adjusted ANI threshold for filtering (default: 95)',
    ),
    'outdir': NextflowParameter(
        type=typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})],
        default=None,
        section_title=None,
        description='The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.',
    ),
}

