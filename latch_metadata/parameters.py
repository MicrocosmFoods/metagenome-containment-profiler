
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
    'samples_list': NextflowParameter(
        type=LatchFile,
        default=None,
        section_title=None,
        description='Path to TSV list of sample accessions to download, first column is the name of the sample, second column is the SRA accession.',
    ),
    'outdir': NextflowParameter(
        type=typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})],
        default=None,
        section_title=None,
        description='The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.',
    ),
}

