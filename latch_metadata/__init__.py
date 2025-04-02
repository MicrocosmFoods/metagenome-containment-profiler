
from latch.types.metadata import (
    NextflowMetadata,
    LatchAuthor,
    NextflowRuntimeResources
)
from latch.types.directory import LatchDir

from .parameters import generated_parameters

NextflowMetadata(
    display_name='metagenome-containment-profiler',
    author=LatchAuthor(
        name="Elizabeth McDaniel",
    ),
    parameters=generated_parameters,
    runtime_resources=NextflowRuntimeResources(
        cpus=10,
        memory=50,
        storage_gib=500,
    ),
    log_dir=LatchDir("latch:///your_log_dir"),
)
