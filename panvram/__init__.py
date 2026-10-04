"""panvram - a pangenome cohort resident on the GPU as refrel3 v1 archives (edits against one decoded reference);
random windows and coordinate slices decoded on the card by a queue kernel, a CPU decoder for the same.

    c = panvram.Cohort.open(dir, device="cuda")      # reference + assemblies resident
    x = c.sample(n, W, generator=g)                  # uint8 [n, W] (bytes or tokens) on the device
    y = c.fetch(assembly, contig, start, length)     # a slice by coordinates
"""
from . import _C
from .cohort import DATASETS, Cohort

FormatError = _C.FormatError
with_cuda = _C.with_cuda
__version__ = "0.1.0"
__all__ = ["Cohort", "DATASETS", "FormatError", "with_cuda"]
