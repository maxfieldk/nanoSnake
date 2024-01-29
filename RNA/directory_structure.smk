import os
import pandas as pd
from pathlib import Path

paths = [
    "envs",
    "report",
    "qc",
    "scripts",
    "qc/{sample}",
    "rawdata/{sample}",
    "logs/{sample}",
    "intermediates/{sample}/methylation",
    "intermediates/{sample}/dorado",
    "intermediates/{sample}/fastqs",
    "intermediates/{sample}/counts",
    "intermediates/{sample}/counts/genome",
    "intermediates/{sample}/counts/transcriptome",
    "intermediates/{sample}/alignments",
    "intermediates/{sample}/alignments/genome",
    "intermediates/{sample}/alignments/transcriptome",
    "results/{sample}",
    "results/plots",
    "results/tables"
    ]
paths = paths + []
for path in paths:
    for sample in samples:
        os.makedirs(path.format(sample=sample), exist_ok = True)
