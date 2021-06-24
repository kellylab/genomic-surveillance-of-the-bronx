# This script gets quality metrics for the outputs.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from covid_bronx.quality import fasta_files, sam_files
from covid_bronx.quality.gaps import *

primer_binding_sites = "data/external/amplicon_binding_sites.csv"

for sample_id in tqdm(fasta_files.keys()):

    try:
        consensus_fasta = fasta_files[sample_id]
        consensus_sam = sam_files[sample_id]
        out = f"data/processed/quality/{sample_id}_gaps"
        df = compute_primer_coverages(consensus_sam, consensus_fasta, primer_binding_sites, out)
        df.to_csv(f"data/processed/quality/{sample_id}.csv")
    except:
        print(f"Could not read files for {sample_id}.")

