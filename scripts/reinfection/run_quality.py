# This script gets quality metrics for the outputs.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from covid_bronx.quality import fasta_files, sam_files
from covid_bronx.quality.gaps import *

primer_binding_sites = "data/external/amplicon_binding_sites.csv"

sample_ids = {
    "reinfection01": "sample_barcode01",
    "reinfection02": "sample_barcode02",
    "reinfection03": "sample_barcode03",
    "reinfection04": "sample_barcode04",
    "reinfection05": "sample_barcode05",
    "reinfection06": "sample_barcode06",
    "reinfection07": "sample_barcode07",
    "reinfection08": "sample_barcode08",
    "reinfection09": "sample_barcode09",
}

fasta_dict = {k: f"data/final/reinfection/output/{v}.consensus.fasta" for k,v in sample_ids.items()}
sam_dict = {k: f"data/final/reinfection/output/{v}.primertrimmed.rg.sorted.bam" for k,v in sample_ids.items()}

for sample_id, fasta_file in tqdm(fasta_dict.items()):


    consensus_fasta = fasta_file
    consensus_sam = sam_dict[sample_id]
    out = f"data/processed/reinfection/{sample_id}_gaps"
    df = compute_primer_coverages(consensus_sam, consensus_fasta, primer_binding_sites, out)
    df.to_csv(f"data/processed/reinfection/{sample_id}.csv")
