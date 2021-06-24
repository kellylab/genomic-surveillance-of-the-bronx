import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from covid_bronx.quality import fasta_files, sam_files
from covid_bronx.quality.gaps import *
import plot_helper
from plot_helper import Plotter

from skbio import Sequence
from pysam import AlignmentFile

plot_helper.SHOW = True
    
primer_binding_sites = "data/external/amplicon_binding_sites.csv"
coverage = {}
gaps = {}
for sample_id in tqdm(fasta_files.keys()):
    consensus_fasta = Sequence.read(fasta_files[sample_id])
    consensus_sam = list(AlignmentFile(sam_files[sample_id]).fetch())
    coverage[sample_id], gaps[sample_id] = compute_coverage(consensus_sam, consensus_fasta)
    
# Positions where frameshifts were detected

frameshifts = {
    'AECOM_080': [18546],
    'AECOM_068': [8296],
    'AECOM_034': list(range(20271,20271+58)),
    'AECOM_123': list(range(27970,27970+37)),
    'AECOM_132': list(range(27970,27970+37)),
    'AECOM_070': list(range(27948,27948+20)),
    'AECOM_109': list(range(27766,27766+20)),
    'AECOM_058': [3746,3747],
    'AECOM_059': [15158,15159] + list(range(19749,19749+77)),
    'AECOM_104': list(range(29575,29575+13)),
}

frameshift_coverages = {k: coverage[k][v] for k,v in frameshifts.items()}
windowed_coverage = {k: coverage[k][list(range(v[0]-10,v[0])) + list(range(v[-1],v[-1]+10))].mean() for k,v in frameshifts.items()}

for k,v in frameshifts.items():
    plt.figure(figsize=(30,10))
    plt.plot(coverage[k], linewidth=2)
    fs = coverage[k]*0
    fs[v] = coverage[k][v]
    plt.plot(fs, linewidth=2)
    plt.xlabel("Genomic Position")
    plt.title(f"Coverage by Position for {k} (Frameshifts Highlighted)")
    plt.ylabel("Coverage")
    plt.savefig(f"data/processed/frameshifts/{k}.pdf")
    plt.clf()

# Compare coverage within and outside of frameshift regions
df = pd.DataFrame({
    "Frameshifts":{k: v.mean() for k,v in frameshift_coverages.items()}, 
    "All": {k: coverage[k].mean() for k in frameshift_coverages.keys()},
    "Length": {k: len(v) for k,v in frameshift_coverages.items()},
    "Window Around Frameshift": {k: v for k,v in windowed_coverage.items()},
    })

df.sort_values('Frameshifts').plot.barh(figsize=(16,10))
plt.title("Coverage in Frameshift Regions vs. Entire Genome")
plt.xlabel("Coverage")
plt.ylabel("Sample ID")
plt.savefig("data/processed/frameshifts/coverages.pdf")
plt.clf()
