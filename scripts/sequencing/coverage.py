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
    plt.clf()
    
completeness = pd.Series({k: 1-np.mean(v) for k,v in gaps.items()})

with Plotter(filename=f"data/processed/sequencing/coverage.pdf", figsize=(12,12)) as ax:
    
    completeness.plot.hist(ax=ax, bins=40, alpha=.8)
    ax.axvline(x=.95, color='green', ls='--')
    num_complete = (completeness >= .95).sum()
    props = {'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': .5}
    ax.text(0.64, .95, f"{num_complete}/{len(completeness)} Genomes with\n>=95% Completeness", transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)
    ax.set_title("Histogram of Completeness Values")
    ax.set_xlabel("Completeness %")
    ax.set_ylabel("Count")



completeness = completeness.sort_values()
completeness.to_csv('data/processed/sequencing/coverage.csv')