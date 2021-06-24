from pysam import VariantFile, AlignmentFile
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from itertools import combinations
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
import skbio
from tqdm import tqdm
from covid_bronx.quality import fasta_files, sam_files, variant_files
from plot_helper import Plotter

coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
coverage_levels = coverage_levels[~coverage_levels.index.isin(["AECOM_127","AECOM_128", "AECOM_129", "AECOM_130","AECOM_131"])] # Drop Nosocomials
passed = coverage_levels[coverage_levels>=.95].index.intersection(sam_files.keys())

num_samples = len(passed)

for sample_id in tqdm(passed):
    # Get reads to assess depth
    sam_filename = sam_files[sample_id]
    fasta_filename = fasta_files[sample_id]
    alignments = AlignmentFile(sam_filename).fetch()
    consensus = list(skbio.read(fasta_filename, format="fasta"))[0]

    coverage = np.zeros(29903)
    for alignment in alignments:
        coverage[alignment.positions] += 1

    # Get variants

meta = []
variants = []
from covid_bronx.variants import parse_nextclade

a,b = parse_nextclade("gisaid/aecom_nextclade.json")
meta.append(a)
variants.append(b)

da = pd.concat(meta)
db = pd.concat(variants)

db['position'] += 1

x = [i for i in range(30000)]                                                                                                                                              
y = [0 for _ in x]                                                                                                                                                         
y_ = [sum(db['position']==i)/num_samples for i in x]

plt.fill_between(x, y_, y, linewidth=1, color='green')
plt.title("Variant Frequency by Position on Genome")
plt.ylabel("Percent of Samples with Variant")
plt.savefig("data/processed/variants_unfiltered.pdf")
plt.show()

coocurrence = pd.crosstab(index=db['seqName'], columns=db['position'])

counts = np.zeros(29801)
for k,v in coocurrence.sum().iteritems():
    counts[k] += v
counts = pd.Series(counts)

from covid_bronx.quality.gaps import *

primer_binding_sites = "data/external/amplicon_binding_sites.csv"
coverage = {}
gaps = {}
for sample_id in tqdm(passed):

    try:
        consensus_fasta = consensus_fasta = Sequence.read(fasta_files[sample_id])
        consensus_sam = list(AlignmentFile(sam_files[sample_id]).fetch())
        
        coverage[sample_id], gaps[sample_id] = compute_coverage(consensus_sam, consensus_fasta)
    except:
        print(f"Could not read files for {sample_id}.")

completeness = pd.Series({k: 1-np.mean(v) for k,v in gaps.items()})

l = 29933
all_coverage = {k: np.zeros(l) for k in coverage.keys()}
all_gaps = {k: np.zeros(l) for k in gaps.keys()}
all_variants = np.zeros(l)

for k, c in coverage.items():
    all_coverage[k][0:len(c)] = c
all_coverage = pd.DataFrame(all_coverage).T

for k, g in gaps.items():
    all_gaps[k][0:len(g)] = g
all_gaps = pd.DataFrame(all_gaps).T

all_variants[counts.index] = counts.values

with Plotter(filename="data/processed/variants_and_gaps.pdf", figsize=(30, 10)) as ax:

    ax.plot(all_coverage.mean(), label="Mean Read Coverage (+- 2 STEM)", color='blue')
    ste = all_coverage.std()/np.sqrt(len(all_coverage.columns))   
    ax.fill_between(all_coverage.columns, all_coverage.mean()+2*ste, all_coverage.mean()-2*ste, alpha=.8, color='blue')
    ax.set_ylabel("Absolute Coverage")

    ag = all_gaps * all_coverage.max().max()
    ax.plot(ag.mean(), label="Mean Gap Frequency (+- 2 STEM)", color='orange')
    ste = ag.std()/np.sqrt(len(ag.columns))
    ax.fill_between(ag.columns, ag.mean()+2*ste, ag.mean()-2*ste, alpha=.8, color='orange')
    ax.fill_between(ag.columns, ag.mean()+2*ste, 0, alpha=.3, color='orange')
    ax.plot(all_variants* all_coverage.max().max()/all_variants.max(), label="Mean Variant Frequency (+- 2 STEM)", color='green')
    ax.set_title("Distribution of Coverage, Gaps, and Variants")
    ax.set_xlabel("Position on Genome")

    secax = ax.secondary_yaxis('right', functions = (lambda x: x / 800, lambda x: x*800))
    secax.set_ylabel('Percent Gap / Variant Frequency')

    ax.legend(loc='upper right')
    
coocurrence.to_csv("data/processed/variants/variant_positions.csv")