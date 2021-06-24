
from skbio import Sequence
from typing import List
import os
import numpy as np
from itertools import count
import click
from pysam import AlignmentFile
import pandas as pd
from typing import Tuple
import matplotlib.pyplot as plt


def match_gaps(consensus: Sequence, reference: Sequence, title:str) -> List[Sequence]:
    
    consensus_gaps = parse_segments(np.where(consensus.values==b'N')[0])
    response = [
        reference[segment]
        for segment in consensus_gaps
    ]
    for sequence, i, gap in zip(response, count(), consensus_gaps):
        start = gap[0]
        end = gap[-1]
        sequence.metadata['description'] = f"Gap {i} for {title}, spanning nucleotides: {start}-{end}"

    return response

def parse_segments(sequence: np.array) -> List[np.array]:
    """
    Parses an array of integers into disjoint lists.
    """
    segments = []
    current = -1
    current_segment = []
    for i in sequence:
        if current >=0 and i > current+1: # New segment
            segments.append(np.array(current_segment))
            current_segment = [i]
        else:
            current_segment.append(i)
        
        current = i

    segments.append(current_segment)

    return segments

def dropna(l):
    return [x for x in l if x]

def overlap(interval1: Tuple, interval2: Tuple) -> Tuple:

    intersection = (max(interval1[0],interval2[0]), min(interval1[1],interval2[1]))
    if intersection[1] >= intersection[0]:
        return intersection

def compute_coverage(consensus_sam, consensus_fasta):

    consensus_gaps = [(x[0], x[-1]) for x in parse_segments(np.where(consensus_fasta.values==b'N')[0])]
    positions = [np.array(x.positions) for x in consensus_sam]
    coverage = np.zeros(len(consensus_fasta)+20)
    gap_sites = np.zeros(len(consensus_fasta)+20)
    for pos in positions:
        coverage[pos] += 1
    for g in consensus_gaps:
        gap_sites[g[0]:g[1]] += 1

    return coverage, gap_sites

def compute_primer_coverages(consensus_sam, consensus_fasta, primer_binding_sites, out):

    consensus_fasta = Sequence.read(consensus_fasta)
    consensus_gaps = [(x[0], x[-1]) for x in parse_segments(np.where(consensus_fasta.values==b'N')[0])]
    consensus_sam = list(AlignmentFile(consensus_sam).fetch())    
    coverage, gap_sites = compute_coverage(consensus_sam, consensus_fasta)
    
    binding_sites = pd.read_csv(primer_binding_sites)
    sites = [(x[1]['Start'], x[1]['End']) for x in binding_sites[['Start','End']].iterrows()]
    primer_names = binding_sites['amplicon'].tolist()

    # Compute the regions of overlap between gaps for each primer.
    overlapping_regions = {
    name: dropna([overlap(gap, site) for gap in consensus_gaps])
    for name, site
    in zip(primer_names, sites)
    }
    # Compute the total coverage of a primer.
    
    gap_coverage = []
    for name, regions in overlapping_regions.items():
        if len(regions):
            mean = []
            for region in regions:
                mean.append(np.sum(coverage[region[0]:region[1]] / (region[1]-region[0])))
            gap_coverage.append({
                "name": name,
                "gap_coverage": np.mean(mean),
            })

    df = pd.DataFrame(gap_coverage).dropna().sort_values("gap_coverage")
    # eo = coverage / coverage.mean()
    eo = coverage
    go = max(eo) * gap_sites / max(gap_sites)
    sequence_coverage = (1 -  sum([x == 'N' for x in str(consensus_fasta)]) / len(consensus_fasta))*100
    fig, ax = plt.subplots(figsize=(20,10))
    ax.set_title("coverage Vs. Gaps for {0} (Coverage={1:.2f}%".format(out.split("/")[-1], sequence_coverage))
    ax.set_xlabel("Position on Reference Genome")
    ax.set_ylabel("Coverage")
    ax.plot(eo, color="blue", label="Reads")
    ax.fill_between([i for i in range(len(go))], go, color="red", label="Gaps", alpha=.5)
    ax.legend()
    plt.savefig(out+".svg")

    return df

