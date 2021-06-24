from pysam import AlignmentFile
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import skbio
import numpy as np
from scipy.optimize import curve_fit
from plot_helper import Plotter

file_dir = "data/processed/reinfection/"
files = [f for f in os.listdir(file_dir) if f.endswith(".csv") and "mean" not in f]
coverage = []
reference_length = 29903

for filename in files:
    barcode = filename.split(".")[0]
    sample_id = int(barcode[-2:])
    alignments = list(AlignmentFile(f"data/final/reinfection/output/sample_barcode{barcode[-2:]}.trimmed.rg.sorted.bam").fetch())
    consensus = list(skbio.read(f"data/final/reinfection/output/sample_barcode{barcode[-2:]}.consensus.fasta", format="fasta"))[0]
    coverage.append({
        "sample_id": sample_id,
        "mean_aligned_coverage": sum([len(x.positions)/29903 for x in alignments]),
        "mean_coverage": sum([len(x.seq)/29903 for x in alignments]),
        "total_num_reads": len(alignments),
        "completeness": 1-consensus.frequencies()['N']/len(consensus),
    })

coverage = pd.DataFrame(coverage)
coverage.to_csv("data/processed/reinfection/mean_coverage_reinfection.csv", index=False)
coverage.index = coverage['sample_id']
coverage = coverage.sort_index()

def sigmoid(x, L, k, x_0):
    x_0 = 0
    L = 1
    return L/(1+np.exp(-k*(x-x_0)))

xdata = np.array(coverage['mean_aligned_coverage'])
ydata = np.array(coverage['completeness'])
popt, pcov = curve_fit(sigmoid, xdata, ydata, p0=[1, .001, 0])

def f(x):
    return sigmoid(x, *popt) - .94

from scipy.optimize import root_scalar
root = root_scalar(f, bracket=[0,4000]).root
with Plotter(filename="data/processed/reinfection/mean_aligned_coverage_vs_completeness.svg", figsize=(12,12)) as ax:
    coverage.plot.scatter("mean_aligned_coverage", "completeness", ax=ax, label="samples")
    xdata.sort()
    ax.plot(xdata, sigmoid(xdata, *popt), label=f"y={popt[0]:.2f}/(1+e^-({popt[1]:.2f}*(x-{popt[2]:.2f}))")
    ax.axvline(root, linestyle='--', color='green', label=f"Minimum mean coverage for success = {root}")
    ax.axhline(.94, linestyle='--', color='green')
    ax.legend()
    ax.set_title("Mean Coverage of Aligned Reads vs. Completeness of Covid Genomes.")

xdata = np.array(coverage['mean_coverage'])
ydata = np.array(coverage['completeness'])
popt, pcov = curve_fit(sigmoid, xdata, ydata, p0=[1, .001, 0])

def f(x):
    return sigmoid(x, *popt) - .94

root = root_scalar(f, bracket=[100,4000]).root
with Plotter(filename="data/processed/reinfection/mean_unaligned_coverage_vs_completeness.svg", figsize=(12,12)) as ax:
    coverage.plot.scatter("mean_coverage", "completeness", ax=ax, label="samples")
    xdata.sort()
    ax.plot(xdata, sigmoid(xdata, *popt), label=f"y={popt[0]:.2f}/(1+e^-({popt[1]:.2f}*(x-{popt[2]:.2f}))")
    ax.axvline(root, linestyle='--', color='green', label=f"Minimum mean coverage for success = {root}")
    ax.axhline(.94, linestyle='--', color='green')
    ax.legend()
    ax.set_title("Mean Coverage of Unaligned Reads vs. Completeness of Covid Genomes.")
    
# How does number of reads correlate with aligned coverage?
# The reads in the final assembly do not fully align with the reference. Hence,
# the number of aligned reads is less than the total number of reads after
# adjusting for length. How does the distribution of this discrepancy vary?

with Plotter(filename="data/processed/reinfection/aligned_vs_unaligned_reads.svg", figsize=(9,9)) as ax:
    coverage.plot.scatter("mean_coverage","mean_aligned_coverage", c="completeness", ax=ax, colormap='Greens')
    ax.set_title("Aligned vs. Unaligned Coverage for Covid Genomes.")