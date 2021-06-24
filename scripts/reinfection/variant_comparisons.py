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
passed = coverage_levels[coverage_levels>=.95].index

reinfection_samples = {
    "data/final/reinfection/output/sample_barcode01" : "AECOM-123",
    "data/final/reinfection/output/sample_barcode02" : "AECOM-124",
    "data/final/reinfection/output/sample_barcode03" : "AECOM-125",
    "data/final/reinfection/output/sample_barcode04" : "AECOM-126",
    "data/final/reinfection/output/sample_barcode05" : "AECOM-127",
    "data/final/reinfection/output/sample_barcode06" : "AECOM-128",
    "data/final/reinfection/output/sample_barcode07" : "AECOM-129",
    "data/final/reinfection/output/sample_barcode08" : "AECOM-130",
    "data/final/reinfection/output/sample_barcode09" : "AECOM-131",
    "data/final/reinfection2b/output/sample_barcode01" : "AECOM-103",
    "data/final/reinfection2b/output/sample_barcode02" : "AECOM-104",
    "data/final/reinfection2b/output/sample_barcode03" : "AECOM-105",
    "data/final/reinfection2b/output/sample_barcode04" : "AECOM-106",
    "data/final/reinfection2b/output/sample_barcode05" : "AECOM-107",
    "data/final/reinfection2b/output/sample_barcode06" : "AECOM-108",
    "data/final/reinfection2b/output/sample_barcode07" : "AECOM-109",
    "data/final/reinfection2b/output/sample_barcode08" : "AECOM-110",
    "data/final/reinfection2b/output/sample_barcode09" : "AECOM-111",
    "data/final/reinfection2b/output/sample_barcode10" : "AECOM-112",
    "data/final/reinfection2b/output/sample_barcode11" : "AECOM-113",
    "data/final/reinfection2b/output/sample_barcode12" : "AECOM-114",
    "data/final/reinfection2b/output/sample_barcode13" : "AECOM-115",
    "data/final/reinfection2b/output/sample_barcode14" : "AECOM-116",
    "data/final/reinfection2b/output/sample_barcode15" : "AECOM-117",
    "data/final/reinfection2b/output/sample_barcode16" : "AECOM-118",
    "data/final/reinfection2b/output/sample_barcode17" : "AECOM-119",
    "data/final/reinfection2b/output/sample_barcode18" : "AECOM-120",
    "data/final/reinfection2b/output/sample_barcode19" : "AECOM-121",
    "data/final/reinfection2b/output/sample_barcode20" : "AECOM-122",
    "data/final/reinfection2b/output/sample_barcode21" : "AECOM-126b",
}

vardf_all = []

all_samples = reinfection_samples
num_samples = len(all_samples)
base = "data/final/reinfection/output/"
# sam_files = {**sam_files, **{x: base + x + '.primertrimmed.rg.sorted.bam' for x in reinfection_samples}}
# fasta_files = {**fasta_files, **{x: base + x + '.consensus.fasta' for x in reinfection_samples}}
# variant_files = {**variant_files, **{x: base + x + '.merged.vcf' for x in reinfection_samples}}

for filename, sample_id in tqdm(all_samples.items()):
    # Get reads to assess depth
    sam_filename = filename + '.primertrimmed.rg.sorted.bam'
    fasta_filename = filename + '.consensus.fasta'
    variant_filename = filename + '.merged.vcf'

    alignments = AlignmentFile(sam_filename).fetch()
    consensus = list(skbio.read(fasta_filename, format="fasta"))[0]

    coverage = np.zeros(29903)
    for alignment in alignments:
        coverage[alignment.positions] += 1

    # Get variants
    variants = VariantFile(variant_filename)
    vardf_all.extend([
        {
            **{
                key: value
                for key, value
                in var.info.items()
                },
            **{
                "sample_id": sample_id,
                "position": var.pos,
                "quality": var.qual,
                "reference": var.ref,
                "alternates": var.alts,
                "depth": coverage[var.pos],
                },
        }
        for var in variants.fetch()
    ])

vardf_all = pd.DataFrame(vardf_all).set_index("sample_id")
vardf_all['score'] = vardf_all['quality'] / vardf_all['depth']

vardf = vardf_all[vardf_all['score']>3] # Filter out variants that are not high quality
x = [i for i in range(30000)]                                                                                                                                              
y = [0 for _ in x]                                                                                                                                                         
y_ = [sum(vardf['position']==i)/num_samples for i in x]

coocurrence = pd.crosstab(index=vardf.index, columns=vardf['position'])
coocoo = pd.crosstab(index=vardf.index, columns=vardf['position'], values=vardf['score'], aggfunc=sum).fillna(0.) # Has scores included

dx = coocurrence.loc[['AECOM-126', 'AECOM-126b']]
da = coocurrence.loc[['AECOM-123', 'AECOM-124']]
db = coocurrence.loc[['AECOM-125', 'AECOM-126']]
x = dx.T[dx.sum() == 1]
a = da.T[da.sum() == 1]
b = db.T[db.sum() == 1]


dxx = coocoo.loc[['AECOM-126', 'AECOM-126b']]
daa = coocoo.loc[['AECOM-123', 'AECOM-124']]
dbb = coocoo.loc[['AECOM-125', 'AECOM-126']] 
xx = dxx.T[dxx.sum() != 0]
aa = daa.T[daa.sum()!=0] 
bb = dbb.T[dbb.sum()!=0]

with Plotter(filename="data/processed/reinfection/sample_126_variants.pdf", figsize=(14,10)) as ax:
    sns.heatmap(x.T, square=True, linewidths=.5, cbar=False, ax=ax)

with Plotter(filename="data/processed/reinfection/sample_123_124_variants.pdf", figsize=(14,10)) as ax:
    sns.heatmap(a.T, square=True, linewidths=.5, cbar=False, ax=ax)

with Plotter(filename="data/processed/reinfection/sample_125_126_variants.pdf", figsize=(14,10)) as ax:
    sns.heatmap(b.T, square=True, linewidths=.5, cbar=False, ax=ax)

with Plotter(filename="data/processed/reinfection/sample_126_variants_numerical.pdf", figsize=(14,14)) as ax:
    sns.heatmap(xx.T, square=True, linewidths=.5, ax=ax)

with Plotter(filename="data/processed/reinfection/sample_123_124_variants_numerical.pdf", figsize=(14,14)) as ax:
    sns.heatmap(aa.T, square=True, linewidths=.5, ax=ax)

with Plotter(filename="data/processed/reinfection/sample_125_126_variants_numerical.pdf", figsize=(14,14)) as ax:
    sns.heatmap(bb.T, square=True, linewidths=.5, ax=ax)

relvar = coocoo.var()/coocoo.mean()