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


vardf_all = []

reinfection_samples = [
    "sample_barcode01",
    "sample_barcode02",
    "sample_barcode03",
    "sample_barcode04",
    "sample_barcode05",
    "sample_barcode06",
    "sample_barcode07",
    "sample_barcode08",
    "sample_barcode09",
]
all_samples = reinfection_samples + list(passed)
num_samples = len(all_samples)
base = "data/final/reinfection/output/"
sam_files = {**sam_files, **{x: base + x +'.primertrimmed.rg.sorted.bam' for x in reinfection_samples}}
fasta_files = {**fasta_files, **{x: base + x + '.consensus.fasta' for x in reinfection_samples}}
variant_files = {**variant_files, **{x: base + x + '.merged.vcf' for x in reinfection_samples}}

for sample_id in tqdm(all_samples):
    # Get reads to assess depth
    sam_filename = sam_files[sample_id]
    fasta_filename = fasta_files[sample_id]
    alignments = AlignmentFile(sam_filename).fetch()
    consensus = list(skbio.read(fasta_filename, format="fasta"))[0]

    coverage = np.zeros(29903)
    for alignment in alignments:
        coverage[alignment.positions] += 1

    # Get variants
    variants = VariantFile(variant_files[sample_id])
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
vardf_all.to_csv("data/processed/reinfection/variants.csv")

vardf = vardf_all[vardf_all['score']>3] # Filter out variants that are not high quality
x = [i for i in range(30000)]                                                                                                                                              
y = [0 for _ in x]                                                                                                                                                         
y_ = [sum(vardf['position']==i)/num_samples for i in x]

plt.fill_between(x, y_, y, linewidth=1, color='green')
plt.title("Variant Frequency by Position on Genome")
plt.ylabel("Percent of Samples with Variant")
plt.savefig("data/processed/reinfection_unfiltered.pdf")
plt.show()

var_groups = {position: df for position, df in vardf.groupby('position')}
var_lengths = {k: len(v) for k,v in var_groups.items()}

rare_variants = {k: v for k,v in var_groups.items() if var_lengths[k] == 1}
uncommon_variants = {k: v for k,v in var_groups.items() if var_lengths[k] > 1 and var_lengths[k] < 20}
common_variants = {k: v for k,v in var_groups.items() if var_lengths[k] >= 20}

rare_alleles = pd.DataFrame([
    {"position": k, "reference": v.iloc[0]['reference'], "alternate": v.iloc[0]['alternates'], "count": len(v)}
    for k,v in rare_variants.items()
    ])

uncommon_alleles = pd.DataFrame([
    {"position": k, "reference": v.iloc[0]['reference'], "alternate": v.iloc[0]['alternates'], "count": len(v)}
    for k,v in uncommon_variants.items()
    ])

common_alleles = pd.DataFrame([
    {"position": k, "reference": v.iloc[0]['reference'], "alternate": v.iloc[0]['alternates'], "count": len(v)}
    for k,v in common_variants.items()
    ])

rare_alleles.to_csv("data/processed/rare_alleles_unfiltered.csv")
uncommon_alleles.to_csv("data/processed/uncommon_alleles_unfiltered.csv")
common_alleles.to_csv("data/processed/common_alleles_unfiltered.csv")

coocurrence = pd.crosstab(index=vardf.index, columns=vardf['position'])
coocurrence.index.name="Sample ID"
sns.clustermap(coocurrence, cbar_pos=None, yticklabels=True, xticklabels=True, figsize=(30,20))
# plt.title("Co-Ocurrence of Variant Alleles in COVID Genomes")
plt.xlabel("Variant Position")
plt.ylabel("Sample ID")
plt.savefig("data/processed/reinfection/allele_coocurrence.pdf")
plt.show()

sns.clustermap(coocurrence.loc[reinfection_samples], cbar_pos=None, yticklabels=True, xticklabels=True, figsize=(30,20))
# plt.title("Co-Ocurrence of Variant Alleles in COVID Genomes")
plt.xlabel("Variant Position")
plt.ylabel("Sample ID")
plt.savefig("data/processed/reinfection/allele_coocurrence_reinfection.pdf")
plt.show()

ordered_coocurrence = coocurrence.copy()
ordered_coocurrence.columns = ordered_coocurrence.columns.sort_values()

sns.clustermap(ordered_coocurrence, col_cluster=False, cbar_pos=None, yticklabels=True, xticklabels=True, figsize=(30,20))
# plt.title("Co-Ocurrence of Variant Alleles in COVID Genomes")
plt.xlabel("Variant Position")
plt.ylabel("Sample ID")
plt.savefig("data/processed/reinfection/allele_coocurrence_sorted.pdf")
plt.show()

sns.clustermap(ordered_coocurrence.loc[reinfection_samples], col_cluster=False, cbar_pos=None, yticklabels=True, xticklabels=True, figsize=(30,20))
# plt.title("Co-Ocurrence of Variant Alleles in COVID Genomes")
plt.xlabel("Variant Position")
plt.ylabel("Sample ID")
plt.savefig("data/processed/reinfection/allele_coocurrence_reinfection_sorted.pdf")
plt.show()


links = coocurrence.to_numpy()
num_samples, num_variants = links.shape
sample_self_links = np.zeros((num_samples, num_samples))
variant_self_links = np.zeros((num_variants, num_variants))

biadj = np.block([[sample_self_links, links], [links.T, variant_self_links]]) 
bigraph = nx.from_numpy_matrix(biadj) 
labels = {
    **{i: sample_num for i, sample_num in enumerate(coocurrence.index)},
    **{i+num_samples: variant_num for i, variant_num in enumerate(coocurrence.columns)},
}
colors = { 
    **{i: 'red' for i, sample_num in enumerate(coocurrence.index)}, 
    **{i+num_samples: 'green' for i, variant_num in enumerate(coocurrence.columns)}, 
    }   


fig, ax = plt.subplots(figsize=(15,15))
nx.draw(bigraph, labels=labels, node_color=colors.values(), pos=nx.bipartite_layout(bigraph, [i for i in range(num_samples)]), ax=ax,)
plt.title("Samples vs. Variants Bigraph")
plt.savefig("data/processed/reinfection/variant_bigraph_unfiltered.pdf")
plt.show()

fig, ax = plt.subplots(figsize=(15,15))
nx.draw(bigraph, labels=labels, node_color=colors.values(), pos=nx.kamada_kawai_layout(bigraph), ax=ax)
plt.title("Samples vs. Variants Links")
plt.savefig("data/processed/reinfection/variant_kamada_kawai_graph_unfiltered.pdf")
plt.show()

# Sample linkage similarity
def plot_dendrogram(model, **kwargs): 
    # Create linkage matrix and then plot the dendrogram 
       # create the counts of samples under each node 
    counts = np.zeros(model.children_.shape[0]) 
    n_samples = len(model.labels_) 
    for i, merge in enumerate(model.children_): 
        current_count = 0 
        for child_idx in merge: 
            if child_idx < n_samples: 
                current_count += 1  # leaf node 
            else: 
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count
    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)
    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)
    return linkage_matrix

# Get similarity measurements
jack = list(nx.jaccard_coefficient(bigraph, combinations([i for i in range(num_samples)],2)))
linkage = np.zeros((num_samples, num_samples))
for x,y,v in jack: 
     linkage[x,y] = v
     linkage[y,x] = v
linkgraph = nx.Graph(linkage)
labels = {i:i for i in range(num_samples)}
fig, ax = plt.subplots(figsize=(15,15))
nx.draw(linkgraph, labels=labels, ax=ax)
plt.title("Samples by Jaccard Similarity")
plt.savefig("data/processed/reinfection/sample_similarity_unfiltered.pdf")
plt.show()

fig, ax = plt.subplots(figsize=(15,15))
clustering = AgglomerativeClustering(distance_threshold=0, n_clusters=None).fit(linkage)
plot_dendrogram(clustering, ax=ax)
plt.title("Samples Linkage Clustering")
plt.savefig("data/processed/reinfection/sample_linkage_unfiltered.pdf")
plt.show()

jack = list(nx.jaccard_coefficient(bigraph, combinations([i+num_samples for i in range(num_variants)],2)))
variant_linkage = np.zeros((num_variants, num_variants))
for x,y,v in jack:
    variant_linkage[x-num_samples,y-num_samples] = v
    variant_linkage[y-num_samples,x-num_samples] = v
variant_linkgraph = nx.Graph(variant_linkage)
labels = {i:j for i,j in enumerate(coocurrence.columns)}
fig, ax = plt.subplots(figsize=(15,15))
nx.draw(variant_linkgraph, labels=labels, ax=ax)
plt.title("Variants by Jaccard Similarity")
plt.savefig("data/processed/reinfection/variant_similarity_unfiltered.pdf")
plt.show()

fig, ax = plt.subplots(figsize=(25,13))
clustering = AgglomerativeClustering(distance_threshold=0, n_clusters=None).fit(variant_linkage)
clustering.labels_ = [labels[i] for i in clustering.labels_]
plot_dendrogram(clustering, labels=clustering.labels_, ax=ax)
plt.title("Variant Linkage Clustering")
plt.savefig("data/processed/reinfection/variant_linkage_unfiltered.pdf")
plt.show()

# Plot a Heatmap showing Variant clusters
counts = pd.Series(index=pd.RangeIndex(0,coocurrence.columns.max()+1)).fillna(0.)
counts[coocurrence.columns] = coocurrence.sum()
with Plotter(filename="data/processed/reinfection/reinfection_by_position.pdf", figsize=(12,12)) as ax:
    counts.T.plot(ax=ax)
    ax.set_title("Variant Counts by Position")
    ax.set_xlabel("Position on Genome")
    ax.set_ylabel("Count")


from covid_bronx.quality.gaps import *

primer_binding_sites = "data/external/amplicon_binding_sites.csv"
coverage = {}
gaps = {}
#for sample_id in tqdm(all_samples):
for sample_id in tqdm(reinfection_samples):

    consensus_fasta = consensus_fasta = Sequence.read(fasta_files[sample_id])
    consensus_sam = list(AlignmentFile(sam_files[sample_id]).fetch())
    
    coverage[sample_id], gaps[sample_id] = compute_coverage(consensus_sam, consensus_fasta)

completeness = pd.Series({k: 1-np.mean(v) for k,v in gaps.items()})

l = 29924
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

with Plotter(filename="data/processed/reinfection/reinfection_and_gaps.pdf", figsize=(30, 10)) as ax:

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

with Plotter(filename="data/processed/reinfection/reinfection_and_gaps.pdf", figsize=(30, 10)) as ax:

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

coocurrence.to_csv("data/processed/reinfection/variant_positions.csv")

private_variant_positions = sorted([28917, 28767, 27679, 27389, 23065, 23064, 21245, 11514, 5180, 10448, 8637, 2784, 13665])
private_variants = vardf[vardf['position'].isin(private_variant_positions)]

private_variants[['reference', 'alternates', 'position']].sort_values('position').drop_duplicates()