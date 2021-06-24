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
from covid_bronx.geography import gen_points_in_gdf_polys, blank_background_choropleth, get_zip_codes_metadata_geo
import geopandas as gpd
from typing import Tuple

coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
passed = coverage_levels[coverage_levels>=.95].index

num_samples = len(passed)

RADIUS_OF_THE_EARTH = 6371 # kilometers

def latlon(lat:float, lon: float) -> Tuple[float, float]:
    """
    Converts lat / lon into x / y coordinates.
    """
    x = RADIUS_OF_THE_EARTH * np.cos(lat) * np.cos(lon)
    y = RADIUS_OF_THE_EARTH * np.cos(lat) * np.sin(lon)

    return (x, y)

vardf = []

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
    variants = VariantFile(variant_files[sample_id])
    vardf.extend([
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

vardf = pd.DataFrame(vardf).set_index("sample_id")
vardf['score'] = vardf['quality'] / vardf['depth']
vardf.to_csv("data/processed/variants/variants.csv")

vardf = vardf[vardf['score']>3] # Filter out variants that are not high quality
x = [i for i in range(30000)]                                                                                                                                              
y = [0 for _ in x]                                                                                                                                                         
y_ = [sum(vardf['position']==i)/num_samples for i in x]

plt.fill_between(x, y_, y, linewidth=1, color='green')
# plt.title("Variant Frequency by Position on Genome")
plt.ylabel("Percent of Samples with Variant")
plt.savefig("data/processed/variants_unfiltered.pdf")
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
plt.savefig("data/processed/allele_coocurrence_unfiltered.pdf")
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

# Get similarity measurements
jack = list(nx.jaccard_coefficient(bigraph, combinations([i for i in range(num_samples)],2)))
linkage = np.zeros((num_samples, num_samples))
for x,y,v in jack: 
     linkage[x,y] = v
     linkage[y,x] = v
linkgraph = nx.Graph(linkage)
labels = {i:i for i in range(num_samples)}

fig, ax = plt.subplots(figsize=(15,15))
clustering = AgglomerativeClustering(distance_threshold=0, n_clusters=None).fit(linkage)

jack = list(nx.jaccard_coefficient(bigraph, combinations([i+num_samples for i in range(num_variants)],2)))
variant_linkage = np.zeros((num_variants, num_variants))
for x,y,v in jack:
    variant_linkage[x-num_samples,y-num_samples] = v
    variant_linkage[y-num_samples,x-num_samples] = v
variant_linkgraph = nx.Graph(variant_linkage)
labels = {i:j for i,j in enumerate(coocurrence.columns)}

fig, ax = plt.subplots(figsize=(25,13))
clustering = AgglomerativeClustering(distance_threshold=0, n_clusters=None).fit(variant_linkage)
clustering.labels_ = [labels[i] for i in clustering.labels_]

# Plot a Heatmap showing Variant clusters
counts = pd.Series(index=pd.RangeIndex(0,coocurrence.columns.max()+1)).fillna(0.)
counts[coocurrence.columns] = coocurrence.sum()

from covid_bronx.quality.gaps import *

primer_binding_sites = "data/external/amplicon_binding_sites.csv"
coverage = {}
gaps = {}
for sample_id in tqdm(passed):

    try:
        consensus_fasta = fasta_files[sample_id]
        consensus_sam = sam_files[sample_id]
        
        coverage[sample_id], gaps[sample_id] = compute_coverage(consensus_sam, consensus_fasta)
    except:
        print(f"Could not read files for {sample_id}.")

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

zip_df, bronx_zip_codes = get_zip_codes_metadata_geo()

from covid_bronx.metadata import preprocess_metadata

metadata = preprocess_metadata()
coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
passed = coverage_levels[coverage_levels>=.95].index.map(lambda x: x.replace("_","-").replace("-0", "-").replace("-0", "-"))

df = pd.read_csv("data/external/pangolin2.csv")
df.index = df['Sequence name'].apply(lambda x: x.split(" ")[0])
metadata[df.columns] = df

clades = pd.read_csv("data/external/nextclade.csv", sep=";")
clades.index = clades['seqName'].apply(lambda x: x.split(" ")[0])
metadata[clades.columns] = clades

zips = metadata.loc[passed]['zip_code'].to_numpy()
zips = np.array(sorted(zips)[2:])

# Get a listing of coordinates by zip code

gdf = gpd.read_file("data/external/ZIP_CODE_040114/ZIP_CODE_040114.geojson")
# dd = pd.read_csv("data/external/us-zip-code-latitude-and-longitude.csv", sep=";", index_col=0)

# lats = dd.loc[zips]['geopoint'].apply(lambda x: float(str(x).split(",")[0]))
# lons = dd.loc[zips]['geopoint'].apply(lambda x: float(str(x).split(",")[1]))

# latlons = [(la, lo) for la, lo in zip(lats, lons)]
latlons = gpd.GeoDataFrame({"ZIPCODE": gdf['ZIPCODE'], 'geometry': gdf['geometry'].centroid}).set_index("ZIPCODE")
# distances = np.array([[np.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2) for x in latlons] for y in latlons])
from scipy.spatial.distance import pdist, squareform 

jaccard_sample_distances = squareform(pdist(coocurrence,metric='jaccard'))

g = sns.jointplot(
    x=pd.Series(jaccard_sample_distances.flatten(), name="Jaccard Distance of SNP Presence"), 
    y=pd.Series(distances.flatten(), name="Normalized Geographic Distance"),
    kind='reg',
)
g.plot_marginals(sns.rugplot, color="r", height=-.15, clip_on=False)
# plt.title("Genomic Distance vs. Geographic Distance")
plt.savefig("data/processed/geography/variant_vs_position.pdf")

from matplotlib.collections import LineCollection

positions = np.vstack([lons, lats])
with Plotter(filename="data/processed/geography/variant_position_graph.pdf", figsize=(40,32)) as ax:
    ax.scatter(positions[0], positions[1])
    non_zero = (np.abs(np.triu(jaccard_sample_distances, k=1)) > 0.02)
    values = np.abs(jaccard_sample_distances[non_zero])
    start_idx, end_idx = np.where(non_zero)
    # a sequence of (*line0*, *line1*, *line2*), where::
    #            linen = (x0, y0), (x1, y1), ... (xm, ym)
    segments = [[positions[:, start], positions[:, stop]]
                for start, stop in zip(start_idx, end_idx)]
    lc = LineCollection(segments,
                        zorder=0, cmap=plt.cm.RdYlGn,
                        norm=plt.Normalize(0, values.mean()),
                        linewidths=(0.25, .5, .75, 1)
                        )
    lc.set_array(values)
    lc.set_linewidths(15 * values)
    ax.add_collection(lc)

coocurrence.index = coocurrence.index.map(lambda x: x.replace("_","-").replace("-0", "-").replace("-0", "-"))
zipgroups = {z: v.index for z,v in metadata.loc[passed].groupby('zip_code')}
zipvectors = pd.DataFrame({z: coocurrence.loc[zipgroups[z]].mean() for z in zipgroups.keys()}).T
zipdistances = squareform(pdist(zipvectors, metric='jaccard'))
zippositions = dd.loc[zipvectors.index]['geopoint']
ziplon = zippositions.map(lambda x: float(str(x).split(",")[0])).drop([6907, 6513])
ziplat = zippositions.map(lambda x: float(str(x).split(",")[1])).drop([6907, 6513])

positions = np.vstack([ziplon, ziplat])
with Plotter(filename="data/processed/geography/variant_position_graph_zip.pdf", figsize=(40,32)) as ax:
    #ax.scatter(positions[0], positions[1])
    gdf.plot(ax=ax)
    ddd.plot(ax=ax, color='red')
    non_zero = (np.abs(np.triu(zipdistances, k=1)) > 0.02)
    values = np.abs(zipdistances[non_zero])
    start_idx, end_idx = np.where(non_zero)
    # a sequence of (*line0*, *line1*, *line2*), where::
    #            linen = (x0, y0), (x1, y1), ... (xm, ym)
    #segments = [[positions[:, start], positions[:, stop]]
    #            for start, stop in zip(start_idx, end_idx)]
    #lc = LineCollection(segments,
    #                    zorder=0, cmap=plt.cm.RdYlGn,
    #                    norm=plt.Normalize(0, values.mean()),
    #                    linewidths=(0.25, .5, .75, 1)
    #                    )
    for i, row in enumerate(zipdistances):
        for j, value in enumerate(row):
            line = LineString()
    lc.set_array(values)
    lc.set_linewidths(15 * values)
    ax.add_collection(lc)
