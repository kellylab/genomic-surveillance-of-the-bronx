""" 
Makes a figure providing an overview of our dataset with a focus on lineages
laid out as follows:

a - variant clustermap + counts histogram
b - variant and lineage associations
c - geographic associations
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from typing import Dict
from tqdm import tqdm
from pysam import VariantFile, AlignmentFile
import skbio
from itertools import count
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

savefile = "figures/figure4.pdf"

# a) Variant clustermap
logger.info("Plotting 5a")

from covid_bronx.quality import fasta_files, sam_files, variant_files

coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
passed = coverage_levels[coverage_levels>=.95].index.intersection(sam_files.keys())

num_samples = len(passed)
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
vardf.index = vardf.index.map(lambda x: x.replace("_0", "_").replace("_0", "_"))

vardf = vardf[vardf['score']>3] # Filter out variants that are not high quality
x = [i for i in range(30000)]                                                                                                                                              
y = [0 for _ in x]                                                                                                                                                         
y_ = [sum(vardf['position']==i)/num_samples for i in x]

var_groups = {position: df for position, df in vardf.groupby('position')}
var_lengths = {k: len(set(v.index)) for k,v in var_groups.items()}

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

from covid_bronx.metadata import preprocess_metadata

metadata = preprocess_metadata()
coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
passed = coverage_levels[coverage_levels>=.95].index.map(lambda x: x.replace("_","-").replace("-0", "-").replace("-0", "-")).intersection(sam_files.keys())

df = pd.read_csv("data/external/Monte_clade_meta.csv")
df = df.loc[df['Sequence name'].dropna().index]
df.index = df['Sequence name'].apply(lambda x: x.split(" ")[0]).apply(lambda x: x.replace("-0","-").replace("-0","-"))
metadata[df.columns] = df

clades = pd.read_csv("data/external/nextclade.csv", sep=";")
clades.index = clades['seqName'].apply(lambda x: x.split(" ")[0])
metadata[clades.columns] = clades

metadata = metadata.loc[passed]

metadata.index = metadata.index.map(lambda x: x.replace("-", "_"))
# metadata.index = metadata.index.map(lambda x: x.replace("0", "") if x[-3] == '0' or x[-4] == 0 else x)

metadata2 = pd.read_csv("data/confidential/sample_metadata.csv")
metadata2.index = metadata2.number.map(lambda x: f"AECOM-{x}")
metadata2 = metadata2.reindex(metadata2.index.tolist() + list(passed.difference(metadata2.index)))
metadata2.index = metadata2.index.map(lambda x: x.replace("-", "_"))
metadata2 = metadata2.loc[metadata2.index.intersection(metadata.index)]

symptoms = set([y.strip() for x in metadata2['COVID sx'] for y in str(x).split(",")])
comorbidities = set([y.strip() for x in metadata2['Pt co-morbidities'] for y in str(x).split(",")])
treatments = set([y.strip() for x in metadata2['Tx'] for y in str(x).split(",")])

symptom_df = pd.DataFrame({symptom: metadata2['COVID sx'].map(lambda x: symptom in str(x)) for symptom in symptoms})
comorbidities_df = pd.DataFrame({como: metadata2['Pt co-morbidities'].map(lambda x: como in str(x)) for como in comorbidities})
treatments_df = pd.DataFrame({tx: metadata2['Tx'].map(lambda x: tx in str(x)) for tx in treatments})

"""

# Skip remaining figures

# b) Select Variant Associations
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

# Test each variant associations using a contingency test

contingencies = {}
for pos, var in tqdm({**uncommon_variants, **common_variants, **rare_variants}.items()):
    data = {}
    data['sex'] = pd.DataFrame([metadata.loc[var.index]['sex'].value_counts(), metadata.loc[metadata.index.difference(var.index)]['sex'].value_counts(),], index=['With', 'Without']).fillna(0.)
    for col in ['age', 'num_side_effects', 'num_comorbidities', 'num_treatments']:
        with_variant = metadata.loc[var.index][col]
        without_variant = metadata.loc[metadata.index.difference(var.index)][col]
        data[col] = {'with': with_variant, 'without': without_variant}
    for col in symptom_df.columns:
        with_variant = symptom_df.loc[var.index][col]
        without_variant = symptom_df.loc[symptom_df.index.difference(var.index)][col]
        d = pd.DataFrame({'with': with_variant.value_counts(), 'without': without_variant.value_counts()}).fillna(0.)
        if True not in d.index:
            d = d.append(pd.Series({'with': None, 'without': None}, name=True))
        if False not in d.index:
            d = d.append(pd.Series({'with': None, 'without': None}, name=False))
        data[col] = pd.DataFrame([d.loc[True], d.loc[False]]) # Order the rows.
    for col in comorbidities_df.columns:
        with_variant = comorbidities_df.loc[var.index][col]
        without_variant = comorbidities_df.loc[comorbidities_df.index.difference(var.index)][col]
        d = pd.DataFrame({'with': with_variant.value_counts(), 'without': without_variant.value_counts()}).fillna(0.)
        if True not in d.index:
            d = d.append(pd.Series({'with': None, 'without': None}, name=True))
        if False not in d.index:
            d = d.append(pd.Series({'with': None, 'without': None}, name=False))
        data[col] = pd.DataFrame([d.loc[True], d.loc[False]]) # Order the rows.
    for col in treatments_df.columns:
        with_variant = treatments_df.loc[var.index][col]
        without_variant = treatments_df.loc[treatments_df.index.difference(var.index)][col]
        d = pd.DataFrame({'with': with_variant.value_counts(), 'without': without_variant.value_counts()}).fillna(0.)
        if True not in d.index:
            d = d.append(pd.Series({'with': None, 'without': None}, name=True))
        if False not in d.index:
            d = d.append(pd.Series({'with': None, 'without': None}, name=False))
        data[col] = pd.DataFrame([d.loc[True], d.loc[False]]) # Order the rows.

    contingencies[pos] = data

from scipy.stats import fisher_exact, chi2_contingency, mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
pvalues = {}
effects = {}
for pos, data in tqdm(contingencies.items()):
    scores = {}
    effect = {}
    oddsratio, pvalue = fisher_exact(data['sex'])
    scores['sex'] = pvalue
    effect['sex'] = oddsratio
    for col in ['age', 'num_side_effects', 'num_comorbidities', 'num_treatments']:
        stat = mannwhitneyu(*data[col].values())
        if min(len(x) for x in data[col].values()) > 5:
            scores[col] = stat.pvalue
            effect[col] = data[col]['with'].mean() / data[col]['without'].mean()
    for col in [*symptom_df.columns] + [*comorbidities_df.columns] + [*treatments_df.columns]:
        try:
            oddsratio, pvalue = fisher_exact(data[col])
        except ValueError:
            pass
        if min(min(data[col].sum(0)), min(data[col].sum(1))) > 4:
            scores[col] = pvalue
            effect[col] = oddsratio
    pvalues[pos] = scores
    effects[pos] = effect

ppe = pd.DataFrame(effects).T.sort_index()
pp = pd.DataFrame(pvalues).T.sort_index()
pps = pp.T.apply(lambda x:pd.Series( fdrcorrection(x.dropna())[1], index=x.dropna().index))
pps = pp.T.loc[((pp<.05).sum() != 0)]
ppee = ppe.T.loc[((pp<.05).sum() != 0)]
ppee = ppee[pps.columns[np.where((pps<.05).sum() != 0)]]
pps = pps[pps.columns[np.where((pps<.05).sum() != 0)]]


# Break it down by category

demographic_columns = {
    'age': 'Age',
    'sex': 'Sex',
    'former smoker': 'Former Smoker',
    'asymptomatic': 'Asymptomatic',
}

symptoms_columns = {
    'chest tightness': 'Chest Tightness',
    'loss of taste/smell': 'Anosmia',
    'postnasal drip': 'Postnasal Drip',
    'AMS': 'Altered Mental Status',
    'mild cough': 'Mild Cough',
}

comorbidities_columns = {
    'mild intellectual disability': 'Mild Intellectual Disability',
    'CAD': 'Coronary Artery Disease',
    'DM2': 'Type 2 Diabetes',
    'dementia': 'Dementia',
    'OA': 'Osteoarthritis',
    'PPM': 'PPM',
    'CKD3': 'Chronic Kidney Disease',
    'schizophrenia': 'Schizophrenia',
}

pps = pd.DataFrame({v: pps.loc[k] for k,v in {**demographic_columns, **symptoms_columns, **comorbidities_columns}.items() if k in pps.index}).T
ppee = pd.DataFrame({v: ppee.loc[k] for k,v in {**demographic_columns, **symptoms_columns, **comorbidities_columns}.items() if k in ppee.index}).T
variable_counts = {}
variable_counts['Age'] = 0 # Fix this later
variable_counts['Sex'] = sum(metadata['sex']=='male')
variable_counts['Asymptomatic'] = symptom_df['asymptomatic'].sum()

for index in ppee.index:
    if index in ['Age', 'Sex', 'Asymptomatic']:
        pass
    else:
        for i,v in symptoms_columns.items():
            if v == index:
                variable_counts[index] = symptom_df[i].sum()
        for i,v in comorbidities_columns.items():
            if v == index:
                variable_counts[index] = comorbidities_df[i].sum()

variable_counts = pd.Series(variable_counts)

log_transformed_effects = (np.log(ppee).replace(np.inf, 10).replace(-np.inf, -10) * (pps<.05).astype(int)).replace(0, np.nan).dropna(axis=1, how='all')

plt.clf()
plt.close()
fig5 = plt.figure(figsize=(24,24))
gs = fig5.add_gridspec(40,40)
"""

# a) Variant clustermap
logger.info("Plotting 5a")
plt.clf()
plt.close()
fig5 = plt.figure(figsize=(24,24))
gs = fig5.add_gridspec(40,40)

ax_a = fig5.add_subplot(gs[2:21,:])

metadata.index = metadata.index.map(lambda x: x.replace("-", "_"))
coocurrence = pd.crosstab(index=vardf.index, columns=vardf['position'])
dates = metadata.loc[coocurrence.index]['collection_date'].sort_values()
coocurrence = coocurrence.loc[dates.index].astype(bool).astype(int)
coocurrence.index.name="Sample ID"

counts = np.zeros(29801)
for k,v in coocurrence.sum().iteritems():
    counts[k] += v
counts = pd.Series(counts)
ax_a1 = fig5.add_subplot(gs[0:2,:])
coocurrence.sum().plot.bar(ax=ax_a1)
ax_a1.set_ylabel("Count")
ax_a1.xaxis.set_ticks([])
ac = coocurrence.copy()
ac.columns = [str(c) + ('------------' if i%2 == 1 else '') for i,c in zip(count(),coocurrence.columns)]
# ac.index = [c.split("_")[1] + ("-----" if i%2 == 1 else '') for i,c in zip(count(),coocurrence.index)]
sns.heatmap(ac, yticklabels=False, xticklabels=True, cbar=False, ax=ax_a)

# Add Secondary axis with y ticks
y_ticks = list(map(lambda x: x.date(), dates.to_list()))
y_0 = y_ticks[0]
reduced_y_ticks = [y_0]
tick = 1/len(y_ticks)
tock = 1
last_tock = tock
reduced_y_tick_positions = [tock]
for y in y_ticks: # Make the ticks show by week
    if (y - y_0).days >= 7:
        reduced_y_ticks.append(y)
        if last_tock - tock >= 2*tick:
            reduced_y_tick_positions.append(tock)
        else: # Add extra spacing if needed
            reduced_y_tick_positions.append(tock - .6* tick)
        y_0 = y
        last_tock = tock
    tock -= tick

reduced_y_ticks.append(y) # Get the last one as well
reduced_y_tick_positions.append(tock)

reduced_y_ticks = reduced_y_ticks
reduced_y_tick_positions = reduced_y_tick_positions

ax_a2 = ax_a.twinx()
ax_a2.set_ylabel("Date")
ax_a2.tick_params(axis='y')
ax_a2.set_yticklabels(reduced_y_ticks)
ax_a2.set_yticks(reduced_y_tick_positions)

from dna_features_viewer import BiopythonTranslator

class CustomTranslator(BiopythonTranslator):
    def compute_feature_fontdict(self, feature):
        return {'size': 12}

ax_a3 = fig5.add_subplot(gs[25:28,:])
graphic_record = CustomTranslator().translate_record("data/external/covid_track.gff")
graphic_record.plot(ax=ax_a3)

plt.savefig(savefile)
plt.show()

grid = sns.clustermap(ac, yticklabels=False, row_cluster = False, figsize=(40,20), cbar_pos=None)
ax = grid.ax_heatmap
ax_a2 = ax.twinx()
ax_a2.set_ylabel("Date")
ax_a2.tick_params(axis='y')
ax_a2.set_yticklabels(reduced_y_ticks)
ax_a2.set_yticks(reduced_y_tick_positions)
ax_a2.yaxis.tick_left()
plt.savefig("data/processed/variants/variants_clustered.pdf")

# b) Associations

"""
cbar_ax = fig5.add_subplot(gs[38:, 14:17])
heatmap_ax = fig5.add_subplot(gs[24:36, 0:14])
y_ax = fig5.add_subplot(gs[24:36, 14:17])
x_ax = fig5.add_subplot(gs[37:, 0:14])
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
# y_ax2 = fig5.add_subplot(gs[0:14, 36:])
y_ax2 = plt.axes([0,0,1,1])
ip = InsetPosition(y_ax, [0, 10/11, 1, 1/11])
y_ax2.set_axes_locator(ip)
sns.heatmap(data=log_transformed_effects, cmap='coolwarm', linewidths=.5, linecolor='grey', xticklabels=False, yticklabels=True, ax=heatmap_ax, cbar_ax=cbar_ax,cbar_kws={'orientation':'horizontal'})
#.set_title("Selected Patient Data and Variant Associations Colored by Effect Size (FDR Corrected).")
pd.Series([var_lengths[k] for k in ppee.columns], index=ppee.columns).plot.bar(ax=x_ax, alpha=.6)
x_ax.invert_yaxis()
x_ax.xaxis.tick_top()
x_ax.xaxis.set_label_position('top')
x_ax.set_ylabel("Count")
x_ax.set_xlabel("Variant")

variable_counts[::-1].plot.barh(ax=y_ax, alpha=.6)
pd.DataFrame({'Age':metadata['age']}, columns=list(reversed(variable_counts.index)))['Age'].plot.box(vert=False, ax=y_ax2)

y_ax2.yaxis.set_ticks([])
y_ax2.xaxis.set_ticks([])
y_ax.yaxis.set_ticks([])
y_ax.set_xlabel("Count")

"""
"""
# c) Geographic Correlations
logger.info("Plotting 5c")

from covid_bronx.geography import gen_points_in_gdf_polys, blank_background_choropleth, get_zip_codes_metadata_geo
import geopandas as gpd
coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
passed = coverage_levels[coverage_levels>=.95].index

num_samples = len(passed)
from covid_bronx.metadata import preprocess_metadata
from matplotlib import colors

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
# Remove extraneous zip codes
gdf = gdf[~(gdf['PO_NAME']=='Staten Island')]

latlons = gpd.GeoDataFrame({"ZIPCODE": gdf['ZIPCODE'], 'geometry': gdf['geometry'].centroid}).set_index("ZIPCODE")

zdf, bzip = get_zip_codes_metadata_geo()
zdf.index = zdf['zip_code']
gdf.index = gdf['ZIPCODE']
gdf[zdf.columns] = zdf

# Compute graph
adjacency = np.random.rand(len(gdf.index), len(gdf.index)) # Sample
adjacency = ((adjacency + adjacency.T) > 1.95).astype(int)
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree

jaccard_sample_distances = squareform(pdist(coocurrence,metric='jaccard'))
ind = coocurrence.index.map(lambda x: x.replace("_", "-"))
jaccard_df = pd.DataFrame(jaccard_sample_distances, index=ind, columns=ind)

sampleid_to_zip = pd.DataFrame(index=gdf.index, columns=metadata.index).fillna(0.)
for z, samples in {k:v.index for k,v in metadata['zip_code'].groupby(metadata['zip_code'])}.items():
    z = str(z)
    if z in sampleid_to_zip.index:
        sampleid_to_zip.loc[z, samples] = 1
sampleid_to_zip = sampleid_to_zip.drop_duplicates()

adjacency = 1/np.zeros((len(gdf.index), len(gdf.index)))
zip_positions = {z: i for i, z in enumerate(gdf.index)}
from itertools import combinations
# Fill in values
for z1, z2 in combinations(sampleid_to_zip.index, 2):
    a = sampleid_to_zip.loc[z1]
    b = sampleid_to_zip.loc[z2]
    dist = jaccard_df.loc[a==1, b==1].to_numpy().flatten().mean()
    if np.isnan(dist):
        dist = np.inf
    i = zip_positions[z1]
    j = zip_positions[z2]
    adjacency[i,j] = 1/dist
    adjacency[j,i] = 1/dist

from shapely.geometry.linestring import LineString

ax_d = fig5.add_subplot(gs[22:, 15:])

gdf.plot(ax=ax_d, zorder=1, alpha=.4)
length = len(gdf)
centroid = gdf.centroid
gg = []
weights = []
mst = minimum_spanning_tree(adjacency).toarray() # Smaller distances are better
for i in range(length):
    for j in range(length):
        weights.append(mst[i,j])
        gg.append( LineString([centroid.iloc[i], centroid.iloc[j]]) )

weights = np.array(weights) * np.random.rand(len(weights)) * 2
gg = gpd.GeoSeries(gg)
gg.plot(linewidth=weights, ax=ax_d, alpha=.9, color='red')
ax_d.set_axis_off()

# fig5.tight_layout(pad=0.5)
plt.savefig(savefile)
plt.show()

# # Also plot the MST
# fig, ax = plt.subplots(figsize=(20,24))
# gdf.plot(ax=ax, zorder=1, alpha=.5)
# gg.plot(linewidth=weights, ax=ax, alpha=.9, color='red')
# plt.savefig("data/processed/geography/variants_choropleth.png")
# plt.show()
# 
# x, y = np.where(mst != 0)  
# z = adjacency[x,y]
# x_zip = [gdf.index[i] for i in x]
# y_zip = [gdf.index[j] for j in y]
# 
# from collections import defaultdict
# mst_centralities = defaultdict(lambda: 0)
# for _ in x_zip:
#     mst_centralities[_] += 1
# 
# mst_centralities = pd.Series(dict(mst_centralities))
# counts = pd.Series({k: len(v) for k,v in metadata.groupby('zip_code')})                                                                      
# counts.index = counts.index.map(lambda x: str(x))
# df = pd.DataFrame([mst_centralities, counts],index=['centrality', 'count']).T.fillna(0.)
# df = df.sort_values('centrality', ascending=False)
# df.plot.bar(figsize=(16,10))
# plt.title("Highest Centrality Zip Codes")
# plt.xlabel("Zip Code")
# plt.ylabel("Centrality / Count")
# plt.savefig("data/processed/geography/zip_code_centralities.png")
# plt.show()

"""