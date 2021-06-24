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
from covid_bronx.variants import parse_nextclade

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

savefile = "figures_final/figure4a.pdf"

# a) Variant clustermap
logger.info("Plotting 5a")

from covid_bronx.quality import fasta_files, sam_files, variant_files

coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
passed = coverage_levels[coverage_levels>=.95].index.intersection(sam_files.keys())

num_samples = len(passed)
vardf = []

# for sample_id in tqdm(passed):
#     # Get reads to assess depth
#     sam_filename = sam_files[sample_id]
#     fasta_filename = fasta_files[sample_id]
#     alignments = AlignmentFile(sam_filename).fetch()
#     consensus = list(skbio.read(fasta_filename, format="fasta"))[0]
# 
#     coverage = np.zeros(29903)
#     for alignment in alignments:
#         coverage[alignment.positions] += 1
# 
#     # Get variants
#     variants = VariantFile(variant_files[sample_id])
#     vardf.extend([
#         {
#             **{
#                 key: value
#                 for key, value
#                 in var.info.items()
#                 },
#             **{
#                 "sample_id": sample_id,
#                 "position": var.pos,
#                 "quality": var.qual,
#                 "reference": var.ref,
#                 "alternates": var.alts,
#                 "depth": coverage[var.pos],
#                 },
#         }
#         for var in variants.fetch()
#     ])
# 
# vardf = pd.DataFrame(vardf).set_index("sample_id")
# vardf['score'] = vardf['quality'] / vardf['depth']
# vardf.index = vardf.index.map(lambda x: x.replace("_0", "_").replace("_0", "_"))
# 
# vardf = vardf[vardf['score']>3] # Filter out variants that are not high quality
# x = [i for i in range(30000)]                                                                                                                                              
# y = [0 for _ in x]                                                                                                                                                         
# y_ = [sum(vardf['position']==i)/num_samples for i in x]
# 
# var_groups = {position: df for position, df in vardf.groupby('position')}
# var_lengths = {k: len(set(v.index)) for k,v in var_groups.items()}


from covid_bronx.metadata import preprocess_metadata

df = pd.read_csv("data/external/Monte_clade_meta.csv")

metadata = preprocess_metadata()
metadata.index = metadata.index.map(lambda x: x.replace("-","_"))
# passed = coverage_levels[coverage_levels>=.95].index.map(lambda x:
# x.replace("_0", "_").replace("_0", "_")).intersection(sam_files.keys())
passed = passed.map(lambda x: x.replace("_0", "_").replace("_0", "_"))


df = df.loc[df['Sequence name'].dropna().index]
df.index = df['Sequence name'].apply(lambda x: x.split(" ")[0]).apply(lambda x: x.replace("-0","-").replace("-0","-").replace("-","_"))
metadata[df.columns] = df

# a) Variant clustermap
logger.info("Plotting 4a")

meta = []
variants = []

for i in tqdm(range(14)):
    a,b = parse_nextclade(f"gisaid/output_{i}.json")
    meta.append(a)
    variants.append(b)

a,b = parse_nextclade("gisaid/aecom_nextclade.json")
meta.append(a)
variants.append(b)

da = pd.concat(meta)
db = pd.concat(variants)
db['country'] = db['seqName'].apply(lambda x: x.split("/")[0])
db['aaMutation'] = (db['aaMutation'].apply(lambda x: x[0] if type(x) is tuple else x))
# Remove nosocomial sequences from analysis
da = da.loc[~da['seqName'].isin([" AECOM_127"," AECOM_128", " AECOM_129", "AECOM_130","AECOM_131"])]
db = db.loc[~db['seqName'].isin([" AECOM_127"," AECOM_128", " AECOM_129", "AECOM_130","AECOM_131"])]
# db['position'] += 1

# Shift AA Changes down
db['aaMutation'] = db['aaMutation'].replace({"":None}).apply(lambda x: f"{x[0]}{int(x[1:-1])+1}{x[-1]}" if type(x) is str else 0)
aecom_variants = db[['AECOM_' in x for x in db['seqName']]][['seqName', 'position','mutation','aaMutation', 'gene']].drop_duplicates()
aecom_variants['seqName'] = aecom_variants['seqName'].apply(lambda x: x.strip())
aecom_variants['gene'] = aecom_variants['gene'].apply(lambda x: x[0] if type(x) is tuple else x)

# Consolidate aaMutations for positions that result in two aa mutations
# Find all pairs where position and id are same
grouped_variants = {k:v for k,v in aecom_variants.groupby(['seqName','position']) if len(v) >1}
for k,v in grouped_variants.items():
    if len(v) == 2:
        total_aa_mutation = " / ".join(v['aaMutation'])
        total_genes = " / ".join(v['gene'])
        v = v.iloc[0]
        v['aaMutation'] = total_aa_mutation
        v['gene'] = total_genes
        grouped_variants[k] = v
    elif len(v) > 2:
        assert False

x = pd.DataFrame(grouped_variants.values())

for k in grouped_variants.keys(): # Remove rows from original df
     aecom_variants = aecom_variants.drop(aecom_variants.loc[(aecom_variants['seqName'] == k[0]) & (aecom_variants['position'] == k[1])].index)
     
aecom_variants = aecom_variants.append(x)

aecom_variants_history = {
    row['mutation']: db[(db['position']==row['position']) & (db['mutation']==row['mutation'])]
    for _, row in tqdm(aecom_variants.iterrows(), total=len(aecom_variants))
}

rare_variants = {k:v for k,v in aecom_variants.groupby('mutation') if len(v) == 1}  
uncommon_variants =  {k:v for k,v in aecom_variants.groupby('mutation') if len(v) > 1 and len(v) < 20}
common_variants =  {k:v for k,v in aecom_variants.groupby('mutation') if len(v) >= 20}

# Consolidate variants which cause multiple AA Changes

# rare_alleles = pd.DataFrame([
#     {"position": v['seqName'],  "reference": v.iloc[0]['reference'], "alternate": v.iloc[0]['alternates'], "count": len(v)}
#     for k,v in rare_variants.items()
#     ])
# 
# uncommon_alleles = pd.DataFrame([
#     {"position": k, "reference": v.iloc[0]['reference'], "alternate": v.iloc[0]['alternates'], "count": len(v)}
#     for k,v in uncommon_variants.items()
#     ])
# 
# common_alleles = pd.DataFrame([
#     {"position": k, "reference": v.iloc[0]['reference'], "alternate": v.iloc[0]['alternates'], "count": len(v)}
#     for k,v in common_variants.items()
#     ])

world_counts_by_position = db['position'].value_counts().sort_index()
aecom_counts_by_position = db[['AECOM_' in x for x in db['seqName']]]['position'].value_counts().sort_index()

# Interpolate missing positions with 0
world_counts_by_position = pd.Series({i: world_counts_by_position.get(i) or 0 for i in range(29904)})
aecom_counts_by_position = pd.Series({i: aecom_counts_by_position.get(i) or 0 for i in range(29904)})

world_frequencies_by_position = world_counts_by_position / len(set(db['seqName']))
aecom_frequencies_by_position = aecom_counts_by_position / len(set(db[['AECOM_' in x for x in db['seqName']]]['seqName']))

counts_in_world = pd.Series({k: len(v) for k,v in aecom_variants_history.items()}).sort_values()[1:]
counts_in_aecom = pd.Series({k: sum(['AECOM_' in x for x in v['seqName']])  for k,v in aecom_variants_history.items()}).sort_values()[1:]

overlap = np.minimum(aecom_frequencies_by_position, world_frequencies_by_position)

metadata.index = metadata.index.map(lambda x: x.replace("-", "_"))
# coocurrence = pd.crosstab(index=vardf.index, columns=vardf['position'])
coocurrence = pd.crosstab(index=aecom_variants.seqName, columns=aecom_variants.position)
coocurrence.index = coocurrence.index.map(lambda x: x.replace("_0","_").replace("_0","_"))
coocurrence = coocurrence.loc[coocurrence.index.intersection(passed)]
dates = metadata.loc[coocurrence.index]['collection_date'].sort_values()
coocurrence = coocurrence.loc[dates.index].astype(bool).astype(int)
coocurrence.index.name="Sample ID"

counts = np.zeros(29801)
for k,v in coocurrence.sum().iteritems():
    counts[k] += v
counts = pd.Series(counts)

ac = coocurrence.copy()
aecom_variants.index = aecom_variants['position']

# ac.index = [c.split("_")[1] + ("-----" if i%2 == 1 else '') for i,c in zip(count(),coocurrence.index)]
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

# ax_a2 = ax_a.twinx()
# ax_a2.set_ylabel("Date")
# ax_a2.tick_params(axis='y')
# ax_a2.set_yticklabels(reduced_y_ticks)
# ax_a2.set_yticks(reduced_y_tick_positions)


waves = pd.DataFrame({'Diff': coocurrence.iloc[0:52].sum() - coocurrence.iloc[52:].sum(), 'Early': coocurrence.iloc[0:52].sum(), 'Late': coocurrence.iloc[52:].sum(), 'Total': coocurrence.sum()}).sort_values(by='Diff')
high = list(waves[waves['Total']>=26].index)
wave_1 = list(waves[(waves['Late']<=1) & (waves['Total']>3)].index)
wave_2 = list(waves[(waves['Early']<=1) & (waves['Total']>3)].index)
wave_both = list(waves[(waves['Early']>1) & (waves['Late']>1) & (waves['Total']>3) & (waves['Total'] < 26)].index)

from matplotlib import colors as mcolors

pallette = {
    'black': np.array(mcolors.hex2color(mcolors.CSS4_COLORS['black'])),
    'white': np.array(mcolors.hex2color(mcolors.CSS4_COLORS['white'])),
    'blue': np.array(mcolors.hex2color(mcolors.CSS4_COLORS['blue'])),
    'red': np.array(mcolors.hex2color(mcolors.CSS4_COLORS['red'])),
    'orange': np.array(mcolors.hex2color(mcolors.CSS4_COLORS['orange'])),
    'green': np.array(mcolors.hex2color(mcolors.CSS4_COLORS['green'])),
}

ac_rgb = ac.copy()

for col in ac.columns.difference(high+wave_1+wave_2+wave_both):
    ac_rgb[col] = ac_rgb[col].map({0: pallette['black'], 1: pallette['white']})

for col in high:
    ac_rgb[col] = ac_rgb[col].map({0: pallette['black'], 1: pallette['green']})

for col in wave_1:
    ac_rgb[col] = ac_rgb[col].map({0: pallette['black'], 1: pallette['blue']})

for col in wave_2:
    ac_rgb[col] = ac_rgb[col].map({0: pallette['black'], 1: pallette['red']})

for col in wave_both:
    ac_rgb[col] = ac_rgb[col].map({0: pallette['black'], 1: pallette['orange']})
important = high + wave_1 + wave_2 + wave_both     
rgb = np.zeros((104,357,3))
rgb = np.zeros((104,27,3))
for i, x in enumerate(ac_rgb.index): 
    for j,y in enumerate(sorted(important)): # ac_rgb.columns): 
        for z in range(3): 
            rgb[i,j,z] = ac_rgb[y][x][z] 

# ac.columns = [c if ac[c].sum() >= 4 else '' for i,c in zip(count(),coocurrence.columns)]
# sns.heatmap(ac,yticklabels=False, xticklabels=True, cbar=False, ax=ax_a)

def expand(positions: np.array, spacing: float = .3, delta: float = 1.1) -> np.array:
    """
    Adds a small amount of spacing between tick positions that are too close.
    """
    deltas = [positions[i+1]-positions[i] for i in range(len(positions)-1)]    
    total = sum(deltas)
    deltas = [d if d > delta else d + spacing for d in deltas]
    new_total = sum(deltas)
    deltas = [d * total / new_total for d in deltas]

    return np.array([positions[0] + sum(deltas[:i]) for i in range(len(deltas)+1)])

plt.clf()
plt.close()
import matplotlib
matplotlib.rcParams.update({'font.size':9})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(56,24))
gs = fig.add_gridspec(20,20)

ax_a = fig.add_subplot(gs[:, :])

# ax_a.imshow(np.swapaxes(rgb,0,1))
ax_a.imshow(rgb)

xticks = {i: f"{aecom_variants['mutation'].drop_duplicates().get(c)} {aecom_variants['aaMutation'].replace(0,'').apply(lambda x: f'({x})').drop_duplicates().get(c)}".replace(" None","").replace("()", "")  for i,c in enumerate(ac.columns) if ac[c].sum() >= 4}
xticks = {i: f"{aecom_variants['mutation'].drop_duplicates().get(c)} {aecom_variants['aaMutation'].replace(0,'').apply(lambda x: f'({x})').drop_duplicates().get(c)}".replace(" None","").replace("()", "")  for i,c in enumerate(sorted(important)) if ac[c].sum() >= 4}

# Nextclade shifts variant positions down 1 so we adjust here
to_replace = list(map(lambda x: x.split(" ")[0][1:-1], (xticks.values())))
replacement = list(map(lambda x: str(int(x.split(" ")[0][1:-1]) + 1), (xticks.values())))
for k, to_rep, repl in zip(xticks.keys(), to_replace, replacement):
    xticks[k] = xticks[k].replace(to_rep, repl)

# There is one position where multiple variants were found (27602)
# xticks[293.0] = 'C27602 A/T'

ax_a.set_xticks(expand(list(xticks.keys()), spacing=.15))
ax_a.set_xticklabels(list(xticks.values()), rotation=90)

# ax_a.set_ylabel("Date")
ax_a.tick_params(axis='x')
ax_a.set_yticklabels(list(reversed(reduced_y_ticks)))
ax_a.set_yticks([i*103 for i in reduced_y_tick_positions])

# plt.savefig(savefile)
plt.show()
plt.clf()

fig, ax_a1 = plt.subplots(figsize=(28,4))

# coocurrence.sum().plot.bar(ax=ax_a1)
aecom_frequencies_by_position[ac.columns].plot.bar(ax=ax_a1,label='AECOM',color='green', alpha=.8, width=1.)
world_frequencies_by_position[ac.columns].plot.bar(ax=ax_a1,label='World',color='blue', alpha=.8, width=1.)
overlap[ac.columns].plot.bar(ax=ax_a1,color='cyan', label='Overlap', width=1.)
ax_a1.legend()
# ax_a1.set_xlabel("Variant Frequency")
ax_a1.xaxis.set_ticks([])
ax_a1.yaxis.set_ticks([])
plt.axis('off')
plt.savefig("figures_final/figure4a_top.pdf", facecolor='black')
plt.savefig("figures_final/figure4a_top.svg", facecolor='black')
plt.show()
plt.clf()

from dna_features_viewer import BiopythonTranslator

class CustomTranslator(BiopythonTranslator):
    def compute_feature_fontdict(self, feature):
        return {'size': 12}

fig, ax_a3 = plt.subplots(figsize=(24,2))

graphic_record = CustomTranslator().translate_record("data/external/covid_track.gff")
graphic_record.plot(ax=ax_a3)

plt.savefig("figures_final/figure4a_bottom.pdf")
plt.savefig("figures_final/figure4a_bottom.svg")
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

from plot_helper import Plotter

boocurrence = coocurrence.copy()
boocurrence.index = dates

with Plotter(filename="data/processed/variants/low_level_variants_ordered.pdf", figsize=(20,20)) as ax:
    coocurrence[coocurrence.sum().sort_values()[:-33].index].cumsum().plot(ax=ax, legend=False)
    ax.set_title("Variants showing up < 4 Times Chronologically")
    ax.set_xlabel("Sample Number (Ordered)")
    ax.set_ylabel("Cumulative Count")

with Plotter(filename="data/processed/variants/low_level_variants_timeline.pdf", figsize=(20,20)) as ax:
    boocurrence[coocurrence.sum().sort_values()[:-33].index].cumsum().plot(ax=ax, legend=False)
    ax.set_title("Variants showing up < 4 Times Temporally")
    ax.set_xlabel("Date")
    ax.set_ylabel("Cumulative Count")

with Plotter(filename="data/processed/variants/high_level_variants_ordered.pdf", figsize=(20,20)) as ax:
    coocurrence[coocurrence.sum().sort_values()[-33:].index].cumsum().plot(ax=ax, legend=False)
    ax.set_title("Variants showing up >= 4 Times Chronologically")
    ax.set_xlabel("Sample Number (Ordered)")
    ax.set_ylabel("Cumulative Count")

with Plotter(filename="data/processed/variants/high_level_variants_timeline.pdf", figsize=(20,20)) as ax:
    boocurrence[coocurrence.sum().sort_values()[-33:].index].cumsum().plot(ax=ax, legend=False)
    ax.set_title("Variants showing up >= 4 Times Temporally")
    ax.set_xlabel("Date")
    ax.set_ylabel("Cumulative Count")

waves = pd.DataFrame({'Diff': coocurrence.iloc[0:52].sum() - coocurrence.iloc[52:].sum(), 'Early': coocurrence.iloc[0:52].sum(), 'Late': coocurrence.iloc[52:].sum(), 'Total': coocurrence.sum()}).sort_values(by='Diff')

high = waves[waves['Total']>=26].index
wave_1 = waves[(waves['Late']<=1) & (waves['Total']>3)].index  
wave_2 = waves[(waves['Early']<=1) & (waves['Total']>3)].index
wave_both = waves[(waves['Early']>1) & (waves['Late']>1) & (waves['Total']>3) & (waves['Total'] < 26)].index

from labellines import labelLine, labelLines
from matplotlib.lines import Line2D

matplotlib.rcParams.update({'font.size':16})
plt.clf()
fig, ax = plt.subplots(figsize=(30,30))

coocurrence[coocurrence.sum().sort_values()[high].index].cumsum().plot(ax=ax, legend=False, color='green', alpha=.8, linewidth=3)
coocurrence[coocurrence.sum().sort_values()[wave_1].index].cumsum().plot(ax=ax, legend=False, color='blue', alpha=.8, linewidth=3)
coocurrence[coocurrence.sum().sort_values()[wave_both].index].cumsum().plot(ax=ax, legend=False, color='orange', alpha=.8, linewidth=3)
coocurrence[coocurrence.sum().sort_values()[wave_2].index].cumsum().plot(ax=ax, legend=False, color='red', alpha=.8, linewidth=3)

xvals = \
    [i/len(high) * 100 for i in range(len(high))] + \
    [i / len(wave_1) * 45 for i in range(len(wave_1))] + \
    [i / len(wave_both) * 45 + 45 for i in range(len(wave_both))] + \
    [i / len(wave_2) * 45 + 45 for i in range(len(wave_2))]


labelLines(plt.gca().get_lines(),zorder=2.5,fontsize=18, align=True, xvals=xvals)
ax.set_title("Variants showing up Over 3 Times Longitudinally")
ax.set_xlabel("Sample ID Ordered By Collection Date")
ax.set_ylabel("Cumulative Count")
custom_lines = [
    Line2D([0], [0], color='green', lw=3),
    Line2D([0], [0], color='blue', lw=3),
    Line2D([0], [0], color='orange', lw=3),
    Line2D([0], [0], color='red', lw=3),
]

ax.legend(custom_lines, ['Common', 'Wave 1', 'Wave 1 + 2', 'Wave 2'])
plt.savefig("figures_final/figure4b.pdf")
plt.savefig("figures_final/figure4b.svg")
plt.show()
plt.clf()

fig, ax = plt.subplots(figsize=(30,30))

boocurrence[coocurrence.sum().sort_values()[high].index].cumsum().plot(ax=ax, legend=False, color='green', alpha=.8, linewidth=3)
boocurrence[coocurrence.sum().sort_values()[wave_1].index].cumsum().plot(ax=ax, legend=False, color='blue', alpha=.8, linewidth=3)
boocurrence[coocurrence.sum().sort_values()[wave_both].index].cumsum().plot(ax=ax, legend=False, color='orange', alpha=.8, linewidth=3)
boocurrence[coocurrence.sum().sort_values()[wave_2].index].cumsum().plot(ax=ax, legend=False, color='red', alpha=.8, linewidth=3)

xvals = \
    [i/len(high) * 100 for i in range(len(high))] + \
    [i / len(wave_1) * 33 for i in range(len(wave_1))] + \
    [i / len(wave_both) * 33 + 33 for i in range(len(wave_both))] + \
    [i / len(wave_2) * 33 + 66 for i in range(len(wave_2))]

xvals = [ i/100 * (boocurrence.index.max()-boocurrence.index.min()) + boocurrence.index.min() for i in xvals]

labelLines(plt.gca().get_lines(),zorder=2.5,fontsize=18, align=True, xvals=xvals)
ax.set_title("Variants showing up Over 3 Times Longitudinally")
ax.set_xlabel("Date")
ax.set_ylabel("Cumulative Count")
custom_lines = [
    Line2D([0], [0], color='green', lw=3),
    Line2D([0], [0], color='blue', lw=3),
    Line2D([0], [0], color='orange', lw=3),
    Line2D([0], [0], color='red', lw=3),
]

ax.legend(custom_lines, ['Common', 'Wave 1', 'Wave 1 + 2', 'Wave 2'])

plt.savefig("figures_final/figure4b_supplemental.pdf")
plt.savefig("figures_final/figure4b_supplemental.svg")
plt.show()

counts_dict = {x:b[b['position']==x] for x in counts.index.map(lambda x: x-1)}
aa = pd.Series({k+1: v.iloc[0]['aaMutation'] if len(v) else None for k,v in counts_dict.items()})
dam = pd.DataFrame({"count": counts, "aa": aa})

# Make table
names = {c: aecom_variants['mutation'].drop_duplicates().get(c) for c in important}
names[27602] = names[27602].iloc[0]
aanames = {c: aecom_variants['aaMutation'].drop_duplicates().get(c) for c in important}
genenames = {c: aecom_variants['gene'].get(c).iloc[0] for c in important}

dx = pd.DataFrame({'nucleotide': names, 'amino acid': aanames, 'gene': genenames})
dx['aecom_frequency'] = aecom_frequencies_by_position.loc[important]
dx['world'] = world_frequencies_by_position.loc[important]
dx['type'] = None
dx['type'].loc[wave_1] = 'wave 1' 
dx['type'].loc[wave_2] = 'wave 2' 
dx['type'].loc[wave_both] = 'wave both' 
dx['type'].loc[high] = 'high' 
dx.index += 1
dx['nucleotide'] = dx['nucleotide'].apply(lambda x: x[0] + str(int(x[1:-1])+1) + x[-1])

aecom_variants.index += 1
aecom_variants['position'] += 1
aecom_variants['mutation'] = aecom_variants['mutation'].apply(lambda x: x[0] + str(int(x[1:-1])+1) + x[-1] if type(x) is str else x)

positions = pd.DataFrame({
    "gene": x.label,
    "start": x.start.position,
    "end": x.end.position,
}
    for x in graphic_record.features
)

positions_genes = {
    k: " / ".join([r['gene'] for _, r in positions.iterrows() if (k>= r['start']) and (k<= r['end'])])
    for k in ac.columns
}

ac.columns += 1

orfs = {x:aecom_variants['gene'].loc[x] for x in ac.columns}
orfs2 = {}
for k, v in orfs.items():
    if type(v) is str:
        orfs2[k] = v
    elif type(v) is float:
        orfs2[k] = ''
    elif type(v) is pd.Series:
        orfs2[k] = v.iloc[0]

for k,v in orfs2.items():
    if v is np.nan:
        orfs2[k] = ''

# Interpolate ORFs that are missing
min_vals = {}
max_vals = {}
for k,v in orfs2.items():
    if v not in min_vals:
        min_vals[v] = k
    max_vals[v] = k

for k,v in orfs2.items():
    if v == "":
        # Try to fill it in
        for gene, mi in min_vals.items():
            if k >= mi:
                if max_vals[gene] >= k:
                    orfs2[k] = gene

orfs2 = pd.Series(orfs2)
orfs2.to_csv("aecom_variant_genes.csv")