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
from plot_helper import Plotter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# a) Variant clustermap

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
passed = coverage_levels[coverage_levels>=.95].index.intersection(sam_files.keys()).map(lambda x: x.replace("_","-").replace("-0", "-").replace("-0", "-"))

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

metadata.index = metadata.index.map(lambda x: x.replace("-", "_"))
coocurrence = pd.crosstab(index=vardf.index, columns=vardf['position'])
dates = metadata.loc[coocurrence.index]['collection_date'].sort_values()
coocurrence = coocurrence.loc[dates.index].astype(bool).astype(int)
coocurrence.index.name="Sample ID"

counts = np.zeros(29801)
for k,v in coocurrence.sum().iteritems():
    counts[k] += v
counts = pd.Series(counts)
ac = coocurrence.copy()

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

# high = [1059, 25563, 241, 3037, 14408, 23403]

#wave_1 = [18457, 18322, 21245, 28917, 28767, 21697, 13914, 23065, 11514, 26885,]

#wave_2 = [27389, 27679, 24138,  2416, 21697, 23064,  1348, 18877, 28882, 28883, 5180, 28887, 28881, 28917, 28767]

# wave_both = [27603, 16616,]

high = waves[waves['Total']>=26].index
wave_1 = waves[(waves['Late']<=1) & (waves['Total']>3)].index  
wave_2 = waves[(waves['Early']<=1) & (waves['Total']>3)].index
wave_both = waves[(waves['Early']>1) & (waves['Late']>1) & (waves['Total']>3) & (waves['Total'] < 26)].index

from labellines import labelLine, labelLines
from matplotlib.lines import Line2D
import matplotlib
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
