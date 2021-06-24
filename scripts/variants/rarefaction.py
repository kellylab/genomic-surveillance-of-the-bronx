# Makes rarefaction plots for a set of genomes.

from covid_bronx.variants import parse_nextclade
import pandas as pd
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import plot_helper
from plot_helper import Plotter

input_file = 'wadsworth/bronx.json'
output_folder = 'wadsworth'

plot_helper.BASE_FILEPATH = output_folder
plot_helper.SHOW = True

df = pd.read_csv("wadsworth/GISIAD_Meta_Bronx.csv")

lineages = df['Lineage']
lineages.index = df['Collection date']
lineages = lineages.sort_index()

date_range = pd.date_range(lineages.index.min(), lineages.index.max())

pdf = pd.DataFrame(index=date_range, columns = set(lineages)).fillna(0)

for d,l in lineages.items():
    try:
        pdf.loc[d,l] += 1
    except:
        pass

cdf = pdf.cumsum()
cdf = cdf[cdf.iloc[-1].sort_values(ascending=False).index]

with Plotter(filename='lineages_gisaid.pdf', figsize=(20,20)) as ax:
    cdf.loc[:, cdf.iloc[-1]>2000].plot(ax=ax, cmap='tab20b')
    ax.set_xlabel("Date")
    ax.set_ylabel("Count")
    ax.set_title("Bronx Rarefaction of Most Common Lineages")

a, b = parse_nextclade(input_file)

meta = []
variants = []

for i in tqdm(range(41)):
    a,b = parse_nextclade(f"gisaid_jan21/output_{i}.json")
    meta.append(a)
    variants.append(b)

da = pd.concat(meta)
db = pd.concat(variants)
db['country'] = db['seqName'].apply(lambda x: x.split("/")[0])
 
dates = {r['Virus name']: r['Collection date'] for _,r in df.iterrows()} 
b['Date'] = b['seqName'].apply(lambda x: x.split("|")[0]).map(dates)

# spike = b[(b['position']>26523)&(b['position']<=27191)]
spike = b[(b['position']>21563)&(b['position']<=25384)]
spike.index = spike['Date']
spike = spike.sort_index()
spike['Variant'] = spike['mutation'] + " " +  spike['aaMutation'].apply(str)
variants = spike['Variant'].dropna()

date_range2 = pd.date_range(spike.index.min(), spike.index.max())

pdf2 = pd.DataFrame(index=date_range, columns = set(variants)).fillna(0)

for d,l in variants.items():
    pdf2.loc[d,l] += 1

cdf2 = pdf2.cumsum()
cdf2 = cdf2[cdf2.iloc[-1].sort_values(ascending=False).index]

with Plotter(filename='spike_variants.pdf', figsize=(20,20)) as ax:
    cdf2.loc[:, (cdf2.iloc[-1]>2) & (cdf2.iloc[-1]<100)].plot(ax=ax, cmap='tab20b')
    ax.set_xlabel("Date")
    ax.set_ylabel("Count")
    ax.set_title("Bronx Rarefaction of Most Common Spike Variants")

def cave_plot(df, ax=None, **kwargs):

    if ax:
        object = ax
    else:
        object = plt

    bottom = 0.
    for col in df.columns:
        top = df[col] + bottom
        object.fill_between(df.index, bottom, top, label=col, **kwargs)
        bottom = top

    return object

ratios = (cdf.T / cdf.T.sum()).T

with Plotter(filename='gisaid_lineages_cave.pdf', figsize=(20,20)) as ax:
    ax = cave_plot(ratios, ax, cmap='tab20b')
    plt.legend(loc='upper left')
    ax.set_title("Ratios of Lineages in World Over Time")
    ax.set_xlabel("Date")
    ax.set_ylabel("Relative Abundance")

ratios2 = (cdf2.T / cdf2.T.sum()).T

with Plotter(filename='variants_cave.pdf', figsize=(20,20)) as ax:
    ax = cave_plot(ratios2.loc[:, (cdf2.iloc[-1]>1)], ax, cmap='tab20b')
    plt.legend(loc='upper left')
    ax.set_title("Ratios of Variants in Bronx Over Time")
    ax.set_xlabel("Date")
    ax.set_ylabel("Relative Abundance")

with Plotter(filename='variants_cave_full.pdf', figsize=(20,20)) as ax:
    ax = cave_plot(ratios2, ax, cmap='tab20b')
    plt.legend(loc='upper left')
    ax.set_title("Ratios of Variants in Bronx Over Time")
    ax.set_xlabel("Date")
    ax.set_ylabel("Relative Abundance")    