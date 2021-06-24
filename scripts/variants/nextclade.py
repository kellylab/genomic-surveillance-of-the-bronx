# Parses NextClade JSON output

import json
from tqdm import tqdm
import pandas as pd
from typing import Tuple, Dict
from collections import defaultdict
import matplotlib.pyplot as plt
import plot_helper
from plot_helper import Plotter
import numpy as np
from covid_bronx.variants import parse_nextclade

plot_helper.BASE_FILEPATH = "data/processed/nextclade"

import matplotlib
matplotlib.rcParams.update({'font.size':9})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

if __name__ == "__main__":

    meta = []
    variants = []

    for i in tqdm(range(14)):
        a,b = parse_nextclade(f"gisaid/output_{i}.json")
        meta.append(a)
        variants.append(b)

    a,b = parse_nextclade("gisaid/aecom_nextclade.json")
    a = a[a['seqName'] != ' AECOM_D2']
    b = b[b['seqName'] != ' AECOM_D2']
    meta.append(a)
    variants.append(b)

    da = pd.concat(meta)
    db = pd.concat(variants)
    db['country'] = db['seqName'].apply(lambda x: x.split("/")[0])

    # Remove nosocomial sequences from analysis
    da = da.loc[~da['seqName'].isin([" AECOM_127"," AECOM_128", " AECOM_129", "AECOM_130","AECOM_131"])]
    db = db.loc[~db['seqName'].isin([" AECOM_127"," AECOM_128", " AECOM_129", "AECOM_130","AECOM_131"])]

    with Plotter(filename="variant_histogram.pdf", figsize=(16,16)) as ax:
        s = pd.Series([np.mean(db['position'].value_counts() == i) for i in range(20)]) * 100
        s.plot.bar(ax=ax)
        ax.set_title("Distribution of Variant Presence Frequencies")
        ax.set_xlabel("Number of Times Variant is Seen in GISAID")
        ax.set_ylabel("Percent of Total Variants")

    db['position'] += 1 # Adjust positions by 1

    aecom_variants = db[['AECOM_' in x for x in db['seqName']]][['position','mutation']].drop_duplicates()
    aecom_variants_history = {
        row['mutation']: db[(db['position']==row['position']) & (db['mutation']==row['mutation'])]
        for _, row in tqdm(aecom_variants.iterrows(), total=len(aecom_variants))
    }
    
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
    
    with Plotter(filename="variant_distribution_vs_world.pdf", figsize=(16,16)) as ax:
        aecom_frequencies_by_position.plot(ax=ax,label='AECOM',color='blue')
        world_frequencies_by_position.plot(ax=ax,label='World',color='green')
        overlap.plot(ax=ax,color='yellow')
        ax.set_title("Distribution of Variant Presence By Position")
        ax.set_xlabel("Nucleotide")
        ax.set_ylabel("Percent of Samples Having Variant (Overlap in Yellow)")
        plt.legend()

    relative_differences = np.abs(aecom_frequencies_by_position-world_frequencies_by_position).sort_values(ascending=False)

    vdf = db[['AECOM_' in x for x in db['seqName']]][['seqName', 'position']]
    vardf = pd.DataFrame(index=set(vdf['seqName']), columns=set(vdf['position'])).sort_index().fillna(0)
    vardf.columns = vardf.columns.sort_values()

    for _, row in tqdm(vdf.iterrows(), total=len(vdf)):
        vardf.loc[row['seqName']][row['position']] = 1

    # The variants at each position are always the same except at position
    # 27602, where we have a C27602A and a C27602T
    print(aecom_variants[aecom_variants['position']==27602])

    block_1 = [
        28917,
        21245,
        28767,
        5180,
        16466,
        # 27679,
        # 27389,
        # 23065,
        # 24138,
    ]

    block_2 = [
        28881,
        28882,
        28883,
    ]
    
    # [
    #     28766,
    #     21244,
    #     28916,
    #     16465,
    #     5179,
    #     23064,
    #     24137,
    #     27388,
    #     27678,
    #     11513,
    #     10447,
    #     13664,
    #     2783,
    #     23063,
    # ]

    country_counts = pd.Series(db['seqName'].apply(lambda x: x.split("/")[0]).value_counts())

    block_1_variants = aecom_variants[aecom_variants['position'].isin(block_1)]['mutation']
    
    block_1_in_world = {k: aecom_variants_history[k] for k in block_1_variants}
    block_1_in_world = {k: v[ ['AECOM_' not in x for x in v['seqName']] ] for k,v in block_1_in_world.items()}
    
    dx = pd.concat(block_1_in_world.values())
    block_1_countries = dx['seqName'].apply(lambda x: x.split("/")[0])
    dx['country'] = dx['seqName'].apply(lambda x: x.split("/")[0])
    block_1_country_counts = block_1_countries.value_counts()
    block_1_country_ratios = block_1_country_counts / country_counts[block_1_country_counts.index]

    with Plotter(filename="block_1_countries.pdf", figsize=(16,16)) as ax:
        s = block_1_country_ratios.sort_values(ascending=False)
        s.plot.barh(ax=ax)
        ax.set_title("Countries Where Block 1 Variants Were Also Found")
        ax.set_ylabel("Country")
        ax.set_xlabel("Frequency")

    with Plotter(filename="Global_Sequences_With_Block_1_Variants.pdf", figsize=(30,16)) as ax:
        db[db['position'].isin(block_1)]['seqName'].value_counts().loc[lambda x: x>1].plot.barh()
        ax.set_title("Sequences Containing Block 1 Variants")
        ax.set_xlabel("Count")
        ax.set_ylabel("Sequence ID")
    
    block_2_variants = aecom_variants[aecom_variants['position'].isin(block_2)]['mutation']
    
    block_2_in_world = {k: aecom_variants_history[k] for k in block_2_variants}
    block_2_in_world = {k: v[ ['AECOM_' not in x for x in v['seqName']] ] for k,v in block_2_in_world.items()}

    dx = pd.concat(block_2_in_world.values())
    block_2_countries = dx['seqName'].apply(lambda x: x.split("/")[0])
    dx['country'] = dx['seqName'].apply(lambda x: x.split("/")[0])

    block_2_country_counts = block_2_countries.value_counts()
    block_2_country_ratios = block_2_country_counts / country_counts[block_2_country_counts.index]

    with Plotter(filename="block_2_countries.pdf", figsize=(16,16)) as ax:
        s = block_2_country_ratios.sort_values(ascending=False)
        s.plot.barh(ax=ax)
        ax.set_title("Countries Where Block 2 Variants Were Also Found")
        ax.set_ylabel("Country")
        ax.set_xlabel("Frequency")

    with Plotter(filename="Global_Sequences_With_Block_2_Variants.pdf", figsize=(20,16)) as ax:
        db[db['position'].isin(block_2)]['seqName'].value_counts().loc[lambda x: x>1].plot.barh()
        ax.set_title("Sequences Containing Block 2 Variants")
        ax.set_xlabel("Count")
        ax.set_ylabel("Sequence ID")