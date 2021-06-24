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
import gffutils

plot_helper.BASE_FILEPATH = "data/processed/frameshifts"

import matplotlib
matplotlib.rcParams.update({'font.size':9})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

if __name__ == "__main__":

    meta = []
    variants = []

    for i in tqdm(range(41)):
        a,b = parse_nextclade(f"gisaid_jan21/output_{i}.json")
        meta.append(a)
        variants.append(b)

    da = pd.concat(meta)
    db = pd.concat(variants)
    db['country'] = db['seqName'].apply(lambda x: x.split("/")[0])

    deletions_distributions = db[db['type'] == 'deletion']['length'].value_counts().sort_index()
    deletions_distributions.index = deletions_distributions.index.map(lambda x: int)
    with Plotter(filename="deletions.pdf", figsize=(20,10)) as ax:
        deletions_distributions[0:25].plot.bar(ax=ax)
        ax.set_title("Number of deletions of a given length in GISAID")
        ax.set_xlabel("Length of Deletion")
        ax.set_ylabel("Count")
        
                
    deletions = db[db['type'] == 'deletion']
    deletions['length'] = deletions['length'].apply(int)

    frameshift_deletions = deletions.iloc[np.where(deletions['length'] % 3)[0]]

    insertions = db[db['type'] == 'insertion']
    insertions['length'] = insertions['mutation'].apply(lambda x: len(x.lstrip("1234567890")))
    frameshift_insertions = insertions.iloc[np.where(insertions['length'] % 3)[0]]

    all_frameshifts = pd.concat([frameshift_insertions, frameshift_deletions])

    l = 29993
    frameshift_insertion_positions = np.zeros(l)
    frameshift_deletion_positions = np.zeros(l)
    frameshift_positions = np.zeros(l)
    for _, row in frameshift_deletions.iterrows():
        frameshift_positions[row['position']:row['position'] + row['length']] += 1
        frameshift_deletion_positions[row['position']:row['position'] + row['length']] += 1
        
    for _, row in frameshift_insertions.iterrows():
        frameshift_positions[row['position']:row['position'] + row['length']] += 1
        frameshift_insertion_positions[row['position']:row['position'] + row['length']] += 1

    aecom_frameshifts = {
        'AECOM_080': [18546],
        'AECOM_068': [8296],
        'AECOM_034': list(range(20271,20271+58)),
        'AECOM_123': list(range(27970,27970+37)),
        'AECOM_132': list(range(27970,27970+37)),
        'AECOM_070': list(range(27948,27948+20)),
        'AECOM_109': list(range(27766,27766+20)),
        'AECOM_058': [3746,3747],
        'AECOM_059': [15158,15159] + list(range(19749,19749+77)),
        'AECOM_104': list(range(29575,29575+13)),
    }

    aecom_positions = np.zeros(l)

    for v in aecom_frameshifts.values():
        aecom_positions[v] +=1

    x = [i for i,_ in enumerate(frameshift_positions)] 

    with Plotter(filename="aecom_vs_world_frameshifts.pdf", figsize=(16,10)) as ax:
        ax.plot(frameshift_positions, label='Global')
        ax.fill_between(x, (aecom_positions != 0).astype(int) * frameshift_positions.max(), color='orange', label='AECOM', alpha=.5)
        ax.set_title("Counts of Frameshifts in GISAID")
        ax.set_xlabel("Position on Genome")
        ax.set_ylabel("Count (AECOM Frameshift Positions Highlighted in Orange)")

    shifts_by_country = pd.DataFrame({
        "Frameshift Ratio": pd.Series({c: len(set(frameshifts[frameshifts['country']==c]['seqName'])) / len(set(db[db['country']==c]['seqName'])) for c in tqdm(set(db['country']))}),
        'count': frameshifts['country'].value_counts()
    }).fillna(0)
    
    with Plotter(filename="frameshift_frequency_by_country.pdf", figsize=(10,20)) as ax:
        shifts_by_country[shifts_by_country['count'] > 1000]['Frameshift Ratio'].plot.barh(ax=ax)
        ax.axvline(x=.1, linestyle='--')
        ax.set_title("Average Frequency of Frameshifts per Sample by Country (With at Least 1000 Samples)")
        ax.set_ylabel("Country")
        ax.set_xlabel("Frequency")

    import seaborn as sns
    with Plotter(filename="frameshift_distributions_by_country.pdf", figsize=(10,20)) as ax:
        #sns.boxplot(shifts_by_country[shifts_by_country['count'] > 200]['Frameshift Ratio'],ax=ax)
        sns.swarmplot(shifts_by_country[shifts_by_country['count'] > 200]['Frameshift Ratio'],ax=ax)
        ax.set_title("Distribution of Frameshift Frequencies Among Countries With at Least 200 Samples")


    top = shifts_by_country[shifts_by_country['count']>200]['Frameshift Ratio']  # Get the ratios for countries with at least 200 samples
    
    with Plotter(filename="frameshift_distplot_by_country.pdf", figsize=(8,8)) as ax:
        sns.distplot(top,ax=ax)
        ax.set_title("Distribution of Frameshift Frequencies Among Countries With at Least 200 Samples")

    with Plotter(filename="frameshift_log_distplot_by_country.pdf", figsize=(9,9)) as ax:
        sns.distplot(np.log(top[top>0]),ax=ax)
        ax.set_title("Distribution of Log-Frequency of Frameshifts Among Countries With at Least 200 Samples")

    from statsmodels.graphics.gofplots import qqplot     
    with Plotter(filename="frameshift_log_qq_by_country.pdf", figsize=(8,8)) as ax:        
        qqplot(np.log(top[top>0]),ax=ax)
        ax.set_title("The Distribution of Frequencies of Frameshift Mutations is Log-Normally Distributed")

    dab = gffutils.create_db("data/external/GCF_009858895.2_ASM985889v3_genomic.gff", "data/external/annotations5.db", keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    df = pd.DataFrame([
        {
            "id": f.id,
            "start": f.start,
            "end": f.end,
            "type": f.featuretype,
            "chromosome": f.chrom,
            "strand": f.strand,
        }
        for f in dab.all_features()
    ])

    gene_names = {
        "gp01": 'ORF1ab',
        "gp02": 'ORF2',
        "gp03": 'ORF3a',
        "gp04": 'Envelope Protein',
        "gp05": 'Membrane Glycoprotein',
        "gp06": 'ORF6',
        "gp07": 'ORF7a',
        "gp08": 'ORF7b',
        "gp09": 'ORF8',
        "gp10": 'Nucleocapsid Phosphoprotein',
        "gp11": 'ORF10',
    }

    with Plotter(filename="world_frameshifts_by_orf.pdf", figsize=(16,10)) as ax:
        ax.plot(frameshift_insertion_positions, label='Insertion Frameshifts',alpha=.5)
        ax.plot(frameshift_deletion_positions, label='Deletion Frameshifts', alpha=.5)
        ax.set_title("Counts of Frameshifts in GISAID")
        ax.set_xlabel("Position on Genome")
        ax.set_ylabel("Log-Count (AECOM Frameshift Positions Highlighted in Orange)")
        for _, row in df[df['type']=='gene'].iterrows():
            ax.fill_between(x, [frameshift_positions.max()*(_ in range(row['start'],row['end'])) for _ in x], label=gene_names[row['id'].split('_')[-1]], alpha=.5)
        ax.set_yscale('log')
        ax.legend()

    # Save Down the sequences which have frameshifts
    frameshifts_by_gene = {}
    for _, row in tqdm(df[df['type'] == 'gene'].iterrows()):
        start = row['start']
        stop = row['end']
        fs_by_gene = all_frameshifts[(all_frameshifts['position']>=start) & (all_frameshifts['position']<=stop)].append(all_frameshifts[(all_frameshifts['position'] + all_frameshifts['length']>=start) & (all_frameshifts['position']+all_frameshifts['length']<=stop)]).drop_duplicates()
        name = row['id'].split('_')[-1]
        frameshifts_by_gene[name] = fs_by_gene
        fs_by_gene.to_csv(f"data/processed/{name}_sequences.csv")
