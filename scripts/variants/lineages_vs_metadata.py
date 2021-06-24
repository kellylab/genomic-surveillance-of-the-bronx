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

num_samples = len(passed)

color_dict={
'A.1': '#d80416',
'B.1': '#305bbe',
'B.1.3': '#1e6322',
# 'B.1.26': '#a1184f',
'B.2/B.2.1': '#f7d93d',
'B.2': '#f5d327',
'B.2.1': '#debe1d',
}

from covid_bronx.metadata import preprocess_metadata
metadata = preprocess_metadata()
coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
passed = coverage_levels[coverage_levels>=.95].index.map(lambda x: x.replace("_","-").replace("-0", "-").replace("-0", "-"))

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
metadata2.index = metadata2.index.map(lambda x: x.replace("-", "_"))
metadata2 = metadata2.loc[metadata.index]

symptoms = set([y.strip() for x in metadata2['COVID sx'] for y in str(x).split(",")])
comorbidities = set([y.strip() for x in metadata2['Pt co-morbidities'] for y in str(x).split(",")])
treatments = set([y.strip() for x in metadata2['Tx'] for y in str(x).split(",")])

symptom_df = pd.DataFrame({symptom: metadata2['COVID sx'].map(lambda x: symptom in str(x)) for symptom in symptoms})
comorbidities_df = pd.DataFrame({como: metadata2['Pt co-morbidities'].map(lambda x: como in str(x)) for como in comorbidities})
treatments_df = pd.DataFrame({tx: metadata2['Tx'].map(lambda x: tx in str(x)) for tx in treatments})

# Test each variant associations using a contingency test
lineage_groups = {k: v for k,v in metadata.groupby('Lineage')}
lineage_lengths = {k: len(v) for k,v in lineage_groups.items()}
contingencies = {}
for pos, var in tqdm(lineage_groups.items()):
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
statsmodels.stats.multitest import fdrcorrection       
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
        if min(min(data[col].sum(0)), min(data[col].sum(1))) > 2:
            scores[col] = pvalue
            effect[col] = oddsratio
    pvalues[pos] = scores
    effects[pos] = effect

ppe = pd.DataFrame(effects).T.sort_index()
pp = pd.DataFrame(pvalues).T.sort_index()
pps = pp.T.apply(lambda x:pd.Series( fdrcorrection(x.dropna())[1], index=x.dropna().index))
pps = pp.T.loc[((pp<.5).sum() != 0)]
ppee = ppe.T.loc[((pp<.5).sum() != 0)]
ppee = ppee[pps.columns[np.where((pps<.5).sum() != 0)]]
pps = pps[pps.columns[np.where((pps<.5).sum() != 0)]]


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

treatments_columns = {
    'haloperidol': 'Haloperidol',
    'nifedipine': 'Nifedipine',
    'aspirin': 'Aspirin',
    'senna': 'Senna',
    'tamsulosin': 'Tamsulosin',
    'iron': 'Iron',
    'flonase': 'Flonase',
    'hydroxyzin': 'Hydroxyzine',
    'lorazepam': 'Lorazepam',
    'diphenhydramine': 'Diphenhydramine',
    'maalox': 'Maalox'
}


pps = pd.DataFrame({v: pps.loc[k] for k,v in {**demographic_columns, **symptoms_columns, **comorbidities_columns, **treatments_columns}.items() if k in pps.index}).T
ppee = pd.DataFrame({v: ppee.loc[k] for k,v in {**demographic_columns, **symptoms_columns, **comorbidities_columns, **treatments_columns}.items() if k in ppee.index}).T

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
        for i,v in treatments_columns.items():
            if v == index:
                variable_counts[index] = treatments_df[i].sum()
variable_counts = pd.Series(variable_counts)

log_transformed_effects = (np.log(ppee).replace(np.inf, 10).replace(-np.inf, -10) * (pps<.05).astype(int)).replace(0, np.nan)

for i in range(2):
    fig = plt.figure(constrained_layout=True, figsize=(14,8)) #, axs = plt.subplots(nrows=3, ncols=2, figsize=(20,24))
    gs = fig.add_gridspec(20, 18)
    cbar_ax = fig.add_subplot(gs[0,:14])
    heatmap_ax = fig.add_subplot(gs[1:16,:14])
    x_ax = fig.add_subplot(gs[16:, :14])
    y_ax = fig.add_subplot(gs[1:16, 14:16])
    sns.heatmap(data=log_transformed_effects, cmap='coolwarm', linewidths=.5, linecolor='grey', xticklabels=False, yticklabels=True, ax=heatmap_ax, cbar_ax=cbar_ax,cbar_kws={'orientation':'horizontal'})
    #.set_title("Selected Patient Data and Variant Associations Colored by Effect Size (FDR Corrected).")
    pd.Series([lineage_lengths[k] for k in ppee.columns], index=ppee.columns).plot.bar(ax=x_ax, alpha=.6)
    x_ax.invert_yaxis()
    x_ax.xaxis.tick_top()
    x_ax.xaxis.set_label_position('top')
    x_ax.set_ylabel("Count")
    x_ax.set_xlabel("Variant")
    if i == 0:
        pd.DataFrame({'Age':metadata['age']}, columns=list(reversed(variable_counts.index))).plot.box(vert=False, ax=y_ax)
    else:
        variable_counts[::-1].plot.barh(ax=y_ax, alpha=.6)
    y_ax.yaxis.set_ticks([])
    y_ax.set_xlabel("Count")
    plt.savefig(f"data/processed/variants/lineages_vs_metadata_{i}.pdf")
    plt.savefig(f"data/processed/variants/lineages_vs_metadata_{i}.png")
    plt.show()