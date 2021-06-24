# This script creates a single figure which depicts the metadata association
# boxplots.

import pandas as pd
import numpy as np
import plot_helper
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("data/external/pangolin.csv")
df.index = df['Sequence name'].apply(lambda x: x.split(" ")[0])

df = pd.read_csv("data/external/Monte_clade_meta.csv")
df = df.loc[df['Sequence name'].dropna().index]
df.index = df['Sequence name'].apply(lambda x: x.split(" ")[0]).apply(lambda x: x.replace("-0","-").replace("-0","-"))

age = pd.read_csv("data/processed/metadata/age.tsv", sep="\t")
index = age['name']
age.index = index
sex = pd.read_csv("data/processed/metadata/sex.tsv", sep="\t")
sex.index = index
num_comorbidities = pd.read_csv("data/processed/metadata/num_comorbidities.tsv", sep="\t")
num_comorbidities.index = index
num_side_effects = pd.read_csv("data/processed/metadata/num_side_effects.tsv", sep="\t")
num_side_effects.index = index
num_treatments = pd.read_csv("data/processed/metadata/num_treatments.tsv", sep="\t")
num_treatments.index = index
zip_code = pd.read_csv("data/processed/metadata/zipcode.tsv", sep="\t")
zip_code.index = index
outcome = pd.read_csv("data/processed/metadata/outcome.tsv", sep="\t")
outcome.index = index

metadata = pd.DataFrame(index=sex['name'])
metadata['sex'] = sex['annotation']
metadata['age'] = age['annotation']
metadata['num_side_effects'] = num_side_effects['annotation']
metadata['num_comorbidities'] = num_comorbidities['annotation']
metadata['num_treatments'] = num_treatments['annotation']
metadata['zip_code'] = zip_code['annotation']
metadata['outcome'] = outcome['annotation']
metadata[df.columns] = df

clades = pd.read_csv("data/external/nextclade.csv", sep=";")
clades.index = clades['seqName'].apply(lambda x: x.split(" ")[0])
metadata[clades.columns] = clades
coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
passed = coverage_levels[coverage_levels>=.95].index.map(lambda x: x.replace("_","-").replace("-0", "-").replace("-0", "-"))
metadata = metadata.loc[passed]
metadata.to_csv("data/processed/metadata.csv")

metadata2 = pd.read_csv("data/confidential/sample_metadata.csv")
metadata2.index = metadata2.number.map(lambda x: f"AECOM-{x}")
metadata2 = metadata2.loc[metadata2.index.intersection(passed)]

fig, axs = plt.subplots(2,5, figsize=(30,20))
for j, column in enumerate(['age', 'num_side_effects', 'num_comorbidities', 'num_treatments', 'Days since last sampling']):
    for i, category in enumerate(['Lineage']):
        ax = axs[i,j]
        # sns.boxplot(x=column, y=category, data=metadata, ax=ax)
        sns.violinplot(x=column, y=category, data=metadata, ax=ax)
        ax.set_xlabel('')
        ax.set_ylabel('')
        # Selectively remove ticks
        if i == 0: # Remove xticks from top row
            ax.set_xticks([])
        if j != 0: # Remove yticks from every column except first one
            ax.set_yticks([])
## Label Axes
# x
for j, title in enumerate(["Age", "Number of Side Effects", "Number of Comorbidities", "Number of Treatments"]):
    axs[1,j].set_xlabel(title)
# y
for i, category in enumerate(['Lineage', "Clade"]):
    axs[i,0].set_ylabel(category)
fig.subplots_adjust(wspace=0, hspace=0)
fig.suptitle("Associations Between SARS-CoV2 Strains and Patient Metadata")
plt.savefig("figures/metadata.png", dpi=300)
plt.clf()
