import pandas as pd
import numpy as np
import plot_helper
from plot_helper import Plotter
import seaborn as sns
from covid_bronx.metadata import preprocess_metadata
import matplotlib.pyplot as plt

plot_helper.SHOW = False

color_dict={
'A.1': '#d80416',
'B.1': '#305bbe',
'B.1.3': '#1e6322',
'B.1.26': '#a1184f',
'B.2/B.2.1': '#f7d93d',
'B.2': '#f5d327',
'B.2.1': '#debe1d',
}

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
metadata.to_csv("data/processed/metadata.csv")

with Plotter(filename="data/processed/metadata/sparsity.pdf", figsize=(12,12)) as ax:
    sns.heatmap(metadata.isna(), ax=ax, cbar=False)
    ax.set_title("Missing Values for COVID Samples Metadata (White = Missing)")

for column in ['age', 'num_side_effects', 'num_comorbidities', 'num_treatments', 'Days since last sampling']:
    for category in ['Lineage']:
        with Plotter(filename=f"data/processed/metadata/boxplots/{column}_vs_{category}.pdf", figsize=(10,10), show=False) as ax:
            sns.boxplot(category, column, data=metadata, ax=ax, palette=color_dict)
            sns.swarmplot(category, column, data=metadata, ax=ax, palette=color_dict)
            ax.set_title(f"{column} vs. {category}")

for col in ['sex', 'Lineage', 'outcome']: # Plot pie charts
    with Plotter(filename=f"data/processed/metadata/{col}_pie.pdf", figsize=(6,6)) as ax:
        metadata[col].value_counts().plot.pie(ax=ax, autopct='%1.1f%%', shadow=True, startangle=0, title=f"Distribution of {col}")

for col in ['age', 'num_side_effects', 'num_comorbidities', 'num_treatments']:
    with Plotter(filename=f"data/processed/metadata/{col}_box.pdf", figsize=(8,6)) as ax:
        sns.boxplot(metadata[col], ax=ax)
        sns.swarmplot(metadata[col], ax=ax, color='green')
        ax.set_title(f"Distribution of {col}.")
    with Plotter(filename=f"data/processed/metadata/{col}_hist.pdf", figsize=(8,6)) as ax:
        metadata[col].plot.hist(bins=30, alpha=.75)
        ax.set_title(f"Distribution of {col}.")

metadata2 = pd.read_csv("data/confidential/sample_metadata.csv")
metadata2.index = metadata2.number.map(lambda x: f"AECOM-{x}")
metadata2 = metadata2.loc[metadata2.index.intersection(passed)]

symptoms = set([y.strip() for x in metadata2['COVID sx'] for y in str(x).split(",")])
comorbidities = set([y.strip() for x in metadata2['Pt co-morbidities'] for y in str(x).split(",")])
treatments = set([y.strip() for x in metadata2['Tx'] for y in str(x).split(",")])

symptom_df = pd.DataFrame({symptom: metadata2['COVID sx'].map(lambda x: symptom in str(x)) for symptom in symptoms})
comorbidities_df = pd.DataFrame({como: metadata2['Pt co-morbidities'].map(lambda x: como in str(x)) for como in comorbidities})
treatments_df = pd.DataFrame({tx: metadata2['Tx'].map(lambda x: tx in str(x)) for tx in treatments})

sns.clustermap(symptom_df.T, yticklabels=True, xticklabels=True, figsize=(20,30))
plt.savefig("data/processed/metadata/symptoms_clustermap.pdf")
sns.clustermap(comorbidities_df.T, yticklabels=True, xticklabels=True, figsize=(20,30))
plt.savefig("data/processed/metadata/comorbidities_clustermap.pdf")
sns.clustermap(treatments_df.T, yticklabels=True, xticklabels=True, figsize=(20,30))
plt.savefig("data/processed/metadata/treatments_clustermap.pdf")

with Plotter(filename="data/processed/metadata/symptoms_bar.pdf", figsize=(30,30)) as ax:
    symptom_df.sum().sort_values().plot.barh(ax=ax)
    ax.set_title("Frequency of Symptoms in Patients")
    ax.set_xlabel("Number of Patients")
    ax.set_ylabel("Symptom")

with Plotter(filename="data/processed/metadata/comorbidities.pdf", figsize=(30,30)) as ax:
    comorbidities_df.sum().sort_values().plot.barh(ax=ax)
    ax.set_title("Frequency of Comorbidities in Patients")
    ax.set_xlabel("Number of Patients")
    ax.set_ylabel("Comorbidity")

with Plotter(filename="data/processed/metadata/treatments_bar.pdf", figsize=(30,30)) as ax:
    treatments_df.sum().sort_values().plot.barh(ax=ax)
    ax.set_title("Frequency of Treatments in Patients")
    ax.set_xlabel("Number of Patients")
    ax.set_ylabel("Treatment")

# Make it for most common
with Plotter(filename="data/processed/metadata/symptoms_bar_most_common.pdf", figsize=(15,10)) as ax:
    s = symptom_df.sum().sort_values()
    s = s[s>1]
    s = s.drop('unknown')
    s.plot.barh(ax=ax)
    ax.set_title("Frequency of Most Common Symptoms in Patients")
    ax.set_xlabel("Number of Patients")
    ax.set_ylabel("Symptom")

with Plotter(filename="data/processed/metadata/comorbidities_bar_most_common.pdf", figsize=(15,10)) as ax:
    c = comorbidities_df.sum().sort_values()
    c = c[c>1]
    c = c.rename(
        {
            'decompensated alcoholic cirrhosis': 'OH cirrhosis',
            })
    c = c.drop('unk')
    c.plot.barh(ax=ax)
    ax.set_title("Frequency of Most Common Comorbidities in Patients")
    ax.set_xlabel("Number of Patients")
    ax.set_ylabel("Comorbidity")

with Plotter(filename="data/processed/metadata/treatments_bar_most_common.pdf", figsize=(15,10)) as ax:
    t = treatments_df.sum().sort_values()
    t = t[t>2]
    t = t.rename({'none prior to admission': 'none'})
    t = t.drop('')
    t = t.drop('unk')
    t.plot.barh(ax=ax)
    ax.set_title("Frequency of Most Common Treatments in Patients")
    ax.set_xlabel("Number of Patients")
    ax.set_ylabel("Treatment")

# Dates
dates = pd.date_range(metadata['collection_date'].min(), metadata['collection_date'].max())
mdates = metadata[['collection_date', 'Lineage']]
mdates.index = mdates['collection_date']
mdates = mdates['Lineage']
date_df = pd.DataFrame(index=dates, columns=['A.1', 'B', 'B.1', 'B.1.3', 'B.2', 'B.2.1']).fillna(0.)
for l in ['A.1', 'B', 'B.1', 'B.1.3', 'B.2', 'B.2.1']:
    for d in mdates[mdates==l].index:
        date_df[l].loc[d] += 1 #mdates[mdates == l].index] += 1
date_df['B.2/B.2.1'] = date_df['B.2'] + date_df['B.2.1']
date_df = date_df.drop(columns=['B.2', 'B.2.1'])

lineages = metadata['Lineage'].value_counts()
cum_date_df = pd.DataFrame(index=dates, columns=['A.1', 'B', 'B.1', 'B.1.3', 'B.2/B.2.1']).fillna(0.)
for i in range(len(date_df)):
    cum_date_df.iloc[i] = date_df.iloc[0:i].sum()

with Plotter(filename=f"data/processed/metadata/lineage_pie.pdf", figsize=(6,6)) as ax:
    lineages.plot.pie(ax=ax, autopct='%1.1f%%', shadow=True, startangle=0, title=f"Distribution of Lineages", colors=[color_dict.get(x, '#333333') for x in cum_date_df.columns])

with Plotter(filename="data/processed/metadata/lineage_frequency.pdf", figsize=(15,10)) as ax:
    date_df.plot(
        ax=ax, 
        color=[color_dict.get(x, '#333333') for x in date_df.columns]
    )
    ax.set_title("Count of Samples by Date by Lineage")
    ax.set_xlabel("Date")
    ax.set_ylabel("Count")

with Plotter(filename="data/processed/metadata/lineage_cumulative.pdf", figsize=(15,10)) as ax:
    cum_date_df.plot(ax=ax,color=[color_dict.get(x, '#333333') for x in cum_date_df.columns])
    ax.set_title("Cumulative Samples by Date by Lineage")
    ax.set_xlabel("Date")
    ax.set_ylabel("Count")

symptom_df['Lineage'] = metadata['Lineage']
for col in ['fever','cough', 'SOB', 'diarrhea', 'nausea', 'hypoxia', 'weakness', 'chills']:
    with Plotter(filename=f"data/processed/metadata/pie/{col}.pdf") as ax:
        dx = symptom_df[symptom_df[col]]['Lineage'].value_counts()
        dx.plot.pie(
            ax=ax,
            title=f"Number of Patients with {col} by Lineage.",
            autopct='%1.1f%%',
            shadow=True,
            startangle=0,
            colors=[color_dict.get(x, '#333333') for x in dx.index]
        )

treatments_df['Lineage'] = metadata['Lineage']
for col in ['atorvastatin', 'aspirin', 'amlodipine', 'tamsulosin', 'albuterol', 'insulin', 'metoprolol', 'hydroxychloroquine', 'iron', 'pantoprazole']:
    with Plotter(filename=f"data/processed/metadata/pie/{col}.pdf") as ax:
        dx = treatments_df[treatments_df[col]]['Lineage'].value_counts()
        dx.plot.pie(
            ax=ax,
            title=f"Number of Patients with {col} by Lineage.",
            autopct='%1.1f%%',
            shadow=True,
            startangle=0,
            colors=[color_dict.get(x, '#333333') for x in dx.index]
        )

comorbidities_df['Lineage'] = metadata['Lineage']
for col in ['HTN', 'DM', 'HLD', 'CKD', 'asthma', 'CAD']:
    with Plotter(filename=f"data/processed/metadata/pie/{col}.pdf") as ax:
        dx = comorbidities_df[comorbidities_df[col]]['Lineage'].value_counts()
        dx.plot.pie(
            ax=ax,
            title=f"Number of Patients with {col} by Lineage.",
            autopct='%1.1f%%',
            shadow=True,
            startangle=0,
            colors=[color_dict.get(x, '#333333') for x in dx.index]
        )
