# This breaks out the sample_metadata.csv file in the metadata folder into
# individual csv files which can be used as input to alignmentviewer
# https://github.com/sanderlab/alignmentviewer

import pandas as pd

df  = pd.read_csv("data/confidential/sample_metadata.csv")

filenames = {
    'Collection Date': 'collection_date', 
    'Location': 'location',
    'Sex': 'sex',
    'Age': 'age',
    'zipcode': 'zipcode',
    'travel/sickcontact': 'sick_contact',
    'Pt status at time of collection': 'status_upon_collection',
    'Pt Outcome': 'outcome',
    'Specimen source': 'source',
    'Pt co-morbidities': 'num_comorbidities', 
    'Tx': "num_treatments", 
    'COVID sx': "num_side_effects",
}

df['id'] = df['number'].apply(lambda x: f'AECOM-{x}')
df['Sex'] = df['Sex'].map({'F': 'female', 'f': 'female', 'M': 'male', 'm': 'male'})
df['Pt Outcome'] = df['Pt Outcome'].apply(lambda x: 'deceased' if 'deceased' in str(x).lower() else ('unknown' if 'unk' in str(x).lower() or len(str(x)) < 2 else 'alive'))

for column in [
    'Collection Date', 'Location', 'Sex', 'Age', 'zipcode', 'travel/sickcontact', 
    'Pt status at time of collection', 'Pt Outcome', 'Specimen source',
    ]:
    meta_df = df[['id', column]]
    meta_df = meta_df.rename(columns = {'id': 'name', column: 'annotation'})
    meta_df.to_csv(f"data/processed/metadata/{filenames[column]}.tsv", sep="\t", index=False)

# Separately break out categories for 'Pt co-morbidities', 'Tx', 'Covid sx', 'Pt
# co-morbidities'
for column in ['Pt co-morbidities', 'Tx', 'COVID sx']:
    meta_df = df[['id', column]]
    meta_df = meta_df.rename(columns = {'id': 'name', column: 'annotation'})
    meta_df['annotation'] = meta_df['annotation'].apply(lambda x: len(str(x).split(",")))
    meta_df.to_csv(f"data/processed/metadata/{filenames[column]}.tsv", sep="\t", index=False)


date = pd.to_datetime(df['Collection Date'])
deltas = date - min(date)
days = deltas.apply(lambda x: x.days)
meta_df['annotation'] = days
meta_df.to_csv("data/processed/metadata/days.tsv", sep="\t", index=False)

# Make a version of the dates csv for treetime
meta_df = df[['Collection Date', 'id']]
meta_df = meta_df.rename(columns={'Collection Date': 'date', 'id': 'name'})
meta_df['date'] = pd.to_datetime(meta_df['date'])
meta_df['date'] = meta_df['date'].apply(lambda x: f"{x.year}.{x.month}.{x.day}")
meta_df.to_csv("data/processed/metadata/treetime_dates.csv", index=False)
