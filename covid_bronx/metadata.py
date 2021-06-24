import pandas as pd

def preprocess_metadata():

    collection_date = pd.read_csv("data/processed/metadata/collection_date.tsv", sep="\t")
    index = collection_date['name']
    collection_date.index = index
    age = pd.read_csv("data/processed/metadata/age.tsv", sep="\t")
    age.index = age['name']
    sex = pd.read_csv("data/processed/metadata/sex.tsv", sep="\t")
    sex.index = sex['name']
    num_comorbidities = pd.read_csv("data/processed/metadata/num_comorbidities.tsv", sep="\t")
    num_comorbidities.index = num_comorbidities['name']
    num_side_effects = pd.read_csv("data/processed/metadata/num_side_effects.tsv", sep="\t")
    num_side_effects.index = num_side_effects['name']
    num_treatments = pd.read_csv("data/processed/metadata/num_treatments.tsv", sep="\t")
    num_treatments.index = num_treatments['name']
    zip_code = pd.read_csv("data/processed/metadata/zipcode.tsv", sep="\t")
    zip_code.index = zip_code['name']
    outcome = pd.read_csv("data/processed/metadata/outcome.tsv", sep="\t")
    outcome.index = outcome['name']

    metadata = pd.DataFrame(index=collection_date['name'])
    metadata['sex'] = sex['annotation']
    metadata['age'] = age['annotation']
    metadata['num_side_effects'] = num_side_effects['annotation']
    metadata['num_comorbidities'] = num_comorbidities['annotation']
    metadata['num_treatments'] = num_treatments['annotation']
    metadata['zip_code'] = zip_code['annotation']
    metadata['outcome'] = outcome['annotation']
    metadata['collection_date'] = pd.to_datetime(collection_date['annotation'])

    metadata.to_csv("data/processed/metadata.csv")

    return metadata

def get_metadata():

    df = pd.read_csv("data/processed/metadata.csv")
    df['collection_date'] = pd.to_datetime(df['collection_date'])
    
    return df