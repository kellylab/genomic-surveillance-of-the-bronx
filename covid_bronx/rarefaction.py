import pandas as pd
from covid_bronx.variants import parse_nextclade
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import plot_helper
from plot_helper import Plotter
import numpy as np

plot_helper.SHOW = False

def to_cdf(df: pd.DataFrame, dates_column: str = None, target_column: str = None):


    if type(df) is pd.Series:
        df = pd.DataFrame(df)
        df.index = pd.to_datetime(df.index, errors='coerce')
        df[dates_column] = df.index
        df = df.dropna()

    date_range = pd.date_range(df[dates_column].min(), df[dates_column].max())
    target = df[target_column]
    pdf = pd.DataFrame(index=date_range, columns = set(target)).fillna(0)
    
    counts = {k: v.value_counts(target_column) for k, v in df.groupby(dates_column)}
    for date, count in counts.items():
        pdf.loc[pd.to_datetime(date), count.index] += count
    
    cdf = pdf.cumsum()
    try:
        cdf = cdf.drop(columns=[np.nan])
    except:
        pass

    cdf = cdf[cdf.iloc[-1].sort_values(ascending=False).index]
    cdf = cdf[cdf.index > '2019-12-01']

    return cdf

def load_lineages(metadata_file: str, sep=None):
    """
    Given a GISAID metadata file, returns a cumulative density of lineages over time.
    """
    if sep:
        df = pd.read_csv(metadata_file, sep=sep)
    elif metadata_file.endswith('.tsv'):
        df = pd.read_csv(metadata_file, sep="\t")
    else:
        try:
            df = pd.read_csv(metadata_file)
        except UnicodeDecodeError:
            df = pd.read_csv(metadata_file, encoding='latin-1', index_col=0)

    if 'Lineage' in df:
        lineages = df['Lineage']
        target_column = 'Lineage'
    elif 'pangolin_lineage' in df:
        lineages = df['pangolin_lineage']
        target_column = 'pangolin_lineage'
    else:
        raise ValueError("Could not figure out which column in metadata file contains the lineages.")
    if 'Collection date' in df:
        lineages.index = df['Collection date']
        dates_column = 'Collection date'
    elif 'date' in df:
        lineages.index = df['date']
        dates_column = 'date'        
    elif 'Collection.date':
        lineages.index = pd.to_datetime(df['Collection.date'])
        dates_column = 'Collection.date'

    else:
        raise ValueError("Could not figure out which column in metadata file contains the dates info.")

    lineages = lineages.sort_index()

    date_range = pd.date_range(lineages.index.min(), lineages.index.max())
    if 'Virus name' in df:
        sample_dates = pd.Series({
            r['Virus name'].lstrip("hCoV-19/"): r[dates_column]
            for _,r in df.iterrows()
        })
    elif 'strain' in df:
        sample_dates = pd.Series({
            r['strain'].lstrip("hCoV-19/"): r[dates_column] 
            for _,r in df.iterrows()
        })
    elif 'Virus.name' in df:
        sample_dates = pd.Series({
            r['Virus.name'].lstrip("hCoV-19/"): r[dates_column] 
            for _,r in df.iterrows()
        })
    else:
        raise ValueError("Could not figure out which column in metadata file contains the strain names.")

    # sample_dates = df['Virus name'].apply(lambda x: x.lstrip("hCoV-19/"))
    # sample_dates.index = df['Collection date']
    lineages.index = lineages.index.rename('')
    cdf = to_cdf(lineages, dates_column=dates_column, target_column=target_column)

    return cdf, sample_dates


def load_nextclade(json_folder: str):
    """
    Loads GISAID variants jsons in json_folder into dataframes
    """
    
    meta = []
    variants = []

    if 'meta.csv' in os.listdir(json_folder) and 'variants.csv' in os.listdir(json_folder):
        da = pd.read_csv(os.path.join(json_folder, "meta.csv"), index_col=0)
        db = pd.read_csv(os.path.join(json_folder, "variants.csv"), index_col=0)
        da = da.drop_duplicates()
        db = db.drop_duplicates()

        return da, db

    for filename in tqdm(os.listdir(json_folder)):
        if filename.endswith(".json"):
            a,b = parse_nextclade(os.path.join(json_folder, filename))
            meta.append(a)
            variants.append(b)
    
    da = pd.concat(meta)
    db = pd.concat(variants)
    db['country'] = db['seqName'].apply(lambda x: x.split("/")[0])
    db.index = db['seqName'].apply(lambda x: x.lstrip("hCoV-19/").split("|")[0])

    da = da.drop_duplicates()
    db = db.drop_duplicates()

    da.to_csv(os.path.join(json_folder, 'meta.csv'))
    db.to_csv(os.path.join(json_folder, 'variants.csv'))

    return da, db

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

def plot_lineages(metadata_file: str, filename: str, title="", number_to_include=20, sep=None):

    cdf, dates = load_lineages(metadata_file, sep=sep)
    number_to_include = min(cdf.shape[1], number_to_include)

    reduced_cdf = cdf.loc[:, cdf.columns[:number_to_include]]

    cdf.to_csv(filename + "lineages.csv")
    reduced_cdf.to_csv(filename + "lineages_reduced.csv")

    with Plotter(filename=filename + '_lineages_rarefaction.pdf', figsize=(20,20)) as ax:
        reduced_cdf.plot(ax=ax, cmap='tab20b')
        ax.set_xlabel("Date")
        ax.set_ylabel("Count")
        ax.set_title(title + " Rarefaction of Most Common Lineages")

    ratios = (reduced_cdf.T / reduced_cdf.T.sum()).T

    with Plotter(filename=filename + '_lineages_cave.pdf', figsize=(20,20)) as ax:
        ax = cave_plot(ratios, ax, cmap='tab20b')
        plt.legend(loc='upper left')
        ax.set_title(title + " Relative Abundance of Lineages Over Time")
        ax.set_xlabel("Date")
        ax.set_ylabel("Relative Abundance")    

    
orfs = {
    "ORF1ab": [266, 21555],
    "Spike": [21563, 25384],
    "ORF3a": [25393, 26220],
    "Envelope": [26245, 26472],
    "Membrane Glycoprotein": [26523, 27191],
    "ORF6": [27202, 27387],
    "ORF7a": [27394, 27759],
    "ORF7b": [27756, 27887],
    "ORF8": [27894, 28259],
    "N": [28274, 29533],
    "ORF10": [29558, 29674],
}

def variants_cdf(metadata_file: str, json_folder: str, sep=None):
    """
    Loads both GISAID variants JSON and metadata and groups variants by date.
    """
    da, db = load_nextclade(json_folder)
    cdf, dates = load_lineages(metadata_file, sep=sep)

    db['Date'] = dates

    return da, db

def plot_variants(metadata_file: str, json_folder: str, filename: str, title="", number_to_include=20, sep=None):

    da, db = variants_cdf(metadata_file, json_folder, sep=sep)
    db = db.dropna(subset=['Date'])
    
    for protein, (start, end) in tqdm(orfs.items()):
        region = db[(db['position']>start)&(db['position']<=end)]
        region.index = region['Date']
        region = region.sort_index()
        region['Variant'] = region['mutation'] + " " +  region['aaMutation'].fillna('').apply(str)

        variants = region['Variant'].dropna()

        date_range = pd.date_range(region.index.min(), region.index.max())

        cdf = to_cdf(variants, dates_column='Dates', target_column='Variant')

        number_to_include = min(cdf.shape[1], number_to_include)
        reduced_cdf = cdf.loc[:, cdf.columns[:number_to_include]]

        cdf.to_csv(filename + f"{protein}_variants.csv")
        reduced_cdf.to_csv(filename + f"{protein}_variants_reduced.csv")
    
        with Plotter(filename=filename + f'{protein}_rarefaction.pdf', figsize=(20, 20), show=False) as ax:
            reduced_cdf.plot(ax=ax, cmap='tab20b')
            ax.set_xlabel("Date")
            ax.set_ylabel("Count")
            ax.set_title(title + f" Rarefaction of Most Common {protein} Variants")

        ratios = (reduced_cdf.T / reduced_cdf.T.sum()).T

        with Plotter(filename=filename + f'{protein}_cave.pdf', figsize=(20,20), show=False) as ax:
            ax = cave_plot(ratios, ax, cmap='tab20b')
            plt.legend(loc='upper left')
            ax.set_title(title + f" Ratios of Variants in {protein} Over Time")
            ax.set_xlabel("Date")
            ax.set_ylabel("Relative Abundance")

    