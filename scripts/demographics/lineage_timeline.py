import pandas as pd
import numpy as np
import plot_helper
from plot_helper import Plotter
import seaborn as sns
from covid_bronx.metadata import get_metadata

metadata = get_metadata()
index = pd.date_range(metadata['collection_date'].min(), metadata['collection_date'].max())
lineages = set(metadata['Lineage'])
ddf = pd.DataFrame([ # Incremental Values
    {
        l: (metadata[metadata['collection_date'] == d]['Lineage']==l).sum()
        for l in lineages
    }
    for d in index
    ],
    index=index
)

cdf = pd.DataFrame([ # Cumulative Values
    {
        l: (metadata[metadata['collection_date'] <= d]['Lineage']==l).sum()
        for l in lineages
    }
    for d in index
    ],
    index=index
)


with Plotter(filename="data/processed/demographics/incremental_lineage.png", figsize=(10,10)) as ax:
    ddf.plot.line(ax=ax)
    ax.set_title("Lineages per Sample Per Day")

with Plotter(filename="data/processed/demographics/cumulative_lineage.png", figsize=(10,10)) as ax:
    cdf.plot.line(ax=ax)
    ax.set_title("Cumulative Distribution of Lineages Over Time")