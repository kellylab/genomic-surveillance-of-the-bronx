""" 
Makes a figure providing an overview of our dataset with a focus on lineages
laid out as follows:

a - Patient metadata
b - Donut plot of our lineage distributions vs the world
c - Timeline of patient sampling vs lineages identified
d - Choropleth of lineages by region
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from typing import Dict

import logging
import matplotlib
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

from covid_bronx import lineage_colors_dict, lineage_colors_dict_rgb
from covid_bronx.quality import fasta_files, sam_files, variant_files

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

savefile_a = "figures_final/figure1a"
savefile_b = "figures_final/figure1b"

months = {
    1: 'Jan',
    2: 'Feb',
    3: 'Mar',
    4: 'Apr',
    5: 'May',
    6: 'Jun',
    7: 'Jul',
    8: 'Aug',
    9: 'Sep',
    10: 'Oct',
    11: 'Nov',
    12: 'Dec',
}

from covid_bronx.metadata import preprocess_metadata
from matplotlib.colors import colorConverter

# a) Timeline of lineages
logger.info("Plotting 1a")
timeline = pd.read_csv("data/external/global_lineages.csv")
from covid_bronx.metadata import get_metadata
metadata = get_metadata()
index = pd.date_range(metadata['collection_date'].min(), metadata['collection_date'].max())
metadata.index = metadata['name']
df = pd.read_csv("data/external/pangolin2.csv")
df.index = df['Sequence name'].apply(lambda x: x.split(" ")[0])
df.index = df.index.map(lambda x: "AECOM-" + str(int(x.split("-")[1])))

metadata[df.columns] = df
lineages_df = pd.read_csv("data/external/Lineages_updated.csv", index_col=0)
lineages = lineages_df['lineage'].dropna()
lineages.index = lineages.index.map(lambda x: x.replace("_", "-"))
metadata['Lineage'] = lineages
metadata = pd.concat([metadata,metadata.loc[['AECOM-126','AECOM-127','AECOM-128','AECOM-129','AECOM-130']]]).drop_duplicates(keep=False)

ddf = pd.DataFrame([ # Incremental Values
    {
        l: (metadata[metadata['collection_date'] == d]['Lineage']==l).sum()
        for l in lineages
    }
    for d in index
    ],
    index=index
)

ddf.index = ddf.index.map(lambda x: months[x.month])
ddmf = pd.DataFrame({k: v.sum(0) for k,v in ddf.groupby(ddf.index)})

cdf = pd.DataFrame([ # Cumulative Values
    {
        l: (metadata[metadata['collection_date'] <= d]['Lineage']==l).sum()
        for l in lineages
    }
    for d in index
    ],
    index=index
)

dd = pd.read_csv("data/external/data-by-day.csv", index_col=0)
dd.index = pd.to_datetime(dd.index)
dd['month'] = dd.index.map(lambda x: months[x.month])

bronx_sampling = ddmf.sum(0)
sampling = pd.read_csv("data/external/sampling.csv", index_col=0) # TODO: Verify this
sampling['date'] = pd.to_datetime(sampling['date'])
sampling['month'] = sampling['date'].apply(lambda x: months[x.month])
# deathsdmf = pd.Series({k:v['Deaths'].sum() for k,v in sampling.groupby('month')})
# casesdmf = pd.Series({k:v['Cases'].sum() for k,v in sampling.groupby('month')})
# hospitalizationdmf = pd.Series({k:v['Hospitalizations'].sum() for k,v in sampling.groupby('month')})

deathsdmf = pd.Series({k:v['DEATH_COUNT'].sum() for k,v in dd.groupby('month')})
casesdmf = pd.Series({k:v['CASE_COUNT'].sum() for k,v in dd.groupby('month')})
hospitalizationdmf = pd.Series({k:v['HOSPITALIZED_COUNT'].sum() for k,v in dd.groupby('month')})

sampling_df = pd.DataFrame({"Sampling": bronx_sampling, "Cases": casesdmf, "Deaths": deathsdmf, "Hospitalizations": hospitalizationdmf}).fillna(0.)

##########################################################

# Start Plotting 
matplotlib.rcParams.update({'font.size': 16})
plt.clf()
plt.close()
fig1a = plt.figure(figsize=(24,24))

from covid_bronx.geography import gen_points_in_gdf_polys, blank_background_choropleth, get_zip_codes_metadata_geo
import geopandas as gpd

metadata = preprocess_metadata()
coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
passed = coverage_levels[coverage_levels>=.95].index.intersection(sam_files.keys()).map(lambda x: x.replace("_","-").replace("-0", "-").replace("-0", "-"))
num_samples = len(passed)

from covid_bronx.metadata import preprocess_metadata
from matplotlib import colors

def colorizer(df: pd.DataFrame, color_dict: Dict) -> pd.Series:
    """
    Given a dataframe where the rows are zip codes and columns are lineages,
    along with a dict explaining what the RGB color values are, returns a series
    linking zip codes to a color output.
    """
    scale_factors = df.sum(1) / max(df.sum(1))
    weights = (df.T / df.sum(1))
    color_series = pd.DataFrame( [np.sum(weights[z][c]*v for c,v in color_dict.items()) for z in weights.columns], index=weights.columns, columns=['r','g','b'])

    return color_series.T

df = pd.read_csv("data/external/lineages_final.csv", index_col=0)
df.index = df['taxon'].apply(lambda x: x.split(" ")[0])
metadata[df.columns] = df

zips = metadata.loc[metadata.index.intersection(passed)]['zip_code'].to_numpy()
zips = np.array(sorted(zips)[2:])

# Get a listing of coordinates by zip code
bronx_zip_codes = [10453, 10457, 10460,	10458, 10467, 10468,10451, 10452, 10456,10454, 10455, 10459, 10474,	10463, 10471,10466, 10469, 10470, 10475,10461, 10462,10464, 10465, 10472, 10473]
gdf = gpd.read_file("data/external/ZIP_CODE_040114/ZIP_CODE_040114.geojson")
gdf.index = gdf['ZIPCODE']
gdf = gdf.loc[list(map(str, bronx_zip_codes))]
# Remove extraneous zip codes

latlons = gpd.GeoDataFrame({"ZIPCODE": gdf['ZIPCODE'], 'geometry': gdf['geometry'].centroid}).set_index("ZIPCODE")

zdf, bzip = get_zip_codes_metadata_geo()
zdf.index = zdf['zip_code'].map(lambda x: str(int(float(x))))
gdf[zdf.columns] = zdf
gdf = gdf.fillna(0.)

geocolor_dict = {k: lineage_colors_dict_rgb[k] for k in ['B.1', 'B.1.3', 'B.1.1']} # {'B.1': np.array([1,0,0]), 'B.1.3': np.array([0,1,0]), 'B.1.1': np.array([0,0,1])}
lineage_colors = colorizer(zdf[['B.1', 'B.1.3', 'B.1.1']], geocolor_dict).to_numpy()
lineage_colors = np.nan_to_num(lineage_colors, 0.)
gdf['lineage_colors'] = pd.Series([colors.to_rgba(lineage_colors[:,i]/256) for i in range(len(lineage_colors.T))], index=zdf.index)
gdf['lineage_colors'] = gdf['lineage_colors'].fillna('#000000')

fig, ax = plt.subplots()
gdf.fillna(0.).plot(column='count', cmap='Purples',ax=ax, legend=True, legend_kwds={'shrink': 0.3})
gdf.boundary.plot(color='black', ax=ax)
ax.set_axis_off()

# Plot hospital locations
from shapely.geometry import Point

hospitals = [Point(-73.846184,40.849010)]
hospitals_df = gpd.GeoDataFrame(geometry=hospitals)
# hospitals_df.plot(ax=ax, markersize=500, color='black', marker='.', label="Collection Site") # We decided not to do this

plt.tight_layout(pad=.3)

plt.savefig(savefile_a + '.pdf')
plt.savefig(savefile_a + '.svg')
plt.clf()

# Plot lineage colored distribution

geocolor_dict = {k: lineage_colors_dict_rgb[k] for k in ['B.1', 'B.1.3', 'B.1.1']} # {'B.1': np.array([1,0,0]), 'B.1.3': np.array([0,1,0]), 'B.1.1': np.array([0,0,1])}
lineage_colors = colorizer(zdf[['B.1', 'B.1.3', 'B.1.1']], geocolor_dict).to_numpy()
lineage_colors = np.nan_to_num(lineage_colors, 0.)
gdf['lineage_colors'] = pd.Series([colors.to_rgba(lineage_colors[:,i]/256) for i in range(len(lineage_colors.T))], index=zdf.index)
gdf['lineage_colors'] = gdf['lineage_colors'].fillna('#000000')
fig, ax = plt.subplots()
gdf.plot(ax=ax, color=gdf['lineage_colors'])
gdf.boundary.plot(color='black', ax=ax)
ax.set_axis_off()
plt.savefig("figures_final/figure1a_lineage.pdf")
plt.savefig("figures_final/figure1a_lineage.svg")
plt.show()
plt.clf()

# Figure 1b. Sampling Density

fig, ax = plt.subplots(figsize=(15,10))
ax_2 = ax.twinx()
sampling_df[['Cases', 'Hospitalizations', 'Deaths']].loc[['Feb','Mar','Apr','May','Jun','Jul','Aug','Sep', 'Oct']].plot(ax=ax, label=True, color=['yellowgreen','orange','red'], linewidth=6)
ax.grid(linestyle='--', linewidth=1)
ax.set_ylim([0,115000])
ax_2.set_ylim([0,80])
ax.set_ylabel("Count of Cases / Hospitalizations / Deaths")
ax.legend()
ax_2.set_ylabel("Count of Genomes Sequenced")
ax.set_xlabel("Month")
ax.set_xticklabels(['Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct'])
sampling_df['Sampling'][['Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct']].plot.bar(ax=ax_2, alpha=.5)
ax_2.grid(linestyle='--', color='blue', alpha=.5, linewidth=1)
ax_2.spines['right'].set_color('blue')
ax_2.yaxis.label.set_color('blue')
ax_2.tick_params(axis='y', colors='blue')

plt.savefig(savefile_b + '.pdf')
plt.savefig(savefile_b + '.svg')
plt.show()
plt.clf()

# Figure 2b. Lineages Over Time
fig, ax = plt.subplots(figsize=(30,15))
cdf_colors = [lineage_colors_dict[k] for k in ['B.1.26', 'B.1', 'B.2', 'B.2.1', 'A.1', 'B.1.3', 'B.1.1.1', 'B.1.1']]
cdf[['A.1', 'B.1',  'B.1.1', 'B.1.1.1', 'B.1.26', 'B.1.3', 'B.2', 'B.2.1',]].plot.line(legend=True, color=cdf_colors, ax=ax, linewidth=6)
ax.set_ylabel("Cumulative Sample Counts by Lineage")

plt.savefig('figures_final/figure2a' + '.pdf')
plt.savefig('figures_final/figure2a' + '.svg')
plt.clf()

# b) Donut Plot showing lineage distributions in world, US, NYS, and Bronx
# ax_q = fig1.add_subplot(gs[0:7, 13:])
import matplotlib
matplotlib.rcParams.update({'font.size':24})
fig, ax_q = plt.subplots(figsize=(30,30))
facecolor = colorConverter.to_rgba('white', alpha=0)
circulo = lambda r: plt.Circle((0,0), r, ec='white', fc=facecolor, lw=2)
logger.info("Plotting 1b")
donut = pd.read_csv("data/external/Donut_churro_plot.csv", index_col=0)
donut_colors = [lineage_colors_dict[k] for k in donut.index]
artist = donut['world'].plot.pie(radius=1, ax=ax_q, colors=donut_colors)
circle_1 = circulo(.8)
ax_q.add_artist(circle_1)
donut['USA'].plot.pie(radius=.8, ax=ax_q, labels=None, colors=donut_colors)
circle_1a = circulo(.6)
ax_q.add_artist(circle_1a)
donut['NYS'].plot.pie(radius=.6, ax=ax_q, labels=None, colors=donut_colors)
circle_2 = circulo(.4)
ax_q.add_artist(circle_2)
donut['Bronx'].plot.pie(radius=.4, ax=ax_q, labels=None, colors=donut_colors)
circle_3 = circulo(.2)
circle_4 = plt.Circle((0,0), .2, color='white')
ax_q.add_artist(circle_3)
ax_q.add_artist(circle_4)
ax_q.set_ylabel('')
plt.savefig("figures_final/figure2b.pdf")
plt.savefig("figures_final/figure2b.svg")
plt.show()

# Plot a triangular legend
fig, ax = plt.subplots()
x = np.array([-1,0])
y = np.array([1,0])
z = np.array([0,1])
x_c = geocolor_dict['B.1']/256
y_c = geocolor_dict['B.1.3']/256
z_c = geocolor_dict['B.1.1']/256

# Do convex combinations of everything
coordinates = []
k = 100
for lambd in np.linspace(0,1,k):
    for mu in np.linspace(0, 1-lambd, int(k*(1-lambd))):
        for w in np.linspace(0, 1-lambd-mu, int(k*(1-mu))):
            combo = lambd*x + mu*y + w*z
            color = colors.to_hex(max(lambd,0)*x_c + max(mu,0)*y_c + max(w,0)*z_c)
            coordinates.append([combo[0], combo[1], color])

coordinates = np.array(coordinates)
xy = coordinates[:, 0:2].astype(float)
ax.scatter(xy[:,0],xy[:,1], c=coordinates[:,2])
ax.text(-1.4,-.1, 'B.1')
ax.text(1.05,-.1, 'B.1.3')
ax.text(-.25,1.1, 'B.1.1')
ax.set_axis_off()
plt.savefig("figures_final/figure1a_lineage_legend.pdf")
plt.savefig("figures_final/figure1a_lineage_legend.svg")
plt.show()
plt.clf()
