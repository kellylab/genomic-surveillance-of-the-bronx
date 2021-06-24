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

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

savefile = "figures/figure1_v2"


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

bronx_sampling = ddmf.sum(0)
sampling = pd.read_csv("data/external/sampling.csv", index_col=0)
sampling['date'] = pd.to_datetime(sampling['date'])
sampling['month'] = sampling['date'].apply(lambda x: months[x.month])
deathsdmf = pd.Series({k:v['Deaths'].sum() for k,v in sampling.groupby('month')})
casesdmf = pd.Series({k:v['Cases'].sum() for k,v in sampling.groupby('month')})
hospitalizationdmf = pd.Series({k:v['Hospitalizations'].sum() for k,v in sampling.groupby('month')})

sampling_df = pd.DataFrame({"Sampling": bronx_sampling, "Cases": casesdmf, "Deaths": deathsdmf, "Hospitalizations": hospitalizationdmf}).fillna(0.)

##########################################################

# Start Plotting 
matplotlib.rcParams.update({'font.size': 16})
plt.clf()
plt.close()
fig1 = plt.figure(figsize=(24,24))
gs = fig1.add_gridspec(20,20)


# a) Sampling Timeline

ax_c = fig1.add_subplot(gs[0:8, 10:])
ax_c2 = ax_c.twinx()
sampling_df[['Cases', 'Deaths', 'Hospitalizations']].loc[['Feb','Mar','Apr','May','Jun','Jul','Aug','Sep']].plot(ax=ax_c, label=True, color=['yellowgreen','red','orange'], linewidth=6)
ax_c.grid(linestyle='--', linewidth=1)
ax_c.set_ylim([0,100000])
ax_c2.set_ylim([0,80])
ax_c.set_ylabel("Count of Cases / Hospitalizations / Deaths")
ax_c.legend()
ax_c2.set_ylabel("Count of Genomes Sequenced")
ax_c.set_xlabel("Month")
ax_c.set_xticklabels(['Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'])
sampling_df['Sampling'][['Feb','Mar','Apr','May','Jun','Jul','Aug','Sep']].plot.bar(ax=ax_c2, alpha=.5)
ax_c2.grid(linestyle='--', color='blue', alpha=.5, linewidth=1)
ax_c2.spines['right'].set_color('blue')
ax_c2.yaxis.label.set_color('blue')
ax_c2.tick_params(axis='y', colors='blue')

# d) Choropleth by Lineage
logger.info("Plotting 1d")

# ax_d = fig1.add_subplot(gs[6:, 8:])

from covid_bronx.geography import gen_points_in_gdf_polys, blank_background_choropleth, get_zip_codes_metadata_geo
import geopandas as gpd

metadata = preprocess_metadata()
coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
passed = coverage_levels[coverage_levels>=.95].index.map(lambda x: x.replace("_","-").replace("-0", "-").replace("-0", "-"))
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

# clades = pd.read_csv("data/external/nextclade.csv", sep=";")
# clades.index = clades['seqName'].apply(lambda x: x.split(" ")[0])
# metadata[clades.columns] = clades

zips = metadata.loc[metadata.index.intersection(passed)]['zip_code'].to_numpy()
zips = np.array(sorted(zips)[2:])

# Get a listing of coordinates by zip code

gdf = gpd.read_file("data/external/ZIP_CODE_040114/ZIP_CODE_040114.geojson")
# Remove extraneous zip codes
gdf = gdf[~(gdf['PO_NAME']=='Staten Island')]

latlons = gpd.GeoDataFrame({"ZIPCODE": gdf['ZIPCODE'], 'geometry': gdf['geometry'].centroid}).set_index("ZIPCODE")

zdf, bzip = get_zip_codes_metadata_geo()
zdf.index = zdf['zip_code'].map(lambda x: str(int(float(x))))
gdf.index = gdf['ZIPCODE']
gdf[zdf.columns] = zdf
gdf = gdf.fillna(0.)
# gdf.plot(ax=ax_d)

geocolor_dict = {k: lineage_colors_dict_rgb[k] for k in ['B.1', 'B.1.3', 'B.1.1']} # {'B.1': np.array([1,0,0]), 'B.1.3': np.array([0,1,0]), 'B.1.1': np.array([0,0,1])}
lineage_colors = colorizer(zdf[['B.1', 'B.1.3', 'B.1.1']], geocolor_dict).to_numpy()
lineage_colors = np.nan_to_num(lineage_colors, 0.)
gdf['lineage_colors'] = pd.Series([colors.to_rgba(lineage_colors[:,i]/256) for i in range(len(lineage_colors.T))], index=zdf.index)
gdf['lineage_colors'] = gdf['lineage_colors'].fillna('#000000')
# gdf.plot(ax=ax_d, color=gdf['lineage_colors'])


# ax_d.set_axis_off()
# # Plot a triangular legend
# x = np.array([-1,0])
# y = np.array([1,0])
# z = np.array([0,1])
# x_c = geocolor_dict['B.1']/256
# y_c = geocolor_dict['B.1.3']/256
# z_c = geocolor_dict['B.1.1']/256
# 
# # Do convex combinations of everything
# coordinates = []
# k = 30
# for lambd in np.linspace(0,1,k):
#     for mu in np.linspace(0, 1-lambd, int(k*(1-lambd))):
#         for w in np.linspace(0, 1-lambd-mu, int(k*(1-mu))):
#             combo = lambd*x + mu*y + w*z
#             color = colors.to_hex(max(lambd,0)*x_c + max(mu,0)*y_c + max(w,0)*z_c)
#             coordinates.append([combo[0], combo[1], color])
# 
# coordinates = np.array(coordinates)
# xy = coordinates[:, 0:2].astype(float)
# ax2 = plt.axes([0,0,1,1])
# ip = InsetPosition(ax_d, [-.15,.7,0.3,0.2])
# ax2.set_axes_locator(ip)
# ax2.scatter(xy[:,0],xy[:,1], c=coordinates[:,2])
# ax2.text(-1.4,-.1, 'B.1')
# ax2.text(1.05,-.1, 'B.1.3')
# ax2.text(-.25,1.1, 'B.1.1')
# ax2.set_axis_off()

ax_3 = fig1.add_subplot(gs[0:10, 0:10])
gdf.fillna(0.).plot(column='count', cmap='pink', ax=ax_3, legend=True, legend_kwds={'shrink': 0.3})
ax_3.set_axis_off()

# Plot hospital locations
from shapely.geometry import Point

hospitals = [Point(-73.846184,40.849010)]
hospitals_df = gpd.GeoDataFrame(geometry=hospitals)
hospitals_df.plot(ax=ax_3, markersize=500, color='black', marker='.', label="Collection Site")

plt.tight_layout(pad=.3)

plt.savefig(savefile + '.pdf')
plt.savefig(savefile + '.svg')

fig, ax = plt.subplots(figsize=(20,10))
cdf_colors = [lineage_colors_dict[k] for k in ['B.1.26', 'B.1', 'B.2', 'B.2.1', 'A.1', 'B.1.3', 'B.1.1.1', 'B.1.1']]
cdf[['B.1.26', 'B.1', 'B.2', 'B.2.1', 'A.1', 'B.1.3', 'B.1.1.1', 'B.1.1']].plot.line(legend=True, color=cdf_colors, ax=ax)
ax.set_ylabel("Cumulative Sample Counts by Lineage")
plt.savefig("figures/figure2b_v2.pdf")
plt.savefig("figures/figure2b_v2.svg")

# b) Donut Plot showing lineage distributions in world, US, NYS, and Bronx
plt.clf()
# ax_q = fig1.add_subplot(gs[0:7, 13:])
fig, ax_q = plt.subplots(figsize=(15,15))
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
plt.savefig("figures/figure2c_v2.pdf")
plt.savefig("figures/figure2c_v2.svg")

