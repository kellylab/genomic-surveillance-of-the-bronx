import folium
import pandas as pd
import json
from covid_bronx.geography import gen_points_in_gdf_polys, blank_background_choropleth, get_zip_codes_metadata_geo
from folium.plugins import MarkerCluster
import geopandas as gpd
from shapely.geometry import Polygon
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

zip_df, bronx_zip_codes = get_zip_codes_metadata_geo()

for column, legend in tqdm({
    "age": "Mean Age",
    "male": "Percent Male",
    "count": "Number of Samples",
    "num_side_effects": "Number of Side Effects",
    "num_comorbidities": "Number of Comorbidities",
    "death_rate": "Death Rate",
    "freq_A.1": "Percent A.1 Lineage",
    "freq_B": "Percent B Lineage",
    "freq_B.1": "Percent B.1 Lineage",
    "freq_B.1.26": "Percent B.1.26 Lineage",
    "freq_B.1.3": "Percent B.1.3 Lineage",
    "freq_B.2": "Percent B.2 Lineage",
    "freq_B.2.1": "Percent B.2.1 Lineage",
}.items()):

    m = folium.Map(
        location=[40.8448, -73.8648],
        tiles='Stamen Toner',
    )

    blank_background_choropleth(m)

    folium.Choropleth(
        bronx_zip_codes,
        zip_df,
        ('zip_code', column),
        key_on='feature.properties.ZIPCODE',
        fill_color="Blues",
        fill_opacity=1,
        line_opacity=1,
        legend_name=legend,
        smooth_factor=0,
    ).add_to(m)

    folium.GeoJson(
        bronx_zip_codes,
        name='NYC Zip Codes'
    ).add_to(m)

    folium.LayerControl().add_to(m)

    m.save(f"data/processed/geography/bronx_{column}.html")

# ax = g.plot()
# random_points.plot(color='black', ax=ax, alpha=.7)
# plt.show()

df = pd.read_csv("data/processed/metadata.csv")
lineages = set(df['Lineage'])

dg = pd.DataFrame({l: (df['Lineage']==l).to_numpy().astype(int) for l in lineages}, index=df['Unnamed: 0'])
dg['zip_code'] = df['zip_code'].to_numpy()
dg2 = pd.DataFrame({k: v.mean() for k,v in dg.groupby('zip_code')}).T[lineages]

sns.clustermap(dg2, linewidth=.75)
plt.savefig("data/processed/geography/zip_code_clustermap.png")
plt.clf()

dg3 = pd.DataFrame({k: v.sum() for k,v in dg.groupby('zip_code')}).T[lineages]
sns.clustermap(dg3, linewidth=.75)
plt.savefig("data/processed/geography/zip_code_sum_clustermap.png")
plt.clf()

multiple = dg3[dg3.sum(1)>1]

maxima = {k: v.idxmax() for k,v in multiple.iterrows()}
maxima_df = pd.DataFrame(
    index = maxima.keys(), 
    columns=['zip_code', 'A.1', 'B', 'B.1', 'B.1.3', 'B.2', 'B.2.1'],
    ).fillna(0.)
maxima_df['zip_code'] = maxima_df.index
for k,v in maxima.items():
    maxima_df.loc[k,v] = 1

m = folium.Map(
    location=[40.8448, -73.8648],
    tiles='Stamen Toner',
)

blank_background_choropleth(m)

for column, color in {
    'A.1': 'Reds',
    'B.1': 'Blues',
    'B.1.3': 'Greens',
    # 'B.1.26': '#a1184f',
    #'B.2/B.2.1': '#f7d93d',
    'B.2': 'RdGy',
    'B.2.1': 'RdGy',
    }.items():
    folium.Choropleth(
        bronx_zip_codes,
        zip_df,
        ('zip_code', column),
        key_on='feature.properties.ZIPCODE',
        fill_color=color,
        fill_opacity=1,
        line_opacity=1,
        legend_name=column,
        name=column,
        smooth_factor=0,
    ).add_to(m)

folium.GeoJson(
    bronx_zip_codes,
    name='NYC Zip Codes'
).add_to(m)

folium.LayerControl().add_to(m)

m.save(f"data/processed/geography/bronx_lineage_points.html")
