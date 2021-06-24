import folium
import pandas as pd
import json
from covid_bronx.geography import gen_points_in_gdf_polys, blank_background_choropleth, get_zip_codes_metadata_geo
from folium.plugins import MarkerCluster
import geopandas as gpd
from shapely.geometry import Polygon
from tqdm import tqdm

zip_df, bronx_zip_codes = get_zip_codes_metadata_geo()

m = folium.Map(
    location=[40.8448, -73.8648],
    tiles='Stamen Toner',
)

blank_background_choropleth(m)

for lineage, color in {
    "B":'Blues',
    "B.1":'Greens',
    "B.1.3":'Reds',
}.items():

    folium.Choropleth(
        bronx_zip_codes,
        zip_df,
        ('zip_code', f"freq_{lineage}"),
        key_on='feature.properties.ZIPCODE',
        fill_color=color,
        fill_opacity=1,
        line_opacity=1,
        legend_name=f"Frequency of {lineage}.",
        smooth_factor=0,
        name=lineage,
    ).add_to(m)

folium.GeoJson(
    bronx_zip_codes,
    name='NYC Zip Codes'
).add_to(m)

folium.LayerControl().add_to(m)

m.save(f"data/processed/geography/bronx_lineages_shaded.html")
