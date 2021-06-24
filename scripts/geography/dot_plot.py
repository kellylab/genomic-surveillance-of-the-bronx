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
column = 'count'
legend = 'Number of Samples'

folium.Choropleth(
    bronx_zip_codes,
    zip_df,
    ('zip_code', column),
    key_on='feature.properties.ZIPCODE',
    fill_color="Blues",
    fill_opacity=1,
    line_opacity=1,
    legend_name=legend,
    smooth_factor=0
).add_to(m)

folium.GeoJson(
    bronx_zip_codes,
    name='NYC Zip Codes'
).add_to(m)

g = gpd.GeoSeries([Polygon(x['geometry']['coordinates'][0]) for x in bronx_zip_codes['features']])
for lineage, color in {
    "A.1":'red',
    "B":'blue',
    "B.1":'green',
    "B.1.3":'cyan',
    "B.2":'magenta',
    "B.2.1":'yellow',
}.items():
    vals = zip_df[lineage]
    random_points = gen_points_in_gdf_polys(geometry=g, values=vals*40)
    for point in random_points:
        folium.CircleMarker(location=[point.coords[0][1], point.coords[0][0]], color=color, radius=8, fill=True, fill_color=color).add_to(m)

folium.LayerControl().add_to(m)

m.save("data/processed/geography/bronx_dotplot_40.html")
