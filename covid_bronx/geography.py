# Credit to Andrew Gaidus
# These functions were taken from a script for analyzing census data
# http://andrewgaidus.com/Dot_Density_County_Maps/
# https://github.com/agaidus/census_data_extraction/blob/master/census_mapper.py

import geopandas as gpd
from shapely.geometry import shape  
import pandas as pd
from shapely.geometry import Point
from numpy.random import RandomState, uniform
import numpy as np
import json
import folium

def gen_random_points_poly(poly, num_points, seed = None):
    """
    Returns a list of N randomly generated points within a polygon. 
    """
    min_x, min_y, max_x, max_y = poly.bounds
    points = []
    i=0
    while len(points) < num_points:
        s=RandomState(seed+i) if seed else RandomState(seed)
        random_point = Point([s.uniform(min_x, max_x), s.uniform(min_y, max_y)])
        if random_point.within(poly):
            points.append(random_point)
        i+=1
    return points


def gen_points_in_gdf_polys(geometry, values, points_per_value = None, seed = None):
    """
    Take a GeoSeries of Polygons along with a Series of values and returns randomly generated points within
    these polygons. Optionally takes a "points_per_value" integer which indicates the number of points that 
    should be generated for each 1 value.
    """
    if points_per_value:
        new_values = (values/points_per_value).astype(int)
    else:
        new_values = values
    new_values = new_values[new_values>0]
    g = gpd.GeoDataFrame(data = {'vals':new_values}, geometry = geometry)
    
    a = g.apply(lambda row: tuple(gen_random_points_poly(row['geometry'], row['vals'], seed)),1)
    b = gpd.GeoSeries(a.apply(pd.Series).stack(), crs = geometry.crs)
    b.name='geometry'
    return b

def blank_background():
    """
    Returns a dataframe which can be used to make a Choropleth which will add a
    blank layer to the map.
    """

    with open("data/external/tl_2017_us_state/tl_2017_us_state.geojson", 'r') as f:
        us_states = json.load(f)

    blank = pd.DataFrame(
        [
            {
                "state": s['properties']['NAME'],
                "value": 1,
            }
            for s in us_states['features']
        ]
    )    

    return blank, us_states

def blank_background_choropleth(m):
    """
    Applied a Choropleth layer which containing a blank layer to map m.
    """

    blank, us_states = blank_background()

    folium.Choropleth(
        us_states,
        blank,
        ('state', 'value'),
        key_on='feature.properties.NAME',
        fill_color='Greys',
        fill_opacity=1,
    ).add_to(m)

    folium.GeoJson(
        us_states,
        name="US States"
    ).add_to(m)

    return m

##

def get_zip_codes_metadata_geo():

    with open("data/external/ZIP_CODE_040114/ZIP_CODE_040114.geojson", 'r') as f:
        bronx_zip_codes = json.load(f)

    from covid_bronx.metadata import preprocess_metadata
    metadata = preprocess_metadata()
    coverage_levels = pd.read_csv("data/processed/sequencing/coverage.csv", index_col=0)['0']
    passed = coverage_levels[coverage_levels>=.95].index.map(lambda x: x.replace("_","-").replace("-0", "-").replace("-0", "-"))

    df = pd.read_csv("data/external/pangolin2.csv")
    df = df.loc[df['Sequence name'].dropna().index]
    df.index = df['Sequence name'].apply(lambda x: x.split(" ")[0]).apply(lambda x: x.replace("-0","-").replace("-0","-"))
    metadata[df.columns] = df
    zip_codes = set([x['properties']['ZIPCODE'] for x in bronx_zip_codes['features']])
    # Remove duplicates
    # zd = {x['properties']['ZIPCODE']: x for x in bronx_zip_codes['features']}
    # bronx_zip_codes['features'] = list(zd.values())

    # Keep only passing samples

    zip_codes_df = {k:v for k,v in metadata.groupby('zip_code')}

    zip_df = pd.DataFrame(
        [
            {
                "zip_code": str(z),
                "age": v['age'].mean(),
                "male": (v['sex']=='male').mean(),
                "female": (v['sex']=='female').mean(),
                "count": len(v.index),
                "num_side_effects": v['num_side_effects'].mean(),
                "num_comorbidities": v['num_comorbidities'].mean(),
                "death_rate": (v['outcome']=='deceased').mean(),
                "B.1": sum(v['Lineage']=='B.1'),
                "A.1": sum(v['Lineage']=='A.1'),
                "B": sum(v['Lineage']=='B'),
                "B.1.3": sum(v['Lineage']=='B.1.3'),
                "B.2": sum(v['Lineage']=='B.2'),
                "B.2.1": sum(v['Lineage']=='B.2.1'),
                "B.1": sum(v['Lineage']=='B.1'),
                "B.1.26": sum(v['Lineage']=='B.1.26'),
                "B.1.1": sum(v['Lineage']=='B.1.1'),
                "freq_A.1": (v['Lineage']=='A.1').sum() / sum(df['Lineage']=='A.1'),
                "freq_B": (v['Lineage']=='B').sum() / sum(df['Lineage']=='B'),
                "freq_B.1.3": (v['Lineage']=='B.1.3').sum() / sum(df['Lineage']=='B.1.3'),
                "freq_B.2": (v['Lineage']=='B.2').sum() / sum(df['Lineage']=='B.2'),
                "freq_B.2.1": (v['Lineage']=='B.2.1').sum() / sum(df['Lineage']=='B.2.1'),
                "freq_B.1": (v['Lineage']=='B.1').sum() / sum(df['Lineage']=='B.1'),
                "freq_B.1.26": (v['Lineage']=='B.1.26').sum() / sum(df['Lineage']=='B.1.26'),
                "freq_B.1.1": (v['Lineage']=='B.1.1').sum() / sum(df['Lineage']=='B.1.1'),
            }
            for z, v in zip_codes_df.items()
        ]
    )

    return zip_df, bronx_zip_codes