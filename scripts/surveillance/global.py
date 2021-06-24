from covid_bronx import rarefaction
from tqdm import tqdm
import multiprocessing

values = {
    'Meta_North_America.csv': ('north_america/north_america', "North America"),
    'Meta_Europe.csv': ('europe/europe', 'Europe'),
    'Meta_South_America.csv': ('south_america/south_america', 'South America'),
    'Meta_Oceania.csv': ('oceania/oceania', 'Oceania'),
    'Meta_Africa.csv': ('africa/africa', 'Africa'),
    'Meta_Asia.csv': ('asia/asia', 'Asia'),
    }

# values = {
#     'NY.csv': ('ny/ny', 'New York State'),
#     'NYC.csv': ('nyc/nyc', 'New York City'),
#     'USA.csv': ('usa/usa', 'USA'),
# }

values = [(a, b[0], b[1]) for a,b in values.items()]

def func(x):
    csv = x[0]
    filename = x[1]
    title = x[2]
    rarefaction.plot_lineages("Bronx_5_12/" + csv, filename="gisaid_plots/" + filename, title=title, number_to_include=60)
    print(f"Done with {title}")

def func2(x):
    csv = x[0]
    filename = x[1]
    title = x[2]
    rarefaction.plot_variants("bronx_5_12/" + csv, json_folder='bronx_5_12/', filename="gisaid_plots/" + filename, title=title, number_to_include=60)
    print(f"Done with {title}")

p = multiprocessing.Pool(6)
p.map(func, values)
p.map(func2, values)