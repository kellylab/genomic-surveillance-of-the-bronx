from ete3 import Tree, ClusterTree, PhyloTree, TreeStyle, NodeStyle, AttrFace, ImgFace, faces
import pandas as pd

continents = {
 'Algeria': 'Africa',
 'Andorra': 'Europe',
 'Argentina': 'South America',
 'Aruba': 'South America',
 'Australia': 'Australasia',
 'Austria': 'Europe',
 'Bahrain': 'Middle East',
 'Bahrein': 'Middle East',
 'Bangladesh': 'Middle East',
 'Belarus': 'Europe',
 'Belgium': "Europe",
 'Belize': 'North America',
 'Benin': 'Africa',
 'BosniaandHerzegovina': 'Europe',
 'Botswana': 'Africa',
 'Brazil': 'South America',
 'Brunei': 'Asia',
 'Bulgaria': 'Europe',
 'Cambodia': 'Asia',
 'Canada': 'North America',
 'Chile': 'South America',
 'Colombia': 'South America',
 'Congo': 'Africa',
 'CostaRica': 'North America',
 'CotedIvoire': 'Africa',
 'Crimea': 'Europe',
 'Croatia': 'Europe',
 'Cuba': 'North America',
 'Curacao': 'North America',
 'Cyprus': 'Europe',
 'CzechRepublic': 'Europe',
 'DRC': 'Africa',
 'Denmark': 'Europe',
 'DominicanRepublic': 'North America',
 'Ecuador': 'South America',
 'Egypt': 'Africa',
 'England': 'Europe',
 'Estonia': 'Europe',
 'Finland': 'Europe',
 'France': 'Europe',
 'Gabon': 'Africa',
 'Gambia': 'Africa',
 'Georgia': 'Europe',
 'Germany': 'Europe',
 'Ghana': 'Africa',
 'Greece': 'Europe',
 'Guadeloupe': 'North America',
 'Guam': 'North America',
 'Guangdong': 'Asia',
 'Guatemala': 'North America',
 'HongKong': 'Asia',
 'Hungary': 'Europe',
 'Iceland': 'Europe',
 'India': 'Asia',
 'Indonesia': 'Asia',
 'Iran': 'Middle East',
 'Iraq': 'Middle East',
 'Ireland': 'Europe',
 'Israel': 'Middle East',
 'Italy': 'Europe',
 'Jamaica': 'North America',
 'Japan': 'Asia',
 'Jordan': 'Middle East',
 'Kazakhstan': 'Asia',
 'Kenya': 'Africa',
 'Kuwait': 'Middle East',
 'Latvia': 'Europe',
 'Lebanon': 'Middle East',
 'Liaoning': 'Asia',
 'Lithuania': 'Europe',
 'Luxembourg': 'Europe',
 'Madagascar': 'Africa',
 'Malaysia': 'Asia',
 'Mali': 'Africa',
 'Malta': 'Africa',
 'Mexico': 'North America',
 'Moldova': 'Europe',
 'Mongolia': 'Asia',
 'Montenegro': 'Africa',
 'Morocco': 'Africa',
 'Myanmar': 'Asia',
 'Nepal': 'Asia',
 'Netherlands': 'Europe',
 'NewZealand': 'Australasia',
 'Nigeria': 'Africa',
 'NorthMacedonia': 'Europe',
 'Norway': 'Europe',
 'Oman': 'Middle East',
 'Pakistan': 'Asia',
 'Palestine': 'Middle East',
 'Panama': 'North America',
 'Peru': 'South America',
 'Philippines': 'Asia',
 'Poland': 'Europe',
 'Portugal': 'Europe',
 'Romania': 'Europe',
 'Russia': 'Europe',
 'Rwanda': 'Africa',
 'SaintBarthelemy': 'North America',
 'SaintMartin': 'North America',
 'SaudiArabia': 'Middle East',
 'Scotland': 'Europe',
 'Senegal': 'Africa',
 'Serbia': 'Europe',
 'SierraLeone': 'Africa',
 'Singapore': 'Asia',
 'Slovakia': 'Europe',
 'Slovenia': 'Europe',
 'SouthAfrica': 'Africa',
 'SouthKorea': 'Asia',
 'Spain': 'Europe',
 'SriLanka': 'Asia',
 'Suriname': 'Africa',
 'Sweden': 'Europe',
 'Switzerland': 'Europe',
 'Taiwan': 'Asia',
 'Thailand': 'Asia',
 'Timor-Leste': 'Africa',
 'Tunisia': 'Africa',
 'Turkey': 'Middle East',
 'USA': 'North America',
 'Uganda': 'Africa',
 'Ukraine': 'Europe',
 'UnitedArabEmirates': 'Middle East',
 'Uruguay': 'South America',
 'Venezuela': 'North America',
 'Vietnam': 'Asia',
 'Wuhan': 'Asia',
 'Zambia': 'Africa',
 'cat': 'Non-Human',
 'env': 'Non-Human',
 'mink': 'Non-Human',
}

for k in continents:
    if continents[k] == 'Middle East':
        continents[k] = 'Asia' # Just use continents

continental_color_palette = {

}

t = PhyloTree("data/external/global_tree.nwk", format=1)

h = 2

base_style = NodeStyle()
base_style['hz_line_width'] = h
base_style['vt_line_width'] = h

styles = {
    'AECOM': NodeStyle(),
    'North America': NodeStyle(),
    'South America': NodeStyle(),
    'Europe': NodeStyle(),
    'Asia': NodeStyle(),
    'Australasia': NodeStyle(),
    'Africa': NodeStyle(),
    'Non-Human': NodeStyle(),
}

colors = {
    'AECOM': '#37dede',
    'North America': '#2945a3',
    'South America': '#31ad50',
    'Europe': '#99408f',
    'Asia': '#8a731a',
    'Australasia': '#6e8a3e',
    'Africa': '#8c5c0e',
    'Non-Human': '#8f071b', 
}

for k in list(styles.keys()):
    styles[k]['fgcolor'] = colors[k]
    styles[k]['bgcolor'] = '#ffffff'
    styles[k]['shape'] = 'circle'
    # styles[k]['draw_descendants'] = False
    styles[k]['size'] = 30
    styles[k]['hz_line_width'] = h
    styles[k]['vt_line_width'] = h

styles['AECOM']['fgcolor'] = colors['North America']
styles['AECOM']['bgcolor'] = colors['AECOM']
countries = []

leaves = [leaf.name for leaf in t]
for leaf in t:
    if "/" in leaf.name:
        if "NC_" in leaf.name:
            assert False
        country = leaf.name.split("/")[0] # Pull out just the country
        leaf.name = country
        countries.append(country)
        

countries = pd.Series(countries)

def get_node_layout(text=False):
    def node_layout(node):
        if node.is_leaf():
            if text:
                name_face = AttrFace("name", fsize=6)
                faces.add_face_to_node(name_face, node, column=0, position='branch-right')
            # node.set_style(nstyle)
            if 'AECOM' in node.name:
                node.set_style(styles['AECOM'])
            elif node.name in continents:                
                node.set_style(styles[continents[node.name]])
            else:
                node.set_style(base_style)
        else:
            node.set_style(base_style)
    
    return node_layout

ts = TreeStyle()
ts.mode = 'c'
ts.arc_start = 90
ts.arc_span = 180 # 184
ts.layout_fn = get_node_layout()
ts.allow_face_overlap = True
t.convert_to_ultrametric()
ts.show_leaf_name = False

from ete3 import TextFace, CircleFace
# Add Legend
for k, c in colors.items():
    ts.legend.add_face(TextFace(k + "     ", fsize=10), column=0)
    ts.legend.add_face(CircleFace(color=c, radius=5), column = 1)

ts.legend_position = 1

ts.mode = 'c'
ts.arc_span=360
ts.arc_start=0
ts.layout_fn = get_node_layout(text=False)
t.render("figures/figure2a_v2.pdf", w=180, units="mm", tree_style=ts)

t.render("figures/figure3.pdf", w=180, units="mm", tree_style=ts)
t.render("figures/figure3_flat.pdf", w=180, units="mm")
t.render("figures/figure3.svg", w=180, units="mm", tree_style=ts)
t.render("figures/figure3.png", w=180, units="mm", tree_style=ts)

ts.layout_fn = get_node_layout(text=True)
ts.arc_start=270
ts.arc_span=180
ts.mode='r'
endemic = t.get_common_ancestor("AECOM_027", "AECOM_118") # Could also try AECOM_049 or AECOM_107
endemic.render("figures/figure3_inset.pdf", w=180, units="mm", tree_style=ts)
t.render("figures/figure3_circulo.pdf", w=240, units="mm", tree_style=ts)


for k in list(styles.keys()):
    styles[k]['bgcolor'] = "#ffffff"
    # styles[k]['draw_descendants'] = False
    styles[k]['hz_line_width'] = h
    styles[k]['vt_line_width'] = h

styles['AECOM']['bgcolor'] = colors['AECOM']
t.render("figures/figure2a_v2_underlay.pdf", w=180, units="mm", tree_style=ts)