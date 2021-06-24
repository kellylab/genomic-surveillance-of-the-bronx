from ete3 import Tree, ClusterTree, PhyloTree, TreeStyle, NodeStyle, AttrFace, ImgFace, faces
from covid_bronx import lineage_colors_dict
import pandas as pd
from matplotlib import cm

from covid_bronx.metadata import get_metadata

colormap = cm.Purples
metadata = get_metadata()

def formatter(x):
    x = x.replace("-", "-00")
    x = x.split("-")
    x[1] = x[1][-3:]
    x = "-".join(x)
    return x

metadata['name'] = metadata['name'].apply(formatter)
delta = (metadata['collection_date'].max() - metadata['collection_date'].min()).days
start = metadata['collection_date'].min()
metadata['days'] = metadata['collection_date'].apply(lambda x: (x-start).days / delta)
metadata['age_color'] = metadata['days'].apply(lambda x: '#%02x%02x%02x' %  tuple(int(i*256) for i in colormap(x)[0:3]))

#t = PhyloTree("data/external/AECOM_complete_tree_raw.nwk")
#t = PhyloTree("data/external/timetree2.nwk")
t = PhyloTree("data/external/final_local_tree.nwk", format=1)
# t.unroot()
# df = pd.read_csv("data/external/pangolin2.csv", index_col=0)
# df = pd.read_csv("data/external/lineages_pangolin_1.csv")
df = pd.read_csv("data/external/lineages_final.csv", index_col=0)
# df.index = df['taxon'].apply(lambda x: str(x).lstrip("_").replace("_", "-"))

def formatter(x):

    y = str(x)
    y = y.replace("_", "-")
    pieces = y.split("-")
    if len(pieces[-1]) == 2:
        pieces[-1] = '0' + pieces[-1]
    return "-".join(pieces)

df.index = df['taxon'].apply(formatter)
lineages = df['lineage'].to_dict()

# df2 = pd.read_csv("data/processed/reinfection/pangolin.csv", index_col=0)
# df2.index = df2.index.map(lambda x: sample_id_metadata_mapping[x.split("/")[0]])
# lineages = {**lineages, **df2['Lineage'].to_dict()}

h = 2

nstyles = {}
for k, color in lineage_colors_dict.items():
    nstyles[k] = NodeStyle()
    nstyles[k]['bgcolor'] = color
    nstyles[k]['size'] = 3
    nstyles[k]['hz_line_width'] = h
    nstyles[k]['vt_line_width'] = h

space = "                 "
for leaf in t:
    if leaf.name.startswith("sample_"):
        leaf.name = leaf.name.lstrip("sample_")
    leaf.name = leaf.name.split("/")[0]

    leaf.name = leaf.name.replace("_", "-")
    # if leaf.name in sample_id_mapping:
    #     leaf.name = sample_id_mapping[leaf.name]

base_style = NodeStyle()
base_style['hz_line_width'] = h
base_style['vt_line_width'] = h

# def node_layout(node):
#     if node.is_leaf():
#         name_face = AttrFace("name", fsize=4)
#         faces.add_face_to_node(name_face, node, column=0, position='branch-right')
#         # node.set_style(nstyle)
#         if node.name in lineages and lineages[node.name] in nstyles:
#             node.set_style(nstyles[lineages[node.name]])
#         if node.name in ['AECOM-123', 'AECOM-124']: # https://commons.wikimedia.org/wiki/File:Full_Star_Yellow.svg
#             image_face = ImgFace("data/external/star.svg", width=15, height=15)
#             faces.add_face_to_node(image_face, node, column=1, position='branch-right')
#         if node.name in ['AECOM-125', 'AECOM-126']: # https://en.wikipedia.org/wiki/File:Diamond_warning_sign_(fluorescent_green).svg
#             image_face = ImgFace("data/external/diamond.svg", width=9, height=9)
#             faces.add_face_to_node(image_face, node, column=1, position='branch-right')
#         # if node.name in ['AECOM-127', 'AECOM-128', 'AECOM-129', 'AECOM-130', 'AECOM-131']: # <div>Icon made from <a href="http://www.onlinewebfonts.com/icon">Icon Fonts</a> is licensed by CC BY 3.0</div>
#         #     image_face = ImgFace("data/external/nose.svg", width=15, height=15)
#         #     faces.add_face_to_node(image_face, node, column=1, position='branch-right')
#             
#     else:
#         node.set_style(base_style)    
# 
# ts = TreeStyle()
# ts.show_leaf_name = False
# ts.mode = 'c'
# # ts.arc_start = 52
# # ts.arc_span = 220 # 184
# ts.layout_fn = node_layout
# ts.allow_face_overlap = True
# # t.convert_to_ultrametric()



# t.render("figures/figure2.pdf", w=180, units="mm", tree_style=ts)
# t.render("figures/figure2.svg", w=180, units="mm", tree_style=ts)

# Color by Age
def node_layout2(node):
    if node.is_leaf():
        name_face = AttrFace("name", fsize=4)
        faces.add_face_to_node(name_face, node, column=0, position='branch-right')
        # node.set_style(nstyle)
        node_style = NodeStyle()
        node_style['hz_line_width'] = h
        node_style['vt_line_width'] = h
        my_color = metadata.loc[metadata['name'] == node.name]['age_color'].iloc[0]
        node_style['bgcolor'] = my_color
        node.set_style(node_style)

        if node.name in ['AECOM-123', 'AECOM-124']: # https://commons.wikimedia.org/wiki/File:Full_Star_Yellow.svg
            image_face = ImgFace("data/external/star.svg", width=15, height=15)
            faces.add_face_to_node(image_face, node, column=1, position='branch-right')
        if node.name in ['AECOM-125', 'AECOM-126', 'AECOM-132']: # https://en.wikipedia.org/wiki/File:Diamond_warning_sign_(fluorescent_green).svg
            image_face = ImgFace("data/external/diamond.svg", width=9, height=9)
            faces.add_face_to_node(image_face, node, column=1, position='branch-right')
        # if node.name in ['AECOM-127', 'AECOM-128', 'AECOM-129', 'AECOM-130', 'AECOM-131']: # <div>Icon made from <a href="http://www.onlinewebfonts.com/icon">Icon Fonts</a> is licensed by CC BY 3.0</div>
        #     image_face = ImgFace("data/external/nose.svg", width=15, height=15)
        #     faces.add_face_to_node(image_face, node, column=1, position='branch-right')
            
    else:
        node.set_style(base_style)

ts2 = TreeStyle()
ts2.show_leaf_name = False
ts2.mode = 'c'
ts2.allow_face_overlap = True
ts2.layout_fn = node_layout2
from ete3 import TextFace, CircleFace
# Add Legend
for k, c in lineage_colors_dict.items():
    ts2.legend.add_face(TextFace(k + "     ", fsize=10), column=0)
    ts2.legend.add_face(CircleFace(color=c, radius=5), column = 1)

ts2.legend_position = 1

t.render("figures/figure2_age.pdf", w=180, units="mm", tree_style=ts2)
t.render("figures/figure2_age.svg", w=180, units="mm", tree_style=ts2)
    
# Subfigures
ss = TreeStyle()
ss.show_leaf_name = False
ss.layout_fn = node_layout
ss.allow_face_overlap = True
n = t.get_common_ancestor("AECOM-123", "AECOM-124")
# n.prune(["AECOM-105", "AECOM-113", "AECOM-090", "AECOM-032", "AECOM-075", "AECOM-106", "AECOM-105"])
n.convert_to_ultrametric()
n.render("figures/figure2b.pdf", w=90, units="mm", tree_style=ss)
n.render("figures/figure2b.svg", w=90, units="mm", tree_style=ss)

n2 = t.get_common_ancestor("AECOM-125", "AECOM-126")
n2.convert_to_ultrametric()

n2.prune([
    "AECOM-070",
    "AECOM-058",
    "AECOM-128",
    "AECOM-127",
    "AECOM-126",
    "AECOM-132",
    "AECOM-116",
    "AECOM-130",
    "AECOM-114",
    "AECOM-109",
    "AECOM-119",
    "AECOM-118",
    "AECOM-125",
    "AECOM-104",
    "AECOM-100",
    "AECOM-120",
    "AECOM-020",
    "AECOM-011",
    "AECOM-026",
    "AECOM-024",
    "AECOM-080",
    "AECOM-042",
    "AECOM-055",
    "AECOM-036",
    "AECOM-099",
    "AECOM-056",
    "AECOM-039",
    "AECOM-018",
    "AECOM-007",
    "AECOM-102",
])
n2.render("figures/figure2c.pdf", w=90, units="mm", tree_style=ss)
n2.render("figures/figure2c.svg", w=90, units="mm", tree_style=ss)