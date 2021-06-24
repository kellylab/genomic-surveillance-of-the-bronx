import numpy as np

lineage_colors_dict = {
    "A": "#987200",
    "A.1": "#987284",
    "A.3": "#9872ff",
    "B.1": "#385639",
    "B.1.1": "#2855d1",
    "B.1.1.1": "#dae7da",
    "B.1.3": "#d13328",
    "B.1.5": "#eac43f",
    "B.1.26": "#eac435",
    "B.2": "#f9b5ac",
    "B.2.1": "#ee7674",
}

lineage_colors_dict_rgb = {
    "A.1": np.array([152,114,132]),
    "B.1": np.array([56,86,57]),
    "B.1.1": np.array([40,85,209]),
    "B.1.1.1": np.array([218,231,218]),
    "B.1.3": np.array([209,51,40]),
    "B.1.26": np.array([234, 196, 53]),
    "B.2": np.array([249,181,172]),
    "B.2.1": np.array([238,118,116]),
}

lineage_colors_dict = {
    "A": "#E29E0C",
    "A.1": "#FC7000",
    "A.3": "#CEE20C",
    "B.1": "#226B01",
    "B.1.1": "#0BD7A0",
    "B.1.1.1": "#099BA3",
    "B.1.3": "#0969AB",
    "B.1.5": "#0B37CE",
    "B.1.26": "#6618D5",
    "B.2": "#B71878",
    "B.2.1": "#E20C1D",
}

lineage_colors_dict = {
    "A": "#E29E0C",
    "A.1": "#FC7000",
    "A.3": "#E2720C",
    "B.1": "#226B01",
    "B.1.1": "#03fcbe", # "#0BD7A0", # "#F8D51D",
    "B.1.1.1": "#099BA3",
    "B.1.3": "#0B37CE",
    "B.1.5": "#8231D0",
    "B.1.26": "#6618D5",
    "B.2": "#FE7AC5",
    "B.2.1": "#E20C1D",
}

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

lineage_colors_dict_rgb = {
    k: np.array(hex_to_rgb(v))
    for k,v in lineage_colors_dict.items()
    }