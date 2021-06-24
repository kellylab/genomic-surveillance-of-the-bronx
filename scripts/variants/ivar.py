# Parses out iVar results

import pandas as pd
from tqdm import tqdm
import os
import plot_helper
from plot_helper import Plotter
import numpy as np
import matplotlib.pyplot as plt

base_dir = "data/external/ivar"

df_dict = {
    run.split(".")[0]: pd.read_csv(
        os.path.join(base_dir, run),
        sep="\t",
    )
    for run in os.listdir(base_dir)
}

df = pd.DataFrame()

for k,v in tqdm(df_dict.items()):
    v['run'] = k
    df = df.append(v)

df = df[df['PASS']]
positions = {p: v for p,v in tqdm(df.groupby('POS'))}
counts = pd.Series({p: len(v) for p,v in positions.items()}).sort_values(ascending=False)

X = {k: i for i, k in enumerate(set(df['POS']))}
Y = {k: i for i, k in enumerate(set(df['run']))}

matrix = np.zeros((len(X),len(Y)))

for _, row in tqdm(df.iterrows()):
    x = row['POS']
    y = row['run']
    matrix[X[x], Y[y]] = 1

dm = pd.DataFrame(matrix, index=X.keys(), columns=Y.keys())

sns.clustermap(dm, figsize=(20,20))

plt.savefig("data/processed/variants/coocurrence.png")
