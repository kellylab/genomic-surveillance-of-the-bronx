from skbio import read
import numpy as np
import pandas as pd

f = list(read("data/external/msa_outputs_alignment.msa", format="fasta"))

values = np.array([x.values for x in f])
positions = [i for i in range(len(values[0]))]
index = [x.metadata['id'] + x.metadata['description'] for x in f]

df = pd.DataFrame(values, index=index, columns=positions)
