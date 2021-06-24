from skbio import read
import os
import numpy as np
from typing import Dict
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from pysam import AlignmentFile, VariantFile
from tqdm import tqdm
from covid_bronx.quality import sam_files, fasta_files, variant_files
import skbio

def count_it(l): 
    l = map(lambda x: x.upper(), l) 
    counts = defaultdict(lambda : 0) 
    for x in l: 
        counts[x] += 1 
    return dict(counts)


df_dict = {}
samples = ['AECOM_123', 'AECOM_124', 'AECOM_125', 'AECOM_126']

for sample_id in tqdm(samples):
    filename = sam_files[sample_id]
    bamfile = AlignmentFile(filename)
    pileup = bamfile.pileup()
    df = pd.DataFrame([{"pos": p.pos, "n": p.n, "counts": count_it(p.get_query_sequences())} for p in pileup])
    df.index = df['pos']
    df.drop(columns=['pos'])
    df_dict[sample_id] = df


vardf_all = []
for sample_id in tqdm(samples):
    # Get reads to assess depth
    sam_filename = sam_files[sample_id]
    fasta_filename = fasta_files[sample_id]
    variant_filename = variant_files[sample_id]

    alignments = AlignmentFile(sam_filename).fetch()
    consensus = list(skbio.read(fasta_filename, format="fasta"))[0]

    coverage = np.zeros(29903)
    for alignment in alignments:
        coverage[alignment.positions] += 1

    # Get variants
    variants = VariantFile(variant_filename)
    vardf_all.extend([
        {
            **{
                key: value
                for key, value
                in var.info.items()
                },
            **{
                "sample_id": sample_id,
                "position": var.pos,
                "quality": var.qual,
                "reference": var.ref,
                "alternates": var.alts,
                "depth": coverage[var.pos],
                },
        }
        for var in variants.fetch()
    ])

vardf_all = pd.DataFrame(vardf_all).set_index("sample_id")
vardf_all['score'] = vardf_all['quality'] / vardf_all['depth']