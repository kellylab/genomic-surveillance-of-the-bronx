# This script combines barcoded reads from different runs which correspond to
# the same sample ID in order to improve coverage.

import pandas as pd
import os
from tqdm import tqdm

df = pd.read_csv("data/external/swabs.csv")

samples = {
    s: df[df['sample_id'] == s]
    for s in set(df['sample_id'])
}

multiple = {k:v for k,v in samples.items() if len(v)>1}

libraries = { # The folder name for each library. Some of them were ran multiple times, and we use the most recent.
    1: '1b',
    2: '2',
    3: '3c',
    4: '4',
    5: '5c',
    6: '6',
    7: '7',
    8: '8',
    9: '9',
}

barcodes = { # The string representation of each barcode number.
    1: "01",
    2: "02",
    3: "03",
    4: "04",
    5: "05",
    6: "06",
    7: "07",
    8: "08",
    9: "09",
    10: "10",
    11: "11",
    12: "12",
    13: "13",
    14: "14",
    15: "15",
    16: "16",
    17: "17",
    18: "18",
    19: "19",
    20: "20",
    21: "21",
    22: "22",
    23: "23",
    24: "24",
}

# Consolidate multiplexing
drive_dir = "/media/saad/Samsung_T5"

for sample_id, run_info in tqdm(samples.items()):

    files = [
        os.path.join(
            drive_dir,
            f"run{libraries[run['library']]}",
            'multiplexed',
            f"barcode{barcodes[run['barcode']]}_barcode{barcodes[run['barcode']]}.fastq",
        )
        for _, run in run_info.iterrows()
    ]
    with open(os.path.join(drive_dir, "consolidated/multiplexed", f"{sample_id}.fastq"), 'w') as f:
        for file in files:
            with open(file, 'r') as g:
                f.write(g.read())

