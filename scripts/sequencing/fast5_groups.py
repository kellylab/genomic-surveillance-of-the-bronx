# This script parses out which fast5 folders to use for each run.

import pandas as pd
from collections import defaultdict
import tqdm as tqdm

df = pd.read_csv('data/external/swabs.csv')
#groups = {k: set(v['sample_id']) for k,v in df.groupby('library')}
groups = {k: set(v['library']) for k,v in df.groupby('sample_id')}

# Invert dictionary

chunks = defaultdict(lambda: [])
for k,v in groups.items():
    chunks[tuple(v)].append(k)

chunks = dict(chunks)

folder_names = {
    1: 'gs://covid-bronx/SARS-2_Bx036-045_pos/SARS-2_Bx036-045_pos/20200506_1754_MN33931_FAN49865_b9eb023c',
    2: 'gs://covid-bronx/SARS-2_Bx001-010/SARS-2_Bx001-010/20200513_1745_MN33931_FAN49886_94626588/fast5',
    3: 'gs://covid-bronx/SARS-2-Bx051-060/SARS-2-Bx051-060/20200521_1715_MN33931_FAN50732_36a676a3/fast5',
    4: 'gs://covid-bronx/SARS-2_Bx036-045_Run4/SARS-2_Bx036-045_Run4/20200604_1555_MN33931_FAN54363_2f319835/fast5',
    5: 'gs://covid-bronx/run5/20200611_1834_MN33931__506766b1/fast5',
    6: 'gs://covid-bronx/run6/SARS-2_AECOM46-51_62-75_run6/20200623_2300_MN33931_FAN51083_4b6dc355/fast5',
    7: 'gs://covid-bronx/run7/SARS-2_AECOM-081-102_run7/20200626_1956_MN33931_FAN50814_dd816694/fast5',
    8: 'gs://covid-bronx/run8/20200702_1925_MN33931_FAN51071_2d4f788d/fast5',
    9: 'gs://covid-bronx/MSample1/S1/20200709_1911_MN33931_FAO16669_5176212e/fast5',
}

libraries = { # The folder name for each library. Some of them were ran multiple times, and we use the most recent.
    1: '1b',
    2: '2',
    3: '3c',
    4: '4',
    5: '5',
    6: '6',
    7: '7',
    8: '8',
    9: '9',
}

command = "gsutil -m cp -r {source} {destination} \n"

commands = [] # Make commands to consolidate fast5 files
for chunk in chunks.keys():
    destination = "gs://covid-bronx/fast5_" + "_".join(map(lambda x: str(x), chunk)) 
    for source in chunk:
        commands.append(command.format(source=folder_names[source], destination=destination))

with open("data/processed/sequencing/fast5_consolidate.sh", 'w') as f:
    f.writelines(commands)

commands = [] # Make commands to consolidate multiplexed reads
for chunk, ids in chunks.items():
    destination = "gs://covid-bronx/multiplexed_" + "_".join(map(lambda x: str(x), chunk)) + "/"
    for id in ids:
        commands.append(
            command.format(
                source=f"gs://covid-bronx/consolidated/multiplexed/{id}.fastq",
                destination=destination,
                )
            )

with open("data/processed/sequencing/multiplexed_consolidate.sh", 'w') as f:
    f.writelines(commands)

for chunk, ids in tqdm(chunks.items()):
    
    try:
        os.mkdir(os.path.join(drive_dir, "sequencing_summaries", "_".join(map(lambda x: str(x), chunk))))
    except:
        pass

    df = pd.DataFrame()
    for c in chunk:
        c_dir = os.path.join(drive_dir,'sequencing_summaries', f'run{libraries[c]}_basecalled_sequencing_summary.txt')
        print(f"Reading {c_dir}...")
        df = df.append(
            pd.read_csv(c_dir, sep="\t")
        )

    save_dir = os.path.join(drive_dir, 'sequencing_summaries', "_".join(map(lambda x: str(x), chunk)), 'sequencing_summary.txt')
    print(f"Saving to {save_dir}...")
    df.to_csv(save_dir, sep="\t", index=False)

commands = []
for chunk in chunks.keys():
    
    destination = "gs://covid-bronx/multiplexed_" + "_".join(map(lambda x: str(x), chunk)) + "/" + "sequencing_summary.txt"
    source = os.path.join(drive_dir, 'sequencing_summaries', "_".join(map(lambda x: str(x), chunk)), 'sequencing_summary.txt')
    commands.append(f"gsutil -m cp {source} {destination}\n")

with open("data/processed/sequencing/sequencing_summaries_consolidate.sh", "w") as f:
    f.writelines(commands)