# This script parses out the read ids for each of the runs.
import os
from tqdm import tqdm

drive_dir = "/media/saad/Samsung_T5"

nums = {int(x[6:9]):x for x in os.listdir(os.path.join(drive_dir, 'consolidated/multiplexed'))}
for run_id, multiplexed in tqdm(nums.items()):
    # Get read_ids for that sample
    with open(os.path.join(drive_dir, 'consolidated/multiplexed', multiplexed), 'r') as f:
        read_ids = [x.split(' ')[0][1:] for x in f.readlines() if 'runid' in x]
    with open(os.path.join(drive_dir, f'consolidated/read_ids/{run_id}'), 'w') as f:
        f.write("\n".join(read_ids))
    