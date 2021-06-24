# Break Down the Gisaid genomes into smaller chunks for processing.

from skbio import read, write
from tqdm import tqdm

fastas = list(tqdm(read("sequences.fasta", format="fasta")))

chunk_size = 10000

current = 0
l = len(fastas)
index = 0

while current < l:

    print(f"Progress: {current / l * 100:.1f}%")
    chunk = fastas[current:min(current+chunk_size,l)]
    with open(f"gisaid/sequences_{index}.fasta", "w") as f:
        for c in chunk:
            c.write(f)
    current += chunk_size
    index += 1
