from covid_bronx.quality import fasta_files, sam_files

lines = []
for sample_id, filename in fasta_files.items():
    with open(filename, 'r') as f:
        lines.append(f.read())

with open("data/processed/sequencing/sequences.fasta", 'w') as f:
    f.write("\n".join(lines))