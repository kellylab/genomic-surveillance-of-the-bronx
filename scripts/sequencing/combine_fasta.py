import re
from covid_bronx.quality import fasta_files, sam_files

# sample_ids = {
#     "reinfection01": "sample_barcode01",
#     "reinfection02": "sample_barcode02",
#     "reinfection03": "sample_barcode03",
#     "reinfection04": "sample_barcode04",
#     "reinfection05": "sample_barcode05",
#     "reinfection06": "sample_barcode06",
#     "reinfection07": "sample_barcode07",
#     "reinfection08": "sample_barcode08",
#     "reinfection09": "sample_barcode09",
# }



# fasta_dict = {k: f"data/final/reinfection/output/{v}.consensus.fasta" for k,v
# in sample_ids.items()}

fasta_dict = {k: v for k,v in fasta_files.items()}
sample_ids = fasta_files
lines = []
for sample_id, filename in fasta_dict.items():
    b = sample_ids[sample_id]
    with open(filename, 'r') as f:
        line = f.read()
        zo = line.split("\n")
        zo[0] = f"> {sample_id}"

        # lines.append(f.read().replace(b, sample_id))
        lines.append("\n".join(zo))

with open("data/processed/reinfection/sequences.fasta", 'w') as f:
    f.write("\n".join(lines))

for sample_id, line in zip(fasta_dict.keys(), lines):
    with open(f"data/processed/sequences/{sample_id}.fasta", "w") as f:
        f.write(line)