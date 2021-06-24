# This script parses out which fast5 folders to use for each run.

import pandas as pd
from collections import defaultdict

command = "gsutil -m cp -R {source} {destination} \n"

commands = []

multiplex_mapping = {
    "gs://covid-bronx/reinfection2/multiplexed/barcode01_barcode01.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-103.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode02_barcode02.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-104.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode03_barcode03.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-105.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode04_barcode04.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-106.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode05_barcode05.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-107.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode06_barcode06.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-108.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode07_barcode07.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-109.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode08_barcode08.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-110.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode09_barcode09.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-112.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode10_barcode10.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-113.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode11_barcode11.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-114.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode12_barcode12.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-115.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode13_barcode13.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-116.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode14_barcode14.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-117.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode15_barcode15.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-118.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode16_barcode16.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-119.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode17_barcode17.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-120.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode18_barcode18.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-121.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode19_barcode19.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-122.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode20_barcode20.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-123.fastq",
    "gs://covid-bronx/reinfection2/multiplexed/barcode21_barcode21.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-D.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode01_barcode01.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-103.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode02_barcode02.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-104.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode03_barcode03.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-105.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode04_barcode04.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-106.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode05_barcode05.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-107.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode06_barcode06.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-108.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode07_barcode07.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-109.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode08_barcode08.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-110.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode09_barcode09.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-112.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode10_barcode10.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-113.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode11_barcode11.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-114.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode12_barcode12.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-115.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode13_barcode13.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-116.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode14_barcode14.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-117.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode15_barcode15.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-118.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode16_barcode16.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-119.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode17_barcode17.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-120.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode18_barcode18.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-121.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode19_barcode19.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-122.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode20_barcode20.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-123.fastq",
    "gs://covid-bronx/reinfection2b/multiplexed/barcode21_barcode21.fastq": "gs://covid-bronx/reinfection2x/multiplexed/AECOM-D.fastq",   
}

# Make Commands to combine fast5
commands.append(
    command.format(source="gs://covid-bronx/SARS-2_AECOM_103-122_D/20201021_2309_MN33931_FAO78801_cf550349/fast5", destination="gs://covid-bronx/reinfection2x/fast5")
)
commands.append(
    command.format(source="gs://covid-bronx/SARS-2_AECOM_103-122_D_run11/SARS-2_AECOM_103-122_D_run11/20201026_2218_MN33931__f3437483/fast5", destination="gs://covid-bronx/reinfection2x/fast5")
)

# Make commands to combine barcodes

for source, target in multiplex_mapping.items():
    commands.append(
        command.format(source=source, destination=target)
    )

df1 = pd.read_csv("data/final/reinfection2/multiplexed/sequencing_summary.txt", sep="\t")
df2 = pd.read_csv("data/final/reinfection2b/multiplexed/sequencing_summary.txt", sep="\t")

df = pd.concat([df1, df2])
s = "data/final/reinfection2x/multiplexed/sequencing_summary.txt"
df.to_csv(s, sep="\t")

commands.append(
    command.format(source=s, destination="gs://covid-bronx/reinfection2x/multiplexed/sequencing_summary.txt")
)
with open("artic/misc/reinfection.sh", "w") as f:
    f.writelines(commands)

