# Functionally annotating variants.
from typing import Tuple, Dict
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
import json

covid_gff = "data/external/covid_track.gff"

def parse_nextclade(json_file: str) -> Tuple[pd.DataFrame, pd.DataFrame]:

    with open(json_file, "r") as f:
        payload = json.load(f)

    rows = []
    variants = []
    for j in tqdm(payload):
        if len(j['errors']):
            continue

        rows.append({
            'seqName': j['seqName'],
            'totalMutations': j['totalMutations'],
            'totalInsertions': j['totalInsertions'],
            'totalMissing': j['totalMissing'],
            'alignmentStart': j['alignmentStart'],
            'alignmentEnd': j['alignmentEnd'],
            'alignmentScore': j['alignmentScore'],
            'totalPcrPrimerChanges': j['totalPcrPrimerChanges'],
            'clade': j['clade'],
            'qcScore': j['qc']['overallScore'],
            'privateMutationsScore': j['qc']['privateMutations']['score'],
            'missingDataScore': j['qc']['missingData']['score'],
            'snpClustersScore': j['qc']['snpClusters']['score'],
            'mixedSitesScore': j['qc']['mixedSites']['score'],
        })

        for v in j['substitutions']:
            var = {
                'seqName': j['seqName'],
                'position': v['pos'],
                'type': 'substitution',
                'mutation': f"{v['refNuc']}{v['pos']}{v['queryNuc']}",
                'gene': "",
                'codon': "",
                'aaMutation': "",
            }

            if 'aaSubstitutions' in v and len(v['aaSubstitutions']):

                for x in v['aaSubstitutions']:
                    var['gene'] = x['gene'],
                    var['codon'] = x['codon'],
                    var['aaMutation'] = f"{x['refAA']}{x['codon']}{x['queryAA']}",
                    variants.append(var.copy()) # This way we get duplicate mutations if they exist
            else:
                variants.append(var)

        for v in j['insertions']:
            var = {
                'seqName': j['seqName'],
                'position': v['pos'],
                'type': 'insertion',
                'mutation': f"{v['pos']}{v['ins']}"
            }

            variants.append(var)

        for v in j['deletions']:
            var = {
                'seqName': j['seqName'],
                'position': v['start'],
                'length': v['length'],
                'type': 'deletion',
            }

            variants.append(var)
    
    return pd.DataFrame(rows), pd.DataFrame(variants)