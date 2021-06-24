from matplotlib.markers import MarkerStyle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pysam import AlignmentFile
from skbio import Sequence
from typing import List
import click
import matplotlib.pyplot as plt

def parse_segments(sequence: np.array) -> List[np.array]:
    """
    Parses an array of integers into disjoint lists.
    """
    segments = []
    current = -1
    current_segment = []
    for i in sequence:
        if current >=0 and i > current+1: # New segment
            segments.append(np.array(current_segment))
            current_segment = [i]
        else:
            current_segment.append(i)
        
        current = i

    segments.append(current_segment)

    return segments

consensus_length = 29903
def compute_probabilities(consensus_sam, consensus_fasta):
    """
    Computes observations of enrichments at each site vs. whether or not the
    site was mapped.
    """
    # Compute 1 or 0 for gap at each site
    consensus_fasta = Sequence.read(consensus_fasta)
    consensus_gaps = [(x[0], x[-1]) for x in parse_segments(np.where(consensus_fasta.values==b'N')[0])]
    gap_sites = np.zeros(consensus_length)
    for g in consensus_gaps:
        gap_sites[g[0]:g[1]] = 1

    # Compute enrichment level at each site    
    consensus_sam = list(AlignmentFile(consensus_sam).fetch())
    positions = [x.positions for x in consensus_sam]
    enrichment = np.zeros(consensus_length)
    for p in positions:
        enrichment[p] += 1

    # Return results as vector
    mapped = enrichment[np.where(gap_sites==0.)[0]] 
    unmapped = enrichment[np.where(gap_sites==1.)[0]] 
    
    return mapped, unmapped, gap_sites, enrichment

@click.argument("out")
@click.argument("consensus-fasta")
@click.argument("consensus-sam")
@click.command()
def _compute_probabilities(consensus_sam, consensus_fasta, out):

    df = compute_probabilities(consensus_sam, consensus_fasta)
    df.to_csv(out + ".csv")
    print(f"Done {out}")

if __name__=="__main__":


    from covid_bronx.quality import sam_files, fasta_files

    stats = []
    
    for key in sam_files.keys():
        mapped, unmapped, gap_sites, enrichment = compute_probabilities(sam_files[key], fasta_files[key])
        stats.append({
            "run": key,
            "mapped": mapped,
            "unmapped": unmapped,
            "gaps": gap_sites,
            "enrichment": enrichment,
        })

    stats = pd.DataFrame(stats)
    stats.to_pickle("data/processed/ml/mapping_probabilities.pickle")

    fig, ax = plt.subplots(figsize=(12,12))
    for _, row in stats.iterrows(): 
        ax.hist(row['mapped'], density=True, bins=20, label=row['run'], alpha=1, histtype="step") 

    ax.set_title("Density of Mapped Region Enrichments by Run")
    ax.set_ylabel("Density")
    ax.set_xlabel("Normalized Enrichment")
    ax.legend()

    plt.savefig("data/processed/ml/mapped_density.svg")

    fig, ax = plt.subplots(figsize=(12,12))
    for _, row in stats.iterrows(): 
        ax.hist(row['unmapped'], density=True, bins=20, label=row['run'], alpha=1, histtype="step") 

    ax.set_title("Density of Mapped Region Enrichments by Run")
    ax.set_ylabel("Density")
    ax.set_xlabel("Normalized Enrichment")
    ax.legend()

    plt.savefig("data/processed/ml/unmapped_density.svg")
    
    all_mapped = np.concatenate(stats['mapped'])
    all_unmapped = np.concatenate(stats['unmapped'])

    fig, ax = plt.subplots(figsize=(12,12))
    ax.hist(all_mapped, density=True, bins=50, color="blue", alpha=.6, label="Mapped")
    ax.hist(all_unmapped, density=True, bins=50, color="red", alpha=.6, label="Unmapped")
    ax.set_title("Density of Mapped Region Enrichments by Run")
    ax.set_ylabel("Density")
    ax.set_xlabel("Normalized Enrichment")
    ax.legend()
    plt.savefig("data/processed/ml/mapped_vs_unmapped.svg")
    plt.show()

    fig, ax = plt.subplots(figsize=(12,12))
    gap_matrix = np.stack(stats['gaps'])
    enrichment_matrix = np.stack(stats['enrichment'])
    ax.fill_between([i for i in range(consensus_length)], enrichment_matrix.max()*np.mean(gap_matrix,0), color='red', label="Gaps", alpha=.6)
    ax.plot(np.mean(np.stack(enrichment_matrix),0), color='blue', label="Enrichment")
    ax.set_title("Mean Enrichment Across Runs") 
    ax.set_xlabel("Position on Reference Genome")
    ax.set_ylabel("Mean Absolute Enrichment")  
    ax.legend()
    plt.savefig("data/processed/ml/mean_enrichment.svg")
    plt.show()

    labels = np.concatenate([np.zeros(len(all_unmapped)),np.ones(len(all_mapped))]).reshape(-1,1)
    data = np.concatenate([all_unmapped, all_mapped]).reshape(-1,1)

    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import ShuffleSplit
    from sklearn.metrics import auc
    from sklearn.metrics import plot_roc_curve, plot_precision_recall_curve
    classifier = LogisticRegression()
    cv = ShuffleSplit(n_splits=10)

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    intercepts = []
    coefficients = []
    fig, ax = plt.subplots()
    for i, (train, test) in enumerate(cv.split(data, labels)):
        classifier.fit(data[train], labels[train])
        intercepts.append(classifier.intercept_)
        coefficients.append(classifier.coef_)
        viz = plot_roc_curve(classifier, data[test], labels[test],
                            name='ROC fold {}'.format(i),
                            alpha=0.3, lw=1, ax=ax)
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)

    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
            label='Chance', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr, color='b',
            label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
            lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
        title="Gap Prediction from Enrichment Level ROC Curve")
    ax.legend(loc="lower right")
    plt.savefig("data/processed/ml/gap_vs_enrichment_roc.svg")
    plt.show()

    fig, ax = plt.subplots()
    for i, (train, test) in enumerate(cv.split(data, labels)):
        classifier.fit(data[train], labels[train])
        intercepts.append(classifier.intercept_)
        coefficients.append(classifier.coef_)
        viz = plot_precision_recall_curve(classifier, data[test], labels[test],
                            name='PR fold {}'.format(i),
                            alpha=0.3, lw=1, ax=ax)
    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
        title="Gap Prediction from Enrichment Level ROC Curve")
    ax.legend(loc="lower right")
    plt.savefig("data/processed/ml/gap_vs_enrichment_prc.svg")
    plt.show()
    
    intercepts = np.concatenate(intercepts)
    coefficients = np.stack(coefficients).reshape(-1)
    d3 = pd.DataFrame([
        {"intercept": i, "coefficient": c}
        for i, c
        in zip(intercepts, coefficients)
    ])

    d3.to_csv("data/processed/ml/gap_cutoffs.csv")

    # Compute probability of assignment as a function of coverage.
    r = range(800)
    x = [i for i in r]
    y = classifier.predict_proba([[i] for i in r])
    plt.scatter(x, y[:,1])
    plt.title("Probability of Gap as a function of enrichment.")
    plt.savefig("data/processed/ml/gap_vs_enrichment_probabilities.svg")

    not_covered = enrichment_matrix[gap_matrix==1]
    covered = enrichment_matrix[gap_matrix==0]
    covered.sort()
    not_covered.sort()
    empirical_coverage_probabilities = np.array([np.sum(covered<i) / (np.sum(covered<i) + np.sum(not_covered<i)) for i in r])
    ec50 = np.where(empirical_coverage_probabilities<=.5)[0][-1]
    plt.axvline(x=list(r)[ec50], label=f"50% Probability at {list(r)[ec50]}x Coverage", marker='.')
    plt.axhline(y=.5, marker='.')
    plt.axhline(y=empirical_coverage_probabilities[20], label=f"20x Coverage Bascalled With {100*empirical_coverage_probabilities[20]:.2f}% Probability", marker='.', color='green')
    plt.axvline(x=20, marker='.', color='green')
    plt.scatter(r, empirical_coverage_probabilities)
    plt.title("Cumulative Density of Basecall Probability as a Function of Coverage")
    plt.xlabel("Coverage")
    plt.ylabel("Cumulative Probability of Basecall if Coverage is <= This Value")
    plt.legend()
    plt.savefig("data/processed/ml/gap_vs_enrichment_cdf.svg", figsize=(40,40))
    plt.show()

    all_enrichments = enrichment_matrix.flatten()
    plt.hist(all_enrichments, cumulative=True, density=True, bins=200, alpha=.5)
    twentyx = np.mean(all_enrichments <= 20)
    plt.title("Cumulative Density of Coverage Levels")
    plt.xlabel("Coverage Level")
    plt.ylabel("Percentage of Positions With Coverage <= This Value")
    plt.axhline(y=twentyx, label=f"{twentyx*100:.2f}% of Samples Have Coverage <= 20x" ,color='green')
    plt.axvline(x=20, color='green')
    plt.legend()
    plt.savefig("data/processed/ml/enrichment_cdf.pdf", figsize=(30,30))
    plt.show()

    above_empirical_coverage_probabilities = np.array([np.sum(covered>i) / (np.sum(covered>i) + np.sum(not_covered>i)) for i in r])
    plt.scatter(r, above_empirical_coverage_probabilities)
    plt.title("Cumulative Density of Basecall Probability as a Function of Coverage")
    plt.xlabel("Coverage")
    plt.ylabel("Cumulative Probability of Basecall if Coverage is >= This Value")
    plt.savefig("data/processed/ml/gap_vs_enrichment_cdf_above.pdf", figsize=(40,40))
    plt.show()
    