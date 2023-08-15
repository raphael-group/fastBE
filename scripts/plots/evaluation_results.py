import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def load_data(csv_path):
    df = pd.read_csv(csv_path)
    return df

def set_alpha(ax, alpha):
    for patch in ax.patches:
        r, g, b, _ = patch.get_facecolor()
        patch.set_facecolor((r, g, b, alpha))

def main():
    parser = argparse.ArgumentParser(description='Plot false positive rate from a CSV file using boxplot.')
    parser.add_argument('csv_path', type=str, help='Path to the CSV file containing the DataFrame.')
    args = parser.parse_args()

    df = load_data(args.csv_path)

    fig, axes = plt.subplots(figsize=(8, 4), nrows=1, ncols=2, sharey=True)

    sns.boxplot(data=df, x='clones', y='false_positive_rate', hue='algorithm', fliersize=0, ax=axes[0])
    set_alpha(axes[0], 0.5)

    sns.stripplot(data=df, x='clones', y='false_positive_rate', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0])
    axes[0].set_xlabel('Number of Clones')
    axes[0].set_ylabel('False Positive Rate')
    axes[0].get_legend().remove()

    sns.boxplot(data=df, x='clones', y='false_negative_rate', hue='algorithm', fliersize=0, ax=axes[1])
    set_alpha(axes[1], 0.5)

    sns.stripplot(data=df, x='clones', y='false_negative_rate', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1])
    axes[1].set_xlabel('Number of Clones')
    axes[1].set_ylabel('False Negative Rate')
    axes[1].legend(title='Algorithm')

    handles, labels = axes[1].get_legend_handles_labels()
    labels = list(map(lambda x: {'pairtree': 'Pairtree', 'allele_minima': 'AlleleMinima'}[x], labels))
    axes[1].legend(handles, labels, title='Algorithm', loc='upper right')

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()

